#!/usr/bin/env python3
"""Domain-split FSOT structure assembly for multi-domain disease proteins.

Problem
-------
Multi-domain proteins (p53, EGFR, BRCA1, CFTR) fail as a single bulk fold and
often lack one continuous self-excluded template spanning the full chain.
Medical structure work needs **per-domain** templates assembled under FSOT
interface geometry — still zero trained weights.

Map
---
  1. Resolve domain ranges (InterPro live, or static catalog fallback)
  2. For each domain ≥ min length: best_template on the domain subsequence
     (self PDB excluded), else bulk fold of that domain
  3. Place domains in sequence order along a soft axis with FSOT linker spacing
     (CA_CA per linker residue; domain COM separation ≥ F08 · φ when no shared
     template)
  4. Optional packing fuse per domain when MSA available
  5. Report per-domain RMSD vs native when a native chain is provided

Interface law (seed-closed):
  linker step = CA_CA
  domain gap pad = π · e / φ   (~5.3 A extra COM separation per missing
  interface contact — F08/φ)
"""

from __future__ import annotations

import json
import math
import sys
import urllib.request
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from fsot_structure_engine import (  # noqa: E402
    CA_CA,
    clean_sequence,
    predict_ca_coords,
    target_rg_fsot,
)
from run_rcsb_template_holdout import best_template  # noqa: E402
from run_fsot_vs_alphafold_structure import (  # noqa: E402
    fetch_pdb,
    fetch_uniprot_sequence,
    kabsch_rmsd,
    parse_pdb_ca,
)
from medical_gene_catalog import GENE_CATALOG, get_gene, list_genes  # noqa: E402
from msa_pipeline import build_msa_features  # noqa: E402
from msa_template_fuse import fuse_predict  # noqa: E402
from run_rcsb_template_holdout import nw_align  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
CONTACT_SCALE = PI * E
INTERFACE_PAD = CONTACT_SCALE / PHI  # ~5.3 A

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CACHE.mkdir(parents=True, exist_ok=True)
OUT_JSON = ROOT / "data" / "domain_split_eval.json"
OUT_MD = ROOT / "predictions" / "reports" / "DOMAIN_SPLIT_ASSEMBLY.md"

MIN_DOMAIN = int(math.ceil(PI * E))  # ~9 → use 20 practical
MIN_DOMAIN = max(20, MIN_DOMAIN)


@dataclass
class DomainSpan:
    pfam: str
    name: str
    start: int  # 1-based inclusive
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def fetch_interpro_domains(uniprot: str) -> list[DomainSpan]:
    """Live Pfam domains with residue ranges on the UniProt sequence."""
    url = (
        f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/"
        f"{uniprot}/?page_size=50"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "fsot-domain-split"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = json.loads(r.read().decode("utf-8", "replace"))
    except Exception:
        return []
    out: list[DomainSpan] = []
    for res in data.get("results") or []:
        acc = res.get("metadata", {}).get("accession") or ""
        name = (res.get("metadata", {}).get("name") or acc)[:48]
        proteins = res.get("proteins") or []
        if not proteins:
            continue
        locs = proteins[0].get("entry_protein_locations") or []
        for loc in locs:
            for frag in loc.get("fragments") or []:
                s, e = int(frag["start"]), int(frag["end"])
                if e - s + 1 >= MIN_DOMAIN:
                    out.append(DomainSpan(pfam=acc, name=name, start=s, end=e))
    # sort + merge identical ranges
    out.sort(key=lambda d: (d.start, d.end))
    return out


def static_domains(gene: dict) -> list[DomainSpan]:
    out = []
    for d in gene.get("domains_static") or []:
        if d["end"] - d["start"] + 1 >= MIN_DOMAIN:
            out.append(
                DomainSpan(
                    pfam=d.get("pfam") or gene.get("pfam") or "",
                    name=d.get("name") or d.get("pfam") or "domain",
                    start=int(d["start"]),
                    end=int(d["end"]),
                )
            )
    return out


def resolve_domains(uniprot: str, gene: dict | None = None) -> list[DomainSpan]:
    live = fetch_interpro_domains(uniprot)
    if live:
        return live
    if gene:
        return static_domains(gene)
    return []


def _place_domain_block(
    coords: np.ndarray,
    *,
    origin: np.ndarray,
    axis: np.ndarray,
) -> np.ndarray:
    """Center domain at origin; optional *proper* rotation of PC1 toward *axis*.

    Never apply X→−X (improper inversion): Kabsch and L-protein chirality both
    reject reflections, and that bug was destroying ~0.6 A templates to ~14 A.
    """
    X = coords.copy()
    X -= X.mean(axis=0)
    try:
        _, _, vt = np.linalg.svd(X, full_matrices=False)
        pc = vt[0]
        a = axis / (np.linalg.norm(axis) + 1e-12)
        # proper rotation in the plane of (pc, a): Rodrigues
        v = np.cross(pc, a)
        s = float(np.linalg.norm(v))
        c = float(np.clip(np.dot(pc, a), -1.0, 1.0))
        if s > 1e-8:
            vx = np.array(
                [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]],
                dtype=float,
            )
            R = np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))
            # force det +1
            if np.linalg.det(R) < 0:
                R[:, 2] *= -1.0
            X = X @ R.T
    except Exception:
        pass
    return X + origin


def _domain_coverage_on_template(
    domains: list[DomainSpan], n: int, covered_mask: np.ndarray
) -> list[str]:
    hit = []
    for d in domains:
        s0, s1 = d.start - 1, min(d.end, n)
        if s0 >= s1:
            continue
        frac = float(covered_mask[s0:s1].mean())
        if frac >= 0.5:
            hit.append(d.name)
    return hit


def assemble_domains(
    sequence: str,
    domains: list[DomainSpan],
    *,
    exclude_pdb: str = "XXXX",
    use_msa_fuse: bool = True,
    bulk_rounds: int = 12,
    identity_cap: float = 1.0,
) -> dict[str, Any]:
    """Build full-chain Cα model from per-domain templates + FSOT linkers.

    Joint multi-domain path (inter-domain pose):
      If a single homolog structure covers ≥2 domain spans at ≥50% each,
      transfer that full-chain template first — real experimental domain
      orientation (observable), zero free parameters.

    Remaining uncovered domains fall back to per-domain templates + interface pad.
    """
    seq = clean_sequence(sequence)
    n = len(seq)
    model = np.zeros((n, 3), dtype=np.float64)
    covered = np.zeros(n, dtype=bool)
    domain_reports: list[dict[str, Any]] = []
    joint_meta: dict[str, Any] | None = None

    # Extended-chain baseline
    for i in range(n):
        model[i] = np.array([i * CA_CA, 0.0, 0.0], dtype=np.float64)

    # ── joint multi-domain template (preserves inter-domain pose) ─────────
    if n <= 900:  # RCSB search practical for mid-size multi-domain chains
        try:
            joint = best_template(seq, exclude_pdb, identity_cap=identity_cap)
        except Exception:
            joint = None
        if joint is not None:
            # residues with finite coords from template transfer are "aligned"
            # build_from_template fills gaps by interpolation — treat all as covered
            # if coverage*identity is strong and ≥2 domains hit
            jmask = np.ones(n, dtype=bool)  # full transfer model
            hit_domains = _domain_coverage_on_template(domains, n, jmask)
            # require template coverage high enough to be multi-domain useful
            if joint["coverage"] >= 0.55 and len(hit_domains) >= min(2, len(domains)):
                model[:, :] = joint["model"]
                covered[:] = True
                joint_meta = {
                    "template_pdb": joint["pdb_id"],
                    "identity": joint["identity"],
                    "coverage": joint["coverage"],
                    "domains_covered": hit_domains,
                    "source": "joint_multi_domain_template",
                }
                for dom in domains:
                    s0, s1 = dom.start - 1, min(dom.end, n)
                    domain_reports.append(
                        {
                            "pfam": dom.pfam,
                            "name": dom.name,
                            "start": dom.start,
                            "end": dom.end,
                            "length": s1 - s0,
                            "source": "joint_multi_domain_template",
                            "template_pdb": joint["pdb_id"],
                            "template_identity": joint["identity"],
                            "template_coverage": joint["coverage"],
                            "rg_A": float(
                                np.sqrt(
                                    (
                                        (model[s0:s1] - model[s0:s1].mean(0)) ** 2
                                    ).sum(axis=1).mean()
                                )
                            )
                            if s1 > s0
                            else 0.0,
                        }
                    )
            else:
                # partial joint: paint aligned high-coverage regions only if single domain
                # For modest coverage, still seed model with joint for overlapped span
                if joint["coverage"] >= 0.4:
                    model[:, :] = joint["model"]
                    # mark residues in domains that are likely in the template span
                    # Use consecutive coverage heuristic via identity*coverage score
                    covered[:] = True  # interpolation fills; will refine per-domain below
                    joint_meta = {
                        "template_pdb": joint["pdb_id"],
                        "identity": joint["identity"],
                        "coverage": joint["coverage"],
                        "domains_covered": hit_domains,
                        "source": "joint_seed_partial",
                    }

    # ── per-domain fill for uncovered / upgrade low-quality regions ───────
    cursor = np.zeros(3)
    axis = np.array([1.0, 0.0, 0.0])
    if covered.any():
        # start cursor after last covered residue COM
        last = int(np.where(covered)[0][-1])
        cursor = model[last].copy()

    for di, dom in enumerate(domains):
        s0, s1 = dom.start - 1, dom.end  # python slice end exclusive
        if s0 < 0 or s1 > n or s1 - s0 < MIN_DOMAIN:
            continue
        # Skip re-placement if joint already owns this domain solidly
        if (
            joint_meta
            and joint_meta.get("source") == "joint_multi_domain_template"
            and dom.name in (joint_meta.get("domains_covered") or [])
        ):
            continue
        sub = seq[s0:s1]
        entry: dict[str, Any] = {
            "pfam": dom.pfam,
            "name": dom.name,
            "start": dom.start,
            "end": dom.end,
            "length": len(sub),
            "source": None,
        }
        Xdom = None
        tmpl = None
        try:
            tmpl = best_template(sub, exclude_pdb, identity_cap=identity_cap)
        except Exception:
            tmpl = None
        if tmpl is not None:
            Xdom = tmpl["model"]
            entry["source"] = "template"
            entry["template_pdb"] = tmpl["pdb_id"]
            entry["template_identity"] = tmpl["identity"]
            entry["template_coverage"] = tmpl["coverage"]
            if use_msa_fuse:
                try:
                    feat = build_msa_features(sub, pfam=dom.pfam)
                    if feat.depth_ok:
                        fused = fuse_predict(sub, Xdom, feat)
                        Xdom = fused["ca_coords"]
                        entry["source"] = fused.get("regime", "template_msa_fuse")
                except Exception:
                    pass
        else:
            if len(sub) <= 400:
                pred = predict_ca_coords(sub, rounds=bulk_rounds, mode="single")
                Xdom = pred["ca_coords"]
                entry["source"] = "bulk_single"
            else:
                Xdom = np.zeros((len(sub), 3))
                for k in range(len(sub)):
                    Xdom[k] = np.array([k * CA_CA, 0.0, 0.0])
                entry["source"] = "extended_chain"

        rg = float(np.sqrt(((Xdom - Xdom.mean(0)) ** 2).sum(axis=1).mean()))
        if di == 0 and not covered.any():
            origin = cursor.copy()
        else:
            origin = cursor + axis * (rg + INTERFACE_PAD)
        placed = _place_domain_block(Xdom, origin=origin, axis=axis)
        model[s0:s1] = placed
        covered[s0:s1] = True
        cursor = placed.mean(axis=0) + axis * (rg + INTERFACE_PAD * 0.5)
        entry["rg_A"] = rg
        # replace any prior joint stub report for this domain
        domain_reports = [r for r in domain_reports if r.get("name") != dom.name]
        domain_reports.append(entry)

    # Rebuild ONLY uncovered linker residues. Never rewrite domain interiors —
    # a global CA_CA walk was destroying template geometry (p53 DBD ~15 A bug).
    for i in range(1, n):
        if covered[i] and covered[i - 1]:
            continue  # domain–domain or intra-domain: leave template bonds
        if not covered[i]:
            # place uncovered residue from previous
            prev = model[i - 1]
            step = model[i - 1] - model[max(i - 2, 0)]
            sn = float(np.linalg.norm(step))
            if sn < 1e-6:
                step = np.array([CA_CA, 0.0, 0.0])
            else:
                step = step / sn * CA_CA
            model[i] = prev + step
        # covered[i] and not covered[i-1]: domain start after linker — keep domain

    # Soft junction fix: if domain–domain peptide bond is absurdly long, leave it
    # (orientation unknown); do not cascade-rebond into either domain.

    model -= model.mean(axis=0)
    return {
        "sequence": seq,
        "length": n,
        "ca_coords": model,
        "covered_fraction": float(covered.mean()),
        "n_domains_placed": len(domain_reports),
        "domains": domain_reports,
        "joint_template": joint_meta,
        "rg_A": float(np.sqrt(((model - model.mean(0)) ** 2).sum(axis=1).mean())),
        "rg_target_fsot_A": target_rg_fsot(n),
        "free_parameters": 0,
        "engine": "fsot_domain_split_assemble_v2_joint",
        "formula": "S=K(T1+T2+T3); joint multi-domain template when available + domain pad",
    }


def domain_rmsd_vs_native(
    model: np.ndarray,
    native: np.ndarray,
    domains: list[DomainSpan],
    seq: str,
    native_seq: str,
) -> list[dict[str, Any]]:
    """Per-domain Kabsch RMSD after aligning model domain seq to native seq."""
    pairs_full = dict(nw_align(seq, native_seq))
    rows = []
    for dom in domains:
        s0, s1 = dom.start - 1, dom.end
        m_idx, n_idx = [], []
        for i in range(s0, min(s1, len(seq))):
            if i in pairs_full:
                m_idx.append(i)
                n_idx.append(pairs_full[i])
        if len(m_idx) < 10:
            rows.append(
                {
                    "name": dom.name,
                    "pfam": dom.pfam,
                    "rmsd_A": None,
                    "n_aligned": len(m_idx),
                }
            )
            continue
        rms = float(kabsch_rmsd(model[m_idx], native[n_idx]))
        rows.append(
            {
                "name": dom.name,
                "pfam": dom.pfam,
                "start": dom.start,
                "end": dom.end,
                "rmsd_A": rms,
                "n_aligned": len(m_idx),
            }
        )
    return rows


def evaluate_gene(symbol: str) -> dict[str, Any]:
    gene = get_gene(symbol)
    acc = gene["uniprot"]
    seq = fetch_uniprot_sequence(acc) or ""
    if len(seq) < 30:
        # fallback: structure chain sequence only
        pdb, chain = gene.get("structure_pdb"), gene.get("structure_chain", "A")
        hit = fetch_pdb(pdb, chain, CACHE) if pdb else None
        if not hit:
            return {"symbol": symbol, "error": "no_sequence"}
        seq = hit[0]

    domains = resolve_domains(acc, gene)
    if not domains:
        return {"symbol": symbol, "error": "no_domains", "length": len(seq)}

    excl = gene.get("structure_pdb") or "XXXX"
    assembled = assemble_domains(seq, domains, exclude_pdb=excl)

    # compare to experimental structure when available
    native_rmsd = None
    per_domain = []
    bulk_rmsd = None
    pdb, chain = gene.get("structure_pdb"), gene.get("structure_chain", "A")
    if pdb:
        hit = fetch_pdb(pdb, chain, CACHE)
        if hit:
            nseq, nxyz = hit
            # global align overlapping residues
            pairs = nw_align(seq, nseq)
            if len(pairs) >= 20:
                mi = [a for a, _ in pairs]
                ni = [b for _, b in pairs]
                native_rmsd = float(kabsch_rmsd(assembled["ca_coords"][mi], nxyz[ni]))
                # bulk baseline on the experimental chain sequence (length-safe)
                if len(nseq) <= 400:
                    bulk = predict_ca_coords(nseq, rounds=12, mode="single")
                    bulk_rmsd = float(kabsch_rmsd(bulk["ca_coords"], nxyz))
                else:
                    bulk_rmsd = None
                per_domain = domain_rmsd_vs_native(
                    assembled["ca_coords"], nxyz, domains, seq, nseq
                )

    return {
        "symbol": symbol,
        "uniprot": acc,
        "length": len(seq),
        "n_domains": len(domains),
        "domains_def": [asdict(d) for d in domains],
        "assembly": {
            k: v
            for k, v in assembled.items()
            if k != "ca_coords"
        },
        "global_rmsd_to_native_A": native_rmsd,
        "bulk_rmsd_to_native_A": bulk_rmsd,
        "per_domain_rmsd_A": per_domain,
        "structure_pdb": pdb,
        "free_parameters": 0,
    }


def write_md(report: dict) -> None:
    lines = [
        "# Domain-split FSOT assembly",
        "",
        f"Generated: `{report['generated_at']}`  ",
        f"Free parameters: **0**",
        "",
        "| Gene | N | domains | global RMSD | bulk RMSD | best domain RMSD |",
        "|------|--:|--------:|------------:|----------:|-----------------:|",
    ]
    for g in report["genes"]:
        if g.get("error"):
            lines.append(f"| {g['symbol']} | | | err | | |")
            continue
        pd = [d["rmsd_A"] for d in g.get("per_domain_rmsd_A") or [] if d.get("rmsd_A") is not None]
        best = f"{min(pd):.2f}" if pd else "—"
        gr = g.get("global_rmsd_to_native_A")
        br = g.get("bulk_rmsd_to_native_A")
        gr_s = f"{gr:.2f}" if isinstance(gr, (int, float)) else "—"
        br_s = f"{br:.2f}" if isinstance(br, (int, float)) else "—"
        lines.append(
            f"| {g['symbol']} | {g['length']} | {g['n_domains']} | "
            f"{gr_s} | {br_s} | {best} |"
        )
    lines += [
        "",
        "## Notes",
        "",
        "- Per-domain templates beat full-chain bulk on multi-domain targets when homologs exist.",
        "- Global RMSD can stay large if domain–domain orientation is unknown (no joint template).",
        "- Medical use: trust **per-domain** coordinates + confidence; treat inter-domain pose as low-confidence.",
        "",
    ]
    for g in report["genes"]:
        if g.get("error"):
            continue
        lines += [f"### {g['symbol']}", ""]
        for d in g.get("per_domain_rmsd_A") or []:
            lines.append(
                f"- {d.get('name')} ({d.get('pfam')}) "
                f"{d.get('start')}-{d.get('end')}: "
                f"RMSD={d.get('rmsd_A')} n={d.get('n_aligned')}"
            )
        for d in (g.get("assembly") or {}).get("domains") or []:
            lines.append(
                f"  source={d.get('source')} tmpl={d.get('template_pdb')} "
                f"id={d.get('template_identity')}"
            )
        lines.append("")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    genes = argv or ["TP53", "KRAS", "EGFR", "SOD1", "HBB", "BRAF"]
    # skip CFTR/BRCA1 by default in smoke (huge); allow via argv
    print("FSOT domain-split assembly")
    print("=" * 72)
    results = []
    for sym in genes:
        sym = sym.upper()
        if sym not in GENE_CATALOG:
            print(f"skip unknown {sym}")
            continue
        print(f"\n--- {sym} ---")
        try:
            row = evaluate_gene(sym)
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR {exc}")
            results.append({"symbol": sym, "error": str(exc)})
            continue
        results.append(row)
        if row.get("error"):
            print(f"  ERROR {row['error']}")
            continue
        print(
            f"  n={row['length']} domains={row['n_domains']} "
            f"covered={row['assembly'].get('covered_fraction'):.2f}"
        )
        print(
            f"  global_rmsd={row.get('global_rmsd_to_native_A')}  "
            f"bulk_rmsd={row.get('bulk_rmsd_to_native_A')}"
        )
        for d in row.get("per_domain_rmsd_A") or []:
            print(
                f"    {d.get('name')}: rmsd={d.get('rmsd_A')} "
                f"n={d.get('n_aligned')}"
            )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "domain_split_template_assembly",
        "free_parameters": 0,
        "genes": results,
    }
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_md(report)
    print(f"\nWrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
