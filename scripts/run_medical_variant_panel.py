#!/usr/bin/env python3
"""Multi-gene medical variant panel — FSOT conservation + trinary DNA path.

Expands the p53-only pilot to a curated disease-gene set (TP53, KRAS, EGFR,
BRAF, CFTR, SOD1, HBB, BRCA1). For each gene:

  1. Fetch UniProt canonical sequence (residue numbers = HGVS protein).
  2. Build Pfam MSA conservation (+ AA frequencies) — domain-aware when the
     variant sits in a known Pfam range.
  3. Score known drivers as conservation × (1 − f_mutant) and report
     genome-of-missense percentile (SIFT-style, zero trained weights).
  4. Route DNA examples through trinary codon layer when provided.

Writes:
  data/medical_variant_panel.json
  predictions/reports/MEDICAL_VARIANT_PANEL.md

Call thresholds (frozen, seed-closed):
  LIKELY DAMAGING  if impact percentile ≥ 75
  uncertain        if 40 ≤ pct < 75
  likely tolerated if pct < 40
  nonsense         → LIKELY DAMAGING (truncation)
  synonymous       → likely benign
"""

from __future__ import annotations

import json
import math
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from medical_gene_catalog import GENE_CATALOG, list_genes  # noqa: E402
from variant_conservation import conservation_profile, AA20  # noqa: E402
from trinary_syntax import dna_to_aa, codon_primary  # noqa: E402
from fsot_structure_engine import clean_sequence  # noqa: E402

OUT_JSON = ROOT / "data" / "medical_variant_panel.json"
OUT_MD = ROOT / "predictions" / "reports" / "MEDICAL_VARIANT_PANEL.md"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "uniprot"
CACHE.mkdir(parents=True, exist_ok=True)

# Frozen call thresholds (not free params — fixed medical policy cutpoints)
PCT_DAMAGING = 75.0
PCT_UNCERTAIN = 40.0


def fetch_uniprot_seq(acc: str) -> str:
    path = CACHE / f"{acc}.fasta"
    if path.exists():
        txt = path.read_text(encoding="utf-8")
        return clean_sequence("".join(l for l in txt.splitlines() if not l.startswith(">")))
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    req = urllib.request.Request(url, headers={"User-Agent": "fsot-genetics-medical"})
    with urllib.request.urlopen(req, timeout=60) as r:
        txt = r.read().decode("utf-8", "replace")
    path.write_text(txt, encoding="utf-8")
    return clean_sequence("".join(l for l in txt.splitlines() if not l.startswith(">")))


def pfam_for_position(gene: dict, pos: int) -> str:
    """Pick the Pfam domain containing pos; else gene default."""
    for d in gene.get("domains_static") or []:
        if d["start"] <= pos <= d["end"] and d.get("pfam"):
            return d["pfam"]
    return gene["pfam"]


def domain_span_for_position(gene: dict, pos: int) -> tuple[int, int] | None:
    """Return 1-based inclusive (start, end) of domain containing pos, if known."""
    for d in gene.get("domains_static") or []:
        if d["start"] <= pos <= d["end"]:
            return int(d["start"]), int(d["end"])
    return None


def shannon_conservation(freq_i: dict) -> float:
    """1 - H/log(20) from observed AA frequencies (query-independent).

    Identity-to-query conservation under-calls sites where the family is
    conserved *as a class* but the query residue is not the modal AA, and
    fails at MSA edge gaps. Shannon conservation is the field-standard
    evolutionary rate proxy (still zero free params).
    """
    if not freq_i:
        return 0.0
    vals = np.array([float(v) for v in freq_i.values() if v > 0], dtype=float)
    if vals.size == 0:
        return 0.0
    vals = vals / vals.sum()
    H = float(-(vals * np.log(vals)).sum())
    return max(0.0, 1.0 - H / math.log(20.0))


def position_conservation(identity_cons: float, freq_i: dict) -> float:
    """Blend identity and Shannon; take the stronger evolutionary signal."""
    return max(float(identity_cons), shannon_conservation(freq_i))


def missense_background(cons: np.ndarray, freq: list[dict], seq: str) -> np.ndarray:
    scores = []
    for i, aa in enumerate(seq):
        fi = freq[i] if i < len(freq) else {}
        c = position_conservation(cons[i], fi)
        for m in AA20:
            if m == aa:
                continue
            fm = fi.get(m, 0.0)
            scores.append(c * (1.0 - fm))
    return np.asarray(scores, dtype=float) if scores else np.zeros(1)


def call_from_pct(pct: float | None, *, kind: str, coverage_ok: bool = True) -> str:
    if kind.startswith("nonsense"):
        return "LIKELY DAMAGING"
    if kind.startswith("synonymous"):
        return "likely benign"
    if not coverage_ok:
        return "insufficient_MSA_coverage"
    if pct is None:
        return "uncertain"
    if pct >= PCT_DAMAGING:
        return "LIKELY DAMAGING"
    if pct >= PCT_UNCERTAIN:
        return "uncertain"
    return "likely tolerated"


def score_missense(
    seq: str,
    pos: int,
    wt: str,
    mut: str,
    cons: np.ndarray,
    freq: list[dict],
    background: np.ndarray,
) -> dict[str, Any]:
    """pos is 1-based UniProt."""
    i = pos - 1
    out: dict[str, Any] = {
        "pos": pos,
        "wt": wt,
        "mut": mut,
        "seq_wt": seq[i] if 0 <= i < len(seq) else None,
        "in_sequence": 0 <= i < len(seq),
    }
    if not out["in_sequence"]:
        out["error"] = "position_out_of_range"
        return out
    if seq[i] != wt:
        out["error"] = f"sequence_mismatch_expected_{wt}_found_{seq[i]}"
        return out
    fi = freq[i] if i < len(freq) else {}
    c = position_conservation(float(cons[i]), fi)
    coverage_ok = bool(fi) and c > 0.0
    if mut == "*":
        out["kind"] = "nonsense"
        out["impact"] = 1.0
        out["impact_percentile"] = 100.0
        out["conservation"] = c
        out["call"] = "LIKELY DAMAGING"
        return out
    fm = fi.get(mut, 0.0)
    impact = c * (1.0 - fm)
    pct = float((background < impact).mean()) * 100.0 if coverage_ok else None
    out.update(
        {
            "kind": "missense",
            "conservation": c,
            "identity_conservation": float(cons[i]),
            "shannon_conservation": shannon_conservation(fi),
            "mutant_freq": float(fm),
            "impact": impact,
            "impact_percentile": pct,
            "msa_coverage_ok": coverage_ok,
            "conservation_percentile": float((cons < cons[i]).mean()) * 100.0,
            "call": call_from_pct(pct, kind="missense", coverage_ok=coverage_ok),
        }
    )
    return out


def score_gene(symbol: str) -> dict[str, Any]:
    gene = GENE_CATALOG[symbol]
    t0 = time.perf_counter()
    seq = fetch_uniprot_seq(gene["uniprot"])
    # cache profiles by (pfam, start, end) on the domain slice — much faster
    # than aligning a multi-domain chain to a single-domain Pfam MSA
    profiles: dict[tuple, tuple] = {}

    def get_profile_at(pos: int):
        pfam = pfam_for_position(gene, pos)
        span = domain_span_for_position(gene, pos)
        if span is None:
            # whole-chain vs default pfam (short proteins / single domain)
            key = (pfam, 1, len(seq))
            s0, s1 = 0, len(seq)
        else:
            key = (pfam, span[0], span[1])
            s0, s1 = span[0] - 1, span[1]
        if key not in profiles:
            sub = seq[s0:s1]
            cons_sub, freq_sub, nrows, pf = conservation_profile(
                sub, gene.get("structure_pdb") or "XXXX", pfam=pfam, max_rows=3000
            )
            # expand to full-chain arrays (zeros outside domain)
            cons = np.zeros(len(seq))
            freq: list[dict] = [{} for _ in range(len(seq))]
            for i, c in enumerate(cons_sub):
                cons[s0 + i] = c
                if i < len(freq_sub):
                    freq[s0 + i] = freq_sub[i]
            profiles[key] = (cons, freq, nrows, pf, s0, s1)
        cons, freq, nrows, pf, s0, s1 = profiles[key]
        # background only over residues inside this domain (fair percentile)
        bg = missense_background(cons[s0:s1], freq[s0:s1], seq[s0:s1])
        return cons, freq, nrows, pf, bg

    # default profile for summary = primary gene pfam domain or whole chain
    cons0, freq0, nrows0, pf0, _bg0 = get_profile_at(
        (gene.get("domains_static") or [{"start": 1}])[0].get("start", 1)
        if gene.get("domains_static")
        else 1
    )

    drivers_out = []
    for d in gene.get("drivers") or []:
        cons, freq, nrows, pf, bg = get_profile_at(d["pos"])
        sc = score_missense(seq, d["pos"], d["wt"], d["mut"], cons, freq, bg)
        sc.update(
            {
                "hgvs_p": d["hgvs_p"],
                "note": d.get("note"),
                "pfam_used": pf,
                "msa_rows": nrows,
                "role": "driver",
            }
        )
        drivers_out.append(sc)

    controls_out = []
    for d in gene.get("controls") or []:
        cons, freq, nrows, pf, bg = get_profile_at(d["pos"])
        sc = score_missense(seq, d["pos"], d["wt"], d["mut"], cons, freq, bg)
        sc.update(
            {
                "hgvs_p": d["hgvs_p"],
                "note": d.get("note"),
                "pfam_used": pf,
                "msa_rows": nrows,
                "role": "control",
            }
        )
        controls_out.append(sc)

    dna_out = []
    for item in gene.get("dna_examples") or []:
        pos, wt_codon, cpos, alt, hgvs = item
        mut = list(wt_codon)
        mut[cpos] = alt
        mut_codon = "".join(mut)
        wt_aa, mut_aa = dna_to_aa(wt_codon), dna_to_aa(mut_codon)
        if mut_aa in ("*", "Stop"):
            kind = "nonsense"
        elif wt_aa == mut_aa:
            kind = "synonymous"
        else:
            kind = "missense"
        entry: dict[str, Any] = {
            "hgvs_c": hgvs,
            "wt_codon": wt_codon,
            "mut_codon": mut_codon,
            "trinary_wt": list(codon_primary(wt_codon)),
            "trinary_mut": list(codon_primary(mut_codon)),
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "pos": pos,
            "kind": kind,
        }
        if kind == "missense":
            cons, freq, nrows, pf, bg = get_profile_at(pos)
            sc = score_missense(seq, pos, wt_aa, mut_aa, cons, freq, bg)
            entry["impact_percentile"] = sc.get("impact_percentile")
            entry["call"] = sc.get("call")
            entry["conservation"] = sc.get("conservation")
        else:
            entry["call"] = call_from_pct(None, kind=kind)
        dna_out.append(entry)

    driver_pcts = [
        d["impact_percentile"]
        for d in drivers_out
        if isinstance(d.get("impact_percentile"), (int, float)) and d.get("kind") == "missense"
    ]
    n_dam = sum(1 for d in drivers_out if d.get("call") == "LIKELY DAMAGING")
    n_ok = sum(1 for d in drivers_out if not d.get("error"))

    return {
        "symbol": symbol,
        "name": gene["name"],
        "uniprot": gene["uniprot"],
        "indication": gene.get("indication"),
        "length": len(seq),
        "default_pfam": pf0,
        "msa_rows_default": nrows0,
        "mean_conservation_default": float(cons0.mean()) if len(cons0) else 0.0,
        "drivers": drivers_out,
        "controls": controls_out,
        "dna_examples": dna_out,
        "n_drivers_scored": n_ok,
        "n_drivers_likely_damaging": n_dam,
        "mean_driver_impact_percentile": float(np.mean(driver_pcts)) if driver_pcts else None,
        "elapsed_s": round(time.perf_counter() - t0, 2),
        "free_parameters": 0,
    }


def write_md(report: dict) -> None:
    lines = [
        "# Medical variant panel — multi-gene FSOT conservation",
        "",
        f"Generated: `{report['generated_at']}`  ",
        f"Free parameters: **0**  ",
        f"Genes: **{report['summary']['n_genes']}**  ",
        f"Drivers scored: **{report['summary']['n_drivers']}**  ",
        f"Drivers called LIKELY DAMAGING: **{report['summary']['n_likely_damaging']}** "
        f"({report['summary']['driver_recall_at_75pct']:.0%})  ",
        f"Mean driver impact percentile: **{report['summary']['mean_driver_percentile']}**",
        "",
        "Thresholds: ≥75 LIKELY DAMAGING · 40–75 uncertain · <40 likely tolerated",
        "",
    ]
    for g in report["genes"]:
        lines += [
            f"## {g['symbol']} — {g['name']}",
            "",
            f"UniProt `{g['uniprot']}` · n={g['length']} · Pfam `{g['default_pfam']}` · "
            f"MSA rows={g['msa_rows_default']} · mean cons={g['mean_conservation_default']:.2f}",
            "",
            f"Indication: {g.get('indication')}",
            "",
            f"Mean driver percentile: **{g.get('mean_driver_impact_percentile')}** · "
            f"LIKELY DAMAGING {g['n_drivers_likely_damaging']}/{g['n_drivers_scored']}",
            "",
            "| HGVS | cons | impact% | call | note |",
            "|------|-----:|--------:|------|------|",
        ]
        for d in g["drivers"]:
            if d.get("error"):
                lines.append(f"| {d.get('hgvs_p')} | err | | | {d['error']} |")
                continue
            cons = d.get("conservation")
            pct = d.get("impact_percentile")
            cons_s = f"{cons:.2f}" if isinstance(cons, (int, float)) else "—"
            pct_s = f"{pct:.0f}" if isinstance(pct, (int, float)) else "—"
            lines.append(
                f"| {d.get('hgvs_p')} | {cons_s} | {pct_s} | {d.get('call')} | {d.get('note','')} |"
            )
        if g.get("dna_examples"):
            lines += ["", "DNA front door:", ""]
            for e in g["dna_examples"]:
                lines.append(
                    f"- `{e['hgvs_c']}` {e['wt_codon']}→{e['mut_codon']} "
                    f"{e['wt_aa']}{e['pos']}{e['mut_aa'] if e['mut_aa']!=e['wt_aa'] else '='} "
                    f"({e['kind']}) → **{e.get('call')}**"
                )
        lines.append("")
    lines += [
        "## Honesty",
        "",
        "- Conservation from real Pfam MSAs; impact = cons × (1 − mutant frequency).",
        "- Not a substitute for ACMG clinical classification; research / triage tool.",
        "- Domain-aware Pfam selection when variant falls in a static domain range.",
        "",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    genes = list_genes()
    if argv:
        genes = [g.upper() for g in argv if g.upper() in GENE_CATALOG]
    print("FSOT multi-gene medical variant panel", flush=True)
    print("=" * 78, flush=True)
    results = []
    for sym in genes:
        print(f"\n--- {sym} ---", flush=True)
        try:
            row = score_gene(sym)
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR {exc}", flush=True)
            results.append({"symbol": sym, "error": str(exc)})
            continue
        results.append(row)
        mp = row.get("mean_driver_impact_percentile")
        print(
            f"  n={row['length']}  pfam={row['default_pfam']}  msa={row['msa_rows_default']}  "
            f"mean_cons={row['mean_conservation_default']:.2f}"
        )
        print(
            f"  drivers LIKELY DAMAGING {row['n_drivers_likely_damaging']}/"
            f"{row['n_drivers_scored']}  mean_pct={mp}"
        )
        for d in row["drivers"]:
            if d.get("error"):
                print(f"    {d.get('hgvs_p')} ERROR {d['error']}")
            else:
                pct = d.get("impact_percentile")
                print(
                    f"    {d.get('hgvs_p'):<10} cons={d.get('conservation',0):.2f}  "
                    f"pct={pct:5.1f}%  {d.get('call')}"
                    if pct is not None
                    else f"    {d.get('hgvs_p')} {d.get('call')}"
                )

    all_drivers = []
    for g in results:
        if g.get("error"):
            continue
        all_drivers.extend(
            d for d in g.get("drivers") or [] if not d.get("error") and d.get("kind") == "missense"
        )
    n_dam = sum(1 for d in all_drivers if d.get("call") == "LIKELY DAMAGING")
    pcts = [d["impact_percentile"] for d in all_drivers if d.get("impact_percentile") is not None]
    summary = {
        "n_genes": sum(1 for g in results if not g.get("error")),
        "n_drivers": len(all_drivers),
        "n_likely_damaging": n_dam,
        "driver_recall_at_75pct": (n_dam / len(all_drivers)) if all_drivers else 0.0,
        "mean_driver_percentile": float(np.mean(pcts)) if pcts else None,
        "threshold_damaging": PCT_DAMAGING,
        "threshold_uncertain": PCT_UNCERTAIN,
    }
    print("\n" + "=" * 78)
    print(
        f"SUMMARY genes={summary['n_genes']} drivers={summary['n_drivers']} "
        f"LIKELY DAMAGING={summary['n_likely_damaging']} "
        f"recall@75={summary['driver_recall_at_75pct']:.0%} "
        f"mean_pct={summary['mean_driver_percentile']}"
    )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "multi_gene_medical_variant_panel",
        "free_parameters": 0,
        "method": "Pfam MSA conservation × (1 - mutant_AA_frequency); domain-aware Pfam",
        "summary": summary,
        "genes": results,
    }
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_md(report)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
