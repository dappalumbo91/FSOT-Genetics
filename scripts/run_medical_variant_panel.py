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
from msa_uniref import best_conservation_profile, shannon_from_freq  # noqa: E402
from trinary_syntax import dna_to_aa, codon_primary  # noqa: E402
from fsot_structure_engine import clean_sequence  # noqa: E402

OUT_JSON = ROOT / "data" / "medical_variant_panel.json"
OUT_MD = ROOT / "predictions" / "reports" / "MEDICAL_VARIANT_PANEL.md"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "uniprot"
CACHE.mkdir(parents=True, exist_ok=True)

# Frozen call thresholds (seed-closed medical policy — not free fit weights)
PCT_DAMAGING = 75.0
PCT_UNCERTAIN = 40.0
# Absolute intolerance (needed for tight UniRef clusters where almost every
# site is conserved → percentiles collapse). Forms from {e, φ}:
#   CONS_HIGH  = 1 - 1/e² ≈ 0.865
#   CONS_MID   = 1 - 1/e  ≈ 0.632
#   FMUT_RARE  = 1/e²     ≈ 0.135
#   FMUT_ABSENT= 1/e³     ≈ 0.050
# 1 - exp(-π/2) ≈ 0.792 — high conservation absolute gate (seed-closed)
CONS_HIGH = 1.0 - math.exp(-math.pi / 2.0)
CONS_MID = 1.0 - 1.0 / math.e
FMUT_RARE = 1.0 / (math.e ** 2)
FMUT_ABSENT = 1.0 / (math.e ** 3)
# Common-allele demotion: population AF is *data* (gnomAD / 1000G), not a
# trained weight. Threshold 1/φ³ ≈ 0.236 — P72R is ~0.46 globally.
_PHI = (1.0 + 5.0 ** 0.5) / 2.0
POP_AF_COMMON = 1.0 / (_PHI ** 3)


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


def call_variant(
    *,
    kind: str,
    coverage_ok: bool,
    cons: float | None,
    f_mut: float | None,
    pct: float | None,
    pop_af: float | None = None,
) -> str:
    """Dual-gate call: absolute evolutionary intolerance OR high percentile.

    Tight UniRef clusters make nearly all sites high-cons, so percentile-only
    gates fail. Absolute gate mirrors SIFT spirit with seed-closed cutpoints.
    Population AF ≥ 1/φ³ demotes a damaging call — common poly is data
    (P72R), not a free specificity weight.
    """
    if kind.startswith("nonsense"):
        return "LIKELY DAMAGING"
    if kind.startswith("synonymous"):
        return "likely benign"
    if not coverage_ok:
        return "insufficient_MSA_coverage"
    c = float(cons or 0.0)
    fm = float(f_mut if f_mut is not None else 1.0)
    # Absolute intolerance
    if c >= CONS_HIGH and fm <= FMUT_ABSENT:
        call = "LIKELY DAMAGING"
    elif c >= CONS_HIGH and fm <= FMUT_RARE:
        call = "LIKELY DAMAGING"
    elif c >= CONS_MID and fm <= FMUT_ABSENT:
        call = "uncertain"
    elif pct is not None and pct >= PCT_DAMAGING:
        call = "LIKELY DAMAGING"
    elif pct is not None and pct >= PCT_UNCERTAIN:
        call = "uncertain"
    elif c >= CONS_MID and fm <= FMUT_RARE:
        call = "uncertain"
    else:
        call = "likely tolerated"
    if (
        pop_af is not None
        and float(pop_af) >= POP_AF_COMMON
        and call == "LIKELY DAMAGING"
    ):
        return "common_polymorphism"
    return call


def score_missense(
    seq: str,
    pos: int,
    wt: str,
    mut: str,
    cons: np.ndarray,
    freq: list[dict],
    background: np.ndarray,
    *,
    pop_af: float | None = None,
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
    call = call_variant(
        kind="missense",
        coverage_ok=coverage_ok,
        cons=c,
        f_mut=fm,
        pct=pct,
        pop_af=pop_af,
    )
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
            "call": call,
            "gates": {
                "cons_high": CONS_HIGH,
                "cons_mid": CONS_MID,
                "fmut_absent": FMUT_ABSENT,
                "fmut_rare": FMUT_RARE,
                "pct_damaging": PCT_DAMAGING,
                "pop_af_common": POP_AF_COMMON,
            },
            "pop_af": pop_af,
        }
    )
    return out


def score_gene(symbol: str) -> dict[str, Any]:
    gene = GENE_CATALOG[symbol]
    t0 = time.perf_counter()
    seq = fetch_uniprot_seq(gene["uniprot"])

    # 1) Protein-specific UniRef MSA (full chain — fixes N-term / edge coverage)
    # 2) Gap-fill remaining columns from domain Pfam MSAs (real observables both)
    cons_u, freq_u, n_u, meta_u = best_conservation_profile(
        seq,
        uniprot=gene["uniprot"],
        pfam=gene.get("pfam"),
        self_pdb=gene.get("structure_pdb"),
    )
    cons = np.array(cons_u, dtype=float)
    freq: list[dict] = list(freq_u)
    # Ensure length match
    if len(cons) < len(seq):
        cons = np.pad(cons, (0, len(seq) - len(cons)))
    while len(freq) < len(seq):
        freq.append({})

    # Gap-fill with domain Pfam where UniRef left zeros
    pfam_fills = 0
    for d in gene.get("domains_static") or []:
        s0, s1 = int(d["start"]) - 1, int(d["end"])
        s0 = max(0, s0)
        s1 = min(len(seq), s1)
        need = [i for i in range(s0, s1) if not freq[i] or cons[i] <= 0]
        if not need:
            continue
        pfam = d.get("pfam") or gene.get("pfam")
        sub = seq[s0:s1]
        try:
            cs, fs, _nr, _pf = conservation_profile(
                sub, gene.get("structure_pdb") or "XXXX", pfam=pfam, max_rows=2500
            )
        except Exception:
            continue
        for j, i in enumerate(range(s0, s1)):
            if j >= len(cs):
                break
            if not freq[i] or cons[i] <= 0:
                # blended local
                fj = fs[j] if j < len(fs) else {}
                cons[i] = max(float(cs[j]), shannon_from_freq(fj))
                freq[i] = fj
                pfam_fills += 1

    bg = missense_background(cons, freq, seq)
    source_label = meta_u.get("chosen", "unknown")
    if pfam_fills:
        source_label = f"{source_label}+pfam_gapfill"
    nrows0 = int(n_u)
    pf0 = meta_u.get("cluster_id") or meta_u.get("pfam") or gene.get("pfam")

    def get_profile_at(pos: int):
        # full-chain UniRef(+gapfill) profile — same arrays for all positions
        return cons, freq, nrows0, pf0, bg

    drivers_out = []
    for d in gene.get("drivers") or []:
        cons_p, freq_p, nrows, pf, bg_p = get_profile_at(d["pos"])
        sc = score_missense(
            seq,
            d["pos"],
            d["wt"],
            d["mut"],
            cons_p,
            freq_p,
            bg_p,
            pop_af=d.get("pop_af"),
        )
        sc.update(
            {
                "hgvs_p": d["hgvs_p"],
                "note": d.get("note"),
                "msa_source": source_label,
                "pfam_used": pf,
                "msa_rows": nrows,
                "role": "driver",
            }
        )
        drivers_out.append(sc)

    controls_out = []
    for d in gene.get("controls") or []:
        cons_p, freq_p, nrows, pf, bg_p = get_profile_at(d["pos"])
        sc = score_missense(
            seq,
            d["pos"],
            d["wt"],
            d["mut"],
            cons_p,
            freq_p,
            bg_p,
            pop_af=d.get("pop_af"),
        )
        sc.update(
            {
                "hgvs_p": d["hgvs_p"],
                "note": d.get("note"),
                "msa_source": source_label,
                "pfam_used": pf,
                "msa_rows": nrows,
                "role": "control",
            }
        )
        controls_out.append(sc)

    context_out = []
    for d in gene.get("context_dependent") or []:
        cons_p, freq_p, nrows, pf, bg_p = get_profile_at(d["pos"])
        sc = score_missense(
            seq,
            d["pos"],
            d["wt"],
            d["mut"],
            cons_p,
            freq_p,
            bg_p,
            pop_af=d.get("pop_af"),
        )
        sc.update(
            {
                "hgvs_p": d["hgvs_p"],
                "note": d.get("note"),
                "msa_source": source_label,
                "pfam_used": pf,
                "msa_rows": nrows,
                "role": "context_dependent",
            }
        )
        context_out.append(sc)

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
            cons_p, freq_p, nrows, pf, bg_p = get_profile_at(pos)
            sc = score_missense(seq, pos, wt_aa, mut_aa, cons_p, freq_p, bg_p)
            entry["impact_percentile"] = sc.get("impact_percentile")
            entry["call"] = sc.get("call")
            entry["conservation"] = sc.get("conservation")
        else:
            entry["call"] = call_variant(
                kind=kind, coverage_ok=True, cons=None, f_mut=None, pct=None
            )
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
        "msa_source": source_label,
        "msa_meta": meta_u,
        "msa_rows_default": nrows0,
        "mean_conservation_default": float(cons.mean()) if len(cons) else 0.0,
        "n_positions_covered": int(sum(1 for f in freq if f)),
        "pfam_gapfill_positions": pfam_fills,
        "drivers": drivers_out,
        "controls": controls_out,
        "context_dependent": context_out,
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
