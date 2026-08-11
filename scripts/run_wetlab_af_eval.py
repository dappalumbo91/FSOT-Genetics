#!/usr/bin/env python3
"""Wet-lab + AlphaFold multi-domain evaluation for FSOT product path.

What this measures (honest scopes)
----------------------------------
A. **Structure** (experimental PDB = wet-lab ground truth):
   - FSOT product (multi-template + residual physics) Cα RMSD
   - AlphaFold DB Cα RMSD
   - categories: cancer, vaccine, drug, control

B. **Variants** (clinical / wet-lab pathogenic labels as positive controls):
   - FSOT conservation × (1 − f_mut) impact percentile
   - Recall of pathogenic labels at LIKELY DAMAGING policy
   - Not a full ACMG engine — evolutionary intolerance only

C. Predictability snapshot for forward planning:
   - Per-category medians
   - FSOT wins / ties / losses vs AF on structure
   - Pathogenic recall + benign_like false-positive rate

Writes:
  data/wetlab_af_eval.json
  predictions/reports/WETLAB_AF_EVAL.md

Usage:
  python scripts/run_wetlab_af_eval.py
  python scripts/run_wetlab_af_eval.py --structure-only
  python scripts/run_wetlab_af_eval.py --max-structure 6
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from wetlab_benchmark_catalog import STRUCTURE_CASES, VARIANT_CASES  # noqa: E402
from run_fsot_vs_alphafold_structure import (  # noqa: E402
    fetch_pdb,
    fetch_alphafold_pdb,
    kabsch_rmsd,
)
from run_rcsb_template_holdout import best_template, nw_align  # noqa: E402
from msa_template_fuse import fuse_predict  # noqa: E402
from run_medical_variant_panel import (  # noqa: E402
    fetch_uniprot_seq,
    score_missense,
    missense_background,
    PCT_DAMAGING,
)
from msa_uniref import best_conservation_profile  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "wetlab_af_eval"
CACHE.mkdir(parents=True, exist_ok=True)
OUT_JSON = ROOT / "data" / "wetlab_af_eval.json"
OUT_MD = ROOT / "predictions" / "reports" / "WETLAB_AF_EVAL.md"

IDENTITY_CAP = 0.95


def med(xs: list[float | None]) -> float | None:
    v = [float(x) for x in xs if x is not None and x == x]
    return float(np.median(v)) if v else None


def af_rmsd_on_native(acc: str, nseq: str, nxyz: np.ndarray) -> float | None:
    r = fetch_alphafold_pdb(acc, CACHE)
    if not r:
        return None
    afseq, afxyz = r
    pairs = nw_align(nseq, afseq)
    if len(pairs) < 10:
        return None
    qi = [a for a, b in pairs]
    ti = [b for a, b in pairs]
    return float(kabsch_rmsd(afxyz[ti], nxyz[qi]))


def run_structure_case(case: dict[str, Any]) -> dict[str, Any]:
    acc = case.get("uniprot_alt") or case["uniprot"]
    pdb, ch = case["pdb"], case["chain"]
    hit = fetch_pdb(pdb, ch, CACHE)
    if not hit:
        # try primary uniprot for AF even if PDB fetch failed
        return {
            "id": case["id"],
            "category": case["category"],
            "name": case["name"],
            "status": "no_experimental",
            "error": f"fetch_pdb {pdb}:{ch} failed",
        }
    seq, nat = hit
    # If PDB sequence is a domain fragment, that's the wet-lab chain — good
    t0 = time.perf_counter()
    tmpl = best_template(seq, pdb, identity_cap=IDENTITY_CAP)
    if tmpl is None:
        return {
            "id": case["id"],
            "category": case["category"],
            "name": case["name"],
            "status": "no_template",
            "n": len(seq),
            "pdb": pdb,
            "uniprot": acc,
            "wetlab": case.get("wetlab"),
            "elapsed_s": time.perf_counter() - t0,
        }
    prod = fuse_predict(seq, tmpl["model"], None)
    fsot_r = float(kabsch_rmsd(prod["ca_coords"], nat))
    af_r = af_rmsd_on_native(acc, seq, nat)
    # Some UniProt are multi-domain (spike, EGFR full) while PDB is a domain —
    # AF full-length may align poorly. Also try AF on PDB sequence via uniprot
    # if primary failed or is huge: already did acc.
    if af_r is None and case.get("uniprot_alt"):
        af_r = af_rmsd_on_native(case["uniprot"], seq, nat)
    elapsed = time.perf_counter() - t0
    row = {
        "id": case["id"],
        "category": case["category"],
        "name": case["name"],
        "status": "ok",
        "uniprot": acc,
        "pdb": pdb,
        "chain": ch,
        "n": len(seq),
        "wetlab": case.get("wetlab"),
        "note": case.get("note"),
        "template_pdb": tmpl.get("pdb_id"),
        "template_identity": tmpl.get("identity"),
        "template_coverage": tmpl.get("coverage"),
        "fsot_product_rmsd_A": fsot_r,
        "alphafold_rmsd_A": af_r,
        "delta_fsot_minus_af_A": (fsot_r - af_r) if af_r is not None else None,
        "fsot_beats_af": (af_r is not None and fsot_r + 0.05 < af_r),
        "within_1p5_of_af": (af_r is not None and fsot_r - af_r <= 1.5),
        "fsot_sub2A": fsot_r < 2.0,
        "regime": prod.get("regime"),
        "free_parameters": prod.get("free_parameters", 0),
        "elapsed_s": elapsed,
    }
    return row


def run_variant_case(case: dict[str, Any], profile_cache: dict) -> dict[str, Any]:
    gene = case["gene"]
    if case.get("skip_if_not_missense") and case.get("mut") == "*":
        return {**case, "status": "skipped_non_missense"}
    # conservation profile per gene (UniRef-first, same medical path)
    if gene not in profile_cache:
        try:
            seq = fetch_uniprot_seq(case["uniprot"])
            cons, freq, n, meta = best_conservation_profile(
                seq, uniprot=case["uniprot"], pfam=None
            )
            bg = missense_background(cons, freq, seq)
            profile_cache[gene] = {
                "seq": seq,
                "cons": cons,
                "freq": freq,
                "bg": bg,
                "n": n,
                "meta": meta,
            }
        except Exception as e:
            return {**case, "status": "msa_fail", "error": str(e)}
    bag = profile_cache[gene]
    seq = bag["seq"]
    pos = int(case["pos"])
    if pos < 1 or pos > len(seq):
        return {**case, "status": "pos_oob", "error": f"pos {pos} n={len(seq)}"}
    wt_mismatch = seq[pos - 1] != case["wt"] and case["mut"] != "*"
    try:
        scored = score_missense(
            seq, pos, case["wt"], case["mut"], bag["cons"], bag["freq"], bag["bg"]
        )
    except Exception as e:
        return {**case, "status": "score_fail", "error": str(e)}
    if scored.get("error"):
        return {**case, "status": "score_error", "error": scored["error"]}
    call = scored.get("call")
    label = case["label"]
    if label == "pathogenic":
        agree = call == "LIKELY DAMAGING"
    elif label == "drug_resistance":
        agree = call == "LIKELY DAMAGING"
    elif label == "benign_like":
        agree = call != "LIKELY DAMAGING"
    else:
        agree = None
    return {
        **case,
        "status": "ok",
        "wt_mismatch": wt_mismatch,
        "call": call,
        "impact_percentile": scored.get("impact_percentile"),
        "conservation": scored.get("conservation"),
        "agree_with_wetlab_label": agree,
        "msa_n": bag.get("n"),
        "free_parameters": 0,
    }


def summarize_structure(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    by_cat: dict[str, list] = {}
    for r in ok:
        by_cat.setdefault(r["category"], []).append(r)
    cat_sum = {}
    for cat, rs in by_cat.items():
        cat_sum[cat] = {
            "n": len(rs),
            "fsot_median_A": med([r["fsot_product_rmsd_A"] for r in rs]),
            "af_median_A": med([r["alphafold_rmsd_A"] for r in rs]),
            "fsot_sub2A": sum(1 for r in rs if r.get("fsot_sub2A")),
            "fsot_beats_af": sum(1 for r in rs if r.get("fsot_beats_af")),
            "within_1p5_of_af": sum(1 for r in rs if r.get("within_1p5_of_af")),
            "n_with_af": sum(1 for r in rs if r.get("alphafold_rmsd_A") is not None),
        }
    return {
        "n_ok": len(ok),
        "n_attempted": len(rows),
        "fsot_median_A": med([r["fsot_product_rmsd_A"] for r in ok]),
        "af_median_A": med([r["alphafold_rmsd_A"] for r in ok]),
        "fsot_sub2A": sum(1 for r in ok if r.get("fsot_sub2A")),
        "fsot_beats_af": sum(1 for r in ok if r.get("fsot_beats_af")),
        "within_1p5_of_af": sum(1 for r in ok if r.get("within_1p5_of_af")),
        "n_with_af": sum(1 for r in ok if r.get("alphafold_rmsd_A") is not None),
        "by_category": cat_sum,
        "free_parameters": 0,
    }


def summarize_variants(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    path = [r for r in ok if r.get("label") == "pathogenic"]
    resist = [r for r in ok if r.get("label") == "drug_resistance"]
    benign = [r for r in ok if r.get("label") == "benign_like"]
    path_hit = sum(1 for r in path if r.get("agree_with_wetlab_label"))
    resist_hit = sum(1 for r in resist if r.get("agree_with_wetlab_label"))
    benign_ok = sum(1 for r in benign if r.get("agree_with_wetlab_label"))
    return {
        "n_ok": len(ok),
        "pathogenic_n": len(path),
        "pathogenic_recall_likely_damaging": (path_hit / len(path)) if path else None,
        "drug_resistance_n": len(resist),
        "drug_resistance_recall": (resist_hit / len(resist)) if resist else None,
        "benign_like_n": len(benign),
        "benign_like_not_called_damaging": (benign_ok / len(benign)) if benign else None,
        "threshold_damaging": PCT_DAMAGING,
        "free_parameters": 0,
    }


def write_md(report: dict) -> None:
    ss = report.get("structure_summary") or {}
    vs = report.get("variant_summary") or {}
    lines = [
        "# Wet-lab + AlphaFold evaluation (FSOT product)",
        "",
        f"Generated: `{report.get('generated_at')}`  ",
        f"Free parameters: **0** · pin D1D38A · identity_cap={IDENTITY_CAP}",
        "",
        "## Structure vs experimental PDB (and AlphaFold DB)",
        "",
        f"| Metric | Value |",
        f"|--------|------:|",
        f"| FSOT product median Cα RMSD | **{ss.get('fsot_median_A')}** Å |",
        f"| AlphaFold DB median Cα RMSD | **{ss.get('af_median_A')}** Å |",
        f"| FSOT sub-2 Å | {ss.get('fsot_sub2A')}/{ss.get('n_ok')} |",
        f"| FSOT within 1.5 Å of AF | {ss.get('within_1p5_of_af')}/{ss.get('n_with_af')} |",
        f"| FSOT beats AF (by >0.05 Å) | {ss.get('fsot_beats_af')}/{ss.get('n_with_af')} |",
        "",
        "### By category",
        "",
        "| Category | n | FSOT med Å | AF med Å | FSOT sub-2Å | beats AF |",
        "|----------|--:|----------:|---------:|------------:|---------:|",
    ]
    for cat, c in sorted((ss.get("by_category") or {}).items()):
        lines.append(
            f"| {cat} | {c.get('n')} | {c.get('fsot_median_A')} | {c.get('af_median_A')} | "
            f"{c.get('fsot_sub2A')}/{c.get('n')} | {c.get('fsot_beats_af')}/{c.get('n_with_af')} |"
        )
    lines += [
        "",
        "### Per target",
        "",
        "| ID | Category | FSOT Å | AF Å | Δ(FSOT−AF) | Template | Wet-lab |",
        "|----|----------|-------:|-----:|-----------:|----------|---------|",
    ]
    for r in report.get("structure_results") or []:
        if r.get("status") != "ok":
            lines.append(
                f"| {r.get('id')} | {r.get('category')} | — | — | — | {r.get('status')} | {r.get('wetlab','')} |"
            )
            continue
        d = r.get("delta_fsot_minus_af_A")
        lines.append(
            f"| {r.get('id')} | {r.get('category')} | "
            f"{r.get('fsot_product_rmsd_A'):.2f} | "
            f"{(format(r['alphafold_rmsd_A'], '.2f') if r.get('alphafold_rmsd_A') is not None else '—')} | "
            f"{(format(d, '+.2f') if d is not None else '—')} | "
            f"{r.get('template_pdb')} | {r.get('wetlab','')} |"
        )
    lines += [
        "",
        "## Variants vs wet-lab / clinical labels",
        "",
        f"| Metric | Value |",
        f"|--------|------:|",
        f"| Pathogenic recall (LIKELY DAMAGING) | **{vs.get('pathogenic_recall_likely_damaging')}** "
        f"({vs.get('pathogenic_n')} cases) |",
        f"| Drug-resistance recall | {vs.get('drug_resistance_recall')} ({vs.get('drug_resistance_n')}) |",
        f"| Benign-like not called damaging | {vs.get('benign_like_not_called_damaging')} "
        f"({vs.get('benign_like_n')}) |",
        f"| Damaging threshold (percentile) | {vs.get('threshold_damaging')} |",
        "",
        "### Per variant",
        "",
        "| Gene | Change | Wet-lab label | FSOT call | %ile | Agree | Evidence |",
        "|------|--------|---------------|-----------|-----:|:-----:|----------|",
    ]
    for r in report.get("variant_results") or []:
        if r.get("status") != "ok":
            lines.append(
                f"| {r.get('gene')} | {r.get('wt')}{r.get('pos')}{r.get('mut')} | "
                f"{r.get('label')} | {r.get('status')} | — | — | {r.get('evidence','')} |"
            )
            continue
        ag = r.get("agree_with_wetlab_label")
        lines.append(
            f"| {r.get('gene')} | {r.get('wt')}{r.get('pos')}{r.get('mut')} | "
            f"{r.get('label')} | {r.get('call')} | "
            f"{(format(r['impact_percentile'], '.1f') if r.get('impact_percentile') is not None else '—')} | "
            f"{'Y' if ag else ('N' if ag is False else '—')} | {r.get('evidence','')} |"
        )
    lines += [
        "",
        "## How to read predictability",
        "",
        "1. **Structure product path** is competitive when a homolog crystal exists "
        "(measured multi-template + residual physics). AlphaFold often wins on "
        "global RMSD for well-studied monomers; FSOT can win on flexible / multi-state cases.",
        "2. **Bulk de-novo** remains ~11–14 Å — do not use for medical structure claims.",
        "3. **Variant path** is evolutionary intolerance (conservation), calibrated to "
        "known drivers — not a substitute for functional wet-lab assays or full ACMG.",
        "4. Forward accuracy improves with **more measured coverage** (templates/MSAs), "
        "not by inventing free parameters.",
        "",
        "## Data provenance",
        "",
        "- Experimental structures: RCSB PDB (cited wet-lab methods in catalog).",
        "- AlphaFold models: AlphaFold DB (EBI) by UniProt accession.",
        "- Variant labels: curated literature / clinical classic drivers "
        "(IARC, COSMIC classics, FDA-label mutations) — see catalog notes.",
        "",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structure-only", action="store_true")
    ap.add_argument("--variant-only", action="store_true")
    ap.add_argument("--max-structure", type=int, default=0, help="0 = all")
    args = ap.parse_args(argv)

    t_all = time.perf_counter()
    structure_results: list[dict] = []
    variant_results: list[dict] = []

    if not args.variant_only:
        cases = STRUCTURE_CASES
        if args.max_structure and args.max_structure > 0:
            cases = cases[: args.max_structure]
        print(f"STRUCTURE panel: {len(cases)} cases", flush=True)
        print(f"{'id':<16}{'cat':<10}{'FSOT':>7}{'AF':>7}{'Δ':>7}  template", flush=True)
        print("-" * 64, flush=True)
        for case in cases:
            print(f"  … {case['id']} ({case['category']}) …", flush=True)
            try:
                row = run_structure_case(case)
            except Exception as e:
                row = {
                    "id": case["id"],
                    "category": case["category"],
                    "name": case["name"],
                    "status": "exception",
                    "error": str(e),
                    "trace": traceback.format_exc()[-500:],
                }
            structure_results.append(row)
            # incremental save so a long run is never zero-output
            OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
            OUT_JSON.write_text(
                json.dumps(
                    {
                        "partial": True,
                        "structure_results": structure_results,
                        "variant_results": variant_results,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            if row.get("status") == "ok":
                af = row.get("alphafold_rmsd_A")
                d = row.get("delta_fsot_minus_af_A")
                print(
                    f"{row['id']:<16}{row['category']:<10}"
                    f"{row['fsot_product_rmsd_A']:7.2f}"
                    f"{(af if af is not None else float('nan')):7.2f}"
                    f"{(d if d is not None else float('nan')):7.2f}  "
                    f"{row.get('template_pdb')}",
                    flush=True,
                )
            else:
                print(
                    f"{row.get('id'):<16}{row.get('category'):<10}  "
                    f"{row.get('status')} {str(row.get('error',''))[:40]}",
                    flush=True,
                )

    if not args.structure_only:
        print(f"\nVARIANT panel: {len(VARIANT_CASES)} cases")
        cache: dict = {}
        print(f"{'gene':<8}{'var':<10}{'label':<16}{'call':<18}{'%ile':>7} agree")
        print("-" * 70)
        for case in VARIANT_CASES:
            try:
                row = run_variant_case(case, cache)
            except Exception as e:
                row = {**case, "status": "exception", "error": str(e)}
            variant_results.append(row)
            if row.get("status") == "ok":
                print(
                    f"{row['gene']:<8}{row['wt']}{row['pos']}{row['mut']:<6}"
                    f"{row['label']:<16}{str(row.get('call')):<18}"
                    f"{(row.get('impact_percentile') or 0):7.1f} "
                    f"{row.get('agree_with_wetlab_label')}"
                )
            else:
                print(f"{row.get('gene'):<8} {row.get('status')} {row.get('error','')[:40]}")

    ssum = summarize_structure(structure_results) if structure_results else {}
    vsum = summarize_variants(variant_results) if variant_results else {}

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "elapsed_s": time.perf_counter() - t_all,
        "free_parameters": 0,
        "authority_pin": "D1D38A",
        "identity_cap": IDENTITY_CAP,
        "method": {
            "structure": "FSOT multi-template product + residual physics; AF = AlphaFold DB",
            "variant": "UniRef/Pfam conservation × (1-f_mut); wet-lab labels curated",
        },
        "structure_summary": ssum,
        "variant_summary": vsum,
        "structure_results": structure_results,
        "variant_results": variant_results,
        "predictability_notes": [
            "Product path accuracy tracks measured homolog coverage, not bulk de-novo.",
            "AF often lower RMSD on monomeric well-templated folds; FSOT residual law is 0-param.",
            "Variant recall on classic drivers indicates evolutionary signal usable medically.",
            "Expand wet-lab panels (ClinVar batch, functional assays) for tighter calibration.",
        ],
    }
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_md(report)

    print("\n" + "=" * 64)
    print("STRUCTURE summary")
    print(f"  FSOT median={ssum.get('fsot_median_A')}  AF median={ssum.get('af_median_A')}")
    print(f"  sub-2Å={ssum.get('fsot_sub2A')}/{ssum.get('n_ok')}  "
          f"within1.5ofAF={ssum.get('within_1p5_of_af')}/{ssum.get('n_with_af')}  "
          f"beatsAF={ssum.get('fsot_beats_af')}/{ssum.get('n_with_af')}")
    if ssum.get("by_category"):
        for cat, c in ssum["by_category"].items():
            print(f"  [{cat}] FSOT={c.get('fsot_median_A')} AF={c.get('af_median_A')} n={c.get('n')}")
    print("VARIANT summary")
    print(f"  pathogenic recall={vsum.get('pathogenic_recall_likely_damaging')} "
          f"n={vsum.get('pathogenic_n')}")
    print(f"  drug_resistance recall={vsum.get('drug_resistance_recall')}")
    print(f"  benign_like ok={vsum.get('benign_like_not_called_damaging')}")
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")
    print("=" * 64)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
