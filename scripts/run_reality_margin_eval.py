#!/usr/bin/env python3
"""Reality-match eval: FSOT product vs experimental wet-lab only.

AlphaFold is NOT the target. Experimental PDB Cα is ground truth.
For each wet-lab case:
  1. Product path (M1 authority + residual physics)
  2. Kabsch Cα RMSD to crystal/NMR/cryo-EM
  3. Error-margin localization (distance bins, hotspots, topology)
  4. Aggregate primary failure modes → FSOT fix queue

Writes:
  data/reality_margin_eval.json
  predictions/reports/REALITY_MARGIN_EVAL.md

Usage:
  python scripts/run_reality_margin_eval.py
  python scripts/run_reality_margin_eval.py --max 8
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from wetlab_benchmark_catalog import STRUCTURE_CASES  # noqa: E402
from run_fsot_vs_alphafold_structure import fetch_pdb, kabsch_rmsd  # noqa: E402
from run_rcsb_template_holdout import best_template  # noqa: E402
from msa_template_fuse import fuse_predict  # noqa: E402
from run_error_margin_log import (  # noqa: E402
    analyze_one,
    ERROR_CATALOG,
    build_fix_queue,
)

CACHE = Path.home() / ".cache" / "fsot-genetics" / "reality_margin"
CACHE.mkdir(parents=True, exist_ok=True)
OUT_JSON = ROOT / "data" / "reality_margin_eval.json"
OUT_MD = ROOT / "predictions" / "reports" / "REALITY_MARGIN_EVAL.md"

# Reality gates (product path vs wet lab — not AF)
REALITY_MEDIAN_TARGET = 2.5  # stretch goal for medical panel
REALITY_SUB2_FRAC = 0.5
PRODUCT_FREEZE_MAX = 1.20  # soft: H2H should stay near freeze


def med(xs: list[float]) -> float | None:
    return float(np.median(xs)) if xs else None


def run_case(case: dict[str, Any]) -> dict[str, Any]:
    pdb, ch = case["pdb"], case["chain"]
    hit = fetch_pdb(pdb, ch, CACHE)
    if not hit:
        return {"id": case["id"], "status": "no_fetch", "category": case["category"]}
    seq, nat = hit
    t0 = time.perf_counter()
    tmpl = best_template(seq, pdb, identity_cap=0.95)
    if tmpl is None:
        return {
            "id": case["id"],
            "status": "no_template",
            "category": case["category"],
            "name": case["name"],
            "n": len(seq),
            "wetlab": case.get("wetlab"),
            "elapsed_s": time.perf_counter() - t0,
        }
    prod = fuse_predict(
        seq, tmpl["model"], None, tertiary_contacts=tmpl.get("tertiary_contacts")
    )
    X = prod["ca_coords"]
    rmsd = float(kabsch_rmsd(X, nat))
    # Error margin on experimental sequence (product model)
    margin = analyze_one(
        case["name"],
        f"{pdb}:{ch}",
        seq,
        nat,
        {
            "sequence": seq,
            "ca_coords": X,
            "engine": prod.get("engine"),
            "long_range_gate": 7,
        },
    )
    return {
        "id": case["id"],
        "status": "ok",
        "category": case["category"],
        "name": case["name"],
        "wetlab": case.get("wetlab"),
        "note": case.get("note"),
        "pdb": pdb,
        "chain": ch,
        "n": len(seq),
        "template_pdb": tmpl.get("pdb_id"),
        "template_mode": tmpl.get("template_mode"),
        "template_identity": tmpl.get("identity"),
        "template_coverage": tmpl.get("coverage"),
        "length_sim": tmpl.get("length_sim")
        or (tmpl.get("score") and None),  # may live on cand only
        "n_candidates": tmpl.get("n_candidates"),
        "expanded_isoform_pool": tmpl.get("expanded_isoform_pool"),
        "fsot_product_rmsd_A": rmsd,
        "sub2A": rmsd < 2.0,
        "sub3A": rmsd < 3.0,
        "primary_error_mode": margin.get("primary_error_mode"),
        "ranked_modes": margin.get("ranked_modes"),
        "per_res_p90_A": margin.get("per_res_p90_A"),
        "contact_mae_A": margin.get("contact_mae_A"),
        "rg_err_A": margin.get("rg_err_A"),
        "hotspots_top5": (margin.get("hotspots") or [])[:5],
        "distance_bins": margin.get("distance_bins"),
        "mode_scores": margin.get("mode_scores"),
        "free_parameters": 0,
        "elapsed_s": time.perf_counter() - t0,
    }


def write_md(report: dict) -> None:
    sm = report.get("summary") or {}
    lines = [
        "# Reality-margin eval — FSOT product vs wet-lab experimental",
        "",
        f"Generated: `{report.get('generated_at')}`  ",
        "**Target = experimental PDB Cα, not AlphaFold.**  ",
        f"Free parameters: **0** · pin D1D38A",
        "",
        "## Summary",
        "",
        f"| Metric | Value |",
        f"|--------|------:|",
        f"| n OK | {sm.get('n_ok')}/{sm.get('n_attempted')} |",
        f"| **Median Cα RMSD vs wet lab** | **{sm.get('median_rmsd_A')}** Å |",
        f"| Sub-2 Å | {sm.get('n_sub2A')}/{sm.get('n_ok')} |",
        f"| Sub-3 Å | {sm.get('n_sub3A')}/{sm.get('n_ok')} |",
        f"| Primary mode (most common) | `{sm.get('dominant_error_mode')}` |",
        "",
        "### By category (median vs experimental)",
        "",
        "| Category | n | Median Å | Sub-2 Å | Dominant mode |",
        "|----------|--:|---------:|--------:|---------------|",
    ]
    for cat, c in sorted((sm.get("by_category") or {}).items()):
        lines.append(
            f"| {cat} | {c.get('n')} | {c.get('median_rmsd_A')} | "
            f"{c.get('n_sub2A')}/{c.get('n')} | `{c.get('dominant_mode')}` |"
        )
    lines += [
        "",
        "### Per target",
        "",
        "| ID | Cat | RMSD Å | Mode | Template | Wet-lab |",
        "|----|-----|-------:|------|----------|---------|",
    ]
    for r in report.get("results") or []:
        if r.get("status") != "ok":
            lines.append(
                f"| {r.get('id')} | {r.get('category')} | — | {r.get('status')} | — | {r.get('wetlab','')} |"
            )
            continue
        lines.append(
            f"| {r.get('id')} | {r.get('category')} | {r['fsot_product_rmsd_A']:.2f} | "
            f"`{r.get('primary_error_mode')}` | {r.get('template_pdb')} "
            f"({r.get('template_mode')}) | {r.get('wetlab','')} |"
        )
    lines += ["", "## Fix queue (from reality modes)", ""]
    for i, item in enumerate(report.get("fix_queue") or [], 1):
        lines.append(f"### {i}. `{item.get('mode')}` (n={item.get('n')})")
        lines.append(f"- {item.get('meaning')}")
        lines.append(f"- FSOT: {item.get('fsot_handle')}")
        lines.append("")
    lines += [
        "## Doctrine",
        "",
        "1. Match **reality** (experimental construct), not AF leaderboard.",
        "2. Localize error → one mechanism → full law only.",
        "3. Product path = measured homolog authority + residual physics.",
        "",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max", type=int, default=0)
    args = ap.parse_args(argv)

    cases = list(STRUCTURE_CASES)
    if args.max and args.max > 0:
        cases = cases[: args.max]

    print(f"REALITY MARGIN — product vs wet lab  n={len(cases)}", flush=True)
    print(f"{'id':<16}{'cat':<10}{'RMSD':>7}  mode  template", flush=True)
    print("-" * 72, flush=True)

    results: list[dict] = []
    for case in cases:
        print(f"  … {case['id']}", flush=True)
        try:
            row = run_case(case)
        except Exception as e:
            row = {
                "id": case["id"],
                "status": "exception",
                "category": case["category"],
                "error": str(e),
            }
        results.append(row)
        if row.get("status") == "ok":
            print(
                f"{row['id']:<16}{row['category']:<10}{row['fsot_product_rmsd_A']:7.2f}  "
                f"{row.get('primary_error_mode')}  {row.get('template_pdb')}",
                flush=True,
            )
        else:
            print(f"{row.get('id'):<16}{row.get('status')} {row.get('error','')}", flush=True)

    ok = [r for r in results if r.get("status") == "ok"]
    rmsds = [r["fsot_product_rmsd_A"] for r in ok]
    modes = [r.get("primary_error_mode") for r in ok if r.get("primary_error_mode")]
    mode_ct = Counter(modes)

    by_cat: dict[str, list] = {}
    for r in ok:
        by_cat.setdefault(r["category"], []).append(r)
    cat_sum = {}
    for cat, rs in by_cat.items():
        mc = Counter(r.get("primary_error_mode") for r in rs)
        cat_sum[cat] = {
            "n": len(rs),
            "median_rmsd_A": med([r["fsot_product_rmsd_A"] for r in rs]),
            "n_sub2A": sum(1 for r in rs if r.get("sub2A")),
            "n_sub3A": sum(1 for r in rs if r.get("sub3A")),
            "dominant_mode": mc.most_common(1)[0][0] if mc else None,
        }

    # Fix queue from margin rows shaped like error_margin_log
    margin_rows = []
    for r in ok:
        margin_rows.append(
            {
                "name": r["name"],
                "primary_error_mode": r.get("primary_error_mode"),
                "mode_scores": r.get("mode_scores") or {},
                "rmsd_A": r["fsot_product_rmsd_A"],
            }
        )
    try:
        fix_q = build_fix_queue(margin_rows)
    except Exception:
        # minimal queue from counts
        fix_q = [
            {
                "mode": m,
                "n": c,
                "meaning": (ERROR_CATALOG.get(m) or {}).get("meaning", m),
                "fsot_handle": (ERROR_CATALOG.get(m) or {}).get("fsot_handle", ""),
            }
            for m, c in mode_ct.most_common()
        ]

    summary = {
        "n_attempted": len(results),
        "n_ok": len(ok),
        "median_rmsd_A": med(rmsds),
        "mean_rmsd_A": float(np.mean(rmsds)) if rmsds else None,
        "n_sub2A": sum(1 for r in ok if r.get("sub2A")),
        "n_sub3A": sum(1 for r in ok if r.get("sub3A")),
        "dominant_error_mode": mode_ct.most_common(1)[0][0] if mode_ct else None,
        "mode_counts": dict(mode_ct),
        "by_category": cat_sum,
        "free_parameters": 0,
        "reality_median_target_A": REALITY_MEDIAN_TARGET,
        "reality_median_ok": (med(rmsds) or 99) <= REALITY_MEDIAN_TARGET,
    }

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "doctrine": "Match experimental wet-lab Cα; localize residual error; FSOT-only fixes",
        "authority_pin": "D1D38A",
        "summary": summary,
        "fix_queue": fix_q,
        "results": results,
    }
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_md(report)

    print("=" * 64, flush=True)
    print(
        f"REALITY median={summary['median_rmsd_A']:.3f} Å  "
        f"sub2={summary['n_sub2A']}/{summary['n_ok']}  "
        f"mode={summary['dominant_error_mode']}",
        flush=True,
    )
    for cat, c in cat_sum.items():
        print(f"  [{cat}] med={c['median_rmsd_A']:.2f} sub2={c['n_sub2A']}/{c['n']}", flush=True)
    print(f"Wrote {OUT_JSON}", flush=True)
    print(f"Wrote {OUT_MD}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
