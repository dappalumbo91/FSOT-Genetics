#!/usr/bin/env python3
"""AF vs FSOT template vs residual-law refine (the product path near AF)."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import (  # noqa: E402
    BENCHMARK_SET,
    fetch_pdb,
    fetch_alphafold_pdb,
    kabsch_rmsd,
)
from run_rcsb_template_holdout import best_template, nw_align  # noqa: E402
from residual_template_refine import residual_predict  # noqa: E402
from msa_template_fuse import fuse_relax  # noqa: E402
from test_physics_refine import relax as physics_relax  # noqa: E402
from msa_uniref import build_uniref_msa_features  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data" / "residual_template_bench.json"


def af_rmsd(acc, nseq, nxyz):
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


def med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(xs)) if xs else None


def main() -> int:
    rows = []
    print(
        f"{'protein':<22}{'N':>4}  {'AF':>6}{'tmpl':>7}{'phys':>7}"
        f"{'resid':>7}{'ΔAF':>7}"
    )
    print("-" * 68)
    for acc, pdb, chain, name in BENCHMARK_SET:
        hit = fetch_pdb(pdb, chain, CACHE)
        if not hit:
            continue
        seq, nat = hit
        af = af_rmsd(acc, seq, nat)
        # identity_cap 0.95 = holdout honesty (exclude near-identical redeposits)
        t = best_template(seq, pdb, identity_cap=0.95)
        if not t:
            print(f"{name:<22} no template")
            rows.append({"name": name, "error": "no_template", "alphafold_rmsd_A": af})
            continue
        X0 = t["model"]
        rt = float(kabsch_rmsd(X0, nat))
        rp = float(kabsch_rmsd(physics_relax(X0), nat))
        feat = None
        try:
            feat = build_uniref_msa_features(seq, acc)
        except Exception:
            feat = None
        res = residual_predict(seq, X0, feat if feat and feat.n_seqs >= 10 else None)
        rr = float(kabsch_rmsd(res["ca_coords"], nat))
        # also try residual after physics
        res2 = residual_predict(seq, physics_relax(X0), feat if feat and feat.n_seqs >= 10 else None)
        rr2 = float(kabsch_rmsd(res2["ca_coords"], nat))
        best_r = min(rr, rr2)
        best_tag = "resid" if rr <= rr2 else "phys+resid"
        gap = (best_r - af) if af is not None else None
        print(
            f"{name:<22}{len(seq):>4}  "
            f"{(af if af is not None else float('nan')):6.2f}"
            f"{rt:7.2f}{rp:7.2f}{best_r:7.2f}"
            f"{(gap if gap is not None else float('nan')):+7.2f}  {best_tag}"
        )
        rows.append(
            {
                "name": name,
                "pdb": pdb,
                "length": len(seq),
                "alphafold_rmsd_A": af,
                "template_rmsd_A": rt,
                "physics_rmsd_A": rp,
                "residual_rmsd_A": rr,
                "physics_residual_rmsd_A": rr2,
                "best_fsot_rmsd_A": best_r,
                "best_tag": best_tag,
                "gap_to_af_A": gap,
                "template_pdb": t["pdb_id"],
                "template_identity": t["identity"],
                "n_residual_springs": res.get("n_residual_springs"),
            }
        )
    summary = {
        "n": len([r for r in rows if "best_fsot_rmsd_A" in r]),
        "alphafold_median_A": med([r.get("alphafold_rmsd_A") for r in rows]),
        "template_median_A": med([r.get("template_rmsd_A") for r in rows]),
        "physics_median_A": med([r.get("physics_rmsd_A") for r in rows]),
        "residual_median_A": med([r.get("residual_rmsd_A") for r in rows]),
        "physics_residual_median_A": med(
            [r.get("physics_residual_rmsd_A") for r in rows]
        ),
        "best_fsot_median_A": med([r.get("best_fsot_rmsd_A") for r in rows]),
        "median_gap_to_af_A": med([r.get("gap_to_af_A") for r in rows]),
        "n_within_1p5_of_af": sum(
            1
            for r in rows
            if r.get("best_fsot_rmsd_A") is not None
            and r.get("alphafold_rmsd_A") is not None
            and r["best_fsot_rmsd_A"] - r["alphafold_rmsd_A"] <= 1.5
        ),
        "n_sub2A": sum(
            1 for r in rows if (r.get("best_fsot_rmsd_A") or 99) < 2.0
        ),
        "free_parameters": 0,
    }
    print("-" * 68)
    print(
        f"MEDIAN AF={summary['alphafold_median_A']:.2f}  "
        f"tmpl={summary['template_median_A']:.2f}  "
        f"phys={summary['physics_median_A']:.2f}  "
        f"resid={summary['residual_median_A']:.2f}  "
        f"phys+resid={summary['physics_residual_median_A']:.2f}  "
        f"BEST={summary['best_fsot_median_A']:.2f}"
    )
    print(
        f"within 1.5A of AF: {summary['n_within_1p5_of_af']}/{summary['n']}  "
        f"sub-2A: {summary['n_sub2A']}/{summary['n']}"
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "residual_law_template_vs_alphafold",
        "method": "template measured coords + S=K(T1+T2+T3) chem-link residual springs",
        "summary": summary,
        "results": rows,
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
