#!/usr/bin/env python3
"""Product path (template + physics/fuse) vs AlphaFold vs bulk."""

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
from template_select_ss import best_template_ss  # noqa: E402
from fsot_structure_engine import predict_ca_coords  # noqa: E402
from msa_template_fuse import fuse_predict  # noqa: E402
from msa_uniref import build_uniref_msa_features  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data" / "product_vs_alphafold.json"


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
    xs = [x for x in xs if x is not None]
    return float(np.median(xs)) if xs else None


def main() -> int:
    rows = []
    print(f"{'protein':<22}{'AF':>6}{'tmpl':>7}{'product':>8}{'bulk':>7}")
    print("-" * 52)
    for acc, pdb, ch, name in BENCHMARK_SET:
        hit = fetch_pdb(pdb, ch, CACHE)
        if not hit:
            continue
        seq, nat = hit
        af = af_rmsd(acc, seq, nat)
        t = best_template(seq, pdb, identity_cap=0.95)
        bulk = predict_ca_coords(
            seq, rounds=24, canonicalize_chirality=True, observer_bulk_dim=25
        )
        rb = float(kabsch_rmsd(bulk["ca_coords"], nat))
        if t:
            feat = None
            try:
                feat = build_uniref_msa_features(seq, acc)
                if feat.n_seqs < 10:
                    feat = None
            except Exception:
                feat = None
            prod = fuse_predict(seq, t["model"], feat)
            rp = float(kabsch_rmsd(prod["ca_coords"], nat))
            rt = float(kabsch_rmsd(t["model"], nat))
            tid = t["pdb_id"]
            regime = prod.get("regime")
        else:
            rp, rt, tid, regime = rb, None, None, "bulk_fallback"
        print(
            f"{name:<22}"
            f"{(af if af is not None else float('nan')):6.2f}"
            f"{(rt if rt is not None else float('nan')):7.2f}"
            f"{rp:8.2f}{rb:7.2f}  {regime}"
        )
        rows.append(
            {
                "name": name,
                "alphafold_rmsd_A": af,
                "fsot_template_rmsd_A": rt,
                "fsot_product_rmsd_A": rp,
                "fsot_bulk_rmsd_A": rb,
                "template_pdb": tid,
                "regime": regime,
            }
        )
    afm = med([r["alphafold_rmsd_A"] for r in rows])
    tm = med([r["fsot_template_rmsd_A"] for r in rows])
    pm = med([r["fsot_product_rmsd_A"] for r in rows])
    bm = med([r["fsot_bulk_rmsd_A"] for r in rows])
    within = sum(
        1
        for r in rows
        if r["fsot_product_rmsd_A"] is not None
        and r["alphafold_rmsd_A"] is not None
        and r["fsot_product_rmsd_A"] - r["alphafold_rmsd_A"] <= 1.5
    )
    sub2 = sum(1 for r in rows if (r["fsot_product_rmsd_A"] or 99) < 2)
    print("-" * 52)
    print(f"MEDIAN AF={afm:.2f} tmpl={tm:.2f} product={pm:.2f} bulk={bm:.2f}")
    print(f"product within 1.5A of AF: {within}/{len(rows)}  sub-2A: {sub2}/{len(rows)}")
    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "n": len(rows),
            "alphafold_median_A": afm,
            "fsot_template_median_A": tm,
            "fsot_product_median_A": pm,
            "fsot_bulk_median_A": bm,
            "product_within_1p5A_of_alphafold": within,
            "product_sub2A": sub2,
            "free_parameters": 0,
        },
        "results": rows,
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
