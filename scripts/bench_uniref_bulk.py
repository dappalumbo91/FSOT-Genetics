#!/usr/bin/env python3
"""Bulk with UniRef MSA bridged into tertiary chem-link (hits/δψ)."""

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
    kabsch_rmsd,
)
from fsot_structure_engine import predict_ca_coords  # noqa: E402
from msa_uniref import build_uniref_msa_features  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data" / "uniref_bulk_bench.json"
PRIOR = 13.65
PFAM = {
    "P69905": "PF00042",
    "P68871": "PF00042",
    "P00918": "PF00194",
    "P00441": "PF00080",
    "P61626": "PF00062",
    "P61823": "PF00074",
    "P0CG47": "PF00240",
    "P01308": "PF00049",
    "P04637": "PF00870",
    "P0DP23": "PF00036",
}


def main() -> int:
    rows = []
    print(f"{'protein':<22}{'N':>4}  {'single':>7}{'uniref':>7}{'d_s':>7}{'d_u':>7}")
    print("-" * 62)
    for acc, pdb, chain, name in BENCHMARK_SET:
        hit = fetch_pdb(pdb, chain, CACHE)
        if not hit:
            continue
        seq, nat = hit
        b = predict_ca_coords(
            seq, rounds=28, canonicalize_chirality=True, observer_bulk_dim=25
        )
        rb = float(kabsch_rmsd(b["ca_coords"], nat))
        ru = None
        try:
            feat = build_uniref_msa_features(seq, acc, pfam=PFAM.get(acc))
            if feat is not None and feat.n_seqs >= 10:
                m = predict_ca_coords(
                    seq,
                    rounds=28,
                    mode="msa",
                    msa_features=feat,
                    canonicalize_chirality=True,
                    observer_bulk_dim=25,
                )
                ru = float(kabsch_rmsd(m["ca_coords"], nat))
        except Exception as exc:  # noqa: BLE001
            print(f"  err {name}: {exc}")
        print(
            f"{name:<22}{len(seq):>4}  {rb:7.2f}"
            f"{(ru if ru is not None else float('nan')):7.2f}"
            f"{rb - PRIOR:+7.2f}"
            f"{((ru - PRIOR) if ru is not None else float('nan')):+7.2f}"
        )
        rows.append(
            {
                "name": name,
                "length": len(seq),
                "single_A": rb,
                "uniref_A": ru,
                "delta_single_vs_prior": rb - PRIOR,
                "delta_uniref_vs_prior": (ru - PRIOR) if ru is not None else None,
            }
        )
    singles = [r["single_A"] for r in rows]
    unis = [r["uniref_A"] for r in rows if r["uniref_A"] is not None]
    med_s = float(np.median(singles))
    med_u = float(np.median(unis)) if unis else None
    print("-" * 62)
    print(
        f"MEDIAN single={med_s:.2f} uniref={med_u} prior={PRIOR} "
        f"d_s={med_s - PRIOR:+.2f} d_u={(med_u - PRIOR) if med_u else None}"
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "prior_median_A": PRIOR,
        "median_single_A": med_s,
        "median_uniref_A": med_u,
        "delta_single_A": med_s - PRIOR,
        "delta_uniref_A": (med_u - PRIOR) if med_u is not None else None,
        "n_uniref_beats_single": sum(
            1
            for r in rows
            if r["uniref_A"] is not None and r["uniref_A"] < r["single_A"] - 0.05
        ),
        "free_parameters": 0,
        "results": rows,
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
