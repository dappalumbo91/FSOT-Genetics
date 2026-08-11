#!/usr/bin/env python3
"""Measure bulk RMSD after chem-link geometry + MSA→tertiary bridge."""

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
from msa_pipeline import build_msa_features  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data" / "chemlink_bulk_bench.json"
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
# Last published bulk median on this set (pre chem-link geometry)
PRIOR_MEDIAN = 13.65


def main() -> int:
    rows = []
    print(f"{'protein':<22}{'N':>4}  {'single':>7}{'msa':>7}  links")
    print("-" * 60)
    for acc, pdb, chain, name in BENCHMARK_SET:
        hit = fetch_pdb(pdb, chain, CACHE)
        if not hit:
            continue
        seq, nat = hit
        b = predict_ca_coords(
            seq, rounds=24, canonicalize_chirality=True, observer_bulk_dim=25
        )
        rb = float(kabsch_rmsd(b["ca_coords"], nat))
        hist = b.get("chem_link_histogram") or {}
        rmsa = None
        pf = PFAM.get(acc)
        if pf:
            try:
                feat = build_msa_features(seq, pfam=pf)
                if feat.n_seqs > 0:
                    m = predict_ca_coords(
                        seq,
                        rounds=24,
                        mode="msa",
                        msa_features=feat,
                        canonicalize_chirality=True,
                        observer_bulk_dim=25,
                    )
                    rmsa = float(kabsch_rmsd(m["ca_coords"], nat))
            except Exception as exc:  # noqa: BLE001
                print(f"  msa fail {name}: {exc}")
        print(
            f"{name:<22}{len(seq):>4}  {rb:7.2f}"
            f"{(rmsa if rmsa is not None else float('nan')):7.2f}  "
            f"ss={hist.get('disulfide_covalent',0)//2} "
            f"salt={hist.get('salt_bridge_electrostatic',0)//2} "
            f"pack={hist.get('hydrophobic_packing',0)//2}"
        )
        rows.append(
            {
                "name": name,
                "length": len(seq),
                "single_A": rb,
                "msa_A": rmsa,
                "chem_link_histogram": hist,
            }
        )
    singles = [r["single_A"] for r in rows]
    msas = [r["msa_A"] for r in rows if r["msa_A"] is not None]
    med_s = float(np.median(singles))
    med_m = float(np.median(msas)) if msas else None
    print("-" * 60)
    print(
        f"MEDIAN single={med_s:.2f}  msa={med_m}  "
        f"prior={PRIOR_MEDIAN:.2f}  Δsingle={med_s - PRIOR_MEDIAN:+.2f}"
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "prior_bulk_median_A": PRIOR_MEDIAN,
        "median_single_A": med_s,
        "median_msa_A": med_m,
        "delta_vs_prior_A": med_s - PRIOR_MEDIAN,
        "free_parameters": 0,
        "results": rows,
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {OUT}")
    # gate check
    import subprocess

    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "verify_cross.py")],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )
    print("verify_cross exit", r.returncode)
    if r.returncode != 0:
        print(r.stdout[-1500:] if r.stdout else r.stderr[-1500:])
    return r.returncode


if __name__ == "__main__":
    raise SystemExit(main())
