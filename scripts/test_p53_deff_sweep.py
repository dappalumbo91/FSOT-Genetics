#!/usr/bin/env python3
"""D_eff-level test for p53: is a DNA-regulatory switch folded at the wrong D_eff?

Sweeps the observer bulk dimension (the D_eff at which the observer solidifies the
fold) for p53 vs non-regulatory controls, scored on the DNA-bound native. If the
regulator prefers a different D_eff than the enzyme/transport controls, that
supports modeling it at a function-specific effective dimension.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import fetch_pdb, kabsch_rmsd  # noqa: E402
from fsot_structure_engine import predict_ca_coords  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
LADDER = [8, 9, 13, 14, 18, 25]  # ChemLink D_eff rungs + base law
TARGETS = [
    ("1TUP", "A", "p53 (DNA on/off regulator)"),
    ("1LZ1", "A", "lysozyme (enzyme)"),
    ("1A3N", "A", "hemoglobin (transport)"),
    ("1UBQ", "A", "ubiquitin (signal tag)"),
]


def main() -> int:
    header = "protein".ljust(30) + "".join(f"D{d:>6}" for d in LADDER) + "   best"
    print(header)
    print("-" * len(header))
    for pdb, chain, name in TARGETS:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            print(f"{name:<30} fetch failed")
            continue
        seq, nxyz = r
        rmsds = []
        for d in LADDER:
            try:
                m = predict_ca_coords(seq, canonicalize_chirality=True, observer_bulk_dim=d)
                rmsds.append(kabsch_rmsd(m["ca_coords"], nxyz))
            except Exception:
                rmsds.append(float("nan"))
        best_d = LADDER[int(np.nanargmin(rmsds))]
        cells = "".join(f"{v:7.2f}" for v in rmsds)
        print(f"{name:<30}{cells}   D{best_d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
