#!/usr/bin/env python3
"""Interrogate the de-novo wall with FSOT: where is the distance information lost?

Compares FSOT's own predicted distance matrix (build_distogram -> proximity_to_distance)
against the native distance matrix, split into contacts (native < 8 A) and non-contacts.
If FSOT predicts contacts well but non-contacts poorly, the wall is precisely the
missing continuous non-contact distances (the many-body medium geometry), quantified.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import fetch_pdb  # noqa: E402
from fsot_structure_engine import build_distogram, proximity_to_distance, LONG_RANGE_GATE  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
TARGETS = [("1UBQ", "A", "Ubiquitin"), ("1LZ1", "A", "Lysozyme"),
           ("7RSA", "A", "RNase A"), ("2C9V", "A", "SOD1")]


def corr(a, b):
    if len(a) < 3 or np.std(a) < 1e-9 or np.std(b) < 1e-9:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def main() -> int:
    gate = int(LONG_RANGE_GATE)
    print(f"{'protein':<12}{'pairs':>7}{'contact r':>11}{'contact MAE':>13}"
          f"{'noncontact r':>14}{'noncontact MAE':>16}")
    print("-" * 73)
    for pdb, chain, name in TARGETS:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            continue
        seq, xyz = r
        M, props, regions, chars, iface = build_distogram(seq)
        if len(chars) != len(xyz):
            m = min(len(chars), len(xyz))
            xyz = xyz[:m]
            M = M[:m, :m]
        n = len(xyz)
        Dn = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
        Df = proximity_to_distance(M, props, regions, iface)[:n, :n]
        ii, jj = np.triu_indices(n, k=gate)
        dn = Dn[ii, jj]
        df = Df[ii, jj]
        cmask = dn < 8.0
        ncmask = ~cmask
        print(f"{name:<12}{len(dn):>7}"
              f"{corr(df[cmask], dn[cmask]):>11.2f}"
              f"{np.mean(np.abs(df[cmask] - dn[cmask])):>13.2f}"
              f"{corr(df[ncmask], dn[ncmask]):>14.2f}"
              f"{np.mean(np.abs(df[ncmask] - dn[ncmask])):>16.2f}")
    print("-" * 73)
    print("r = Pearson correlation of FSOT-predicted vs native distance; MAE in A.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
