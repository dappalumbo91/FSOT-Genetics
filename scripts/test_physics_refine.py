#!/usr/bin/env python3
"""Close the gap: FSOT physics relaxation of the template model.

The single-template model's residual error concentrates in gap/loop regions where
linear interpolation leaves stretched bonds and steric clashes. This applies a
light, zero-parameter energy minimization: idealize CA-CA bonds (3.8 A), relieve
clashes (CA-CA non-bonded floor), anchored weakly to the template so the correct
topology is preserved. No rescaling, no trained weights. Scored on the native.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import BENCHMARK_SET, fetch_pdb, kabsch_rmsd  # noqa: E402
from run_rcsb_template_holdout import best_template  # noqa: E402
from fsot_structure_engine import CA_CA  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CLASH = 4.0   # CA-CA non-bonded floor (geometry constant)


def relax(X0, iters=120, lr=0.08, w_anchor=0.05):
    X = X0.copy()
    n = len(X)
    for _ in range(iters):
        G = w_anchor * (X - X0)
        # bond term: |CA_i - CA_{i+1}| -> 3.8
        d = X[1:] - X[:-1]
        L = np.linalg.norm(d, axis=1) + 1e-9
        f = ((L - CA_CA) / L)[:, None] * d
        G[:-1] -= f
        G[1:] += f
        # clash term: push apart non-bonded CA closer than CLASH
        diff = X[:, None, :] - X[None, :, :]
        D = np.linalg.norm(diff, axis=2)
        np.fill_diagonal(D, 1e9)
        for k in (1,):
            idx = np.arange(n - k)
            D[idx, idx + k] = 1e9
            D[idx + k, idx] = 1e9
        mask = D < CLASH
        if mask.any():
            coef = np.where(mask, (D - CLASH) / (D + 1e-9), 0.0)
            G += np.einsum("ij,ijk->ik", coef, diff)
        X = X - lr * G
    return X


def main() -> int:
    print(f"{'protein':<22}{'template':>10}{'+relax':>9}{'delta':>8}  tmpl")
    print("-" * 55)
    base, refined = [], []
    for acc, pdb, chain, name in BENCHMARK_SET:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            continue
        nseq, nxyz = r
        t = best_template(nseq, pdb)
        if not t:
            print(f"{name:<22}{'n/a':>10}")
            continue
        X0 = t["model"]
        r0 = kabsch_rmsd(X0, nxyz)
        r1 = kabsch_rmsd(relax(X0), nxyz)
        base.append(r0)
        refined.append(r1)
        print(f"{name:<22}{r0:>10.2f}{r1:>9.2f}{r1 - r0:>+8.2f}  {t['pdb_id']} id={t['identity']:.2f}")
    print("-" * 55)
    print(f"{'MEDIAN':<22}{np.median(base):>10.2f}{np.median(refined):>9.2f}"
          f"{np.median(refined) - np.median(base):>+8.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
