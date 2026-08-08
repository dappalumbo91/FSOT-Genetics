#!/usr/bin/env python3
"""Close the AlphaFold gap: multi-template consensus averaging.

Single-template transfer carries one crystal's idiosyncrasies (~1.2 A vs AF ~0.5).
Averaging several real homolog structures (superposed) cancels per-template noise -
ensemble averaging, still zero trained weights. Compares best-single vs structural
medoid vs consensus-average, scored on the experimental native.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import BENCHMARK_SET, fetch_pdb, kabsch_rmsd  # noqa: E402
from test_structural_selection import candidates, medoid  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"


def superpose(mc, ref):
    """Rotate centered coords mc onto centered ref (Kabsch)."""
    H = mc.T @ ref
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return mc @ R.T


def consensus(models, weights=None):
    ref = models[0] - models[0].mean(0)
    stack = [ref]
    for m in models[1:]:
        stack.append(superpose(m - m.mean(0), ref))
    if weights is None:
        return np.mean(stack, axis=0)
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    return np.tensordot(w, np.asarray(stack), axes=(0, 0))


def main() -> int:
    print(f"{'protein':<22}{'#tmpl':>6}{'single':>9}{'medoid':>9}{'consensus':>11}{'wtd-cons':>10}")
    print("-" * 67)
    singles, medoids, consensuses, weighteds = [], [], [], []
    for acc, pdb, chain, name in BENCHMARK_SET:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            continue
        nseq, nxyz = r
        cand = candidates(nseq, pdb, cap_hi=0.95)  # match strict head-to-head cap
        if not cand:
            print(f"{name:<22}{0:>6}   no template")
            continue
        models = [c[3] for c in cand]
        gi = int(np.argmax([c[0] * c[1] for c in cand]))
        mi = medoid(models)
        s_rmsd = kabsch_rmsd(models[gi], nxyz)
        m_rmsd = kabsch_rmsd(models[mi], nxyz)
        c_rmsd = kabsch_rmsd(consensus(models), nxyz)
        wts = [(c[0] * c[1]) ** 8 for c in cand]  # standout templates dominate
        w_rmsd = kabsch_rmsd(consensus(models, wts), nxyz)
        singles.append(s_rmsd)
        medoids.append(m_rmsd)
        consensuses.append(c_rmsd)
        weighteds.append(w_rmsd)
        print(f"{name:<22}{len(cand):>6}{s_rmsd:>9.2f}{m_rmsd:>9.2f}{c_rmsd:>11.2f}{w_rmsd:>10.2f}")
    print("-" * 67)
    print(f"{'MEDIAN':<22}{'':>6}{np.median(singles):>9.2f}"
          f"{np.median(medoids):>9.2f}{np.median(consensuses):>11.2f}{np.median(weighteds):>10.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
