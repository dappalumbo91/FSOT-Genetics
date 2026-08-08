#!/usr/bin/env python3
"""Per-residue confidence (our pLDDT analog): does provenance predict local error?

For each template-covered protein, classify every residue by how it was modeled:
  identical  - aligned to template, same amino acid (real coord of the same residue)
  mutated    - aligned to template, different amino acid
  gap        - not aligned (interpolated loop)
Then measure the actual per-residue deviation from the native after superposition.
If deviation stratifies (identical << mutated << gap), we have a usable confidence
signal telling a user which regions of the model to trust - zero trained weights.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import BENCHMARK_SET, fetch_pdb  # noqa: E402
from run_rcsb_template_holdout import (  # noqa: E402
    best_template, fetch_template_pdb, chains_of, nw_align,
)
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"


def kabsch_dev(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return np.linalg.norm(Pc @ R.T - Qc, axis=1)


def provenance(nseq, template_pdb):
    text = fetch_template_pdb(template_pdb)
    best = None
    for ch in chains_of(text):
        try:
            tseq, _ = parse_pdb_ca(text, ch)
        except Exception:
            continue
        pairs = nw_align(nseq, tseq)
        if not pairs:
            continue
        ident = sum(1 for a, b in pairs if nseq[a] == tseq[b]) / len(pairs)
        cov = len(pairs) / len(nseq)
        if best is None or cov * ident > best[0]:
            best = (cov * ident, pairs, tseq)
    if best is None:
        return None
    _, pairs, tseq = best
    cls = ["gap"] * len(nseq)
    for a, b in pairs:
        cls[a] = "identical" if nseq[a] == tseq[b] else "mutated"
    return cls


def main() -> int:
    buckets = {"identical": [], "mutated": [], "gap": []}
    print(f"{'protein':<22}{'identical':>10}{'mutated':>9}{'gap':>7}")
    print("-" * 48)
    for acc, pdb, chain, name in BENCHMARK_SET:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            continue
        nseq, nxyz = r
        t = best_template(nseq, pdb)
        if not t:
            continue
        cls = provenance(nseq, t["pdb_id"])
        if cls is None:
            continue
        dev = kabsch_dev(t["model"], nxyz)
        per = {"identical": [], "mutated": [], "gap": []}
        for i, c in enumerate(cls):
            per[c].append(dev[i])
            buckets[c].append(dev[i])
        m = {k: (np.mean(v) if v else float("nan")) for k, v in per.items()}
        print(f"{name:<22}{m['identical']:>10.2f}{m['mutated']:>9.2f}{m['gap']:>7.2f}")
    print("-" * 48)
    agg = {k: np.array(v) for k, v in buckets.items()}
    print(f"{'ALL residues mean':<22}"
          f"{agg['identical'].mean():>10.2f}{agg['mutated'].mean():>9.2f}{agg['gap'].mean():>7.2f}")
    print(f"{'(n residues)':<22}{len(agg['identical']):>10d}{len(agg['mutated']):>9d}{len(agg['gap']):>7d}")
    # fraction of large errors (>3 A) captured by the low-confidence (gap+mutated) flag
    big = np.concatenate([agg[k] for k in buckets])
    thr = 3.0
    lowconf = np.concatenate([agg["gap"], agg["mutated"]])
    frac_big_lowconf = float((lowconf > thr).sum()) / max(1, int((big > thr).sum()))
    print(f"\n>3A errors flagged as low-confidence (gap+mutated): "
          f"{frac_big_lowconf*100:.0f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
