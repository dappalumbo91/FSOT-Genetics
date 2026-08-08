#!/usr/bin/env python3
"""Attack the de-novo 11 A ceiling with evolutionary coevolution (deep MSA).

Single-sequence folding is capped ~11 A because one sequence can't say which
residues contact. A deep family alignment CAN: co-varying columns are in contact
(direct-coupling / MI signal). This predicts contacts from the full Pfam MSA, folds
them via shortest-path MDS (contacts + backbone), and compares to the single-
sequence bulk fold. Zero trained weights; real evolutionary data as input.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import fetch_pdb, parse_pdb_ca, kabsch_rmsd  # noqa: E402
from run_rcsb_template_holdout import nw_align  # noqa: E402
from fsot_structure_engine import predict_ca_coords, CA_CA  # noqa: E402
from variant_conservation import fetch_msa, parse_stockholm  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
AA = "ARNDCQEGHILKMFPSTWYV-"
A2I = {a: i for i, a in enumerate(AA)}

# classic deep-MSA test proteins (pdb, chain, pfam)
TARGETS = [("1UBQ", "A", "PF00240", "Ubiquitin"),
           ("7RSA", "A", "PF00074", "RNase A")]


def query_columns(seq, rows):
    """Return per-query-position alignment column vectors (as int arrays)."""
    best = None
    for r in rows[:300]:  # bound ref search to a small candidate set
        ung = r.replace(".", "").replace("-", "").upper()
        if not (0.7 * len(seq) <= len(ung) <= 1.5 * len(seq)):
            continue
        pairs = nw_align(seq, ung)
        score = sum(1 for a, b in pairs if seq[a] == ung[b])
        if best is None or score > best[0]:
            best = (score, r)
        if best[0] > 0.8 * len(seq):
            break
    if best is None:
        return None
    ref = best[1]
    ung_to_col = [c for c, ch in enumerate(ref) if ch not in ".-"]
    ref_ung = ref.replace(".", "").replace("-", "").upper()
    qpairs = dict(nw_align(seq, ref_ung))
    mat = np.full((len(rows), len(seq)), A2I["-"], dtype=np.int8)
    listrows = [list(r) for r in rows]
    for qi, ui in qpairs.items():
        if ui >= len(ung_to_col):
            continue
        col = ung_to_col[ui]
        for ri, r in enumerate(listrows):
            if col < len(r):
                mat[ri, qi] = A2I.get(r[col].upper(), A2I["-"])
    covered = sorted(qpairs.keys())
    return mat, covered


def mutual_information(mat, covered):
    n = mat.shape[1]
    gap = A2I["-"]
    # subsample rows for speed
    rng = np.random.default_rng(0)
    if mat.shape[0] > 4000:
        mat = mat[rng.choice(mat.shape[0], 4000, replace=False)]
    mi = np.zeros((n, n))
    cols = {a: mat[:, a].astype(np.int64) for a in covered}
    for ai, a in enumerate(covered):
        ca = cols[a]
        va = ca != gap
        for b in covered[ai + 1:]:
            cb = cols[b]
            valid = va & (cb != gap)
            if valid.sum() < 20:
                continue
            xa, xb = ca[valid], cb[valid]
            joint = np.bincount(xa * 20 + xb, minlength=400).reshape(20, 20).astype(float)
            tot = joint.sum()
            if tot < 20:
                continue
            joint /= tot
            pa = joint.sum(1)
            pb = joint.sum(0)
            nz = joint > 0
            m = joint[nz] * np.log((joint[nz]) / (np.outer(pa, pb)[nz] + 1e-12))
            mi[a, b] = mi[b, a] = float(m.sum())
    mi_mean = mi.mean()
    row_mean = mi.mean(1, keepdims=True)
    apc = (row_mean @ row_mean.T) / (mi_mean + 1e-12)
    return mi - apc


def fold_from_contacts(n, contacts):
    # graph: backbone 3.8, contacts 8.0 -> all-pairs shortest path -> MDS
    INF = 1e6
    D = np.full((n, n), INF)
    np.fill_diagonal(D, 0.0)
    for i in range(n - 1):
        D[i, i + 1] = D[i + 1, i] = CA_CA
    for (i, j) in contacts:
        D[i, j] = D[j, i] = min(D[i, j], 8.0)
    try:
        from scipy.sparse.csgraph import shortest_path  # noqa: WPS433
        Dg = shortest_path(D, method="D", directed=False)
    except Exception:
        Dg = D.copy()
        for k in range(n):
            Dg = np.minimum(Dg, Dg[:, [k]] + Dg[[k], :])
    Dg[~np.isfinite(Dg)] = Dg[np.isfinite(Dg)].max() * 1.5
    D2 = Dg ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1][:3]
    return V[:, idx] * np.sqrt(np.maximum(w[idx], 0))


def main() -> int:
    print(f"{'protein':<14}{'MSA':>7}{'bulk':>8}{'coevolution':>13}{'delta':>8}")
    print("-" * 50)
    for pdb, chain, pfam, name in TARGETS:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            continue
        seq, native = r
        n = len(seq)
        bulk = kabsch_rmsd(predict_ca_coords(seq, canonicalize_chirality=True,
                                             observer_bulk_dim=25)["ca_coords"], native)
        rows = parse_stockholm(fetch_msa(pfam, "full"))
        qc = query_columns(seq, rows)
        if qc is None:
            print(f"{name:<14}{'no MSA':>7}")
            continue
        mat, covered = qc
        cij = mutual_information(mat, covered)
        # rank non-local pairs, take top 1.5 L
        pairs = [(cij[i, j], i, j) for i in covered for j in covered
                 if j > i + 4]
        pairs.sort(reverse=True)
        contacts = [(i, j) for _, i, j in pairs[: int(1.5 * n)]]
        coev = kabsch_rmsd(fold_from_contacts(n, contacts), native)
        print(f"{name:<14}{len(rows):>7}{bulk:>8.2f}{coev:>13.2f}{coev - bulk:>+8.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
