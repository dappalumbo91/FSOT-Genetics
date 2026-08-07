#!/usr/bin/env python3
"""Cofactor-constrained fold: do metal coordination nodes tighten the structure?

Takes the zero-parameter bulk fold, then forces each metal's co-coordinating
residues to coordination-shell geometry (weighted SMACOF anchored to the bulk
fold, with cofactor pairs overridden to a metal-site CA-CA scale). Ca/Mg nodes
are bioelectrical/electrostatic (Electromagnetism D9); Zn/Cu/Fe are dative
(Atomic_Physics D7). Coordinating residues come from the real structure (data in,
zero trained weights). Scored vs the metal-bound native.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import kabsch_rmsd  # noqa: E402
from fsot_structure_engine import predict_ca_coords, CA_CA  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "T", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}
METALS = {"ZN", "FE", "CU", "CA", "MG", "MN", "CO", "NI"}
SITE_CA = 8.0  # coordination-shell CA-CA scale (geometry constant, not fitted)


def parse(text, chain):
    seq, xyz, nums, allatom = [], [], [], []
    for line in text.splitlines():
        tag = line[:6].strip()
        if tag not in ("ATOM", "HETATM"):
            continue
        res = line[17:20].strip()
        try:
            p = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        if tag == "ATOM" and line[21] == chain and res in AA3:
            allatom.append((line[22:26].strip(), p))
            if line[12:16].strip() == "CA":
                seq.append(AA3[res])
                xyz.append(p)
                nums.append(line[22:26].strip())
        elif tag == "HETATM" and res in METALS:
            allatom.append(("METAL:" + res, p))
    return "".join(seq), np.array(xyz), nums, allatom


def metal_nodes(text, chain, nums):
    seq, xyz, _n, allatom = parse(text, chain)
    num_to_idx = {nm: i for i, nm in enumerate(nums)}
    metals = [(r.split(":")[1], p) for r, p in allatom if r.startswith("METAL:")]
    prot = [(nm, p) for nm, p in allatom if not nm.startswith("METAL:")]
    nodes = []
    for mres, mp in metals:
        idxs = set()
        for nm, p in prot:
            if np.linalg.norm(p - mp) <= 3.0 and nm in num_to_idx:
                idxs.add(num_to_idx[nm])
        if len(idxs) >= 2:
            nodes.append(sorted(idxs))
    return nodes


def smacof(X0, pairs, target, iters=60, w_cof=25.0, w_bb=12.0):
    n = len(X0)
    D0 = np.linalg.norm(X0[:, None, :] - X0[None, :, :], axis=2)
    W = np.ones((n, n))
    Dt = D0.copy()
    for i in range(n - 1):
        W[i, i + 1] = W[i + 1, i] = w_bb
        Dt[i, i + 1] = Dt[i + 1, i] = CA_CA
    for (i, j) in pairs:
        W[i, j] = W[j, i] = w_cof
        Dt[i, j] = Dt[j, i] = target
    np.fill_diagonal(W, 0.0)
    V = -W.copy()
    np.fill_diagonal(V, W.sum(1))
    Vp = np.linalg.pinv(V)
    X = X0.copy()
    for _ in range(iters):
        Dx = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        np.fill_diagonal(Dx, 1.0)
        B = -W * Dt / Dx
        np.fill_diagonal(B, 0.0)
        np.fill_diagonal(B, -B.sum(1))
        X = Vp @ (B @ X)
    return X


def main() -> int:
    targets = [("1CLL", "A", "Calmodulin (4x Ca)"), ("1CA2", "A", "Carbonic anhydrase (Zn)"),
               ("2C9V", "A", "SOD1 (Cu+Zn)"), ("1TUP", "A", "p53 (Zn)")]
    print(f"{'protein':<26}{'bulk':>8}{'+cofactor':>11}{'delta':>8}  nodes")
    print("-" * 62)
    for pdb, chain, name in targets:
        fp = CACHE / f"{pdb}.pdb"
        if not fp.exists():
            print(f"{name:<26} not cached")
            continue
        text = fp.read_text(encoding="utf-8", errors="replace")
        seq, native, nums, _ = parse(text, chain)
        if len(seq) < 20:
            print(f"{name:<26} parse fail")
            continue
        X0 = predict_ca_coords(seq, canonicalize_chirality=True, observer_bulk_dim=25)["ca_coords"]
        rmsd0 = kabsch_rmsd(X0, native)
        nodes = metal_nodes(text, chain, nums)
        pairs = [(i, j) for node in nodes for a, i in enumerate(node) for j in node[a + 1:]]
        if not pairs:
            print(f"{name:<26}{rmsd0:>8.2f}{'--':>11}{'':>8}  no node")
            continue
        Xc = smacof(X0, pairs, SITE_CA)
        rmsd1 = kabsch_rmsd(Xc, native)
        print(f"{name:<26}{rmsd0:>8.2f}{rmsd1:>11.2f}{rmsd1 - rmsd0:>+8.2f}  "
              f"{len(nodes)} node(s), {len(pairs)} pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
