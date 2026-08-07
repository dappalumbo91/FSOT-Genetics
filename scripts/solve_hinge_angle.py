#!/usr/bin/env python3
"""Solve the interdomain hinge angle of 1DUH with the FSOT compactness law.

Real data gives two rigid domains (template transfer, each ~5 A near-native).
The relative bend is the free/context variable. Instead of inheriting the
template's (wrong-context) angle, we let FSOT's zero-parameter radius-of-gyration
prediction choose the angle, then measure RMSD to the native crystal.

Hinge axis + pivot are taken from the TEMPLATE geometry (native never used to set
the angle); native is only the scoring target.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_rna_template_probe import _pdb, parse_rna_c1, rna_chains  # noqa: E402
from run_rcsb_template_holdout import nw_align  # noqa: E402
from run_fsot_vs_alphafold_structure import kabsch_rmsd  # noqa: E402

import fsot_compute as fc  # noqa: E402

PI = float(fc.PI)


def rg(coords):
    c = coords - coords.mean(0)
    return float(np.sqrt((c * c).sum(1).mean()))


def pca_axis(coords):
    c = coords - coords.mean(0)
    _, _, vh = np.linalg.svd(c, full_matrices=False)
    return vh[0]


def rodrigues(v, axis, angle, pivot):
    k = axis / (np.linalg.norm(axis) + 1e-12)
    p = v - pivot
    return (p * np.cos(angle)
            + np.cross(k, p) * np.sin(angle)
            + k * (k @ p.T)[:, None] * (1 - np.cos(angle))) + pivot


def template_hinge(T, w=4):
    """Sharpest backbone turn in the template = the hinge residue (native-free)."""
    best = (1e9, len(T) // 2)
    for m in range(w, len(T) - w):
        a = T[m] - T[m - w]
        b = T[m + w] - T[m]
        cos = (a @ b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12)
        if cos < best[0]:
            best = (cos, m)
    return best[1]


def main():
    query, template = "1DUH", "2PXE"
    q = _pdb(query)
    qs, Q = parse_rna_c1(q, rna_chains(q)[0])

    t = _pdb(template)
    best = None
    for c in rna_chains(t):
        ts, TX = parse_rna_c1(t, c)
        if len(ts) < 15:
            continue
        pairs = nw_align(qs, ts)
        if len(pairs) < 15:
            continue
        ident = sum(1 for a, b in pairs if qs[a] == ts[b]) / len(pairs)
        if best is None or ident > best[0]:
            best = (ident, pairs, TX)
    ident, pairs, TX = best
    qi = np.array([a for a, b in pairs])
    ti = np.array([b for a, b in pairs])
    Qa = Q[qi]           # native, aligned order (scoring target only)
    T = TX[ti].copy()    # template coords transferred to aligned query order

    k = template_hinge(T)
    A, B = np.arange(0, k), np.arange(k, len(T))
    uA, uB = pca_axis(T[A]), pca_axis(T[B])
    axis = np.cross(uA, uB)
    pivot = T[k]

    fsot_rg = PI * len(Qa) ** (1.0 / PI)
    native_rg = rg(Qa)

    sweep = []
    for dth in np.linspace(-np.pi, np.pi, 361):
        M = T.copy()
        M[B] = rodrigues(T[B], axis, dth, pivot)
        sweep.append((dth, kabsch_rmsd(M, Qa), rg(M)))
    sweep = np.array(sweep)

    # FSOT stacking principle: coaxial stack = maximal base-stacking continuity.
    # Angle that makes domain-B helix axis colinear with domain-A axis (growth term).
    phi = np.arccos(np.clip(uA @ uB / (np.linalg.norm(uA) * np.linalg.norm(uB)), -1, 1))
    stack_cands = [-phi, np.pi - phi, phi - np.pi, phi]

    def eval_angle(dth):
        M = T.copy()
        M[B] = rodrigues(T[B], axis, dth, pivot)
        return kabsch_rmsd(M, Qa), rg(M)

    # native-free pick: the colinear direction that CONTINUES the backbone 5'->3'
    dA = T[k - 1] - T[max(k - 4, 0)]
    dA = dA / (np.linalg.norm(dA) + 1e-12)

    def continuation(dth):
        M = T.copy()
        M[B] = rodrigues(T[B], axis, dth, pivot)
        dB = M[min(k + 3, len(T) - 1)] - M[k]
        return float(dA @ (dB / (np.linalg.norm(dB) + 1e-12)))

    templ_rmsd = float(sweep[np.argmin(np.abs(sweep[:, 0])), 1])
    oi = int(np.argmin(sweep[:, 1]))
    fi = int(np.argmin(np.abs(sweep[:, 2] - fsot_rg)))
    stack_dth = max(stack_cands, key=continuation)
    stack_rmsd, stack_rg = eval_angle(stack_dth)

    print(f"template={template} id={ident:.2f} n_aligned={len(Qa)} hinge@res={k}")
    print(f"native Rg={native_rg:.2f}  FSOT globular Rg={fsot_rg:.2f}")
    print(f"[template angle dth=0]          RMSD {templ_rmsd:.2f} A")
    print(f"[FSOT-Rg globular angle]        dth={sweep[fi,0]*180/np.pi:+.0f} deg  "
          f"Rg={sweep[fi,2]:.2f}  RMSD {sweep[fi,1]:.2f} A")
    print(f"[FSOT coaxial-stack angle]      dth={stack_dth*180/np.pi:+.0f} deg  "
          f"Rg={stack_rg:.2f}  RMSD {stack_rmsd:.2f} A")
    print(f"[oracle best-angle ceiling]     dth={sweep[oi,0]*180/np.pi:+.0f} deg  "
          f"Rg={sweep[oi,2]:.2f}  RMSD {sweep[oi,1]:.2f} A")


if __name__ == "__main__":
    raise SystemExit(main())
