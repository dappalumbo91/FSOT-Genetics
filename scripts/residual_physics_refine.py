#!/usr/bin/env python3
"""Physics refine with archive residual law — correct FSOT application.

Archive (fsot_api_predict_lib):
  computed = measured * (1 + |S| * factor)
  S = domain_scalar(named_domain), factor = P_NEW (or domain factor)

Structure:
  measured = template Cα (homolog PDB — real observable)
  each force channel is a named ChemLink / pin domain:
    bond   → Physical_Chemistry  (backbone geometry)
    clash  → Chemistry           (steric / local chemistry)
    anchor → Biochemistry        (fold observation; hold measured)
  residual_r = 1 + |S_domain| * P_NEW
  force_channel *= residual_r

No new free parameters. No invented contact graphs. No medoid tricks.
If product median rises above GATE, do not ship.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from fsot_structure_engine import CA_CA  # noqa: E402
from full_scalar_law import P_NEW, residual_scale  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)

# Seed clash floor (same spirit as e+φ ≈ 4.34; match proven 4.0 physics path closely)
CLASH = E + PHI  # ~4.34 A
# base lr / iters from proven physics refine (test_physics_refine)
ITERS = 120
LR = 0.08
W_ANCHOR0 = 0.05


def _S(name: str) -> float:
    return abs(float(fc.domain_scalar(name)))


def residual_physics_relax(
    X0: np.ndarray,
    *,
    iters: int = ITERS,
    lr: float = LR,
) -> np.ndarray:
    """Template-anchored bond/clash refine with per-domain residual weights."""
    # Named domains — pin table only
    r_bond = residual_scale(_S("Physical_Chemistry"))   # backbone
    r_clash = residual_scale(_S("Chemistry"))            # local steric/chem
    r_anchor = residual_scale(_S("Biochemistry"))        # measured fold observation

    w_anchor = W_ANCHOR0 * r_anchor
    X = X0.copy()
    n = len(X)
    for _ in range(iters):
        G = w_anchor * (X - X0)
        # bonds — Physical_Chemistry residual
        d = X[1:] - X[:-1]
        L = np.linalg.norm(d, axis=1) + 1e-9
        f = ((L - CA_CA) / L)[:, None] * d * r_bond
        G[:-1] -= f
        G[1:] += f
        # clashes — Chemistry residual
        diff = X[:, None, :] - X[None, :, :]
        D = np.linalg.norm(diff, axis=2)
        np.fill_diagonal(D, 1e9)
        idx = np.arange(n - 1)
        D[idx, idx + 1] = 1e9
        D[idx + 1, idx] = 1e9
        mask = D < CLASH
        if mask.any():
            coef = np.where(mask, (D - CLASH) / (D + 1e-9), 0.0) * r_clash
            G += np.einsum("ij,ijk->ik", coef, diff)
        X = X - lr * G
    return X - X.mean(axis=0)


def residuals_report() -> dict[str, float]:
    return {
        "Physical_Chemistry_S": _S("Physical_Chemistry"),
        "Physical_Chemistry_residual": residual_scale(_S("Physical_Chemistry")),
        "Chemistry_S": _S("Chemistry"),
        "Chemistry_residual": residual_scale(_S("Chemistry")),
        "Biochemistry_S": _S("Biochemistry"),
        "Biochemistry_residual": residual_scale(_S("Biochemistry")),
        "P_NEW": float(P_NEW),
    }
