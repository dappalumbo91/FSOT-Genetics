#!/usr/bin/env python3
"""FSOT fusion: template geometry + MSA coevolution polish (zero free params).

Why this exists
---------------
De-novo + MSA contacts alone cannot break the ~11 A information ceiling
(contacts underdetermine full distance structure). The *template* regime
already places the correct topology (~1.2 A). Evolutionary coevolution is
then useful as a *local packing polish* and confidence channel — not as a
from-scratch fold engine.

Map (all seed-closed):
  - CA–CA bond idealization at CA_CA (geometry constant)
  - Non-bonded clash floor at φ·e ≈ 4.4 A (seed product, not free)
  - Soft coevolution clamps: high MI−APC pairs pulled toward F08 contact
    scale π·e ≈ 8.54 A, weighted by √(c_i c_j) and evo_amp (F09 family)
  - Anchor to template with weight 1/φ so topology cannot drift

No trained weights. MSA is data input; amplitudes are domain/seed scalars.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from fsot_structure_engine import CA_CA  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402
from msa_pipeline import EVO_AMP, MsaFeatures, conservation_confidence  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)

CONTACT_SCALE = PI * E  # F08 ~8.54 A
CLASH_FLOOR = E + PHI  # ~4.34 A CA non-bonded floor (seed sum)
GATE = max(7, int(math.ceil(PI * math.sqrt(E))))
ANCHOR_W = 1.0 / PHI  # strong template memory
LR = 1.0 / (PI * E)  # ~0.117
ITERS = int(round((E ** PI) * PHI * PHI))  # ~61 soft-clamp steps

# Archive residual law: force_channel *= (1 + |S_domain| · P_NEW)
# Named domains from pin table only — zero free parameters.
_R_BOND = residual_scale(abs(float(fc.domain_scalar("Physical_Chemistry"))))  # backbone
_R_CLASH = residual_scale(abs(float(fc.domain_scalar("Chemistry"))))  # local steric
_R_ANCHOR = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))  # measured fold


def top_coevolution_pairs(
    features: MsaFeatures, *, gate: int = GATE, top_n: int | None = None
) -> list[tuple[int, int, float]]:
    """Top long-range coevolution pairs (i,j,score), score already ≥0."""
    C = features.coevolution
    n = C.shape[0]
    if top_n is None:
        top_n = n  # classic top-L
    c = features.conservation
    pairs: list[tuple[int, int, float]] = []
    for i in range(n):
        for j in range(i + gate, n):
            raw = float(C[i, j])
            if raw <= 0:
                continue
            w = math.sqrt(max(c[i], 0.0) * max(c[j], 0.0))
            if w <= 0:
                continue
            pairs.append((i, j, raw * w))
    pairs.sort(key=lambda t: t[2], reverse=True)
    return pairs[: max(top_n, 1)]


def template_coevolution_agreement(
    model: np.ndarray, features: MsaFeatures, *, contact_A: float | None = None
) -> float:
    """Fraction of top-L coevolution pairs that are spatial contacts in the model.

    Used as a template-ranking bonus (data agreement), not a free parameter.
    """
    if contact_A is None:
        contact_A = CONTACT_SCALE
    pairs = top_coevolution_pairs(features, top_n=len(features.sequence))
    if not pairs:
        return 0.0
    hits = 0
    for i, j, _ in pairs:
        d = float(np.linalg.norm(model[i] - model[j]))
        if d < contact_A * PHI:  # soft contact envelope
            hits += 1
    return hits / len(pairs)


def fuse_relax(
    X0: np.ndarray,
    features: MsaFeatures | None = None,
    *,
    iters: int = ITERS,
    lr: float = LR,
    anchor_w: float = ANCHOR_W,
) -> np.ndarray:
    """Template-anchored physics with residual-weighted force channels.

    measured = template Cα (homolog authority)
    residual_r = 1 + |S_domain| · P_NEW  (archive law)
      bond   ← Physical_Chemistry
      clash  ← Chemistry
      anchor ← Biochemistry  (holds measured fold; residual strengthens fidelity)
    MSA coevolution is data polish only — not a free residual invent.
    """
    X = X0.copy()
    n = len(X)
    # residual law on named domains (precomputed pin scalars)
    w_anchor = anchor_w * _R_ANCHOR
    r_bond = _R_BOND
    r_clash = _R_CLASH
    evo_pairs: list[tuple[int, int, float]] = []
    if features is not None and features.depth_ok and features.coevolution.max() > 0:
        raw_pairs = top_coevolution_pairs(features, top_n=n)
        # Normalize weights by top mean so evo_amp sets scale
        if raw_pairs:
            mean_top = float(np.mean([p[2] for p in raw_pairs])) + 1e-12
            evo_pairs = [(i, j, EVO_AMP * (s / mean_top)) for i, j, s in raw_pairs]

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
        mask = D < CLASH_FLOOR
        if mask.any():
            coef = np.where(mask, (D - CLASH_FLOOR) / (D + 1e-9), 0.0) * r_clash
            G += np.einsum("ij,ijk->ik", coef, diff)
        # coevolution packing polish ONLY — never rewire topology.
        # Only act on pairs already near contact (clash < d < φ·F08).
        # Distant coevolving pairs are reported in confidence, not forced.
        d_hi = CONTACT_SCALE * PHI  # ~13.8 A soft envelope
        for i, j, w in evo_pairs:
            vec = X[j] - X[i]
            dist = float(np.linalg.norm(vec) + 1e-9)
            if dist <= CLASH_FLOOR or dist >= d_hi:
                continue
            # gentle spring toward F08 contact scale
            pull = (w / PHI) * (dist - CONTACT_SCALE) / dist
            force = pull * vec
            G[i] -= force
            G[j] += force
        X = X - lr * G
    return X - X.mean(axis=0)


def fused_confidence(
    features: MsaFeatures | None,
    provenance: list[str] | None = None,
) -> np.ndarray:
    """Per-residue medical confidence in [0,1].

    Fuse evolutionary confidence with template provenance when available:
      identical → 1, mutated → 1/φ, gap → 1/φ²
    Final: 1 - (1 - evo)·(1 - prov)  (noisy-OR of support channels)
    """
    n = 0
    if features is not None:
        n = len(features.sequence)
        evo = conservation_confidence(features)
    else:
        evo = None
    if provenance is not None:
        n = max(n, len(provenance))
        prov = np.zeros(n)
        rank = {"identical": 1.0, "mutated": 1.0 / PHI, "gap": 1.0 / (PHI * PHI)}
        for i, tag in enumerate(provenance):
            prov[i] = rank.get(tag, 1.0 / (PHI * PHI))
    else:
        prov = None
    if evo is None and prov is None:
        return np.zeros(max(n, 1))
    if evo is None:
        return prov
    if prov is None:
        return evo
    m = min(len(evo), len(prov))
    return 1.0 - (1.0 - evo[:m]) * (1.0 - prov[:m])


def select_regime(
    has_template: bool,
    features: MsaFeatures | None,
) -> str:
    """Deployable regime picker (medical default policy).

    template_available → template_physics (+ optional packing fuse)
    else if deep MSA   → bulk_msa
    else               → bulk_single
    """
    if has_template:
        if features is not None and features.depth_ok:
            return "template_msa_fuse"
        return "template_physics"
    if features is not None and features.depth_ok:
        return "bulk_msa"
    return "bulk_single"


def fuse_predict(
    sequence: str,
    template_model: np.ndarray,
    features: MsaFeatures | None,
) -> dict[str, Any]:
    """Product path near AF: measured template + physics (+ optional MSA packing).

    Energy = bond + clash + template fidelity ‖X−X0‖² so we never drift off the
    measured homolog (archive spirit: measured authority stays primary).
    """
    X0 = template_model
    X_phys = fuse_relax(X0, None)
    X_fuse = fuse_relax(X0, features)

    def _energy(X: np.ndarray) -> float:
        n = len(X)
        bonds = np.linalg.norm(X[1:] - X[:-1], axis=1)
        e = float(((bonds - CA_CA) ** 2).sum())
        # fidelity to measured template (primary authority)
        e += float(ANCHOR_W * ((X - X0) ** 2).sum())
        for i in range(0, n, max(1, n // 40)):
            for j in range(i + 2, n, max(1, n // 40)):
                d = float(np.linalg.norm(X[i] - X[j]))
                if d < CLASH_FLOOR:
                    e += (CLASH_FLOOR - d) ** 2
        return e

    cands = [
        ("template_raw", X0 - X0.mean(0), _energy(X0 - X0.mean(0))),
        ("template_physics", X_phys, _energy(X_phys)),
        ("template_msa_fuse", X_fuse, _energy(X_fuse)),
    ]
    chosen, X, e_best = min(cands, key=lambda c: c[2])
    conf = fused_confidence(features, provenance=None)
    return {
        "ca_coords": X,
        "confidence": conf,
        "regime": chosen,
        "energy_best": e_best,
        "energies": {c[0]: c[2] for c in cands},
        "n_evo_clamps": int(
            len(top_coevolution_pairs(features, top_n=len(sequence)))
            if features is not None and features.depth_ok
            else 0
        ),
        "free_parameters": 0,
        "engine": "fsot_template_product_v5_residual",
        "contact_scale_A": CONTACT_SCALE,
        "clash_floor_A": CLASH_FLOOR,
        "iters": ITERS,
        "residual": {
            "Physical_Chemistry": _R_BOND,
            "Chemistry": _R_CLASH,
            "Biochemistry": _R_ANCHOR,
        },
        "formula": (
            "measured_homolog_Cα + residual-weighted physics; "
            "force*=(1+|S|·P_NEW); S=K(T1+T2+T3) pin domains"
        ),
    }
