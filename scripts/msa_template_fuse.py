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
# Lean ChemLink: tertiaryBiochem D=13, disulfide Atomic_Physics D=7
_R_TERT = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))
_R_SS = residual_scale(abs(float(fc.domain_scalar("Atomic_Physics"))))
# F13 gate at Biochemistry D_eff=13
_TERT_GATE = max(7, int(math.ceil(float(fc.ETA_EFF) * 13.0)))


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


def _contacts_from_measured(
    X0: np.ndarray,
    sequence: str | None,
    extra: list[tuple[int, int, float]] | None,
) -> tuple[list[tuple[int, int, float, float]], int]:
    """Build residual-tagged springs from *measured* homolog geometry.

    Returns (i, j, d_measured, residual_r) and count of SS links.
    Backbone sep≤2 never included (Lean: no residual on backbone).
    """
    n = len(X0)
    springs: list[tuple[int, int, float, float]] = []
    n_ss = 0
    # Consensus / extra measured tertiary (already filtered)
    raw_t: list[tuple[int, int, float]] = []
    if extra:
        for i, j, d0 in extra:
            if abs(j - i) <= 2:
                continue
            raw_t.append((i, j, float(d0)))
    else:
        for i in range(n):
            for j in range(i + _TERT_GATE, n):
                d0 = float(np.linalg.norm(X0[i] - X0[j]))
                if d0 < CONTACT_SCALE:
                    raw_t.append((i, j, d0))
    # Top-L measured contacts only (CASP-style data cap — not a free graph)
    raw_t.sort(key=lambda t: t[2])
    for i, j, d0 in raw_t[: max(n, 1)]:
        springs.append((i, j, d0, _R_TERT))
    # Disulfide: Cys–Cys close in measured map — Atomic_Physics residual
    if sequence and len(sequence) == n:
        cys = [i for i, a in enumerate(sequence.upper()) if a == "C"]
        for a in range(len(cys)):
            for b in range(a + 1, len(cys)):
                i, j = cys[a], cys[b]
                if abs(j - i) < 3:
                    continue
                d0 = float(np.linalg.norm(X0[i] - X0[j]))
                if d0 < CONTACT_SCALE:  # measured SS-like Cα span
                    springs.append((i, j, d0, _R_SS))
                    n_ss += 1
    return springs, n_ss


def fuse_relax(
    X0: np.ndarray,
    features: MsaFeatures | None = None,
    *,
    iters: int = ITERS,
    lr: float = LR,
    anchor_w: float = ANCHOR_W,
    sequence: str | None = None,
    tertiary_contacts: list[tuple[int, int, float]] | None = None,
) -> np.ndarray:
    """Template-anchored physics with residual-weighted ChemLink channels.

    measured = homolog Cα (real observable)
    residual_r = 1 + |S_domain| · P_NEW  (Lean ChemLink D_eff)
      bond     ← Physical_Chemistry  (sep=1; geometry)
      clash    ← Chemistry
      tertiary ← Biochemistry D=13 on *measured* long-range contacts only
      disulfide← Atomic_Physics on measured Cys–Cys close pairs
    Backbone is not residual-invented. No false contact graphs.
    """
    X = X0.copy()
    n = len(X)
    w_anchor = anchor_w * _R_ANCHOR
    r_bond = _R_BOND
    r_clash = _R_CLASH
    springs, _n_ss = _contacts_from_measured(X0, sequence, tertiary_contacts)
    evo_pairs: list[tuple[int, int, float]] = []
    if features is not None and features.depth_ok and features.coevolution.max() > 0:
        raw_pairs = top_coevolution_pairs(features, top_n=n)
        if raw_pairs:
            mean_top = float(np.mean([p[2] for p in raw_pairs])) + 1e-12
            evo_pairs = [(i, j, EVO_AMP * (s / mean_top)) for i, j, s in raw_pairs]

    for _ in range(iters):
        G = w_anchor * (X - X0)
        # bonds — Physical_Chemistry residual (geometry channel)
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
        # Measured tertiary / SS — residual at ChemLink interface, target = measured d
        for i, j, d0, r_link in springs:
            vec = X[j] - X[i]
            dist = float(np.linalg.norm(vec) + 1e-9)
            pull = (r_link / PHI) * (dist - d0) / dist
            force = pull * vec
            G[i] -= force
            G[j] += force
        # MSA polish: data only, already-near-contact envelope
        d_hi = CONTACT_SCALE * PHI
        for i, j, w in evo_pairs:
            vec = X[j] - X[i]
            dist = float(np.linalg.norm(vec) + 1e-9)
            if dist <= CLASH_FLOOR or dist >= d_hi:
                continue
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


def _termini_bond_stress(X: np.ndarray) -> tuple[float, float]:
    """Mean |L−CA_CA| in termini vs core — residual interface diagnosis."""
    n = len(X)
    if n < 8:
        return 0.0, 0.0
    bonds = np.linalg.norm(X[1:] - X[:-1], axis=1)
    err = np.abs(bonds - CA_CA)
    band = max(3, n // 10)
    term_idx = list(range(0, min(band, len(err)))) + list(
        range(max(0, len(err) - band), len(err))
    )
    core_idx = list(range(band, max(band + 1, len(err) - band)))
    if not core_idx:
        core_idx = list(range(len(err)))
    t = float(np.mean(err[term_idx])) if term_idx else 0.0
    c = float(np.mean(err[core_idx])) if core_idx else 1e-9
    return t, c


def fuse_predict(
    sequence: str,
    template_model: np.ndarray,
    features: MsaFeatures | None,
    *,
    tertiary_contacts: list[tuple[int, int, float]] | None = None,
) -> dict[str, Any]:
    """Product path: measured template + residual-weighted ChemLink physics.

    Energy uses residual_r on bond/clash/anchor (named domains). Soft termini
    only when termini bond stress exceeds core × r_bond (Physical_Chemistry
    residual-at-interface). Tertiary springs use measured homolog distances
    at Biochemistry D=13 (Lean ChemLink.tertiaryBiochem).
    """
    from run_rcsb_template_holdout import soft_flexible_termini  # noqa: WPS433

    X0 = template_model
    X_phys = fuse_relax(
        X0, None, sequence=sequence, tertiary_contacts=tertiary_contacts
    )
    X_fuse = fuse_relax(
        X0, features, sequence=sequence, tertiary_contacts=tertiary_contacts
    )

    def _energy(X: np.ndarray) -> float:
        # Residual-weighted energy — same interfaces as fuse_relax
        n = len(X)
        bonds = np.linalg.norm(X[1:] - X[:-1], axis=1)
        e = float(_R_BOND * ((bonds - CA_CA) ** 2).sum())
        e += float(ANCHOR_W * _R_ANCHOR * ((X - X0) ** 2).sum())
        for i in range(0, n, max(1, n // 40)):
            for j in range(i + 2, n, max(1, n // 40)):
                d = float(np.linalg.norm(X[i] - X[j]))
                if d < CLASH_FLOOR:
                    e += float(_R_CLASH * (CLASH_FLOOR - d) ** 2)
        return e

    cands = [
        ("template_raw", X0 - X0.mean(0), _energy(X0 - X0.mean(0))),
        ("template_physics", X_phys, _energy(X_phys)),
        ("template_msa_fuse", X_fuse, _energy(X_fuse)),
    ]
    # Soft termini candidates when residual bond stress localizes to ends
    soft_applied = False
    for label, Xb, _eb in list(cands):
        t_stress, c_stress = _termini_bond_stress(Xb)
        if t_stress > c_stress * _R_BOND and c_stress > 0:
            n = len(Xb)
            n_term = max(1, n // 20)
            c_term = max(2, n // 15)
            Xs = soft_flexible_termini(Xb, n_term=n_term, c_term=c_term)
            # re-physics after termini rebuild (bond residual again)
            Xs = fuse_relax(Xs, None if "msa" not in label else features)
            cands.append((label + "_soft_termini", Xs, _energy(Xs)))
            soft_applied = True

    chosen, X, e_best = min(cands, key=lambda c: c[2])
    conf = fused_confidence(features, provenance=None)
    return {
        "ca_coords": X,
        "confidence": conf,
        "regime": chosen,
        "energy_best": e_best,
        "energies": {c[0]: c[2] for c in cands},
        "soft_termini_considered": soft_applied,
        "n_evo_clamps": int(
            len(top_coevolution_pairs(features, top_n=len(sequence)))
            if features is not None and features.depth_ok
            else 0
        ),
        "free_parameters": 0,
        "engine": "fsot_template_product_v6_multisystem",
        "contact_scale_A": CONTACT_SCALE,
        "clash_floor_A": CLASH_FLOOR,
        "iters": ITERS,
        "residual": {
            "Physical_Chemistry": _R_BOND,
            "Chemistry": _R_CLASH,
            "Biochemistry": _R_ANCHOR,
        },
        "formula": (
            "measured_homolog_Cα + residual-weighted multi-system physics; "
            "termini soft if bond residual localizes to ends; "
            "S=K(T1+T2+T3) pin domains"
        ),
    }
