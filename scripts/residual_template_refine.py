#!/usr/bin/env python3
"""Archive residual law applied to structure (template = measured authority).

Archive pattern (fsot_api_predict_lib):
  computed = measured * (1 + |S| * P_NEW)

Structure analog:
  measured = homolog Cα coordinates (real PDB observable)
  S        = full scalar at each ChemLink D_eff for that pair
  residual = 1 + |S| * P_NEW
  refine   = spring to chem_link geometric target, force ∝ residual
  anchor   = template memory weight 1/φ (topology preserved)

This is the same multi-system math as the green domain panels: pin → named
domain → bridge measured data → residual → measure. Zero free parameters.
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
from fsot_structure_engine import CA_CA, clean_sequence  # noqa: E402
from full_scalar_law import (  # noqa: E402
    P_NEW,
    pair_full_scalar,
    chem_link_target_distance,
    residual_scale,
)
from trinary_syntax import aa_opcode  # noqa: E402
from fsot_structure_engine import SsPropensity  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
GATE = 7
ANCHOR = 1.0 / PHI
LR = 1.0 / (PI * E * PHI)  # ~0.072 — stable step
ITERS = int(round((E ** PI) * PHI))  # ~38
CLASH = E + PHI
# residual weight clip: keep 1+|S|P_NEW in seed-bounded range
W_MAX = PHI * E  # ~4.4


def residual_refine(
    sequence: str,
    X0: np.ndarray,
    *,
    iters: int = ITERS,
    lr: float = LR,
    anchor_w: float = ANCHOR,
    evo_cons: np.ndarray | None = None,
    evo_coev: np.ndarray | None = None,
) -> dict[str, Any]:
    """Template-anchored residual refine under full chem-link S."""
    seq = clean_sequence(sequence)
    n = len(seq)
    assert X0.shape == (n, 3)
    X = X0.copy()
    props = [SsPropensity.from_amino_acid(c) for c in seq]
    ops = [aa_opcode(c) for c in seq]
    spins = [op.spin() for op in ops]
    charges = [op.charge() for op in ops]
    branches = [op.branch for op in ops]
    aros = [op.aromatic for op in ops]
    if evo_cons is None:
        evo_cons = np.zeros(n)
    if evo_coev is None:
        evo_coev = np.zeros((n, n))

    # Precompute per-pair residual + target once (S from sequence + evo, not coords)
    targets: list[tuple[int, int, float, float, str]] = []
    s_abs = 0.0
    n_obs = 0
    for i in range(n):
        for j in range(i + 1, n):
            sep = j - i
            fs = pair_full_scalar(
                sep,
                spins[i],
                spins[j],
                charges[i],
                charges[j],
                branch_i=branches[i],
                branch_j=branches[j],
                aro_i=aros[i],
                aro_j=aros[j],
                long_range_gate=GATE,
                chain_len=n,
                aa1=seq[i],
                aa2=seq[j],
                p_alpha_i=props[i].p_alpha,
                p_alpha_j=props[j].p_alpha,
                p_beta_i=props[i].p_beta,
                p_beta_j=props[j].p_beta,
                evo_cons_i=float(evo_cons[i]),
                evo_cons_j=float(evo_cons[j]),
                evo_coev=float(evo_coev[i, j]),
            )
            link = str(fs["chem_link"])
            res = float(fs["residual"])
            s_abs += abs(float(fs["S"]))
            if fs["observed"]:
                n_obs += 1
            # SS / salt / hbond only — residual springs at those ChemLink D_eff
            if link not in (
                "disulfide_covalent",
                "salt_bridge_electrostatic",
                "hbond_secondary",
            ):
                continue
            d_tgt = chem_link_target_distance(
                link,
                sep,
                res,
                p_alpha_i=props[i].p_alpha,
                p_alpha_j=props[j].p_alpha,
            )
            if d_tgt is None:
                continue
            w = min(res, W_MAX)
            targets.append((i, j, float(d_tgt), w, link))

    for _ in range(iters):
        G = anchor_w * (X - X0)
        # bonds
        dvec = X[1:] - X[:-1]
        L = np.linalg.norm(dvec, axis=1) + 1e-9
        f = ((L - CA_CA) / L)[:, None] * dvec
        G[:-1] -= f
        G[1:] += f
        # clash (local only — O(n) window, not dense O(n²) every iter)
        for i in range(n):
            for j in range(i + 2, min(i + int(PI * E) + 2, n)):
                vec = X[j] - X[i]
                dist = float(np.linalg.norm(vec) + 1e-9)
                if dist < CLASH:
                    pull = (dist - CLASH) / dist
                    force = pull * vec
                    G[i] -= force
                    G[j] += force
        # residual chem-link springs
        for i, j, d_tgt, w, _link in targets:
            vec = X[j] - X[i]
            dist = float(np.linalg.norm(vec) + 1e-9)
            pull = (w / (PHI * E)) * (dist - d_tgt) / dist
            force = pull * vec
            # clip force magnitude (seed Å/step scale)
            fn = float(np.linalg.norm(force) + 1e-12)
            cap = PHI * E
            if fn > cap:
                force = force * (cap / fn)
            G[i] -= force
            G[j] += force
        # clip total gradient
        gn = np.linalg.norm(G, axis=1, keepdims=True) + 1e-12
        gcap = PHI * E
        G = np.where(gn > gcap, G * (gcap / gn), G)
        X = X - lr * G
        if not np.all(np.isfinite(X)):
            X = X0.copy()
            break

    X = X - X.mean(axis=0)
    if not np.all(np.isfinite(X)):
        X = X0 - X0.mean(axis=0)
    return {
        "ca_coords": X,
        "n_residual_springs": len(targets),
        "mean_abs_S_pairs": s_abs / max(n * (n - 1), 1),
        "observed_pair_fraction": n_obs / max(n * (n - 1), 1),
        "engine": "fsot_residual_template_refine_v1",
        "free_parameters": 0,
        "formula": "X_tmpl + residual(S_chem_link)·springs; S=K(T1+T2+T3)",
        "anchor_w": anchor_w,
        "iters": iters,
    }


def residual_predict(
    sequence: str,
    template_model: np.ndarray,
    features: Any | None = None,
) -> dict[str, Any]:
    evo_cons = None
    evo_coev = None
    if features is not None:
        try:
            from msa_pipeline import conservation_confidence  # noqa: WPS433

            evo_cons = conservation_confidence(features)
            evo_coev = getattr(features, "coevolution", None)
        except Exception:
            pass
    out = residual_refine(
        sequence,
        template_model,
        evo_cons=evo_cons,
        evo_coev=evo_coev,
    )
    out["structure_mode"] = "template_residual"
    return out
