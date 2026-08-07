#!/usr/bin/env python3
"""Long-range contact ranking — Haskell Contact.hs parity (evidence tags only).

Zero free parameters: each channel is a seed/SMILES closed form.
Used for (1) D-matrix top-L clamps and (2) rich diagnostics vs native contacts.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Any

import numpy as np

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "vendor"))
import fsot_compute as fc  # noqa: E402

PHI = float(fc.PHI)
E = float(fc.E)
PI = float(fc.PI)

from smiles_aa_chem import (  # noqa: E402
    f18_disulfide_gate,
    formal_charge,
    hydrophobicity_kd,
    polarizability_aa,
)


@dataclass
class Region:
    kind: str
    start: int
    end: int


def helix_heptad_multiplier(i: int, j: int, start_i: int, start_j: int) -> float:
    mi = (i - start_i) % 7
    mj = (j - start_j) % 7
    if mi in (0, 3) and mj in (0, 3):
        return PHI
    return 1.0 / PHI


def beta_register_multiplier(
    i: int, j: int, sa: int, ea: int, sb: int, eb: int
) -> float:
    j_anti = eb - (i - sa)
    j_par = sb + (i - sa)
    off = min(abs(j - j_anti), abs(j - j_par))
    return PHI ** (-off / PI)


@dataclass
class EvidenceBreakdown:
    m_ij: float
    hydrophobic: float
    salt: float
    disulfide: float
    register: float
    polarizability: float
    total: float

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


def evidence_breakdown(
    aa1: str,
    aa2: str,
    m_ij: float,
    sep: int,
    *,
    region_same: bool = False,
    register_mult: float = 1.0,
) -> EvidenceBreakdown:
    """One score = sum of evidence channels (matches Haskell ContactScore)."""
    a, b = aa1.upper(), aa2.upper()
    h1, h2 = hydrophobicity_kd(a), hydrophobicity_kd(b)
    hydro = 0.0
    if h1 > 0 and h2 > 0:
        hydro = PHI * math.sqrt(h1 * h2)
    q1, q2 = formal_charge(a), formal_charge(b)
    salt = 0.0
    if q1 * q2 < 0:
        salt = E * abs(q1 * q2) * PHI
    disul = 0.0
    if a == "C" and b == "C":
        disul = (PHI ** 3) * f18_disulfide_gate(sep)
    reg = register_mult * PHI if region_same else 0.0
    pol = math.sqrt(max(polarizability_aa(a) * polarizability_aa(b), 0.0)) / (PHI * E)
    total = float(m_ij) + hydro + salt + disul + reg + pol
    return EvidenceBreakdown(
        m_ij=float(m_ij),
        hydrophobic=hydro,
        salt=salt,
        disulfide=disul,
        register=reg,
        polarizability=pol,
        total=total,
    )


def rank_long_range_contacts(
    sequence: str,
    M: np.ndarray,
    regions: list[Region],
    gate: int,
) -> list[dict[str, Any]]:
    """Rank all |i-j|≥gate pairs high-score first; full evidence for diagnostics."""
    n = M.shape[0]
    seq = sequence.upper()
    assert len(seq) == n
    # residue → region
    rmap: list[int | None] = [None] * n
    for ri, r in enumerate(regions):
        for i in range(r.start, r.end + 1):
            if i < n:
                rmap[i] = ri

    ranked: list[dict[str, Any]] = []
    for i in range(n):
        for j in range(i + gate, n):
            region_same = False
            reg_mult = 1.0
            ri, rj = rmap[i], rmap[j]
            if (
                ri is not None
                and rj is not None
                and ri != rj
                and regions[ri].kind == regions[rj].kind
                and regions[ri].kind != "C"
            ):
                region_same = True
                ra, rb = regions[ri], regions[rj]
                if ra.kind == "H":
                    reg_mult = helix_heptad_multiplier(i, j, ra.start, rb.start)
                else:
                    reg_mult = beta_register_multiplier(
                        i, j, ra.start, ra.end, rb.start, rb.end
                    )
            ev = evidence_breakdown(
                seq[i],
                seq[j],
                float(M[i, j]),
                j - i,
                region_same=region_same,
                register_mult=reg_mult,
            )
            ranked.append(
                {
                    "i": i,
                    "j": j,
                    "sep": j - i,
                    "aa_i": seq[i],
                    "aa_j": seq[j],
                    "score": ev.total,
                    "evidence": ev.as_dict(),
                    "region_same": region_same,
                }
            )
    ranked.sort(key=lambda r: -r["score"])
    return ranked


def top_l_precision_vs_native(
    ranked: list[dict[str, Any]],
    De: np.ndarray,
    contact_cutoff: float,
    L: int,
    gate: int = 7,
) -> dict[str, Any]:
    """How well evidence ranking recovers experimental contacts (pre-fold metric)."""
    natives_lr: set[tuple[int, int]] = set()
    n = De.shape[0]
    for i in range(n):
        for j in range(i + gate, n):
            if De[i, j] < contact_cutoff:
                natives_lr.add((i, j))
    top = ranked[: max(L, 1)]
    top_pairs = {(r["i"], r["j"]) for r in top}
    hits = top_pairs & natives_lr
    rank_of = {(r["i"], r["j"]): k for k, r in enumerate(ranked)}
    hit_ev = {
        "hydrophobic": 0.0,
        "salt": 0.0,
        "disulfide": 0.0,
        "register": 0.0,
        "m_ij": 0.0,
        "polarizability": 0.0,
    }
    miss_natives = []
    for r in ranked:
        key = (r["i"], r["j"])
        if key in natives_lr and key not in top_pairs:
            miss_natives.append(
                {
                    "i": r["i"],
                    "j": r["j"],
                    "aa": r["aa_i"] + r["aa_j"],
                    "score": r["score"],
                    "rank": rank_of.get(key, -1),
                    "evidence": r["evidence"],
                }
            )
    for r in top:
        if (r["i"], r["j"]) in hits:
            for k in hit_ev:
                hit_ev[k] += float(r["evidence"].get(k, 0.0))
    nh = max(len(hits), 1)
    for k in hit_ev:
        hit_ev[k] /= nh
    return {
        "L": L,
        "n_native_long_range": len(natives_lr),
        "n_top_hits": len(hits),
        "top_L_precision_evidence": (len(hits) / max(L, 1)),
        "top_L_recall_evidence": (len(hits) / max(len(natives_lr), 1)),
        "mean_evidence_on_hits": hit_ev,
        "missed_natives_sample": sorted(miss_natives, key=lambda x: x["rank"])[:15],
        "top_predicted_sample": [
            {
                "i": r["i"],
                "j": r["j"],
                "aa": r["aa_i"] + r["aa_j"],
                "score": r["score"],
                "native": (r["i"], r["j"]) in natives_lr,
                "evidence": r["evidence"],
            }
            for r in top[:12]
        ],
    }
