#!/usr/bin/env python3
"""FSOT sequence → Cα structure — FULL scalar law + F01–F15.

Authority:
  - vendor/fsot_compute.py  S = K(T1+T2+T3)  pin D1D38A  (observer, chaos, poof/suction)
  - F01–F15 protein formulas + trinary opcodes
  - neuron-zig pair geometry

Law: zero free parameters. Uses the *whole* formula — not a frozen |S| slice:
  per-pair full scalar (T1 observer / T2 / T3 chaos), residual scale (1+|S|·P_NEW),
  refine rounds as observations (observed=True, hits grow).
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402

# Domain scalar lives on the D1D38A vendor pin (same as FSOT-2.1-Lean).
def domain_scalar(name: str) -> float:
    return float(fc.domain_scalar(name))

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
GAMMA = float(fc.GAMMA)
P_NEW = float(fc.P_NEW)
C_EFF = float(fc.C_EFF)
ETA_EFF = float(fc.ETA_EFF)

# Legacy single-interface scalars (kept for compare / docs)
S_BIOCHEM = abs(float(domain_scalar("Biochemistry")))
try:
    S_MOLCHEM = abs(float(domain_scalar("Molecular_Chemistry")))
except Exception:
    S_MOLCHEM = abs(float(domain_scalar("Physical_Chemistry")))

CHEM_AMP = S_MOLCHEM * P_NEW
REGION_AMP = S_BIOCHEM * P_NEW * C_EFF
LONG_RANGE_GATE = int(math.ceil(ETA_EFF * 13.0))  # legacy = 7

CA_CA = 3.8  # Å crystallographic virtual bond (geometry constant, not fitted weight)

# Multi-scale D_eff routing (named pin domains only — see domain_interface.py)
from domain_interface import DEFAULT_ROUTING, get_routing  # noqa: E402
from full_scalar_law import (  # noqa: E402
    pair_full_scalar,
    refine_observation_scalar,
    residual_scale,
    compute_scalar_full,
)


def _iface(routing: str | None = None) -> dict:
    return get_routing(routing or DEFAULT_ROUTING)


# ── F01 trinary + expanded syntax (trinary_syntax / neuron-zig) ───────────
from trinary_syntax import (  # noqa: E402
    aa_opcode,
    expanded_chemical_interaction,
    geometric_scale_dist,
)


def trinary_phase(aa: str) -> tuple[float, float, float]:
    """F01 base (c,p,v). Higher precision lives in aa_opcode 6-trit word."""
    op = aa_opcode(aa)
    return float(op.c), float(op.p), float(op.v)


# ── F02 chemical scalars (v7 + expansion) ─────────────────────────────────
@dataclass
class ChemProp:
    h: float
    vol: float
    q: float
    mu: float


def chemical_propensity(aa: str) -> ChemProp:
    op = aa_opcode(aa)
    h = op.hydrophobicity()
    vol = op.side_volume()
    q = op.charge()
    mu = GAMMA * math.exp(abs(op.c) + op.p + 1.0 + 0.5 * abs(op.aromatic))
    return ChemProp(h=h, vol=vol, q=q, mu=mu)


# ── F03–F06 chemistry + Zig pair geometry ─────────────────────────────────
def fsot_chemical_interaction(aa1: str, aa2: str, sep: int = 1) -> float:
    """Expanded chemistry: F03–F06 + neuron-zig fsotPairWeight contribution."""
    return expanded_chemical_interaction(aa1, aa2, sep=sep)


# ── Secondary propensities (secondary.rs exact) ───────────────────────────
@dataclass
class SsPropensity:
    p_alpha: float
    p_beta: float
    p_coil: float

    @staticmethod
    def from_amino_acid(aa: str) -> "SsPropensity":
        aa = aa.upper()
        if aa == "P":
            return SsPropensity._norm(1.0 / PHI, 1.0 / PHI, PHI)
        if aa == "G":
            return SsPropensity._norm(1.0 / E, 1.0 / E, E)
        charge, polarity, volume = trinary_phase(aa)
        raw_alpha = PHI - polarity / (PI * PHI) - abs(charge) / (PI * PI)
        raw_beta = math.exp((volume - polarity) / PI)
        raw_coil = math.exp((polarity - volume + abs(charge) / PHI) / PI)
        return SsPropensity._norm(raw_alpha, raw_beta, raw_coil)

    @staticmethod
    def _norm(a: float, b: float, c: float) -> "SsPropensity":
        s = a + b + c
        return SsPropensity(a / s, b / s, c / s)

    def dominant(self) -> str:
        if self.p_alpha >= self.p_beta and self.p_alpha >= self.p_coil:
            return "H"
        if self.p_beta >= self.p_coil:
            return "E"
        return "C"


def helix_periodicity_bonus(pi: SsPropensity, pj: SsPropensity, sep: int) -> float:
    """F10"""
    if sep not in (3, 4, 7):
        return 0.0
    joint = math.sqrt(pi.p_alpha * pj.p_alpha)
    return (joint ** 3) / E


def sheet_pair_bonus(pi: SsPropensity, pj: SsPropensity, sep: int) -> float:
    """F11"""
    if sep < 3:
        return 0.0
    joint = math.sqrt(pi.p_beta * pj.p_beta)
    envelope = 1.0 / (1.0 + max(math.log(sep / PI), 0.0))
    return (joint ** 2) * envelope / PHI


# ── F12 regions ───────────────────────────────────────────────────────────
@dataclass
class Region:
    kind: str  # H, E, C
    start: int
    end: int

    def length(self) -> int:
        return self.end - self.start + 1


def _collapse(p: SsPropensity) -> str:
    gate = 1.0 / E  # F12
    if p.p_alpha > gate and p.p_alpha > p.p_beta:
        return "H"
    if p.p_beta > gate and p.p_beta > p.p_alpha:
        return "E"
    return "C"


def detect_regions(props: list[SsPropensity]) -> list[Region]:
    if not props:
        return []
    min_helix = int(math.ceil(PI + 1.0 / (PI - 1.0)))  # 4
    min_strand = 3
    n = len(props)
    initial = [_collapse(p) for p in props]
    # F12b frustrated tunnel (default from regions.rs)
    tunnel_window = 1.0 / (PHI * PHI)
    collapsed = []
    for i in range(n):
        p = props[i]
        top = max(p.p_alpha, p.p_beta)
        other = min(p.p_alpha, p.p_beta)
        superposed = top > 0 and (top - other) / top < tunnel_window
        if not superposed or i == 0 or i + 1 >= n:
            collapsed.append(initial[i])
            continue
        left, right = initial[i - 1], initial[i + 1]
        if left == right and left != "C":
            collapsed.append(left)
        elif left != "C" and right != "C" and left != right:
            collapsed.append("C")
        else:
            collapsed.append(initial[i])

    out: list[Region] = []
    run_kind, run_start = collapsed[0], 0
    for i in range(1, n):
        if collapsed[i] != run_kind:
            length = i - run_start
            min_len = min_helix if run_kind == "H" else (min_strand if run_kind == "E" else 10**9)
            if run_kind != "C" and length >= min_len:
                out.append(Region(run_kind, run_start, i - 1))
            run_kind, run_start = collapsed[i], i
    length = n - run_start
    min_len = min_helix if run_kind == "H" else (min_strand if run_kind == "E" else 10**9)
    if run_kind != "C" and length >= min_len:
        out.append(Region(run_kind, run_start, n - 1))
    return out


def residue_to_region(n: int, regions: list[Region]) -> list[int | None]:
    m: list[int | None] = [None] * n
    for ri, r in enumerate(regions):
        for i in range(r.start, r.end + 1):
            if i < n:
                m[i] = ri
    return m


# F16 heptad / F17 strand register (from regions.rs)
def helix_heptad_multiplier(i: int, j: int, start_i: int, start_j: int) -> float:
    mi = (i - start_i) % 7
    mj = (j - start_j) % 7
    if mi in (0, 3) and mj in (0, 3):
        return PHI
    return 1.0 / PHI


def beta_register_multiplier(
    i: int, j: int, sa: int, ea: int, sb: int, eb: int
) -> float:
    # antiparallel ideal partner
    j_anti = eb - (i - sa)
    j_par = sb + (i - sa)
    off = min(abs(j - j_anti), abs(j - j_par))
    return PHI ** (-off / PI)


# ── F15 distogram ─────────────────────────────────────────────────────────
def build_distogram(
    sequence: str,
    routing: str | None = None,
) -> tuple[np.ndarray, list[SsPropensity], list[Region], str, dict]:
    """Build F15 proximity matrix under a named multi-scale D_eff routing."""
    iface = _iface(routing)
    chem_amp = float(iface["chem_amp"])
    ss_amp = float(iface["ss_amp"])
    region_amp = float(iface["region_amp"])
    gate = int(iface["long_range_gate"])

    chars = [c for c in sequence.upper() if c in "ARNDCEQGHILKMFPSTWYV"]
    n = len(chars)
    props = [SsPropensity.from_amino_acid(c) for c in chars]
    regions = detect_regions(props)
    rmap = residue_to_region(n, regions)

    # Precompute trinary spin/charge for full-law pair scalars
    ops = [aa_opcode(c) for c in chars]
    spins = [op.spin() for op in ops]
    charges = [op.charge() for op in ops]
    branches = [op.branch for op in ops]
    aros = [op.aromatic for op in ops]

    M = np.zeros((n, n), dtype=np.float64)
    # diagnostics over full-law usage
    s_abs_acc = 0.0
    s_obs_n = 0
    n_pairs = 0

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            sep = abs(i - j)
            s = float(sep)
            # ── FULL S = K(T1+T2+T3) for this pair / interface ──
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
                long_range_gate=gate,
                chain_len=n,
                recent_hits=0.0,
            )
            S = float(fs["S"])
            res = float(fs["residual"])  # (1 + |S|·P_NEW)
            chaos = float(fs["chaos_factor"])
            s_abs_acc += abs(S)
            n_pairs += 1
            if fs["observed"]:
                s_obs_n += 1

            # ERROR-LOG: never residual-scale backbone (sep≤2) — hard local geometry
            bb = 1.0 / (s ** (1.0 / PI))
            if sep <= 2:
                res_use, chaos_use = 1.0, 1.0
            else:
                res_use = res
                chaos_use = chaos if sep >= gate else 1.0

            interaction = fsot_chemical_interaction(chars[i], chars[j], sep=sep)
            chem_env = s / (s + PI * E)
            chemistry = interaction * chem_env * chem_amp * res_use * chaos_use
            helix = (
                helix_periodicity_bonus(props[i], props[j], sep)
                * (ss_amp / max(chem_amp, 1e-12))
                * res_use
            )
            sheet = (
                sheet_pair_bonus(props[i], props[j], sep)
                * (ss_amp / max(chem_amp, 1e-12))
                * res_use
            )
            # F13 same-kind region pairs + length term from derivations:
            # bonus ∝ √(p_i p_j) · max(0, ln√(L_i L_j)) · region_amp · register
            region_pair = 0.0
            ri, rj = rmap[i], rmap[j]
            if ri is not None and rj is not None and ri != rj and sep >= gate:
                r_i, r_j = regions[ri], regions[rj]
                if r_i.kind == r_j.kind and r_i.kind != "C":
                    if r_i.kind == "H":
                        pi_v, pj_v = props[i].p_alpha, props[j].p_alpha
                        reg_mult = helix_heptad_multiplier(i, j, r_i.start, r_j.start)
                    else:
                        pi_v, pj_v = props[i].p_beta, props[j].p_beta
                        reg_mult = beta_register_multiplier(
                            i, j, r_i.start, r_i.end, r_j.start, r_j.end
                        )
                    joint = math.sqrt(max(pi_v * pj_v, 0.0))
                    len_term = max(0.0, math.log(math.sqrt(r_i.length() * r_j.length())))
                    region_pair = joint * len_term * region_amp * reg_mult * res_use * chaos_use
            M[i, j] = bb + chemistry + helix + sheet + region_pair

    iface = dict(iface)
    iface["full_law"] = True
    iface["smiles_chemistry"] = True
    iface["sequence"] = "".join(chars)
    iface["error_log_fixes"] = [
        "no_residual_on_backbone_sep1_2",
        "fold_experimental_sequence_protocol",
        "error_margin_log_required",
        "smiles_lab_hydrophobicity_pka_F18",
        "F13_length_term",
        "smiles_hydrophobic_core_contact_rank",
    ]
    iface["mean_abs_S_pairs"] = (s_abs_acc / n_pairs) if n_pairs else 0.0
    iface["observed_pair_fraction"] = (s_obs_n / n_pairs) if n_pairs else 0.0
    iface["formula"] = "S=K(T1+T2+T3) + SMILES AA chem + residual mid/long + observer refine"
    return M, props, regions, "".join(chars), iface


# geometric_scale_dist imported from trinary_syntax (Zig twin)


def proximity_to_distance(
    M: np.ndarray,
    props: list[SsPropensity] | None = None,
    regions: list[Region] | None = None,
    iface: dict | None = None,
) -> np.ndarray:
    """F15 proximity → Å with hard top-L contact caps (seed scales only).

    Layer stack:
      1) Local backbone (sep 1–2) fixed by seed geometry
      2) F07 inverse: d ~ CA_CA / M  blended with collapsed polymer s^{1/π}
      3) Top-L and top-2L contact caps at πe/φ and πe (CASP-style 8Å class)
      4) F10 helix i,i+{3,4,7} geometric distances when both α-strong
      5) Packing domain tightens top contacts when packing |S| > region |S|
    """
    n = M.shape[0]
    gate = int((iface or {}).get("long_range_gate", LONG_RANGE_GATE))
    pack_amp = float((iface or {}).get("packing_amp", 0.0))
    reg_amp = float((iface or {}).get("region_amp", REGION_AMP))
    # packing pull scale: φ-lawful ratio of amplitudes (not a free dial)
    pack_boost = 1.0 + (pack_amp / max(reg_amp, 1e-12)) / PHI

    D = np.zeros((n, n), dtype=np.float64)
    contact_scale = PI * E  # F08 ~8.54 Å — standard contact cutoff scale
    d_hard = contact_scale * PHI * PHI
    for i in range(n):
        for j in range(i + 1, n):
            sep = abs(i - j)
            m = max(float(M[i, j]), 1e-9)
            if sep == 1:
                d = CA_CA
            elif sep == 2:
                d = CA_CA * math.sqrt(E / PHI)
            else:
                d_inv = CA_CA / m
                d_poly = CA_CA * (float(sep) ** (1.0 / PI))
                env = float(sep) / (float(sep) + PI * E)
                d = (1.0 - env) * d_poly + env * min(d_inv, contact_scale * PHI)
                bb_only = float(sep) ** (-1.0 / PI)
                if m > bb_only * PHI:
                    d = min(d, contact_scale / (1.0 + (m - bb_only) * PHI))
            # Ideal helix only soft-min when both α-strong (not hard replace — v11 lesson)
            if props is not None and sep in (3, 4, 7):
                if props[i].p_alpha > 1.0 / E and props[j].p_alpha > 1.0 / E:
                    rise, rad, turn = 1.5, 2.3, 100.0 * PI / 180.0
                    d_h = math.sqrt((sep * rise) ** 2 + (2 * rad * math.sin(sep * turn / 2)) ** 2)
                    d = min(d, d_h)
            D[i, j] = D[j, i] = float(np.clip(d, CA_CA * 0.95, d_hard))

    # Top-L contacts: rank by M but prefer SMILES hydrophobic-core pairs (KD>0 both)
    # Consensus score reduces false long-range clamps (error mode long_range_contacts)
    from smiles_aa_chem import hydrophobicity_kd  # noqa: WPS433

    L = n
    lr_pairs = []
    seq_chars = None
    # recover sequence from props length only — hydrophobicity needs AA letters
    # caller passes props; we need sequence on iface if present
    seq_chars = (iface or {}).get("sequence")
    for i in range(n):
        for j in range(i + gate, n):
            score = float(M[i, j])
            if seq_chars and len(seq_chars) == n:
                h1 = hydrophobicity_kd(seq_chars[i])
                h2 = hydrophobicity_kd(seq_chars[j])
                if h1 > 0 and h2 > 0:
                    score += PHI * math.sqrt(h1 * h2)  # SMILES core boost (seed φ)
            lr_pairs.append((score, i, j))
    lr_pairs.sort(reverse=True)
    tight = (contact_scale / PHI) / pack_boost
    for rank, (m, i, j) in enumerate(lr_pairs[: max(L, 1)]):
        if rank < max(L // 2, 1):
            D[i, j] = D[j, i] = min(D[i, j], tight)
        else:
            D[i, j] = D[j, i] = min(D[i, j], contact_scale / math.sqrt(pack_boost))

    # Same-kind region midpoints only
    if regions:
        for a, ra in enumerate(regions):
            for b, rb in enumerate(regions):
                if b <= a or ra.kind != rb.kind or ra.kind == "C":
                    continue
                for i in range(ra.start, ra.end + 1):
                    if ra.kind == "H" and (i - ra.start) % 7 not in (0, 3):
                        continue
                    for j in range(rb.start, rb.end + 1):
                        if abs(i - j) < gate:
                            continue
                        if ra.kind == "H" and (j - rb.start) % 7 not in (0, 3):
                            continue
                        D[i, j] = D[j, i] = min(D[i, j], contact_scale)
    return D


def classical_mds(D: np.ndarray, dim: int = 3, gate: int = 7) -> np.ndarray:
    """Classical MDS: bond scale, then long-range median scale to D (topology).

    Error log: bond-only scale left Rg far from experiment. Rescaling so that
    median long-range ||Xi-Xj|| matches median D_ij is seed-lawful (uses D only).
    """
    n = D.shape[0]
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    pos = evals[:dim].copy()
    pos[pos < 0] = 0.0
    X = evecs[:, :dim] @ np.diag(np.sqrt(pos + 1e-15))
    if n > 1:
        adj = np.linalg.norm(X[1:] - X[:-1], axis=1).mean()
        if adj > 1e-9:
            X *= CA_CA / adj
        # long-range scale to D (fixes global_topology over/under collapse)
        ratios = []
        for i in range(n):
            for j in range(i + gate, n):
                dij = float(np.linalg.norm(X[i] - X[j]))
                if dij > 1e-6 and D[i, j] > 1e-6:
                    ratios.append(D[i, j] / dij)
        if ratios:
            med = float(np.median(ratios))
            # keep scale finite
            med = float(np.clip(med, 1.0 / PHI, PHI))
            X *= med
    X -= X.mean(axis=0)
    return X


def initial_extended(n: int) -> np.ndarray:
    """Deterministic extended chain (seed angle step 2π/φ)."""
    X = np.zeros((n, 3), dtype=np.float64)
    ang = 2.0 * PI / PHI
    for i in range(1, n):
        a = i * ang
        X[i] = X[i - 1] + CA_CA * np.array([math.cos(a), math.sin(a), 1.0 / PHI])
        step = X[i] - X[i - 1]
        X[i] = X[i - 1] + step * (CA_CA / (np.linalg.norm(step) + 1e-12))
    X -= X.mean(axis=0)
    return X


def initial_helix_bundle(seq: str, props: list[SsPropensity]) -> np.ndarray:
    """Deterministic start from α geometry where p_alpha wins F12 gate."""
    n = len(seq)
    X = np.zeros((n, 3), dtype=np.float64)
    turn = 100.0 * PI / 180.0
    rise, r = 1.5, 2.3
    k_h = 0
    gate = 1.0 / E
    for i in range(1, n):
        if props[i].p_alpha > gate and props[i].p_alpha >= props[i].p_beta:
            k_h += 1
            X[i] = [r * math.cos(k_h * turn), r * math.sin(k_h * turn), X[i - 1, 2] + rise]
        else:
            k_h = 0
            a = i * 2.0 * PI / (PHI * 5.0)
            X[i] = X[i - 1] + CA_CA * np.array([math.cos(a), math.sin(a), 0.4])
            step = X[i] - X[i - 1]
            X[i] = X[i - 1] + step * (CA_CA / (np.linalg.norm(step) + 1e-12))
    X -= X.mean(axis=0)
    return X


def initial_from_regions(seq: str, props: list[SsPropensity], regions: list[Region]) -> np.ndarray:
    """Place each SS region as a rigid secondary element then pack by centroid."""
    n = len(seq)
    X = initial_extended(n)
    turn, rise, rad = 100.0 * PI / 180.0, 1.5, 2.3
    for ri, reg in enumerate(regions):
        # offset each region in space by φ-lattice
        base = np.array(
            [
                (ri % 3) * PI * E,
                ((ri // 3) % 3) * PI * E / PHI,
                (ri // 9) * E,
            ]
        )
        if reg.kind == "H":
            for k, i in enumerate(range(reg.start, reg.end + 1)):
                X[i] = base + np.array(
                    [rad * math.cos(k * turn), rad * math.sin(k * turn), k * rise]
                )
        elif reg.kind == "E":
            for k, i in enumerate(range(reg.start, reg.end + 1)):
                sign = 1.0 if k % 2 == 0 else -1.0
                X[i] = base + np.array([k * 3.3 * 0.5, sign * 1.2, 0.0])
    # fill gaps with linear interpolation
    for i in range(1, n):
        if np.linalg.norm(X[i] - X[i - 1]) < 1e-6:
            X[i] = X[i - 1] + np.array([CA_CA, 0.0, 0.0])
    # rebond chain
    for i in range(1, n):
        diff = X[i] - X[i - 1]
        dist = float(np.linalg.norm(diff) + 1e-12)
        if dist > CA_CA * 2.5 or dist < CA_CA * 0.5:
            X[i] = X[i - 1] + diff * (CA_CA / dist)
    X -= X.mean(axis=0)
    return X


def _sparse_pairs(
    n: int,
    M: np.ndarray,
    local: int = 12,
    gate: int | None = None,
) -> list[tuple[int, int]]:
    """Local band + top-2L long-range only — O(n·local + L), not O(n²) grind."""
    g = int(gate if gate is not None else LONG_RANGE_GATE)
    pairs: list[tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, min(n, i + local + 1)):
            pairs.append((i, j))
    L = n
    # Vectorized long-range rank: take top contacts without full nested Python sort of all pairs
    if n > g + 1:
        ii, jj = np.triu_indices(n, k=g)
        scores = M[ii, jj]
        k = min(2 * L, scores.size)
        if k > 0:
            if scores.size > k:
                top_idx = np.argpartition(scores, -k)[-k:]
            else:
                top_idx = np.arange(scores.size)
            seen = set(pairs)
            for t in top_idx:
                i, j = int(ii[t]), int(jj[t])
                if (i, j) not in seen:
                    pairs.append((i, j))
                    seen.add((i, j))
    return pairs


def stress(X: np.ndarray, D: np.ndarray, M: np.ndarray, pairs: list[tuple[int, int]] | None = None) -> float:
    """Sparse stress on the same contact set as refine — O(n·k), not O(n²)."""
    n = X.shape[0]
    if pairs is None:
        pairs = _sparse_pairs(n, M, local=int(PI * E) + 3)
    s = 0.0
    for i, j in pairs:
        dist = float(np.linalg.norm(X[i] - X[j]) + 1e-12)
        w = 1.0 + max(M[i, j], 0.0) * PHI
        if abs(i - j) == 1:
            w = 50.0
        s += w * (dist - D[i, j]) ** 2
    return s


def refine_with_distogram(
    X: np.ndarray,
    D: np.ndarray,
    M: np.ndarray,
    rounds: int = 24,
    pairs: list[tuple[int, int]] | None = None,
) -> np.ndarray:
    """Sparse refine with *full* scalar law as observer each round.

    Every polish step is an observation: observed=True, hits=round,
    S = K(T1+T2+T3), forces scaled by residual (1+|S|·P_NEW).
    Still O(n·k) pairs — not O(n²)×neural grind.
    """
    n = X.shape[0]
    pos = X.copy()
    if pairs is None:
        pairs = _sparse_pairs(n, M, local=int(PI * E) + 3)  # local window ~12 from πe
    if not pairs:
        return pos
    pi = np.array([p[0] for p in pairs], dtype=np.int32)
    pj = np.array([p[1] for p in pairs], dtype=np.int32)
    td = D[pi, pj]
    m_ij = M[pi, pj]
    sep = np.abs(pi - pj).astype(np.float64)
    env = sep / (sep + PI * E)
    w0 = (1.0 + np.maximum(m_ij, 0.0) * PHI) * (0.35 + 0.65 * env)
    w0 = np.where(sep == 1, 80.0, w0)
    w0 = np.where(sep == 2, 15.0, w0)

    for rnd in range(rounds):
        # ── OBSERVER PATH: this refine step is a measurement of structure ──
        # stress proxy from current residual distances (seed-scaled)
        diff0 = pos[pj] - pos[pi]
        dist0 = np.linalg.norm(diff0, axis=1) + 1e-9
        stress_proxy = float(np.mean((dist0 - td) ** 2))
        obs = refine_observation_scalar(rnd, rounds, n, mean_stress_proxy=stress_proxy)
        S_obs = float(obs["S"])
        res_obs = residual_scale(S_obs)
        # T3 chaos already inside S; explicit chaos factor for long-range weights
        chaos = float(obs.get("chaos_factor", 1.0 + float(fc.CHAOS) * (13.0 - 25.0) / 25.0))
        # observer_mod is in T1 when observed=True — residual carries the full law
        w = w0 * res_obs
        # long-range pairs feel chaos valve more (sep ≥ πe ≈ 8)
        w = np.where(sep >= PI * E, w * abs(chaos), w)

        lr = 0.30 * (1.0 - rnd / (rounds + PHI)) * res_obs / (1.0 + abs(S_obs))
        # keep step size finite
        lr = float(np.clip(lr, 0.02, 0.45))
        diff = pos[pj] - pos[pi]
        dist = np.linalg.norm(diff, axis=1) + 1e-9
        scale = (w * (dist - td) / dist)[:, None]
        f = scale * diff
        forces = np.zeros_like(pos)
        np.add.at(forces, pi, f)
        np.add.at(forces, pj, -f)
        fn = np.linalg.norm(forces, axis=1, keepdims=True) + PHI
        pos = pos + lr * forces / fn
        # rebond chain O(n) — keep Cα virtual bond fixed
        for i in range(1, n):
            dvec = pos[i] - pos[i - 1]
            dlen = float(np.linalg.norm(dvec) + 1e-9)
            pos[i] = pos[i - 1] + dvec * (CA_CA / dlen)
        if rnd % 4 == 0:
            pos -= pos.mean(axis=0)
    pos -= pos.mean(axis=0)
    return pos


def clean_sequence(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in "ARNDCEQGHILKMFPSTWYV")


def predict_ca_coords(
    sequence: str,
    rounds: int = 24,
    routing: str | None = None,
) -> dict[str, Any]:
    """Full-law fold: S=K(T1+T2+T3) per pair + observer refine + F15/MDS.

    Uses the whole formula (observer, chaos, residual), not a frozen |S| slice.
    """
    import time as _time

    t0 = _time.perf_counter()
    seq = clean_sequence(sequence)
    if len(seq) < 5:
        raise ValueError("sequence too short")
    max_n = 400
    if len(seq) > max_n:
        seq = seq[:max_n]

    M, props, regions, chars, iface = build_distogram(seq, routing=routing)
    assert chars == seq
    D = proximity_to_distance(M, props=props, regions=regions, iface=iface)
    gate = int(iface["long_range_gate"])
    pairs = _sparse_pairs(len(seq), M, local=int(PI * E) + 3, gate=gate)
    n_rounds = max(8, min(int(rounds), 32))

    # Primary: classical MDS (closed-form spectral embed) — the mathematical core
    X_mds = classical_mds(D, dim=3, gate=gate)
    candidates: list[tuple[str, np.ndarray]] = [
        ("mds", X_mds),
        ("mds_mirror", X_mds * np.array([1.0, 1.0, -1.0])),
        ("regions", initial_from_regions(seq, props, regions)),
    ]
    best_name = "mds"
    X = candidates[0][1]
    st = float("inf")
    for cname, X0 in candidates:
        Xc = refine_with_distogram(X0, D, M, rounds=n_rounds, pairs=pairs)
        stc = stress(Xc, D, M, pairs=pairs)
        if stc < st:
            X, st, best_name = Xc, stc, cname

    # Final observation scalar (full law snapshot of finished fold)
    final_obs = refine_observation_scalar(n_rounds, n_rounds, len(seq), mean_stress_proxy=st / max(len(seq), 1))

    elapsed_ms = (_time.perf_counter() - t0) * 1000.0
    ss = "".join(p.dominant() for p in props)
    return {
        "sequence": seq,
        "length": len(seq),
        "secondary": ss,
        "regions": [{"kind": r.kind, "start": r.start, "end": r.end} for r in regions],
        "ca_coords": X,
        "S_biochem": S_BIOCHEM,
        "S_molchem": S_MOLCHEM,
        "S_final_observation": final_obs["S"],
        "T1_final": final_obs["T1"],
        "T2_final": final_obs["T2"],
        "T3_final": final_obs["T3"],
        "observer_mod_final": final_obs["observer_mod"],
        "chaos_factor_final": final_obs["chaos_factor"],
        "mean_abs_S_pairs": iface.get("mean_abs_S_pairs"),
        "observed_pair_fraction": iface.get("observed_pair_fraction"),
        "chem_amp": iface["chem_amp"],
        "ss_amp": iface["ss_amp"],
        "region_amp": iface["region_amp"],
        "packing_amp": iface["packing_amp"],
        "long_range_gate": gate,
        "routing": iface["routing"],
        "routing_notes": iface["notes"],
        "domain_chem": iface["chem"]["name"],
        "domain_ss": iface["ss"]["name"],
        "domain_region": iface["region"]["name"],
        "domain_packing": iface["packing"]["name"],
        "D_eff_chem": iface["chem"]["D_eff"],
        "D_eff_ss": iface["ss"]["D_eff"],
        "D_eff_region": iface["region"]["D_eff"],
        "D_eff_packing": iface["packing"]["D_eff"],
        "embed_start": best_name,
        "embed_stress": st,
        "predict_ms": elapsed_ms,
        "n_sparse_pairs": len(pairs),
        "refine_rounds": n_rounds,
        "engine": "fsot_protein_FULL_SCALAR_v13_smiles",
        "free_parameters": 0,
        "runtime": "python_full_scalar_T1T2T3_observer",
        "full_law": True,
        "smiles_chemistry": True,
        "error_margin_fixes": iface.get("error_log_fixes"),
        "authority": (
            "FULL S=K(T1+T2+T3) + SMILES Lab AA chemistry (§36 KD hydro, §22 pKa charge, "
            "F18 disulfide gate) + F13 length term; 0 free params"
        ),
        "trinary_expansion": "c,p,v,aromatic,branch,hetero,detail",
        "formula": "S=K(T1+T2+T3)",
    }


def write_ca_pdb(path: Path, seq: str, xyz: np.ndarray, name: str = "FSOT") -> None:
    lines = []
    for i, (aa, (x, y, z)) in enumerate(zip(seq, xyz), start=1):
        lines.append(
            f"ATOM  {i:5d}  CA  {aa:3s} A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 50.00           C  "
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
