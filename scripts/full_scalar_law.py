#!/usr/bin/env python3
"""Full FSOT scalar law for the genetics fold path.

Uses the *entire* formula S = K(T1 + T2 + T3), including:
  - T1 observer-modulated base + observed quirk
  - T2 scale/amplitude
  - T3 valve / chaos / poof / suction / acoustic / phase

Zero free parameters: seeds + pin DomainConfig + trinary-derived δψ only.
Hot path is float64 twin of vendor compute_scalar (same algebra, fast).
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402

# ── Seeds as float (pin values) ───────────────────────────────────────────
PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
GAMMA = float(fc.GAMMA)
ALPHA = float(fc.ALPHA)
PSI_CON = float(fc.PSI_CON)
ETA_EFF = float(fc.ETA_EFF)
BETA = float(fc.BETA)
THETA_S = float(fc.THETA_S)
POOF = float(fc.POOF)
SUCTION = float(fc.SUCTION)
CHAOS = float(fc.CHAOS)
C_EFF = float(fc.C_EFF)
A_BLEED = float(fc.A_BLEED)
A_IN = float(fc.A_IN)
P_VAR = float(fc.P_VAR)
B_IN = float(fc.B_IN)
P_NEW = float(fc.P_NEW)
C_FACTOR = float(fc.C_FACTOR)
K = float(fc.K)


def compute_scalar_full(
    *,
    N: float = 1.0,
    P: float = 1.0,
    D_eff: float = 25.0,
    delta_psi: float = 1.0,
    delta_theta: float = 1.0,
    recent_hits: float = 0.0,
    rho: float = 1.0,
    scale: float = 1.0,
    amplitude: float = 1.0,
    trend_bias: float = 0.0,
    observed: bool = False,
) -> dict[str, float]:
    """Full S = K(T1+T2+T3) with term breakdown. Matches vendor algebra."""
    Nn = max(N, 1e-6)
    D = max(D_eff, 1e-6)
    dp = delta_psi
    dt = delta_theta
    hits = recent_hits

    # ── T1: Observer-Modulated Base ──
    growth = math.exp(ALPHA * (1.0 - hits / Nn) * GAMMA / PHI)
    base = (
        (N * P / math.sqrt(D))
        * math.cos((PSI_CON + dp) / ETA_EFF)
        * math.exp(-ALPHA * hits / Nn + rho + B_IN * dp)
        * (1.0 + growth * C_EFF)
    )
    T1 = base * (1.0 + P_NEW * math.log(D / 25.0))
    observer_mod = 1.0
    if observed:
        observer_mod = math.exp(C_FACTOR * P_VAR) * math.cos(dp + P_VAR)
        T1 = T1 * observer_mod

    # ── T2: Linear Modulation ──
    T2 = scale * amplitude + trend_bias

    # ── T3: Valve-Acoustic-Phase ──
    valve = (
        BETA
        * math.cos(dp)
        * (N * P / math.sqrt(D))
        * (1.0 + CHAOS * (D - 25.0) / 25.0)
        * (1.0 + POOF * math.cos(THETA_S + PI) + SUCTION * math.sin(THETA_S))
    )
    acoustic = 1.0 + (A_BLEED * math.sin(dt) ** 2) / PHI + (A_IN * math.cos(dt) ** 2) / PHI
    phase = 1.0 + B_IN * P_VAR
    T3 = valve * acoustic * phase

    S = K * (T1 + T2 + T3)
    return {
        "S": S,
        "T1": T1,
        "T2": T2,
        "T3": T3,
        "observer_mod": observer_mod,
        "chaos_factor": 1.0 + CHAOS * (D - 25.0) / 25.0,
        "D_eff": D,
        "delta_psi": dp,
        "delta_theta": dt,
        "recent_hits": hits,
        "observed": 1.0 if observed else 0.0,
    }


def residual_scale(S: float, factor: float | None = None) -> float:
    """FSOT residual law factor: (1 + |S| · f). Default f = P_NEW (pin)."""
    f = P_NEW if factor is None else factor
    return 1.0 + abs(S) * f


@dataclass(frozen=True)
class ScaleDomain:
    name: str
    D_eff: int
    delta_psi0: float
    delta_theta0: float
    observed_default: bool


def _domain_row(name: str) -> ScaleDomain:
    d = fc.DOMAINS[name]
    return ScaleDomain(
        name=name,
        D_eff=int(d.D_eff),
        delta_psi0=float(d.delta_psi),
        # pin table often has delta_theta=0 — match vendor domain_scalar exactly
        delta_theta0=float(d.delta_theta),
        observed_default=bool(d.observed),
    )


# Named pin domains for *chemical connection class* (not free continuous D)
SCALE_BACKBONE = _domain_row("Physical_Chemistry")     # virtual bond / local covalent geometry
SCALE_COVALENT = _domain_row("Atomic_Physics")         # disulfide / covalent chemistry
SCALE_HBOND = _domain_row("Chemistry")                 # H-bond secondary structure
SCALE_MOLECULE = _domain_row("Molecular_Chemistry")    # general residue chemistry
SCALE_ELECTRO = _domain_row("Electromagnetism")        # salt bridges / charge interaction
SCALE_TERTIARY = _domain_row("Biochemistry")           # tertiary fold observation
SCALE_PACK = _domain_row("Condensed_Matter")           # hydrophobic packing / dense core
SCALE_BIO = _domain_row("Biology")                     # optional cellular context

# Alias for older callers
SCALE_LOCAL = SCALE_MOLECULE
SCALE_SS = SCALE_HBOND


def scale_for_sep(sep: int, long_range_gate: int) -> ScaleDomain:
    """Legacy sep-only route (fallback). Prefer chem_link_domain()."""
    if sep <= 2:
        return SCALE_BACKBONE
    if sep in (3, 4, 7):
        return SCALE_HBOND
    if sep < long_range_gate:
        return SCALE_MOLECULE
    return SCALE_TERTIARY


def chem_link_domain(
    sep: int,
    long_range_gate: int,
    *,
    aa1: str = "",
    aa2: str = "",
    p_alpha_i: float = 0.0,
    p_alpha_j: float = 0.0,
    p_beta_i: float = 0.0,
    p_beta_j: float = 0.0,
) -> tuple[ScaleDomain, str]:
    """Route pair to D_eff by *connecting chemical system*, not sep alone.

    Doctrine: interacting chemical structures live at different FSOT interfaces.
    Named domains only (pin table). Returns (domain, link_class label).
    """
    a1, a2 = (aa1 or "?").upper(), (aa2 or "?").upper()

    # 1) Backbone geometry — Physical_Chemistry (not an observation of fold)
    if sep <= 2:
        return SCALE_BACKBONE, "backbone_covalent_geometry"

    # 2) Disulfide — Atomic_Physics covalent scale + observed
    if a1 == "C" and a2 == "C" and sep >= 3:
        return SCALE_COVALENT, "disulfide_covalent"

    # SMILES charges / hydrophobicity when available
    q1 = q2 = 0.0
    h1 = h2 = 0.0
    try:
        from smiles_aa_chem import formal_charge, hydrophobicity_kd  # noqa: WPS433

        if a1 in "ARNDCQEGHILKMFPSTWYV" and a2 in "ARNDCQEGHILKMFPSTWYV":
            q1, q2 = formal_charge(a1), formal_charge(a2)
            h1, h2 = hydrophobicity_kd(a1), hydrophobicity_kd(a2)
    except Exception:
        pass

    # 3) Salt bridge — Electromagnetism (charge–charge), observed
    if q1 * q2 < -1e-9 and sep >= 3:
        return SCALE_ELECTRO, "salt_bridge_electrostatic"

    # 4) Hydrophobic packing core — Condensed_Matter (dense soft matter)
    if h1 > 0.0 and h2 > 0.0 and sep >= long_range_gate:
        return SCALE_PACK, "hydrophobic_packing"

    # 5) H-bond secondary — Chemistry (H-bond is chemical, not packing)
    gate_a = 1.0 / E
    helix_pair = (
        sep in (3, 4, 7)
        and p_alpha_i > gate_a
        and p_alpha_j > gate_a
    )
    sheet_pair = sep >= 3 and p_beta_i > gate_a and p_beta_j > gate_a
    if helix_pair or sheet_pair:
        return SCALE_HBOND, "hbond_secondary"

    # 6) Mid-range molecular chemistry
    if sep < long_range_gate:
        return SCALE_MOLECULE, "molecular_sidechain"

    # 7) Long-range tertiary fold — Biochemistry (observed macromolecule)
    return SCALE_TERTIARY, "tertiary_biochem"


def delta_psi_from_trinary(
    dpsi0: float,
    spin_i: float,
    spin_j: float,
    charge_i: float,
    charge_j: float,
) -> float:
    """Local observer/phase offset from AA trinary (seed-only modulation)."""
    # spin mean pulls phase; charge difference pulls via B_IN (already seed)
    spin_m = 0.5 * (spin_i + spin_j)
    charge_m = 0.5 * (charge_i + charge_j)
    return dpsi0 * (1.0 + spin_m / PHI) + B_IN * charge_m / PHI


def delta_theta_from_trinary(dtheta0: float, branch_i: int, branch_j: int, aro_i: int, aro_j: int) -> float:
    """Acoustic angle from structural trits."""
    topo = 0.5 * (branch_i + branch_j + aro_i + aro_j) / 2.0
    return dtheta0 * (1.0 + topo / (PHI * PHI))


@lru_cache(maxsize=65536)
def _cached_S(
    D_eff: int,
    dpsi_q: int,
    dtheta_q: int,
    hits_q: int,
    observed: int,
    N_q: int,
) -> float:
    """Quantized cache for O(n²) pair path — same algebra, no free dials."""
    dpsi = dpsi_q / 1000.0
    dtheta = dtheta_q / 1000.0
    hits = float(hits_q)
    N = max(float(N_q), 1.0)
    out = compute_scalar_full(
        N=N,
        P=1.0,
        D_eff=float(D_eff),
        delta_psi=dpsi,
        delta_theta=dtheta,
        recent_hits=hits,
        observed=bool(observed),
    )
    return out["S"]


def pair_full_scalar(
    sep: int,
    spin_i: float,
    spin_j: float,
    charge_i: float,
    charge_j: float,
    *,
    branch_i: int = 0,
    branch_j: int = 0,
    aro_i: int = 0,
    aro_j: int = 0,
    long_range_gate: int = 7,
    chain_len: int = 1,
    recent_hits: float = 0.0,
    force_observed: bool | None = None,
    aa1: str = "",
    aa2: str = "",
    p_alpha_i: float = 0.0,
    p_alpha_j: float = 0.0,
    p_beta_i: float = 0.0,
    p_beta_j: float = 0.0,
) -> dict[str, float]:
    """Full-law S for one residue pair at the *chemically correct* D_eff interface."""
    sc, link = chem_link_domain(
        sep,
        long_range_gate,
        aa1=aa1,
        aa2=aa2,
        p_alpha_i=p_alpha_i,
        p_alpha_j=p_alpha_j,
        p_beta_i=p_beta_i,
        p_beta_j=p_beta_j,
    )
    dpsi = delta_psi_from_trinary(sc.delta_psi0, spin_i, spin_j, charge_i, charge_j)
    dtheta = delta_theta_from_trinary(sc.delta_theta0, branch_i, branch_j, aro_i, aro_j)

    # Observer: backbone unobserved; connecting chemical systems that are
    # measurements of structure (pack, salt, disulfide, tertiary) observed=True
    if force_observed is not None:
        observed = force_observed
    elif link == "backbone_covalent_geometry":
        observed = False
    elif link in (
        "disulfide_covalent",
        "salt_bridge_electrostatic",
        "hydrophobic_packing",
        "tertiary_biochem",
        "hbond_secondary",
    ):
        observed = True
    else:
        observed = bool(sc.observed_default)

    N = max(float(chain_len), 1.0)
    S = _cached_S(
        sc.D_eff,
        int(round(dpsi * 1000)),
        int(round(dtheta * 1000)),
        int(round(recent_hits)),
        1 if observed else 0,
        int(min(N, 400)),
    )
    chaos_factor = 1.0 + CHAOS * (float(sc.D_eff) - 25.0) / 25.0
    return {
        "S": S,
        "residual": residual_scale(S),
        "domain": sc.name,
        "D_eff": float(sc.D_eff),
        "delta_psi": dpsi,
        "delta_theta": dtheta,
        "observed": 1.0 if observed else 0.0,
        "chaos_factor": chaos_factor,
        "recent_hits": recent_hits,
        "chem_link": link,
        "delta_psi0_domain": sc.delta_psi0,
    }


def refine_observation_scalar(
    round_idx: int,
    n_rounds: int,
    chain_len: int,
    mean_stress_proxy: float = 0.0,
) -> dict[str, float]:
    """Each refine round is an *observation* of the structure (observer path on).

    hits grow with round; δψ from residual stress proxy (seed-scaled).
    Domain: Biochemistry (folding as biochemical observation).
    """
    sc = SCALE_TERTIARY
    # hits ∈ [0, n_rounds] — information accumulated by iterative observation
    hits = float(round_idx)
    # δψ base from domain; stress pulls phase (seed B_IN)
    dpsi = sc.delta_psi0 + B_IN * math.tanh(mean_stress_proxy / PHI) / PHI
    dtheta = sc.delta_theta0 * (1.0 + float(round_idx) / (float(n_rounds) + PHI))
    return compute_scalar_full(
        N=max(float(chain_len), 1.0),
        P=1.0,
        D_eff=float(sc.D_eff),
        delta_psi=dpsi,
        delta_theta=dtheta,
        recent_hits=hits,
        observed=True,  # THE observer path — structure is being measured
        scale=1.0,
        amplitude=1.0,
    )


def parity_check_domain(name: str = "Biochemistry", tol: float = 1e-5) -> dict[str, Any]:
    """Confirm float twin matches vendor domain_scalar for pin domains."""
    d = fc.DOMAINS[name]
    ref = float(fc.domain_scalar(name))
    ours = compute_scalar_full(
        N=1.0,
        P=1.0,
        D_eff=float(d.D_eff),
        delta_psi=float(d.delta_psi),
        delta_theta=float(d.delta_theta),
        recent_hits=float(d.hits),
        observed=bool(d.observed),
    )["S"]
    err = abs(ours - ref)
    return {
        "domain": name,
        "vendor_S": ref,
        "full_law_S": ours,
        "abs_err": err,
        "ok": err < tol or err / max(abs(ref), 1e-12) < 1e-4,
    }


def main() -> int:
    print("Full scalar law — parity vs vendor domain_scalar")
    for name in (
        "Physical_Chemistry",
        "Chemistry",
        "Molecular_Chemistry",
        "Biochemistry",
        "Condensed_Matter",
        "Biology",
    ):
        r = parity_check_domain(name)
        status = "OK" if r["ok"] else "FAIL"
        print(
            f"  {status} {name:22s} vendor={r['vendor_S']:+.8f}  "
            f"full={r['full_law_S']:+.8f}  err={r['abs_err']:.2e}"
        )
    # observer on vs off
    off = compute_scalar_full(D_eff=13, delta_psi=0.35, recent_hits=1, observed=False)
    on = compute_scalar_full(D_eff=13, delta_psi=0.35, recent_hits=1, observed=True)
    print(f"\n  observer OFF S={off['S']:+.8f}  ON S={on['S']:+.8f}  "
          f"mod={on['observer_mod']:.6f}")
    print(f"  chaos_factor@D13={on['chaos_factor']:.6f}  T1={on['T1']:.6f} T2={on['T2']:.6f} T3={on['T3']:.6e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
