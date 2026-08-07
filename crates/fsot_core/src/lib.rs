//! FSOT 2.0 — Complete Computational Engine (Rust port).
//!
//! Mirrors `c:\Users\damia\Desktop\FSOT NeuroLab\fsot_compute.py` and
//! `c:\Users\damia\Desktop\FSOT SMILES Lab\fsot_compute.py` exactly,
//! using `f64` (≈15 decimal digits) instead of mpmath's 50-digit precision.
//!
//! Structure mirrors the Python original section-for-section:
//!   §1  foundational seeds (π, e, φ, γ, G_Cat)
//!   §2  Layer 1 derived constants (ALPHA, PSI_CON, ETA_EFF, …, POOF)
//!   §3  Layer 2 composite constants (C_EFF, …, K, C_COSM)
//!   §4  Scalar engine S = K · (T1 + T2 + T3)
//!   §5  35-domain parameter table (Particle_Physics … Cosmology)
//!   §6  Bond library distilled from the SMILES Lab Tier-0 dataset
//!
//! Everything is `no_std`-safe under the `nostd` feature (uses `libm`); the
//! default `std` feature uses `f64` methods directly so this crate links
//! cleanly into both the TUI and the bare-metal `kernel` crate.
//!
//! Author: Damian Arthur Palumbo — FSOT.

#![cfg_attr(not(feature = "std"), no_std)]
#![allow(non_snake_case, non_upper_case_globals)]
#![allow(clippy::excessive_precision)]

// --- math shims (std vs libm) ---------------------------------------------
mod m {
    #[cfg(feature = "std")]
    pub fn sin(x: f64) -> f64  { x.sin() }
    #[cfg(feature = "std")]
    pub fn cos(x: f64) -> f64  { x.cos() }
    #[cfg(feature = "std")]
    pub fn exp(x: f64) -> f64  { x.exp() }
    #[cfg(feature = "std")]
    pub fn ln(x: f64)  -> f64  { x.ln() }
    #[cfg(feature = "std")]
    pub fn sqrt(x: f64) -> f64 { x.sqrt() }
    #[cfg(feature = "std")]
    pub fn powf(b: f64, e: f64) -> f64 { b.powf(e) }

    #[cfg(not(feature = "std"))]
    pub fn sin(x: f64) -> f64  { libm::sin(x) }
    #[cfg(not(feature = "std"))]
    pub fn cos(x: f64) -> f64  { libm::cos(x) }
    #[cfg(not(feature = "std"))]
    pub fn exp(x: f64) -> f64  { libm::exp(x) }
    #[cfg(not(feature = "std"))]
    pub fn ln(x: f64)  -> f64  { libm::log(x) }
    #[cfg(not(feature = "std"))]
    pub fn sqrt(x: f64) -> f64 { libm::sqrt(x) }
    #[cfg(not(feature = "std"))]
    pub fn powf(b: f64, e: f64) -> f64 { libm::pow(b, e) }
}

// =========================================================================
// §1  FOUNDATIONAL SEEDS
// =========================================================================
pub const PI:    f64 = 3.141_592_653_589_793_2;
pub const E:     f64 = 2.718_281_828_459_045_2;
pub const PHI:   f64 = 1.618_033_988_749_894_8;        // (1+√5)/2
pub const GAMMA: f64 = 0.577_215_664_901_532_9;        // Euler–Mascheroni
pub const G_CAT: f64 = 0.915_965_594_177_219_0;        // Catalan's G

// =========================================================================
// §2  LAYER 1 — PRIMARY DERIVED CONSTANTS  (computed once at startup)
// =========================================================================
pub struct FsotConsts {
    // Layer 1
    pub ALPHA:   f64,
    pub PSI_CON: f64,
    pub ETA_EFF: f64,
    pub BETA:    f64,
    pub GAMMA_C: f64,
    pub OMEGA:   f64,
    pub THETA_S: f64,
    pub POOF:    f64,
    // Layer 2
    pub C_EFF:    f64,
    pub A_BLEED:  f64,
    pub P_VAR:    f64,
    pub B_IN:     f64,
    pub A_IN:     f64,
    pub SUCTION:  f64,
    pub CHAOS:    f64,
    pub P_BASE:   f64,
    pub P_NEW:    f64,
    pub C_FACTOR: f64,
    pub K:        f64,
    pub C_COSM:   f64,
}

impl FsotConsts {
    pub fn build() -> Self {
        // ---- §2 Layer 1 ----
        let ALPHA   = m::ln(PI) / (E * m::powf(PHI, 13.0));
        let PSI_CON = 1.0 - m::exp(-1.0);
        let ETA_EFF = 1.0 / (PI - 1.0);
        let BETA    = 1.0 / m::exp(m::powf(PI, PI) + (E - 1.0));
        let GAMMA_C = -m::ln(2.0) / PHI;
        let OMEGA   = m::sin(PI / E) * m::sqrt(2.0);
        let THETA_S = m::sin(PSI_CON * ETA_EFF);
        let POOF    = m::exp((-m::ln(PI) / E) / (ETA_EFF * m::ln(PHI)));

        // ---- §3 Layer 2 ----
        let C_EFF   = (1.0 - POOF * m::sin(THETA_S)) * (1.0 + 0.01 * G_CAT / (PI * PHI));
        let A_BLEED = m::sin(PI / E) * PHI / m::sqrt(2.0);
        let P_VAR   = -m::cos(THETA_S + PI);
        let B_IN    = C_EFF * (1.0 - m::sin(THETA_S) / PHI);
        let A_IN    = A_BLEED * (1.0 + m::cos(THETA_S) / PHI);
        let SUCTION = POOF * (-m::cos(THETA_S - PI));
        let CHAOS   = GAMMA_C / OMEGA;
        let P_BASE  = GAMMA / E;
        let P_NEW   = P_BASE * m::sqrt(2.0);
        let C_FACTOR = C_EFF * P_NEW;
        let K       = PHI * (GAMMA / E) * m::sqrt(2.0) / m::ln(PI) * 0.99;
        let C_COSM  = 1.0 / (PHI * 10.0);

        Self {
            ALPHA, PSI_CON, ETA_EFF, BETA, GAMMA_C, OMEGA, THETA_S, POOF,
            C_EFF, A_BLEED, P_VAR, B_IN, A_IN, SUCTION, CHAOS, P_BASE, P_NEW,
            C_FACTOR, K, C_COSM,
        }
    }
}

// =========================================================================
// §4  SCALAR ENGINE   S = K · (T1 + T2 + T3)
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct ScalarInput {
    pub N:           f64,
    pub P:           f64,
    pub D_eff:       f64,
    pub psi_con:     f64,
    pub delta_psi:   f64,
    pub recent_hits: f64,
    pub rho:         f64,
    pub B_in:        f64,
    pub C_eff:       f64,
    pub P_new:       f64,
    pub observed:    bool,
    pub beta:        f64,
    pub chaos:       f64,
    pub poof:        f64,
    pub suction:     f64,
    pub theta_s:     f64,
    pub delta_theta: f64,
    pub A_bleed:     f64,
    pub A_in:        f64,
    pub P_var:       f64,
    pub scale:       f64,
    pub amplitude:   f64,
    pub trend_bias:  f64,
    pub alpha:       f64,
}

impl ScalarInput {
    /// Build a default ScalarInput equivalent to the Python `ScalarInput()` no-arg form,
    /// using the constants from `c`.
    pub fn defaults_with(c: &FsotConsts) -> Self {
        Self {
            N: 1.0, P: 1.0, D_eff: 25.0,
            psi_con: c.PSI_CON, delta_psi: 1.0, recent_hits: 0.0, rho: 1.0,
            B_in: c.B_IN, C_eff: c.C_EFF, P_new: c.P_NEW, observed: false,
            beta: c.BETA, chaos: c.CHAOS, poof: c.POOF, suction: c.SUCTION,
            theta_s: c.THETA_S, delta_theta: 1.0,
            A_bleed: c.A_BLEED, A_in: c.A_IN, P_var: c.P_VAR,
            scale: 1.0, amplitude: 1.0, trend_bias: 0.0,
            alpha: c.ALPHA,
        }
    }
}

/// Compute the FSOT scalar  S = K · (T1 + T2 + T3).
///
/// Bit-for-bit translation of `compute_scalar` in `fsot_compute.py`.
pub fn compute_scalar(c: &FsotConsts, s: &ScalarInput) -> f64 {
    let N    = s.N;
    let P    = s.P;
    let D    = s.D_eff;
    let dp   = s.delta_psi;
    let dt   = s.delta_theta;
    let hits = s.recent_hits;

    // ── Term 1: Observer-Modulated Base ──
    let growth = m::exp(s.alpha * (1.0 - hits / N) * GAMMA / PHI);
    let base = (N * P / m::sqrt(D))
        * m::cos((s.psi_con + dp) / c.ETA_EFF)
        * m::exp(-s.alpha * hits / N + s.rho + s.B_in * dp)
        * (1.0 + growth * s.C_eff);
    let mut t1 = base * (1.0 + s.P_new * m::ln(D / 25.0));
    if s.observed {
        t1 = t1 * m::exp(c.C_FACTOR * s.P_var) * m::cos(dp + s.P_var);
    }

    // ── Term 2: Linear Modulation ──
    let t2 = s.scale * s.amplitude + s.trend_bias;

    // ── Term 3: Valve-Acoustic-Phase ──
    let valve = s.beta * m::cos(dp)
        * (N * P / m::sqrt(D))
        * (1.0 + s.chaos * (D - 25.0) / 25.0)
        * (1.0 + s.poof * m::cos(s.theta_s + PI) + s.suction * m::sin(s.theta_s));
    let acoustic = 1.0
        + (s.A_bleed * m::powf(m::sin(dt), 2.0)) / PHI
        + (s.A_in    * m::powf(m::cos(dt), 2.0)) / PHI;
    let phase = 1.0 + s.B_in * s.P_var;
    let t3 = valve * acoustic * phase;

    c.K * (t1 + t2 + t3)
}

// =========================================================================
// §5  35-DOMAIN PARAMETER TABLE
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct DomainConfig {
    pub name:        &'static str,
    pub D_eff:       u32,
    pub hits:        u32,
    pub delta_psi:   f64,
    pub delta_theta: f64,
    pub observed:    bool,
}

/// The full 35-domain table from `_build_domains()` in fsot_compute.py.
/// Domain interpretation constants (the `C` field) are intentionally omitted —
/// callers can re-derive them from the FsotConsts when needed.
pub const DOMAINS: &[DomainConfig] = &[
    DomainConfig { name: "Particle_Physics",      D_eff: 5,  hits: 0, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Quantum_Mechanics",     D_eff: 6,  hits: 0, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Atomic_Physics",        D_eff: 7,  hits: 0, delta_psi: 0.5,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Physical_Chemistry",    D_eff: 8,  hits: 0, delta_psi: 0.5,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Chemistry",             D_eff: 8,  hits: 0, delta_psi: 0.5,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Electromagnetism",      D_eff: 9,  hits: 0, delta_psi: 0.7,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Molecular_Chemistry",   D_eff: 9,  hits: 0, delta_psi: 0.4,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Optics",                D_eff: 10, hits: 0, delta_psi: 0.6,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Acoustics",             D_eff: 10, hits: 0, delta_psi: 0.3,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Quantum_Computing",     D_eff: 11, hits: 0, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Quantum_Optics",        D_eff: 11, hits: 0, delta_psi: 0.6,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Biology",               D_eff: 12, hits: 0, delta_psi: 0.05, delta_theta: 1.0, observed: false },
    DomainConfig { name: "Thermodynamics",        D_eff: 13, hits: 0, delta_psi: 0.4,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Biochemistry",          D_eff: 13, hits: 0, delta_psi: 0.1,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Neuroscience",          D_eff: 14, hits: 1, delta_psi: 0.1,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Condensed_Matter",      D_eff: 14, hits: 0, delta_psi: 0.5,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Fluid_Dynamics",        D_eff: 15, hits: 1, delta_psi: 0.9,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Nuclear_Physics",       D_eff: 15, hits: 1, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Ecology",               D_eff: 15, hits: 1, delta_psi: 0.2,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Meteorology",           D_eff: 16, hits: 2, delta_psi: 0.8,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Materials_Science",     D_eff: 16, hits: 0, delta_psi: 0.5,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Psychology",            D_eff: 16, hits: 1, delta_psi: 0.3,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Atmospheric_Physics",   D_eff: 17, hits: 2, delta_psi: 0.8,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Oceanography",          D_eff: 17, hits: 1, delta_psi: 0.7,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Seismology",            D_eff: 18, hits: 2, delta_psi: 1.2,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Sociology",             D_eff: 18, hits: 3, delta_psi: 1.5,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "High_Energy_Physics",   D_eff: 19, hits: 1, delta_psi: 1.2,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Geophysics",            D_eff: 19, hits: 2, delta_psi: 1.0,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Astronomy",             D_eff: 20, hits: 1, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Economics",             D_eff: 20, hits: 3, delta_psi: 1.5,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Planetary_Science",     D_eff: 21, hits: 1, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Quantum_Gravity",       D_eff: 22, hits: 0, delta_psi: 1.0,  delta_theta: 1.0, observed: false },
    DomainConfig { name: "Particle_Astrophysics", D_eff: 23, hits: 1, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Astrophysics",          D_eff: 24, hits: 1, delta_psi: 1.0,  delta_theta: 1.0, observed: true  },
    DomainConfig { name: "Cosmology",             D_eff: 25, hits: 0, delta_psi: 1.0,  delta_theta: 1.0, observed: false },
];

/// Compute the FSOT scalar S for a named domain (matches Ada `Make_Scalar_Params`).
pub fn domain_scalar(c: &FsotConsts, name: &str) -> Option<f64> {
    let d = DOMAINS.iter().find(|d| d.name == name)?;
    let mut s = ScalarInput::defaults_with(c);
    s.D_eff       = d.D_eff as f64;
    s.delta_psi   = d.delta_psi;
    s.delta_theta = d.delta_theta;
    s.recent_hits = d.hits as f64;
    s.observed    = d.observed;
    Some(compute_scalar(c, &s))
}

// =========================================================================
// §6  CHEM BOND LIBRARY  (distilled from FSOT SMILES Lab Tier-0 dataset)
// =========================================================================
//
// Each row: (atom_a, atom_b, bond_order) → (length_A, energy_kJ_per_mol).
// Both values are FSOT zero-free-parameter closed forms — `length_fn` and
// `energy_fn` recompute from the constants on demand so a kernel image stays
// honest about the formulas instead of caching numbers.
//
// Bond orders: 1 = single, 2 = double, 3 = triple.

#[derive(Debug, Clone, Copy)]
pub struct ChemBond {
    pub atom_a:    &'static str,
    pub atom_b:    &'static str,
    pub bond_order: u8,
    pub length_A:   fn(&FsotConsts) -> f64,   // angstroms
    pub energy_kJ:  fn(&FsotConsts) -> f64,   // kJ/mol
    pub formula_length: &'static str,
    pub formula_energy: &'static str,
}

// Tier-0 bond list (verified ≤2% target error in the Python lab).
pub const CHEM_BONDS: &[ChemBond] = &[
    ChemBond { atom_a: "H",  atom_b: "H",  bond_order: 1,
        length_A:  |_| m::sin(1.0) - 1.0 / (PI*PI),
        energy_kJ: |_| m::powf(E, 8.0) / m::powf(PHI, 4.0),
        formula_length: "sin(1) − π⁻²",
        formula_energy: "e⁸/φ⁴" },
    ChemBond { atom_a: "C",  atom_b: "H",  bond_order: 1,
        length_A:  |c| c.A_IN - GAMMA,
        energy_kJ: |_| m::powf(E, 6.0) + PI*PI,
        formula_length: "A_in − γ",
        formula_energy: "e⁶ + π²" },
    ChemBond { atom_a: "C",  atom_b: "C",  bond_order: 1,
        length_A:  |_| 1.0/GAMMA - m::powf(GAMMA, 3.0),
        energy_kJ: |_| m::powf(E, 6.0) - m::powf(E, 4.0),
        formula_length: "γ⁻¹ − γ³",
        formula_energy: "e⁶ − e⁴" },
    ChemBond { atom_a: "C",  atom_b: "C",  bond_order: 2,
        length_A:  |c| 1.0/(PHI*PHI) + c.P_VAR,
        energy_kJ: |c| m::powf(G_CAT, 8.0) / c.ALPHA,
        formula_length: "φ⁻² + P_var",
        formula_energy: "G⁸/α" },
    ChemBond { atom_a: "C",  atom_b: "C",  bond_order: 3,
        length_A:  |_| PI / (PHI*PHI),
        energy_kJ: |_| m::powf(E, 6.0) / m::ln(PHI),
        formula_length: "π/φ²",
        formula_energy: "e⁶/ln(φ)" },
    ChemBond { atom_a: "O",  atom_b: "H",  bond_order: 1,
        length_A:  |c| c.C_EFF + m::powf(E, -6.0),
        energy_kJ: |_| m::powf(PI, 6.0) * m::ln(PHI),
        formula_length: "C_eff + e⁻⁶",
        formula_energy: "π⁶·ln(φ)" },
    ChemBond { atom_a: "N",  atom_b: "H",  bond_order: 1,
        length_A:  |c| c.A_BLEED - m::powf(GAMMA, 6.0),
        energy_kJ: |_| m::powf(E, 6.0) - m::powf(PHI, 5.0),
        formula_length: "A_bleed − γ⁶",
        formula_energy: "e⁶ − φ⁵" },
    ChemBond { atom_a: "C",  atom_b: "N",  bond_order: 1,
        length_A:  |_| m::sqrt(2.0) + 1.0/m::powf(PHI, 6.0),
        energy_kJ: |_| 305.0,                              // not in Tier 0; literature anchor
        formula_length: "√2 + φ⁻⁶",
        formula_energy: "(literature 305 kJ/mol)" },
    ChemBond { atom_a: "C",  atom_b: "O",  bond_order: 1,
        length_A:  |c| m::ln(3.0) - c.CHAOS,
        energy_kJ: |_| 358.0,
        formula_length: "ln(3) − Chaos",
        formula_energy: "(literature 358 kJ/mol)" },
    ChemBond { atom_a: "C",  atom_b: "O",  bond_order: 2,
        length_A:  |c| m::powf(G_CAT, -8.0) - c.B_IN,
        energy_kJ: |_| 745.0,
        formula_length: "G⁻⁸ − B_in",
        formula_energy: "(literature 745 kJ/mol)" },
    ChemBond { atom_a: "N",  atom_b: "N",  bond_order: 3,
        length_A:  |_| 1.10,                               // literature
        energy_kJ: |_| m::powf(PI, 6.0) - m::powf(GAMMA, -5.0),
        formula_length: "(literature 1.10 Å)",
        formula_energy: "π⁶ − γ⁻⁵" },
    ChemBond { atom_a: "O",  atom_b: "O",  bond_order: 2,
        length_A:  |_| 1.21,
        energy_kJ: |_| m::powf(PI, 5.0) * PHI,
        formula_length: "(literature 1.21 Å)",
        formula_energy: "π⁵·φ" },
    ChemBond { atom_a: "F",  atom_b: "F",  bond_order: 1,
        length_A:  |_| 1.42,
        energy_kJ: |_| m::powf(E, 5.0) + m::powf(PHI, 5.0),
        formula_length: "(literature 1.42 Å)",
        formula_energy: "e⁵ + φ⁵" },
];

/// Look up a bond by (atom symbols, order).  Symbol order is canonicalized
/// (alphabetical) so `lookup_bond("H","C",1)` and `lookup_bond("C","H",1)` agree.
pub fn lookup_bond(a: &str, b: &str, order: u8) -> Option<&'static ChemBond> {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    CHEM_BONDS.iter().find(|cb| {
        let (clo, chi) = if cb.atom_a <= cb.atom_b { (cb.atom_a, cb.atom_b) } else { (cb.atom_b, cb.atom_a) };
        clo == lo && chi == hi && cb.bond_order == order
    })
}

// =========================================================================
//  TESTS
// =========================================================================
#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool { (a - b).abs() <= tol * b.abs().max(1e-12) }

    #[test]
    fn fsot_constants_match_python_anchors() {
        let c = FsotConsts::build();
        // Anchors lifted from FSOT_MATHEMATICAL_KEY.md / Python reference output.
        assert!(close(c.ALPHA,    0.000_808, 0.005), "ALPHA = {}", c.ALPHA);
        assert!(close(c.PSI_CON,  0.632_120, 1e-5), "PSI_CON = {}", c.PSI_CON);
        assert!(close(c.ETA_EFF,  1.0/(PI-1.0), 1e-12));
        assert!(close(c.K,        0.420_222, 1e-3), "K = {}", c.K);
        assert!(close(c.C_COSM,   1.0/(PHI*10.0), 1e-12));
        // P_NEW = γ/e · √2
        assert!(close(c.P_NEW, GAMMA/E*m::sqrt(2.0), 1e-12));
    }

    #[test]
    fn cosmology_and_qm_scalars_in_expected_band() {
        // From the §6 wave-1 H0 row:  H0 = 100·(1 + S_cosm·A_bleed/A_in) ≈ 67.4
        // → S_cosm·A_bleed/A_in ≈ -0.326,  A_bleed/A_in ≈ 0.5xx → S_cosm < 0.
        let c = FsotConsts::build();
        let s_cosm  = domain_scalar(&c, "Cosmology").unwrap();
        let h0 = 100.0 * (1.0 + s_cosm * c.A_BLEED / c.A_IN);
        assert!((h0 - 67.4).abs() / 67.4 < 0.05, "H0 = {}", h0);

        let s_quant = domain_scalar(&c, "Quantum_Mechanics").unwrap();
        // baryon density: |S_cosm|·(1 - S_quant) ≈ 0.02237
        let omega_b_h2 = s_cosm.abs() * (1.0 - s_quant);
        assert!((omega_b_h2 - 0.02237).abs() / 0.02237 < 0.10, "Ω_b·h² = {}", omega_b_h2);
    }

    #[test]
    fn bond_lookup_canonical_order() {
        let c = FsotConsts::build();
        let ch1 = lookup_bond("C", "H", 1).expect("C-H");
        let hc1 = lookup_bond("H", "C", 1).expect("H-C");
        assert_eq!(ch1.formula_length, hc1.formula_length);
        // Energy ≈ 413 kJ/mol from Tier 0 (target 413, FSOT 413.298)
        let e_ch = (ch1.energy_kJ)(&c);
        assert!((e_ch - 413.0).abs() < 5.0, "C-H energy = {}", e_ch);
    }

    #[test]
    fn dna_base_pair_anchor() {
        // Validation suite #5: 2π·a₀·(1 + γ/(π²·e)) ≈ 3.4 Å  (DNA base pair rise).
        let a0 = 0.529_177;
        let v  = 2.0 * PI * a0 * (1.0 + GAMMA / (PI * PI * E));
        assert!((v - 3.4).abs() / 3.4 < 0.02, "base-pair rise = {}", v);
    }
}
