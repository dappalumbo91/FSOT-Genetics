use crate::secondary::trinary_phase;
use fsot_core::{PI, E, PHI, GAMMA};

/// Models the chemical interaction parameters for an amino acid.
///
/// All four scalars are DERIVED from the trinary phase vector
/// (charge, polarity, volume) ∈ {-1, 0, +1}³ using only the FSOT
/// seeds {π, e, φ, γ}. No lookup table; no fitted decimals.
///
/// AAs that share a trinary signature (e.g. I, L, M, F, W all map to
/// [0,-1,+1]) intentionally share these scalars — that degeneracy is
/// the trinary collapse, not a defect. Side-chain features that go
/// beyond [charge, polarity, volume] (aromaticity, branching, sulfur,
/// indole H-bonding) must be added as separate higher-order layers
/// rather than smuggled in here.
#[derive(Debug, Clone, Copy)]
pub struct ChemicalPropensity {
    pub hydrophobicity_fsot: f64,
    pub volume_fsot:         f64,
    pub charge:              f64,
    pub dipole_moment:       f64,
}

impl ChemicalPropensity {
    pub fn from_amino_acid(aa: char) -> Self {
        let (c, p, v) = trinary_phase(aa);

        // ── Hydrophobicity ─────────────────────────────────────────
        //   h = φ^(-p) · exp(v / π)
        // Polarity sets the order of magnitude; volume modulates it.
        //   p=-1, v=+1  (large nonpolar: I/L/M/F/W):  φ·e^(1/π) ≈ 2.218
        //   p=-1, v= 0  (medium nonpolar: V/P):        φ          ≈ 1.618
        //   p=-1, v=-1  (small nonpolar: A/G):         φ·e^(-1/π) ≈ 1.181
        //   p= 0, v= 0  (neutral pivot):               1.000
        //   p=+1, v=*   (polar, e.g. T/S/Y/N/Q):       1/φ·e^(v/π) ≈ 0.45..0.85
        // Centered downstream at 1.0 (the natural FSOT pivot, p=0,v=0).
        let h = PHI.powf(-p) * (v / PI).exp();

        // ── Side-chain volume ─────────────────────────────────────
        //   V = π·e · φ^v
        // Anchored at π·e ≈ 8.54 (the FSOT contact scale), one φ-step
        // per trinary volume bin. Bin ratios 1 : φ : φ² match Zamyatnin
        // small : medium : large side-chain volumes within ~10%.
        let vol = PI * E * PHI.powf(v);

        // ── Net side-chain charge ─────────────────────────────────
        // Direct readout of the charge trit. (His's physiological +0.5
        // is unrepresentable in the trinary — currently lumped at +1.)
        let q = c;

        // ── Dipole moment ─────────────────────────────────────────
        //   μ = γ · e^(|c| + p + 1)
        // γ is the FSOT seed for the dipole scale. Exponent runs from
        // 0 (purely hydrophobic neutral: μ=γ) to 3 (charged polar like
        // K/R/D/E: μ=γ·e³≈11.6). Strictly positive, monotone in both
        // |charge| and polarity.
        let mu = GAMMA * (c.abs() + p + 1.0).exp();

        Self {
            hydrophobicity_fsot: h,
            volume_fsot:         vol,
            charge:              q,
            dipole_moment:       mu,
        }
    }
}

/// Evaluates the pairwise interaction scalar between two amino acid positions.
/// Positive => attractive (contact favorable); Negative => repulsive.
/// All terms derived from FSOT seeds {π, e, φ, γ} — no free parameters.
pub fn fsot_chemical_interaction(aa1: char, aa2: char) -> f64 {
    let c1 = aa1.to_ascii_uppercase();
    let c2 = aa2.to_ascii_uppercase();

    // ── Covalent disulfide bridge (Cys-Cys) ───────────────────────
    // Largest force in the field at φ^6 ≈ 17.94.
    // (Sequence-separation gating is F18, a future layer.)
    if c1 == 'C' && c2 == 'C' {
        return PHI.powi(6);
    }

    let p1 = ChemicalPropensity::from_amino_acid(c1);
    let p2 = ChemicalPropensity::from_amino_acid(c2);

    if p1.volume_fsot == 0.0 || p2.volume_fsot == 0.0 {
        return 0.0;
    }

    // ── Hydrophobic burial ────────────────────────────────────────
    // Centered at the FSOT neutral pivot h=1 (corresponds to p=0,v=0).
    // h_norm = (h - 1) / φ keeps magnitudes O(1).
    // Both-above-pivot AND both-below-pivot → positive (attractive in
    // aqueous phase); mixed → negative (repulsive).
    let h1 = (p1.hydrophobicity_fsot - 1.0) / PHI;
    let h2 = (p2.hydrophobicity_fsot - 1.0) / PHI;
    let hydrophobic_term = h1 * h2;

    // ── Electrostatic (Coulomb sign) ──────────────────────────────
    // Opposite charges attract: -q1*q2 > 0 when signs differ.
    // Scaled by E (FSOT seed). Magnitude e for a K-D salt bridge.
    let electrostatic_term = -p1.charge * p2.charge * E;

    // ── Dipole alignment (geometric mean, sign-blind) ─────────────
    // Damped by γ·π·e² so it never dominates hydrophobic or salt-bridge.
    let dipole_term = (p1.dipole_moment * p2.dipole_moment).sqrt()
        / (GAMMA * PI * E * E);

    hydrophobic_term + electrostatic_term + dipole_term
}

