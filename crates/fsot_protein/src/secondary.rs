//! Secondary structure propensity — α-helix, β-sheet, coil — derived from
//! the same trinary phase vector [Charge, Polarity, Volume] used in the
//! 20-amino-acid map.
//!
//! Each residue maps to a triple of FSOT-derived propensities:
//!     [P_α, P_β, P_coil]   (all >0, normalized so they sum to 1)
//!
//! The propensities are computed from the trinary state via FSOT seeds
//! {π, e, φ, γ} — no fitted coefficients. The physics:
//!
//!   * α-helix likes:  small, nonpolar, neutral residues that can H-bond
//!                     along the backbone (i ↔ i+4). Pro & Gly disrupt it.
//!   * β-sheet likes:  large, hydrophobic, β-branched side chains (V, I, F, T)
//!                     that pack against neighboring strands.
//!   * Coil:           polar, charged, flexible, or helix-breakers.
//!
//! These map naturally onto the [Charge, Polarity, Volume] phase signature.

use fsot_core::{PI, E, PHI};

/// Secondary-structure propensity triple for one amino acid.
#[derive(Debug, Clone, Copy)]
pub struct SsPropensity {
    pub p_alpha: f64,
    pub p_beta:  f64,
    pub p_coil:  f64,
}

impl SsPropensity {
    /// Maps an amino acid character to its (α, β, coil) propensity triple.
    /// Values are derived from the trinary phase vector via FSOT scalars.
    pub fn from_amino_acid(aa: char) -> Self {
        let aa = aa.to_ascii_uppercase();

        // Special cases: helix breakers and end-cappers
        // Proline kinks the backbone; Glycine has no Cβ → flexible.
        if aa == 'P' {
            return Self::normalize(1.0 / PHI, 1.0 / PHI, PHI);  // strongly coil-favoring
        }
        if aa == 'G' {
            return Self::normalize(1.0 / E, 1.0 / E, E);        // moderately coil
        }

        // Trinary phase signature: [Charge, Polarity, Volume] ∈ {-1, 0, +1}
        let (charge, polarity, volume) = trinary_phase(aa);

        // α-helix: high baseline (the i,i+4 backbone H-bond does not care
        // about side-chain polarity — almost any residue except G and P
        // can sit in a helix). Small additive penalties for polar/charged
        // side chains that occasionally H-bond with the backbone.
        //   raw_α = φ − polarity/(π·φ) − |charge|/π²
        //   Range: 1.320 (polar+charged) … 1.815 (nonpolar+neutral)
        let raw_alpha = PHI
                      - polarity / (PI * PHI)
                      - charge.abs() / (PI * PI);

        // β-sheet: extended-backbone steric preference. Peaks at large,
        // nonpolar side chains (v=+1, p=-1).
        //   raw_β = exp((volume − polarity)/π)
        //   Range: exp(-2/π) ≈ 0.530 … exp(+2/π) ≈ 1.886
        let raw_beta = ((volume - polarity) / PI).exp();

        // Coil: small or polar side chains favor random-coil ensembles.
        //   raw_coil = exp((polarity − volume + |charge|/φ)/π)
        //   Range: exp(-2/π) ≈ 0.530 … exp((2 + 1/φ)/π) ≈ 2.30
        let raw_coil = ((polarity - volume + charge.abs() / PHI) / PI).exp();

        // All three strictly positive for trinary inputs in {-1,0,+1}³.
        Self::normalize(raw_alpha, raw_beta, raw_coil)
    }

    fn normalize(a: f64, b: f64, c: f64) -> Self {
        let s = a + b + c;
        Self { p_alpha: a / s, p_beta: b / s, p_coil: c / s }
    }

    /// Returns the most likely secondary structure as a single character:
    /// 'H' = helix, 'E' = sheet, 'C' = coil.
    pub fn dominant(&self) -> char {
        if self.p_alpha >= self.p_beta && self.p_alpha >= self.p_coil { 'H' }
        else if self.p_beta >= self.p_coil { 'E' }
        else { 'C' }
    }
}

/// Returns the trinary phase signature [Charge, Polarity, Volume] for an AA.
/// Mirrors `get_trinary_phase` in the AA map generator — single source of truth.
pub fn trinary_phase(aa: char) -> (f64, f64, f64) {
    let t: [i32; 3] = match aa.to_ascii_uppercase() {
        'A' => [ 0, -1, -1],
        'R' => [ 1,  1,  1],
        'N' => [ 0,  1,  0],
        'D' => [-1,  1,  0],
        'C' => [ 0,  0, -1],
        'Q' => [ 0,  1,  1],
        'E' => [-1,  1,  1],
        'G' => [ 0, -1, -1],
        'H' => [ 1,  1,  1],
        'I' => [ 0, -1,  1],
        'L' => [ 0, -1,  1],
        'K' => [ 1,  1,  1],
        'M' => [ 0, -1,  1],
        'F' => [ 0, -1,  1],
        'P' => [ 0, -1,  0],
        'S' => [ 0,  1, -1],
        'T' => [ 0,  1,  0],
        'W' => [ 0, -1,  1],
        'Y' => [ 0,  1,  1],
        'V' => [ 0, -1,  0],
        _   => [ 0,  0,  0],
    };
    (t[0] as f64, t[1] as f64, t[2] as f64)
}

/// Helical-pair bonus: residues at i, i+4 (also i+3) along an α-helix
/// sit ~5.4 Å / ~5.0 Å apart in 3D, well within contact range.
/// Returns a positive proximity bonus when both residues prefer α and
/// the sequence separation matches the helical period (3, 4, or 7).
pub fn helix_periodicity_bonus(p_i: &SsPropensity, p_j: &SsPropensity, sep: usize) -> f64 {
    let helical_period = matches!(sep, 3 | 4 | 7);
    if !helical_period { return 0.0; }
    // Geometric mean of both α-propensities, cubed for nonlinear amplification.
    // Cube self-attenuates: 0.55^3 = 0.166 (strong helix pair) vs 0.20^3 = 0.008
    // (weak pair). No artificial floor needed — physics handles selection.
    let joint = (p_i.p_alpha * p_j.p_alpha).sqrt();
    joint.powi(3) / E
}

pub fn sheet_pair_bonus(p_i: &SsPropensity, p_j: &SsPropensity, sep: usize) -> f64 {
    if sep < 3 { return 0.0; }
    let joint = (p_i.p_beta * p_j.p_beta).sqrt();
    // Heavy-tailed envelope: peaks at short sep, decays logarithmically
    // so long-range β-pairs (ubiquitin β1-β5 at sep≈60) still register.
    let s = sep as f64;
    let envelope = 1.0 / (1.0 + (s / PI).ln().max(0.0));
    joint.powi(2) * envelope / PHI
}
