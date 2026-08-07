//! Secondary-structure **region** detection — the trinary phase
//! collapses from per-residue propensities into contiguous runs.
//!
//! Why: residue-level helix/sheet bonuses can only score local periodicity
//! (i↔i+4) and short-range pairs. They cannot see **two separate helices
//! folding against each other** or **two β-strands pairing across the
//! whole protein**. Those are the contacts that AlphaFold gets right and
//! that the long-range precision metric measures.
//!
//! FSOT framing: each residue lives in a trinary state {α, β, coil}.
//! A *region* is a coherent run of ≥3 residues collapsed to the same
//! state — the minimum FSOT trinary triad. Once collapsed, the region
//! becomes a higher-order phase that couples to other regions via the
//! biochemistry-domain scalar.
//!
//! All thresholds derive from {π, e, φ, γ}:
//!   * α-strong / β-strong gate:  p > 1/e ≈ 0.368
//!     (uniform-prior 1/3 plus a single nat of FSOT coherence)
//!   * Minimum run length:        3  (3 = 3¹, the trinary triad)
//!   * Cross-region gap floor:    ⌈η_eff · D_biochem⌉ = ⌈0.467·13⌉ = 7
//!     (residues closer than this in sequence are already covered by
//!     local periodicity, so cross-region bonus is suppressed)

use crate::secondary::SsPropensity;
use fsot_core::{E, PI, PHI};

/// Which trinary state a residue or region has collapsed to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionKind {
    Helix,
    Strand,
    Coil,
}

#[derive(Debug, Clone, Copy)]
pub struct Region {
    pub kind:  RegionKind,
    pub start: usize, // inclusive
    pub end:   usize, // inclusive
}

impl Region {
    #[inline]
    pub fn contains(&self, i: usize) -> bool { i >= self.start && i <= self.end }
    #[inline]
    pub fn length(&self) -> usize { self.end - self.start + 1 }
}

/// Collapse each residue's (p_α, p_β, p_coil) triple to a definite kind:
///   α if p_α > 1/e and p_α > p_β
///   β if p_β > 1/e and p_β > p_α
///   else coil
#[inline]
fn collapse(p: &SsPropensity) -> RegionKind {
    // Gate = 1/e ≈ 0.368 — the natural information-theoretic
    // threshold: any class whose probability exceeds the uniform
    // entropy contribution 1/e is judged statistically supported.
    let gate = 1.0 / E;
    let alpha_wins = p.p_alpha > gate && p.p_alpha > p.p_beta;
    let beta_wins  = p.p_beta  > gate && p.p_beta  > p.p_alpha;
    if alpha_wins { RegionKind::Helix }
    else if beta_wins { RegionKind::Strand }
    else { RegionKind::Coil }
}

/// Scan the propensity sequence and emit all non-coil regions.
///
/// Minimum region length is kind-dependent (both derived from {π, e, γ}):
///   * Helix:  ⌈π + η_eff⌉ = ⌈3.14 + 0.467⌉ = 4 residues — one full
///     α-helix turn (3.6 residues per turn rounded up).
///   * Strand: 3 residues — the trinary triad minimum (3 = 3¹).
pub fn detect_regions(props: &[SsPropensity]) -> Vec<Region> {
    if props.is_empty() { return Vec::new(); }
    let min_helix:  usize = (PI + 1.0/(PI - 1.0)).ceil() as usize; // = 4
    let min_strand: usize = 3;

    // ── F12b: Quantum-tunneling bridge collapse ──
    // A residue whose two strongest propensities lie within 1/φ² ≈ 0.382
    // (the golden tie window) is in superposition between α and β.
    //
    // Mode selection (env var FSOT_TUNNEL_MODE, default "bridge"):
    //   "bridge"     — original v10a; left==right==k != coil → tunnel.
    //   "frustrated" — refuses to tunnel when neighbors are α/β-opposite;
    //                  the superposed residue stays as coil at the phase
    //                  boundary. Preserves topology where two ordered
    //                  phases meet.
    //   "coherence"  — like bridge, but also requires the destination
    //                  phase to extend ⌈φ²⌉ = 3 residues on the supporting
    //                  side. Stops short fragments from absorbing pivots.
    let tunnel_window = 1.0 / (PHI * PHI);
    let coherence_depth: usize = (PHI * PHI).ceil() as usize; // 3
    let mode = std::env::var("FSOT_TUNNEL_MODE").unwrap_or_else(|_| "frustrated".to_string());
    let initial: Vec<RegionKind> = props.iter().map(collapse).collect();
    let n = props.len();
    let collapsed: Vec<RegionKind> = (0..n).map(|i| {
        let p = &props[i];
        let top   = p.p_alpha.max(p.p_beta);
        let other = p.p_alpha.min(p.p_beta);
        let superposed = top > 0.0 && (top - other) / top < tunnel_window;
        if !superposed { return initial[i]; }
        if i == 0 || i + 1 >= n { return initial[i]; }
        let left  = initial[i - 1];
        let right = initial[i + 1];
        match mode.as_str() {
            "frustrated" => {
                // If both neighbors agree on a non-coil kind → tunnel.
                // If they disagree on TWO non-coil kinds (α vs β) →
                // refuse to collapse; the residue stays superposed
                // (rendered as Coil so it breaks the run cleanly).
                if left == right && left != RegionKind::Coil {
                    left
                } else if left != RegionKind::Coil
                       && right != RegionKind::Coil
                       && left != right {
                    RegionKind::Coil
                } else {
                    initial[i]
                }
            }
            "coherence" => {
                // Like bridge, but require the destination phase to
                // extend at least coherence_depth on the supporting side.
                if left == right && left != RegionKind::Coil {
                    let need = coherence_depth;
                    let mut depth_l = 0usize;
                    let mut k = i;
                    while k > 0 && initial[k - 1] == left && depth_l < need {
                        depth_l += 1;
                        k -= 1;
                    }
                    let mut depth_r = 0usize;
                    let mut k = i;
                    while k + 1 < n && initial[k + 1] == right && depth_r < need {
                        depth_r += 1;
                        k += 1;
                    }
                    if depth_l >= need && depth_r >= need { left } else { initial[i] }
                } else {
                    initial[i]
                }
            }
            _ /* bridge */ => {
                if left == right && left != RegionKind::Coil { left }
                else { initial[i] }
            }
        }
    }).collect();

    let mut out = Vec::new();
    let mut run_kind = collapsed[0];
    let mut run_start = 0usize;
    for i in 1..n {
        let k = collapsed[i];
        if k != run_kind {
            let len = i - run_start;
            let min_len = match run_kind {
                RegionKind::Helix  => min_helix,
                RegionKind::Strand => min_strand,
                RegionKind::Coil   => usize::MAX,
            };
            if run_kind != RegionKind::Coil && len >= min_len {
                out.push(Region { kind: run_kind, start: run_start, end: i - 1 });
            }
            run_kind = k;
            run_start = i;
        }
    }
    // close trailing run
    let last = props.len() - 1;
    let len = last + 1 - run_start;
    let min_len = match run_kind {
        RegionKind::Helix  => min_helix,
        RegionKind::Strand => min_strand,
        RegionKind::Coil   => usize::MAX,
    };
    if run_kind != RegionKind::Coil && len >= min_len {
        out.push(Region { kind: run_kind, start: run_start, end: last });
    }
    out
}

/// Map each residue index to the region it belongs to, if any.
/// Returns Vec<Option<region_idx>> with len = props.len().
pub fn residue_to_region(n: usize, regions: &[Region]) -> Vec<Option<usize>> {
    let mut map = vec![None; n];
    for (ri, r) in regions.iter().enumerate() {
        for i in r.start..=r.end {
            if i < n { map[i] = Some(ri); }
        }
    }
    map
}

/// F17 — β-strand-pair register multiplier (trinary D × R × P).
///
/// Tests both antiparallel and parallel pairing geometries for two
/// β-regions A and B, picks the one with smaller offset, and returns
/// the φ-graded coupling from the 27-cell trinary map (see
/// `beta_pair_trinary_map.txt`).
///
/// All coefficients in {π, e, φ, γ}. No free parameters.
#[inline]
pub fn beta_register_multiplier(
    i: usize, j: usize,
    a_start: usize, a_end: usize,
    b_start: usize, b_end: usize,
) -> f64 {
    // Local strand positions
    let a = i as isize - a_start as isize;
    let b = j as isize - b_start as isize;
    let b_len = (b_end - b_start) as isize;       // = (e_B - s_B)

    // Two candidate registers — pick the one with smaller |Δ|.
    let delta_ap: isize = (a + b) - b_len;        // antiparallel ideal: a + b == b_len
    let delta_p:  isize = b - a;                  // parallel ideal: b == a

    let (delta, dir_trit, pleat_parity): (isize, isize, isize) =
        if delta_ap.abs() <= delta_p.abs() {
            // antiparallel pleat: same face iff (a + (b_len - b)) even
            let p = (a + (b_len - b)).rem_euclid(2);
            (delta_ap, -1, p)
        } else {
            // parallel pleat: same face iff (a + b) even
            let p = (a + b).rem_euclid(2);
            (delta_p,  1, p)
        };

    // Cutoff for "no register fit" (D = 0): |Δ| greater than half the
    // shorter strand. Half-length is the natural FSOT scale at which
    // pairing decays past coherence. Floor at 2 so very short regions
    // still register-test.
    let half_short = (b_len.max(0) as f64 * 0.5).max(2.0);
    if (delta.abs() as f64) > half_short {
        return 1.0 / (PHI * PHI);   // cell [0, *, *]
    }

    // Register trit R from sign(Δ); |Δ|=0 → R=0 (in register).
    let r_trit: isize = if delta == 0 { 0 } else if delta > 0 { 1 } else { -1 };

    // Pleat trit P:
    //   In-register (R=0) collapses to P=0 (the H-bond pair).
    //   Otherwise P = +1 (same face) when parity even, -1 (opposite) when odd.
    let p_trit: isize = if r_trit == 0 {
        0
    } else if pleat_parity == 0 {
        1
    } else {
        -1
    };

    // Trit-cell multiplier:  φ^((1-|R|) + (1-|P|) - 1)
    let exponent = (1 - r_trit.abs()) + (1 - p_trit.abs()) - 1;
    let cell_mult = PHI.powi(exponent as i32);

    // Continuous decay for |Δ| > 1: φ^(-(|Δ|-1)/π)
    let decay = if delta.abs() > 1 {
        PHI.powf(-((delta.abs() - 1) as f64) / PI)
    } else { 1.0 };

    // Used multiplicatively in F13; consumes the D trit through dir_trit
    // selection. dir_trit is informational here (caller doesn't need it).
    let _ = dir_trit;
    cell_mult * decay
}

/// F16 — α-helix heptad register multiplier (trinary H_i × H_j).
///
/// Each residue's local position mod 7 collapses to a face trit:
///   {0, 3}    → +1 (packing — a/d face the partner helix)
///   {4, 6}    →  0 (ionic — e/g, salt-bridge edge)
///   {1, 2, 5} → -1 (outward — b/c/f, solvent-exposed)
///
/// Pair multiplier: φ^(H_i + H_j). See helix_heptad_trinary_map.txt
/// for the full 9-cell trinary map.
///
/// All coefficients in {π, e, φ, γ}. No free parameters.
#[inline]
pub fn helix_heptad_multiplier(
    i: usize, j: usize,
    a_start: usize,
    b_start: usize,
) -> f64 {
    #[inline]
    fn face_trit(local_pos: usize) -> i32 {
        match local_pos % 7 {
            0 | 3       =>  1,   // packing (a, d)
            4 | 6       =>  0,   // ionic   (e, g)
            _           => -1,   // outward (b, c, f)
        }
    }
    let a = i - a_start;
    let b = j - b_start;
    let exponent = face_trit(a) + face_trit(b);
    PHI.powi(exponent)
}
