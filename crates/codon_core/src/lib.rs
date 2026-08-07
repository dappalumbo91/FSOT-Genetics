//! FSOT Codon Trinary Core — `#![no_std]`-clean.
//!
//! This module is the **bare-metal-ready** core of the genetics layer.
//! It holds:
//!   * the 64-codon → trinary map (primary purine/pyrimidine and secondary A-T axes)
//!   * codon → amino-acid (standard genetic code) translation
//!   * Watson–Crick base-pair bond classification (H-bond count)
//!   * FSOT-derived ΔG closed forms for base-pair and stacking energies
//!
//! Everything here is `no_std`-safe: no allocation, no `String`, no `std::fs`.
//! The `genetics.rs` module wraps these primitives with std I/O for the TUI
//! and CLI; the future `kernel/` crate will link this file directly.
//!
//! Author: Damian Arthur Palumbo — FSOT.

#![allow(non_snake_case)]
#![allow(dead_code)]

// --- FSOT seeds (re-declared so this file has no `use` cross-module deps) ---
pub const PI:    f64 = 3.141_592_653_589_793_2;
pub const E:     f64 = 2.718_281_828_459_045_2;
pub const PHI:   f64 = 1.618_033_988_749_894_8;
pub const GAMMA: f64 = 0.577_215_664_901_532_9;

/// Trit type for codon encodings: -1, 0, +1.
pub type Trit = i8;

/// One codon = 3 nucleotides → 3 trits per axis (primary, secondary).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CodonTrinary {
    pub primary:   [Trit; 3], // A,G = +1 ; C,T = -1
    pub secondary: [Trit; 3], // A = +1 ; T = -1 ; G,C = 0
}

/// Encode one nucleotide on the primary axis (purine/pyrimidine).
#[inline]
pub const fn nt_primary(b: u8) -> Trit {
    match b {
        b'A' | b'a' | b'G' | b'g' =>  1,
        b'C' | b'c' | b'T' | b't' | b'U' | b'u' => -1,
        _ => 0, // N or any ambiguous base
    }
}

/// Encode one nucleotide on the secondary axis (A–T balance).
#[inline]
pub const fn nt_secondary(b: u8) -> Trit {
    match b {
        b'A' | b'a' =>  1,
        b'T' | b't' | b'U' | b'u' => -1,
        b'G' | b'g' | b'C' | b'c' =>  0,
        _ => 0,
    }
}

/// Encode a 3-nt codon to its dual-axis trinary form.
/// Returns `None` if any nucleotide is unknown (N, gap, etc.).
#[inline]
pub fn encode_codon(c: [u8; 3]) -> CodonTrinary {
    CodonTrinary {
        primary:   [nt_primary(c[0]),   nt_primary(c[1]),   nt_primary(c[2])],
        secondary: [nt_secondary(c[0]), nt_secondary(c[1]), nt_secondary(c[2])],
    }
}

/// Pack a CodonTrinary into a single u16 (3 trits × 2 axes = 6 trits).
/// Trit mapping: -1 → 0, 0 → 1, +1 → 2.  Layout: [s0|s1|s2|p0|p1|p2], little-end.
pub fn pack_codon(t: CodonTrinary) -> u16 {
    let to_u = |x: Trit| -> u16 { (x + 1) as u16 };
    let p = to_u(t.primary[0]) + to_u(t.primary[1]) * 3 + to_u(t.primary[2]) * 9;
    let s = to_u(t.secondary[0]) + to_u(t.secondary[1]) * 3 + to_u(t.secondary[2]) * 9;
    p + s * 27 // 0..729 (fits u16)
}

/// Watson–Crick base-pair: returns the complementary nucleotide.
#[inline]
pub const fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' | b'U' | b'u' => b'A',
        b'G' | b'g' => b'C',
        b'C' | b'c' => b'G',
        _ => b'N',
    }
}

/// Hydrogen-bond count for a Watson–Crick pair.  G≡C = 3, A=T = 2, else 0.
#[inline]
pub const fn pair_hbonds(a: u8, b: u8) -> u8 {
    match (a, b) {
        (b'G', b'C') | (b'C', b'G') | (b'g', b'c') | (b'c', b'g') => 3,
        (b'A', b'T') | (b'T', b'A') | (b'a', b't') | (b't', b'a') => 2,
        (b'A', b'U') | (b'U', b'A') | (b'a', b'u') | (b'u', b'a') => 2,
        _ => 0,
    }
}

/// FSOT zero-free-parameter base-pair binding ΔG (kJ/mol, magnitude).
/// G≡C  : 3 H-bonds + π-stack → φ·(π + e)   ≈ 9.479 kJ/mol
/// A=T  : 2 H-bonds + π-stack → φ·π          ≈ 5.083 kJ/mol
/// Anchor: nearest-neighbor SantaLucia 1998 means (GC ≈ 9.5, AT ≈ 5.0).
#[inline]
pub fn pair_dG_kJ_per_mol(a: u8, b: u8) -> f64 {
    match pair_hbonds(a, b) {
        3 => PHI * (PI + E),
        2 => PHI * PI,
        _ => 0.0,
    }
}

/// FSOT zero-free-parameter stacking energy (kJ/mol) between two adjacent base pairs.
/// Closed form: γ · (1 + cos(π·δH/3))  where δH = |H₁ − H₂| (H-bond count diff).
/// Equal stacks (GC/GC or AT/AT): γ·(1+1)·... ≈ 1.155 kJ/mol; mixed: smaller.
#[inline]
pub fn stack_dG_kJ_per_mol(pair_a: (u8, u8), pair_b: (u8, u8)) -> f64 {
    let h1 = pair_hbonds(pair_a.0, pair_a.1) as f64;
    let h2 = pair_hbonds(pair_b.0, pair_b.1) as f64;
    let dh = (h1 - h2).abs();
    GAMMA * (1.0 + (PI * dh / 3.0).cos())
}

/// Standard genetic code (DNA, IUPAC).  Returns single-letter amino acid or
/// `*` for stop codons, `X` for unknown / ambiguous.
pub fn codon_to_aa(c: [u8; 3]) -> u8 {
    let to_upper = |b: u8| if (b'a'..=b'z').contains(&b) { b - 32 } else { b };
    let u = [to_upper(c[0]), to_upper(c[1]), to_upper(c[2])];
    match &u {
        b"TTT"|b"TTC" => b'F',
        b"TTA"|b"TTG"|b"CTT"|b"CTC"|b"CTA"|b"CTG" => b'L',
        b"ATT"|b"ATC"|b"ATA" => b'I',
        b"ATG" => b'M',                                  // also START
        b"GTT"|b"GTC"|b"GTA"|b"GTG" => b'V',
        b"TCT"|b"TCC"|b"TCA"|b"TCG"|b"AGT"|b"AGC" => b'S',
        b"CCT"|b"CCC"|b"CCA"|b"CCG" => b'P',
        b"ACT"|b"ACC"|b"ACA"|b"ACG" => b'T',
        b"GCT"|b"GCC"|b"GCA"|b"GCG" => b'A',
        b"TAT"|b"TAC" => b'Y',
        b"TAA"|b"TAG"|b"TGA" => b'*',                    // STOP
        b"CAT"|b"CAC" => b'H',
        b"CAA"|b"CAG" => b'Q',
        b"AAT"|b"AAC" => b'N',
        b"AAA"|b"AAG" => b'K',
        b"GAT"|b"GAC" => b'D',
        b"GAA"|b"GAG" => b'E',
        b"TGT"|b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT"|b"CGC"|b"CGA"|b"CGG"|b"AGA"|b"AGG" => b'R',
        b"GGT"|b"GGC"|b"GGA"|b"GGG" => b'G',
        _ => b'X',
    }
}

/// Is this codon a translation START?  (ATG only — methionine.)
#[inline]
pub fn is_start(c: [u8; 3]) -> bool { codon_to_aa(c) == b'M' && c[0].to_ascii_uppercase() == b'A' }

/// Is this codon a STOP?  (TAA / TAG / TGA.)
#[inline]
pub fn is_stop(c: [u8; 3]) -> bool { codon_to_aa(c) == b'*' }

/// Single-pass running observable accumulator for a FASTA / nucleotide stream.
/// All counts are u64 (handles a full human genome ≈ 3.1 Gbp without overflow).
#[derive(Debug, Default, Clone, Copy)]
pub struct GenomeObservables {
    pub a: u64, pub c: u64, pub g: u64, pub t: u64, pub n: u64,
    pub codons_total: u64,
    pub codons_start: u64,
    pub codons_stop:  u64,
    /// Cumulative primary-axis sign (Σ trits) — should drift toward 0 in healthy DNA.
    pub primary_sum:   i64,
    /// Cumulative secondary-axis sign — Chargaff: A=T ⇒ this drifts toward 0.
    pub secondary_sum: i64,
}

impl GenomeObservables {
    pub const fn new() -> Self { Self { a:0, c:0, g:0, t:0, n:0, codons_total:0, codons_start:0, codons_stop:0, primary_sum:0, secondary_sum:0 } }

    /// Accumulate a single nucleotide byte.
    #[inline]
    pub fn push_nt(&mut self, b: u8) {
        match b {
            b'A'|b'a' => { self.a += 1; self.primary_sum += 1; self.secondary_sum += 1; }
            b'C'|b'c' => { self.c += 1; self.primary_sum -= 1; }
            b'G'|b'g' => { self.g += 1; self.primary_sum += 1; }
            b'T'|b't'|b'U'|b'u' => { self.t += 1; self.primary_sum -= 1; self.secondary_sum -= 1; }
            _ => { self.n += 1; }
        }
    }

    /// Accumulate one codon (after coverage already counted at the nt level).
    #[inline]
    pub fn push_codon(&mut self, c: [u8; 3]) {
        self.codons_total += 1;
        if is_start(c) { self.codons_start += 1; }
        if is_stop(c)  { self.codons_stop  += 1; }
    }

    /// Total counted nucleotides (excluding N).
    pub fn callable(&self) -> u64 { self.a + self.c + self.g + self.t }

    /// GC content as fraction (callable-only denominator).
    pub fn gc_frac(&self) -> f64 {
        let n = self.callable() as f64;
        if n == 0.0 { 0.0 } else { (self.g + self.c) as f64 / n }
    }

    /// Chargaff parity ratios (A/T and G/C).  Both should be ≈ 1 genome-wide.
    pub fn chargaff(&self) -> (f64, f64) {
        let at = if self.t == 0 { 0.0 } else { self.a as f64 / self.t as f64 };
        let gc = if self.c == 0 { 0.0 } else { self.g as f64 / self.c as f64 };
        (at, gc)
    }

    /// FSOT super/spin density observables (per FSOT_BIOLOGICAL_BRAIN derivations).
    /// ρ_super = secondary_sum / N  (A-T axis density)
    /// ρ_spin  = primary_sum   / N  (purine-pyrimidine axis density)
    pub fn fsot_densities(&self) -> (f64, f64) {
        let n = self.callable() as f64;
        if n == 0.0 { (0.0, 0.0) } else { (self.secondary_sum as f64 / n, self.primary_sum as f64 / n) }
    }
}

// ============================================================================
//  Tests (still no_std-clean except for the `#[test]` attribute — runs in std builds.)
// ============================================================================
#[cfg(test)]
mod codon_core_tests {
    use super::*;

    #[test]
    fn nt_axes_match_user_table() {
        // From 64_codon_trinary_map.txt — spot-check a few rows.
        let aaa = encode_codon(*b"AAA");
        assert_eq!(aaa.primary,   [1, 1, 1]);
        assert_eq!(aaa.secondary, [1, 1, 1]);
        let ttt = encode_codon(*b"TTT");
        assert_eq!(ttt.primary,   [-1, -1, -1]);
        assert_eq!(ttt.secondary, [-1, -1, -1]);
        let gca = encode_codon(*b"GCA");
        assert_eq!(gca.primary,   [1, -1, 1]);
        assert_eq!(gca.secondary, [0, 0, 1]);
        let cgt = encode_codon(*b"CGT");
        assert_eq!(cgt.primary,   [-1, 1, -1]);
        assert_eq!(cgt.secondary, [0, 0, -1]);
    }

    #[test]
    fn pack_is_unique_for_64_codons() {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut seen = [false; 729];
        let mut n = 0;
        for &x in &bases { for &y in &bases { for &z in &bases {
            let t = encode_codon([x,y,z]);
            let p = pack_codon(t) as usize;
            assert!(p < 729);
            assert!(!seen[p], "duplicate pack for {}{}{}", x as char, y as char, z as char);
            seen[p] = true;
            n += 1;
        } } }
        assert_eq!(n, 64);
    }

    #[test]
    fn standard_code_anchors() {
        assert_eq!(codon_to_aa(*b"ATG"), b'M');
        assert_eq!(codon_to_aa(*b"TAA"), b'*');
        assert_eq!(codon_to_aa(*b"TAG"), b'*');
        assert_eq!(codon_to_aa(*b"TGA"), b'*');
        assert_eq!(codon_to_aa(*b"GCT"), b'A');
        assert_eq!(codon_to_aa(*b"GGG"), b'G');
        assert!(is_start(*b"ATG"));
        assert!(is_stop(*b"TGA"));
    }

    #[test]
    fn hbond_counts_textbook() {
        assert_eq!(pair_hbonds(b'G', b'C'), 3);
        assert_eq!(pair_hbonds(b'A', b'T'), 2);
        assert_eq!(pair_hbonds(b'A', b'G'), 0);
    }

    #[test]
    fn pair_dG_within_5pct_of_santalucia() {
        // SantaLucia 1998 mean nearest-neighbor: GC ≈ 9.5 kJ/mol, AT ≈ 5.0 kJ/mol.
        let gc = pair_dG_kJ_per_mol(b'G', b'C');
        let at = pair_dG_kJ_per_mol(b'A', b'T');
        assert!((gc - 9.5).abs() / 9.5 < 0.05, "GC ΔG off: {}", gc);
        assert!((at - 5.0).abs() / 5.0 < 0.05, "AT ΔG off: {}", at);
    }

    #[test]
    fn observables_chargaff_balanced_sample() {
        let dna = b"ACGTACGTACGT"; // exactly balanced
        let mut o = GenomeObservables::new();
        for &b in dna { o.push_nt(b); }
        let (at, gc) = o.chargaff();
        assert!((at - 1.0).abs() < 1e-12);
        assert!((gc - 1.0).abs() < 1e-12);
        assert!((o.gc_frac() - 0.5).abs() < 1e-12);
        let (rs, rp) = o.fsot_densities();
        assert!(rs.abs() < 1e-12 && rp.abs() < 1e-12);
    }
}
