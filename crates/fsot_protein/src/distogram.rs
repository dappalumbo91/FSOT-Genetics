use crate::chemical::fsot_chemical_interaction;
use crate::regions::{beta_register_multiplier, detect_regions, helix_heptad_multiplier, residue_to_region, RegionKind};
use crate::secondary::{SsPropensity, helix_periodicity_bonus, sheet_pair_bonus};
use fsot_core::{domain_scalar, FsotConsts, PI, E};
use std::sync::OnceLock;

/// Cached FSOT canonical constants + biochemistry-domain scalar.
/// Built once on first use; pure function of {π, e, φ, γ} so the value
/// is deterministic across runs.
struct FsotProteinConsts {
    c: FsotConsts,
    s_biochem:  f64,   // domain_scalar("Biochemistry"), D_eff = 13
    s_molchem:  f64,   // domain_scalar("Molecular_Chemistry"), D_eff = 9
}

fn fsot() -> &'static FsotProteinConsts {
    static CACHE: OnceLock<FsotProteinConsts> = OnceLock::new();
    CACHE.get_or_init(|| {
        let c = FsotConsts::build();
        let s_biochem = domain_scalar(&c, "Biochemistry").unwrap_or(0.0);
        let s_molchem = domain_scalar(&c, "Molecular_Chemistry").unwrap_or(0.0);
        FsotProteinConsts { c, s_biochem, s_molchem }
    })
}

/// Represents an N x N interaction matrix for a sequence of amino acids.
/// Cell (i,j) = predicted spatial PROXIMITY (larger ⇒ closer in 3D).
#[derive(Debug)]
pub struct Distogram {
    pub size: usize,
    pub matrix: Vec<Vec<f64>>,
}

impl Distogram {
    /// Generates a new distogram from a raw amino-acid character sequence.
    ///
    /// Five FSOT-derived layers — all amplitudes from the canonical
    /// scalar set {`P_NEW`, `C_EFF`, `ETA_EFF`, `S_biochem`, `S_molchem`}
    /// which themselves derive from {π, e, φ, γ}. No fitted parameters.
    ///
    ///   L1  Backbone — Flory random walk, proximity ∝ 1/√|i-j|.
    ///   L2  Chemistry — hydrophobic + electrostatic + dipole + S-S,
    ///       gated by the molecular-chemistry envelope.
    ///   L3  Helix periodicity — local i↔i+{3,4,7} α-pair bonus.
    ///   L4  Sheet pairing — short-range β-pair bonus.
    ///   L5  Region-pair contact — long-range helix-helix and strand-
    ///       strand bonus once trinary phases have collapsed into
    ///       coherent runs (the FSOT triadic-collapse mechanism).
    pub fn new(sequence: &str) -> Self {
        let chars: Vec<char> = sequence.chars().filter(|c| !c.is_whitespace()).collect();
        let size = chars.len();
        let fk = fsot();

        let props: Vec<SsPropensity> = chars.iter()
            .map(|&c| SsPropensity::from_amino_acid(c))
            .collect();

        // ── L5 pre-pass: detect collapsed α/β regions ──────────────
        let regions   = detect_regions(&props);
        let res_to_rg = residue_to_region(size, &regions);

        // FSOT amplitudes for each layer (zero free parameters).
        // Chemistry amplitude: |S_molchem| · P_NEW. No max()-guard — the
        // value of S_molchem is deterministic; if it's small, chemistry
        // is honestly small.
        let chem_amp   = fk.s_molchem.abs() * fk.c.P_NEW;
        // Region-pair amplitude: |S_biochem| · P_NEW · C_EFF.
        let region_amp = fk.s_biochem.abs() * fk.c.P_NEW * fk.c.C_EFF;
        // Long-range gate: ⌈η_eff · D_biochem⌉  (η_eff ≈ 0.467, D = 13 → 7).
        let long_range_gate = (fk.c.ETA_EFF * 13.0).ceil() as usize;

        let mut matrix = vec![vec![0.0; size]; size];

        for i in 0..size {
            for j in 0..size {
                if i == j { continue; }

                let sep = (i as isize - j as isize).unsigned_abs() as usize;
                let s   = sep as f64;

                // L1: backbone proximity. Folded proteins live in the
                // collapsed-globule regime, not the random-walk regime,
                // so the Flory ½ exponent is wrong physics AND a free
                // parameter. FSOT-native replacement: 1/π ≈ 0.318, the
                // collapsed-globule exponent expressible in seeds.
                let backbone_proximity = 1.0 / s.powf(1.0 / PI);

                // L2: chemistry (envelope peaks at large s, normalized by π·e)
                let interaction = fsot_chemical_interaction(chars[i], chars[j]);
                let chem_env    = s / (s + PI * E);
                let chemistry   = interaction * chem_env * chem_amp;

                // L3/L4: local α and β
                let helix = helix_periodicity_bonus(&props[i], &props[j], sep);
                let sheet = sheet_pair_bonus       (&props[i], &props[j], sep);

                // L5: cross-region trinary-phase contact bonus.
                // Active only when both residues belong to non-coil regions
                // of the same kind AND those regions are different AND the
                // sequence gap exceeds the long-range gate.
                let region_pair = match (res_to_rg[i], res_to_rg[j]) {
                    (Some(ri), Some(rj)) if ri != rj && sep >= long_range_gate => {
                        let r_i = &regions[ri];
                        let r_j = &regions[rj];
                        if r_i.kind == r_j.kind && r_i.kind != RegionKind::Coil {
                            // Joint propensity of the residue pair in their
                            // respective regions, geometric-mean style.
                            let (pi_v, pj_v) = match r_i.kind {
                                RegionKind::Helix  => (props[i].p_alpha, props[j].p_alpha),
                                RegionKind::Strand => (props[i].p_beta,  props[j].p_beta),
                                RegionKind::Coil   => (0.0, 0.0),
                            };
                            let joint = (pi_v * pj_v).sqrt();
                            // F16/F17 — register-aware coupling multipliers.
                            let reg_mult = match r_i.kind {
                                RegionKind::Strand => beta_register_multiplier(
                                    i, j,
                                    r_i.start, r_i.end,
                                    r_j.start, r_j.end,
                                ),
                                RegionKind::Helix => helix_heptad_multiplier(
                                    i, j, r_i.start, r_j.start,
                                ),
                                _ => 1.0,
                            };
                            // F13 — region-pair contact. Magnitude must
                            // stay SECONDARY to chemistry (which is the
                            // direct/local force). No length factor — that
                            // was a free heuristic that made F13 dominate.
                            // The trinary register multiplier IS the
                            // physical localization; joint propensity IS
                            // the per-residue weight.
                            joint * region_amp * reg_mult
                        } else { 0.0 }
                    }
                    _ => 0.0,
                };

                matrix[i][j] = backbone_proximity + chemistry + helix + sheet + region_pair;
            }
        }

        Self { size, matrix }
    }
}

