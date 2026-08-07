//! Phase 1 audit binary — print the canonical FSOT scalars and the
//! detected α/β regions per test protein. Used to confirm the
//! formula layer is actually firing the way the math expects.

use fsot_core::{domain_scalar, FsotConsts, PI, E, PHI, GAMMA};
use fsot_protein::regions::{detect_regions, RegionKind};
use fsot_protein::secondary::{trinary_phase, SsPropensity};
use std::fs;

fn main() {
    let c = FsotConsts::build();
    let s_biochem = domain_scalar(&c, "Biochemistry").unwrap_or(0.0);
    let s_molchem = domain_scalar(&c, "Molecular_Chemistry").unwrap_or(0.0);

    println!("══════════════════════════════════════════════════════════════════");
    println!("  FSOT Canonical Constants (audit)");
    println!("══════════════════════════════════════════════════════════════════");
    println!("  π                  = {:.12}", PI);
    println!("  e                  = {:.12}", E);
    println!("  φ                  = {:.12}", PHI);
    println!("  γ                  = {:.12}", GAMMA);
    println!("  P_NEW  (γ/e)·√2    = {:.12}", c.P_NEW);
    println!("  C_EFF              = {:.12}", c.C_EFF);
    println!("  ETA_EFF 1/(π-1)    = {:.12}", c.ETA_EFF);
    println!("  S_biochem          = {:+.12}", s_biochem);
    println!("  S_molchem          = {:+.12}", s_molchem);
    println!("  chem_amp           = {:+.12}", s_molchem.abs() * c.P_NEW);
    println!("  region_amp         = {:+.12}", s_biochem.abs() * c.P_NEW * c.C_EFF);
    println!("  long_range_gate    = {}", (c.ETA_EFF * 13.0).ceil() as usize);
    println!();

    println!("══════════════════════════════════════════════════════════════════");
    println!("  Trinary phase table (F01)");
    println!("══════════════════════════════════════════════════════════════════");
    println!("  AA  (c, p, v)        p_α      p_β      p_coil    collapse");
    let aas = "ACDEFGHIKLMNPQRSTVWY";
    for ch in aas.chars() {
        let (cc, pp, vv) = trinary_phase(ch);
        let s = SsPropensity::from_amino_acid(ch);
        let kind = if s.p_alpha > 1.0 / E && s.p_alpha > s.p_beta { "H" }
                   else if s.p_beta > 1.0 / E && s.p_beta > s.p_alpha { "E" }
                   else { "C" };
        println!("  {}   ({:+.0}, {:+.0}, {:+.0})    {:.4}   {:.4}   {:.4}    {}",
                 ch, cc, pp, vv, s.p_alpha, s.p_beta, s.p_coil, kind);
    }
    println!();

    let test = [
        ("1UBQ", "data/1ubq.fasta"),
        ("1CRN", "data/1crn.fasta"),
        ("1VII", "data/1vii.fasta"),
        ("2GB1", "data/2gb1.fasta"),
        ("1ENH", "data/1enh.fasta"),
    ];

    println!("══════════════════════════════════════════════════════════════════");
    println!("  Detected regions per protein (F12)");
    println!("══════════════════════════════════════════════════════════════════");
    for (name, path) in &test {
        let raw = fs::read_to_string(path).unwrap_or_default();
        let seq: String = raw.lines()
            .filter(|l| !l.starts_with('>'))
            .collect::<String>()
            .chars()
            .filter(|c| c.is_ascii_alphabetic())
            .collect();
        let props: Vec<SsPropensity> = seq.chars()
            .map(SsPropensity::from_amino_acid)
            .collect();
        let regs = detect_regions(&props);

        let helices: Vec<_> = regs.iter().filter(|r| r.kind == RegionKind::Helix).collect();
        let strands: Vec<_> = regs.iter().filter(|r| r.kind == RegionKind::Strand).collect();
        println!("\n  {} (len {})  helices: {}, strands: {}",
                 name, seq.len(), helices.len(), strands.len());
        for r in &regs {
            let kind = match r.kind {
                RegionKind::Helix  => "H",
                RegionKind::Strand => "E",
                RegionKind::Coil   => "C",
            };
            let frag: String = seq.chars().skip(r.start).take(r.end - r.start + 1).collect();
            println!("    {} [{:3}-{:3}] L={:2}  {}", kind, r.start, r.end, r.length(), frag);
        }
    }
}
