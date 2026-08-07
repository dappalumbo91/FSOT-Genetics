use std::fs::File;
use std::io::Write;
use fsot_protein::chemical::ChemicalPropensity;

fn get_trinary_phase(aa: char) -> [i32; 3] {
    // Maps each amino acid to a [-1, 0, 1] trinary state vector [Charge, Polarity, Volume]
    // matching the 64-codon structural format perfectly.
    match aa {
        'A' => [ 0, -1, -1], // Neutral, Nonpolar, Small
        'R' => [ 1,  1,  1], // Positive, Polar, Large
        'N' => [ 0,  1,  0], // Neutral, Polar, Medium
        'D' => [-1,  1,  0], // Negative, Polar, Medium
        'C' => [ 0,  0, -1], // Neutral, Special/Amphi, Small
        'Q' => [ 0,  1,  1], // Neutral, Polar, Large
        'E' => [-1,  1,  1], // Negative, Polar, Large
        'G' => [ 0, -1, -1], // Neutral, Nonpolar, Small
        'H' => [ 1,  1,  1], // Positive, Polar, Large
        'I' => [ 0, -1,  1], // Neutral, Nonpolar, Large
        'L' => [ 0, -1,  1], // Neutral, Nonpolar, Large
        'K' => [ 1,  1,  1], // Positive, Polar, Large
        'M' => [ 0, -1,  1], // Neutral, Nonpolar, Large
        'F' => [ 0, -1,  1], // Neutral, Nonpolar, Large/Aromatic
        'P' => [ 0, -1,  0], // Neutral, Nonpolar, Medium
        'S' => [ 0,  1, -1], // Neutral, Polar, Small
        'T' => [ 0,  1,  0], // Neutral, Polar, Medium
        'W' => [ 0, -1,  1], // Neutral, Nonpolar, Large/Aromatic
        'Y' => [ 0,  1,  1], // Neutral, Polar, Large/Aromatic
        'V' => [ 0, -1,  0], // Neutral, Nonpolar, Medium
        _   => [ 0,  0,  0],
    }
}

fn main() {
    let amino_acids = [
        ('A', "Alanine"),
        ('R', "Arginine"),
        ('N', "Asparagine"),
        ('D', "Aspartic Acid"),
        ('C', "Cysteine"),
        ('Q', "Glutamine"),
        ('E', "Glutamic Acid"),
        ('G', "Glycine"),
        ('H', "Histidine"),
        ('I', "Isoleucine"),
        ('L', "Leucine"),
        ('K', "Lysine"),
        ('M', "Methionine"),
        ('F', "Phenylalanine"),
        ('P', "Proline"),
        ('S', "Serine"),
        ('T', "Threonine"),
        ('W', "Tryptophan"),
        ('Y', "Tyrosine"),
        ('V', "Valine"),
    ];

    let mut out_file = File::create("20_amino_acid_fsot_map.txt").expect("Failed to create map file");
    
    writeln!(out_file, "FSOT 20-AMINO-ACID UNIFIED STRUCTURAL MAPPING").unwrap();
    writeln!(out_file, "==========================================================================================================").unwrap();
    writeln!(out_file, " AA | TRINARY PHASE [C,P,V] | NAME              | HYDROPHOBICITY (φ) | VOLUME (π, e) | DIPOLE MOMENT (γ)").unwrap();
    writeln!(out_file, "----------------------------------------------------------------------------------------------------------").unwrap();
    
    for (aa, name) in amino_acids.iter() {
        let prop = ChemicalPropensity::from_amino_acid(*aa);
        let trinary = get_trinary_phase(*aa);
        let trinary_str = format!("[{:2}, {:2}, {:2}]", trinary[0], trinary[1], trinary[2]);
        writeln!(
            out_file,
            " {:2} | {:19} | {:<17} | {:>18.4} | {:>13.4} | {:>18.4}",
            aa, trinary_str, name, prop.hydrophobicity_fsot, prop.volume_fsot, prop.dipole_moment
        ).unwrap();
    }
    
    writeln!(out_file, "==========================================================================================================").unwrap();
    println!("Successfully generated comprehensive 20_amino_acid_fsot_map.txt in root.");
}
