use std::fs;
use std::path::Path;
use fsot_protein::distogram::Distogram;

fn main() {
    println!("══════════════════════════════════════════════════════════════════");
    println!("  FSOT Protein Evaluator — 1UBQ Benchmark");
    println!("══════════════════════════════════════════════════════════════════");

    // 1. Load the FASTA file
    let fasta_path = "data/1UBQ.fasta";
    let fasta_raw = fs::read_to_string(Path::new(fasta_path)).expect("Failed to read FASTA");
    
    // Parse out the sequence (ignoring headers)
    let mut sequence = String::new();
    for line in fasta_raw.lines() {
        if !line.starts_with('>') {
            sequence.push_str(line.trim());
        }
    }
    
    println!("Protein Sequence ({} AA):", sequence.len());
    println!("{}", sequence);
    println!("──────────────────────────────────────────────────────────────────");

    // 2. Generate the Distance topology map
    let start_time = std::time::Instant::now();
    let distogram = Distogram::new(&sequence);
    let elapsed = start_time.elapsed();

    println!("Generated N x N ({} x {}) Distogram Matrix in {:?}", distogram.size, distogram.size, elapsed);
    
    // 3. Print a small snippet of the distogram to verify interactions
    println!("Top-left 5x5 Interaction Tensor Snapshot:");
    for i in 0..5.min(distogram.size) {
        let mut row_str = String::new();
        for j in 0..5.min(distogram.size) {
            row_str.push_str(&format!("{:>8.3} ", distogram.matrix[i][j]));
        }
        println!("  [ {}]", row_str);
    }
    
    // Calculate naive interaction sum just to have an early diagnostic scalar
    let mut energy_sum = 0.0;
    for row in &distogram.matrix {
        for &val in row {
            energy_sum += val;
        }
    }
    
    println!("──────────────────────────────────────────────────────────────────");
    println!("Total FSOT Inter-residue Bond Tensor Energy: {:.4}", energy_sum);
    println!("(Target -> PDB atomic cross-validation step next)");
    println!("══════════════════════════════════════════════════════════════════");
}
