use std::fs;
use std::path::Path;
use fsot_protein::distogram::Distogram;

/// Represents a 3D coordinate from a PDB file.
#[derive(Debug, Clone, Copy)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    pub fn distance(self, other: &Point3D) -> f64 {
        ((self.x - other.x).powi(2)
            + (self.y - other.y).powi(2)
            + (self.z - other.z).powi(2))
        .sqrt()
    }
}

/// Parses a PDB file and returns the CA (Alpha-Carbon) backbone sequence and X,Y,Z coordinates.
fn parse_pdb_ca_trace(pdb_path: &str) -> (String, Vec<Point3D>) {
    let raw = fs::read_to_string(Path::new(pdb_path)).expect("Failed to read PDB file");
    let mut coords = Vec::new();
    let mut sequence = String::new();

    // Standard amino acid 3-to-1 letter mapping
    let get_1_letter = |triple: &str| -> char {
        match triple {
            "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D',
            "CYS" => 'C', "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G',
            "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
            "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S',
            "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
            _ => 'X',
        }
    };

    for line in raw.lines() {
        // Find ATOM records representing the Alpha Carbon (CA)
        if line.starts_with("ATOM") && line.len() >= 54 {
            let atom_name = line[12..16].trim();
            if atom_name == "CA" {
                let res_name = line[17..20].trim();
                
                // Extract coordinates (standard PDB format specific indices)
                let x_str = line[30..38].trim();
                let y_str = line[38..46].trim();
                let z_str = line[46..54].trim();

                if let (Ok(x), Ok(y), Ok(z)) = (x_str.parse::<f64>(), y_str.parse::<f64>(), z_str.parse::<f64>()) {
                    coords.push(Point3D { x, y, z });
                    sequence.push(get_1_letter(res_name));
                }
            }
        }
    }

    (sequence, coords)
}

/// Calculates the Pearson Correlation Coefficient between two vectors
fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len());
    let n = x.len() as f64;
    let sum_x: f64 = x.iter().sum();
    let sum_y: f64 = y.iter().sum();
    let sum_x_sq: f64 = x.iter().map(|v| v * v).sum();
    let sum_y_sq: f64 = y.iter().map(|v| v * v).sum();
    let sum_xy: f64 = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum();

    let numerator = n * sum_xy - sum_x * sum_y;
    let denominator = ((n * sum_x_sq - sum_x * sum_x) * (n * sum_y_sq - sum_y * sum_y)).sqrt();

    if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

/// Spearman rank correlation — measures monotonic relationship, not magnitude.
/// Far more informative than Pearson for contact-map quality.
fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
    fn rank(values: &[f64]) -> Vec<f64> {
        let n = values.len();
        let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let mut ranks = vec![0.0; n];
        let mut i = 0;
        while i < n {
            let mut j = i;
            while j + 1 < n && indexed[j + 1].1 == indexed[i].1 {
                j += 1;
            }
            let avg_rank = (i + j) as f64 / 2.0 + 1.0;
            for k in i..=j {
                ranks[indexed[k].0] = avg_rank;
            }
            i = j + 1;
        }
        ranks
    }
    pearson_correlation(&rank(x), &rank(y))
}

/// Contact-map precision: of the top-N predicted contacts (highest predicted
/// proximity), what fraction are actually in contact (< 8 Å Cα-Cα) in reality?
/// This is the CASP-style metric AlphaFold reports.
fn contact_precision(pred: &[f64], actual_dist: &[f64], top_n: usize) -> f64 {
    let mut indices: Vec<usize> = (0..pred.len()).collect();
    indices.sort_by(|&a, &b| pred[b].partial_cmp(&pred[a]).unwrap());
    let top = &indices[..top_n.min(indices.len())];
    let hits = top.iter().filter(|&&i| actual_dist[i] < 8.0).count();
    hits as f64 / top.len() as f64
}

fn benchmark_protein(name: &str, pdb_path: &str) {
    println!("──────────────────────────────────────────────────────────────────");
    println!("  RUNNING BENCHMARK: {} ", name);
    println!("──────────────────────────────────────────────────────────────────");

    // 1. Parse observed reality from PDB
    let (sequence, true_coords) = parse_pdb_ca_trace(pdb_path);
    let size = sequence.len();
    println!("Extracted Sequence ({} AA):", size);
    println!("{}", sequence);

    if size == 0 {
        println!("Error: No CA coordinates found in PDB.");
        return;
    }

    // 2. Generate actual observed distance matrix
    let mut actual_distances = vec![vec![0.0; size]; size];
    for i in 0..size {
        for j in 0..size {
            actual_distances[i][j] = true_coords[i].distance(&true_coords[j]);
        }
    }

    // 3. Generate FSOT-derived predictive mathematical matrix
    let start_time = std::time::Instant::now();
    let distogram = Distogram::new(&sequence);
    let elapsed = start_time.elapsed();
    println!("Generated FSOT Distogram Matrix in {:?}", elapsed);

    // 4. Flatten the upper-triangular non-diagonal matrices to compare sets.
    //    Track long-range pairs (|i-j| >= 6) separately — those are the
    //    folding contacts AlphaFold is actually scored on. Short-range
    //    pairs are trivially predicted by the backbone term.
    let mut fsot_vals = Vec::new();
    let mut true_vals = Vec::new();
    let mut fsot_dist_vals = Vec::new();
    let mut true_dist_vals = Vec::new();
    let mut fsot_lr = Vec::new();
    let mut true_dist_lr = Vec::new();

    for i in 0..size {
        for j in (i + 1)..size {
            let pred_proximity = distogram.matrix[i][j];
            let actual_dist = actual_distances[i][j];
            let actual_proximity = if actual_dist > 0.1 { 1.0 / actual_dist } else { 10.0 };

            fsot_vals.push(pred_proximity);
            true_vals.push(actual_proximity);
            fsot_dist_vals.push(pred_proximity);
            true_dist_vals.push(actual_dist);

            if j - i >= 6 {
                fsot_lr.push(pred_proximity);
                true_dist_lr.push(actual_dist);
            }
        }
    }

    // 5. Evaluate Correlation suite
    let pearson_all = pearson_correlation(&fsot_vals, &true_vals).abs() * 100.0;
    let spearman_all = spearman_correlation(&fsot_vals, &true_vals).abs() * 100.0;
    let top_contacts = size; // CASP "top-L" precision
    let prec_all = contact_precision(&fsot_vals, &true_dist_vals, top_contacts) * 100.0;
    let prec_lr = if !fsot_lr.is_empty() {
        contact_precision(&fsot_lr, &true_dist_lr, size / 2) * 100.0
    } else {
        0.0
    };

    println!("Total pairwise bonds evaluated: {} (long-range: {})", fsot_vals.len(), fsot_lr.len());
    println!("Pearson Correlation       : {:>6.2}%", pearson_all);
    println!("Spearman Rank Correlation : {:>6.2}%   <- contact ordering quality", spearman_all);
    println!("Top-L Contact Precision   : {:>6.2}%   <- CASP-style", prec_all);
    println!("Long-Range Contact Prec.  : {:>6.2}%   <- the real folding test (|i-j|≥6)", prec_lr);
    println!("──────────────────────────────────────────────────────────────────");
}

fn main() {
    println!("══════════════════════════════════════════════════════════════════");
    println!("  FSOT Protein Evaluation Benchmark vs Target PDB Physical Reality");
    println!("══════════════════════════════════════════════════════════════════");

    benchmark_protein("1UBQ (Ubiquitin, α/β fold)",      "data/1UBQ.pdb");
    benchmark_protein("1CRN (Crambin, S-S stabilized)",  "data/1CRN.pdb");
    benchmark_protein("1VII (Villin headpiece, all-α)",  "data/1VII.pdb");
    benchmark_protein("2GB1 (Protein G B1, α+β)",        "data/2GB1.pdb");
    benchmark_protein("1ENH (Engrailed, 3-helix bundle)","data/1ENH.pdb");

    println!("══════════════════════════════════════════════════════════════════");
}
