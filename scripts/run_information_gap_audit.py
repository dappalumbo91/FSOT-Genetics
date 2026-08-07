#!/usr/bin/env python3
"""Audit information loss between protein sequence and FSOT C-alpha coordinates."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import (  # noqa: E402
    SsPropensity,
    build_distogram,
    classical_mds,
    detect_regions_cooperative,
    predict_ca_coords,
    proximity_to_distance,
    stress,
)
from run_fsot_vs_alphafold_structure import (  # noqa: E402
    kabsch_rmsd,
    parse_pdb_ca,
)

SAMPLES = ("1UBQ", "1CRN", "1VII", "2GB1", "1ENH")
OUT = ROOT / "data" / "information_gap_audit.json"


def distance_matrix(coords: np.ndarray) -> np.ndarray:
    return np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)


def negative_eigen_mass(distances: np.ndarray) -> float:
    n = distances.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (distances**2) @ centering
    eigenvalues = np.linalg.eigvalsh(gram)
    total = float(np.abs(eigenvalues).sum())
    return float(np.abs(eigenvalues[eigenvalues < 0]).sum() / total) if total else 0.0


def region_edge_counts(
    sequence: str,
    native: np.ndarray,
    chemistry: np.ndarray,
    helix: np.ndarray,
    sheet: np.ndarray,
) -> dict[str, int | float]:
    props = [SsPropensity.from_expanded_amino_acid(aa) for aa in sequence]
    regions = detect_regions_cooperative(props)
    predicted = 0
    native_edges = 0
    scored_edges: list[float] = []
    for left_index, left in enumerate(regions):
        for right in regions[left_index + 1 :]:
            if left.kind != right.kind:
                continue
            pairs = [
                (i, j)
                for i in range(left.start, left.end + 1)
                for j in range(right.start, right.end + 1)
                if abs(i - j) >= 7
            ]
            if not pairs:
                continue
            predicted += 1
            channel = sheet if left.kind == "E" else helix
            scored_edges.append(
                float(np.mean([chemistry[i, j] + channel[i, j] for i, j in pairs]))
            )
            if min(float(np.linalg.norm(native[i] - native[j])) for i, j in pairs) < 8.0:
                native_edges += 1
    return {
        "candidate_edges": predicted,
        "native_edges": native_edges,
        "precision": native_edges / predicted if predicted else 0.0,
        "mean_edge_score": float(np.mean(scored_edges)) if scored_edges else 0.0,
    }


def audit_sample(pdb_id: str) -> dict:
    pdb_path = ROOT / "data" / "pdb_samples" / f"{pdb_id}.pdb"
    sequence, native = parse_pdb_ca(
        pdb_path.read_text(encoding="utf-8", errors="replace"), "A"
    )
    matrix, props, regions, _, interface = build_distogram(
        sequence, collect_channels=True
    )
    predicted_distances = proximity_to_distance(matrix, props, regions, interface)
    native_distances = distance_matrix(native)
    predicted_coords = classical_mds(predicted_distances)
    oracle_coords = classical_mds(native_distances)
    mirror = np.array([1.0, 1.0, -1.0])
    upper = np.triu_indices(len(sequence), 1)
    channels = interface["channels"]
    direct_stress = stress(
        predicted_coords,
        predicted_distances,
        matrix,
        topology_weight=True,
    )
    mirror_stress = stress(
        predicted_coords * mirror,
        predicted_distances,
        matrix,
        topology_weight=True,
    )
    baseline_fold = predict_ca_coords(sequence)
    cooperative_fold = predict_ca_coords(sequence, cooperative_regions=True)
    oriented_fold = predict_ca_coords(
        sequence,
        cooperative_regions=True,
        canonicalize_chirality=True,
    )
    cooperative_distances = distance_matrix(cooperative_fold["ca_coords"])
    oriented_distances = distance_matrix(oriented_fold["ca_coords"])
    return {
        "pdb_id": pdb_id,
        "length": len(sequence),
        "oracle_mds_rmsd_A": kabsch_rmsd(oracle_coords, native),
        "oracle_mds_best_enantiomer_rmsd_A": min(
            kabsch_rmsd(oracle_coords, native),
            kabsch_rmsd(oracle_coords * mirror, native),
        ),
        "predicted_mds_rmsd_A": kabsch_rmsd(predicted_coords, native),
        "baseline_fold_rmsd_A": kabsch_rmsd(baseline_fold["ca_coords"], native),
        "cooperative_f12c_fold_rmsd_A": kabsch_rmsd(
            cooperative_fold["ca_coords"], native
        ),
        "oriented_f12c_fold_rmsd_A": kabsch_rmsd(
            oriented_fold["ca_coords"], native
        ),
        "chirality_reflected": oriented_fold["chirality_reflected"],
        "chirality_pair_distance_max_delta_A": float(
            np.max(np.abs(cooperative_distances - oriented_distances))
        ),
        "predicted_distance_mae_A": float(
            np.mean(np.abs(predicted_distances[upper] - native_distances[upper]))
        ),
        "predicted_distance_pearson": float(
            np.corrcoef(predicted_distances[upper], native_distances[upper])[0, 1]
        ),
        "negative_eigen_mass_fraction": negative_eigen_mass(predicted_distances),
        "mirror_stress_delta": abs(direct_stress - mirror_stress),
        "cooperative_region_edges": region_edge_counts(
            sequence,
            native,
            channels["chemistry"],
            channels["helix"],
            channels["sheet"],
        ),
    }


def main() -> None:
    samples = [audit_sample(pdb_id) for pdb_id in SAMPLES]
    aggregate = {
        "median_predicted_mds_rmsd_A": float(
            np.median([sample["predicted_mds_rmsd_A"] for sample in samples])
        ),
        "median_baseline_fold_rmsd_A": float(
            np.median([sample["baseline_fold_rmsd_A"] for sample in samples])
        ),
        "median_cooperative_f12c_fold_rmsd_A": float(
            np.median([sample["cooperative_f12c_fold_rmsd_A"] for sample in samples])
        ),
        "median_oriented_f12c_fold_rmsd_A": float(
            np.median([sample["oriented_f12c_fold_rmsd_A"] for sample in samples])
        ),
        "max_chirality_pair_distance_delta_A": max(
            sample["chirality_pair_distance_max_delta_A"] for sample in samples
        ),
        "median_predicted_distance_mae_A": float(
            np.median([sample["predicted_distance_mae_A"] for sample in samples])
        ),
        "median_predicted_distance_pearson": float(
            np.median([sample["predicted_distance_pearson"] for sample in samples])
        ),
        "median_negative_eigen_mass_fraction": float(
            np.median([sample["negative_eigen_mass_fraction"] for sample in samples])
        ),
        "total_region_candidate_edges": sum(
            sample["cooperative_region_edges"]["candidate_edges"] for sample in samples
        ),
        "total_region_native_edges": sum(
            sample["cooperative_region_edges"]["native_edges"] for sample in samples
        ),
        "all_mirror_stress_deltas_zero": all(
            sample["mirror_stress_delta"] == 0.0 for sample in samples
        ),
    }
    result = {
        "purpose": "Localize missing information before further FSOT tuning",
        "production_fold_modified": False,
        "samples": samples,
        "aggregate": aggregate,
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(aggregate, indent=2))


if __name__ == "__main__":
    main()