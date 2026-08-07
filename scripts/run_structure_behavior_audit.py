#!/usr/bin/env python3
"""Measure how FSOT distance and reconstruction layers behave on real structures."""

from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from contact_rank import rank_long_range_contacts  # noqa: E402
from fsot_structure_engine import (  # noqa: E402
    build_distogram,
    canonicalize_l_amino_acid_handedness,
    predict_ca_coords,
    proximity_to_distance,
    refine_with_distogram,
)
from run_fsot_vs_alphafold_structure import kabsch_rmsd, parse_pdb_ca  # noqa: E402
from run_rcsb_holdout import fetch_pdb, sha256_bytes  # noqa: E402

MANIFEST = ROOT / "data" / "rcsb_oriented_backbone_manifest.json"
F19_EVAL = ROOT / "data" / "rcsb_oriented_backbone_eval.json"
OUTPUT = ROOT / "data" / "structure_behavior_audit.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_oriented_backbone"
SEPARATION_BINS = (
    ("local_1_2", 1, 2),
    ("mid_3_6", 3, 6),
    ("medium_7_11", 7, 11),
    ("long_12_23", 12, 23),
    ("long_24_plus", 24, 10**9),
)


def distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    return np.linalg.norm(
        coordinates[:, None, :] - coordinates[None, :, :], axis=2
    )


def radius_of_gyration(coordinates: np.ndarray) -> float:
    centered = coordinates - coordinates.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))))


def raw_classical_mds(distances: np.ndarray) -> np.ndarray:
    n = len(distances)
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (distances**2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1][:3]
    coordinates = eigenvectors[:, order] @ np.diag(
        np.sqrt(np.maximum(eigenvalues[order], 0.0))
    )
    return coordinates - coordinates.mean(axis=0)


def pdb_context(text: str, chain: str, target_ca: np.ndarray) -> dict:
    other_atoms: list[list[float]] = []
    hetero_atoms: list[list[float]] = []
    hetero_names: Counter[str] = Counter()
    chains: Counter[str] = Counter()
    for line in text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        except ValueError:
            continue
        element = line[76:78].strip()
        if element == "H":
            continue
        if line.startswith("ATOM"):
            atom_chain = line[21].strip()
            chains[atom_chain] += 1
            if atom_chain and atom_chain != chain:
                other_atoms.append(xyz)
        else:
            name = line[17:20].strip()
            if name not in ("HOH", "DOD"):
                hetero_names[name] += 1
                hetero_atoms.append(xyz)

    def residue_fraction_near(points: list[list[float]], cutoff: float) -> float:
        if not points:
            return 0.0
        context = np.asarray(points, dtype=np.float64)
        nearest = np.min(
            np.linalg.norm(target_ca[:, None, :] - context[None, :, :], axis=2),
            axis=1,
        )
        return float(np.mean(nearest < cutoff))

    return {
        "protein_chain_count": len([name for name in chains if name]),
        "protein_chains": dict(chains),
        "hetero_components": dict(hetero_names),
        "residue_fraction_near_other_chain_6A": residue_fraction_near(
            other_atoms, 6.0
        ),
        "residue_fraction_near_nonwater_hetero_6A": residue_fraction_near(
            hetero_atoms, 6.0
        ),
    }


def median_summary(values: list[float]) -> dict:
    array = np.asarray(values, dtype=np.float64)
    return {
        "median": float(np.median(array)),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
    }


def main() -> int:
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    f19_eval = json.loads(F19_EVAL.read_text(encoding="utf-8"))
    baseline_by_id = {
        (row["pdb_id"], row["chain"]): row["rmsd_A"]["baseline"]
        for row in f19_eval["results"]
    }
    rows = []
    pooled_target_residuals = {name: [] for name, _, _ in SEPARATION_BINS}
    pooled_realized_residuals = {name: [] for name, _, _ in SEPARATION_BINS}
    tight_hits = tight_lr24_hits = tight_count = 0
    loose_hits = loose_lr24_hits = loose_count = 0
    native_lr_count = covered_native_lr = 0
    false_constraint_residuals: list[float] = []
    diagnostic_rmsd: list[float] = []

    for entry in manifest["entries"]:
        pdb_id, chain = str(entry["pdb_id"]), str(entry["chain"])
        text, pdb_sha256 = fetch_pdb(pdb_id, CACHE)
        sequence, native = parse_pdb_ca(text, chain)
        proximity, props, regions, _, interface = build_distogram(sequence)
        target = proximity_to_distance(proximity, props, regions, interface)
        production = predict_ca_coords(sequence)["ca_coords"]
        native_distances = distance_matrix(native)
        realized_distances = distance_matrix(production)
        length = len(sequence)
        gate = int(interface["long_range_gate"])
        ranked = rank_long_range_contacts(sequence, proximity, regions, gate)
        tight = {(row["i"], row["j"]) for row in ranked[: max(length // 2, 1)]}
        loose = {
            (row["i"], row["j"])
            for row in ranked[max(length // 2, 1) : length]
        }
        native_lr = {
            (i, j)
            for i in range(length)
            for j in range(i + 24, length)
            if native_distances[i, j] < 8.0
        }
        native_gate = {
            (i, j)
            for i in range(length)
            for j in range(i + gate, length)
            if native_distances[i, j] < 8.0
        }
        tight_hits += len(tight & native_gate)
        tight_lr24_hits += len(tight & native_lr)
        tight_count += len(tight)
        loose_hits += len(loose & native_gate)
        loose_lr24_hits += len(loose & native_lr)
        loose_count += len(loose)
        native_lr_count += len(native_lr)
        covered_native_lr += len(native_lr & (tight | loose))
        false_constraints = (tight | loose) - native_gate
        false_constraint_residuals.extend(
            target[i, j] - native_distances[i, j] for i, j in false_constraints
        )

        bin_results = {}
        for name, minimum, maximum in SEPARATION_BINS:
            ii, jj = np.triu_indices(length, minimum)
            selected = jj - ii <= maximum
            ii, jj = ii[selected], jj[selected]
            target_residual = target[ii, jj] - native_distances[ii, jj]
            realized_residual = realized_distances[ii, jj] - native_distances[ii, jj]
            pooled_target_residuals[name].extend(target_residual)
            pooled_realized_residuals[name].extend(realized_residual)
            bin_results[name] = {
                "target_signed_median_A": float(np.median(target_residual)),
                "target_mae_A": float(np.mean(np.abs(target_residual))),
                "realized_signed_median_A": float(np.median(realized_residual)),
                "realized_mae_A": float(np.mean(np.abs(realized_residual))),
            }

        no_rank_target = proximity_to_distance(
            proximity,
            props,
            regions,
            interface,
            apply_ranked_contacts=False,
        )
        diagnostic = refine_with_distogram(
            raw_classical_mds(no_rank_target), no_rank_target, proximity
        )
        diagnostic, reflected = canonicalize_l_amino_acid_handedness(
            diagnostic, regions
        )
        diagnostic_value = kabsch_rmsd(diagnostic, native)
        diagnostic_rmsd.append(diagnostic_value)

        centered = native - native.mean(axis=0)
        shape_eigenvalues = np.sort(
            np.linalg.eigvalsh(centered.T @ centered / length)
        )[::-1]
        possible_lr = max((length - 24) * (length - 23) // 2, 1)
        context = pdb_context(text, chain, native)
        rows.append(
            {
                "pdb_id": pdb_id,
                "chain": chain,
                "length": length,
                "pdb_sha256": pdb_sha256,
                "baseline_rmsd_A": baseline_by_id[(pdb_id, chain)],
                "diagnostic_no_rank_local_rebond_rmsd_A": diagnostic_value,
                "diagnostic_chirality_reflected": reflected,
                "native_radius_of_gyration_A": radius_of_gyration(native),
                "native_shape_axis_variance_ratio": float(
                    shape_eigenvalues[0] / max(shape_eigenvalues[-1], 1e-12)
                ),
                "native_long_range_contact_density": len(native_lr) / possible_lr,
                "tight_ranked_contact_precision": len(tight & native_gate)
                / max(len(tight), 1),
                "tight_ranked_lr24_precision": len(tight & native_lr)
                / max(len(tight), 1),
                "top_L_native_long_range_recall": len(native_lr & (tight | loose))
                / max(len(native_lr), 1),
                "context": context,
                "residual_by_separation": bin_results,
            }
        )

    aggregate = {
        "tight_ranked_contact_precision": tight_hits / tight_count,
        "tight_ranked_lr24_precision": tight_lr24_hits / tight_count,
        "loose_ranked_contact_precision": loose_hits / loose_count,
        "loose_ranked_lr24_precision": loose_lr24_hits / loose_count,
        "top_L_native_long_range_recall": covered_native_lr / native_lr_count,
        "false_constraint_signed_median_A": float(
            np.median(false_constraint_residuals)
        ),
        "baseline_rmsd_A": median_summary(
            [row["baseline_rmsd_A"] for row in rows]
        ),
        "diagnostic_no_rank_local_rebond_rmsd_A": median_summary(diagnostic_rmsd),
        "native_radius_of_gyration_A": median_summary(
            [row["native_radius_of_gyration_A"] for row in rows]
        ),
        "context_counts": {
            "with_other_protein_chain": sum(
                row["context"]["protein_chain_count"] > 1 for row in rows
            ),
            "with_nonwater_hetero_component": sum(
                bool(row["context"]["hetero_components"]) for row in rows
            ),
            "with_either_context": sum(
                row["context"]["protein_chain_count"] > 1
                or bool(row["context"]["hetero_components"])
                for row in rows
            ),
        },
        "residual_by_separation": {
            name: {
                "target_signed_median_A": float(np.median(pooled_target_residuals[name])),
                "target_absolute_median_A": float(
                    np.median(np.abs(pooled_target_residuals[name]))
                ),
                "realized_signed_median_A": float(
                    np.median(pooled_realized_residuals[name])
                ),
                "realized_absolute_median_A": float(
                    np.median(np.abs(pooled_realized_residuals[name]))
                ),
            }
            for name, _, _ in SEPARATION_BINS
        },
    }
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "Observed behavior audit; no production formula selection",
        "manifest_sha256": sha256_bytes(manifest_bytes),
        "aggregate": aggregate,
        "results": rows,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(aggregate, indent=2))
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())