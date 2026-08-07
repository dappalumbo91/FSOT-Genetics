#!/usr/bin/env python3
"""Measure the intrinsic dimensionality of the FSOT proximity object.

Central finding: an FSOT target matrix is not a 3-D spatial object. Its
participation dimension tracks the pinned D_eff ladder and it carries genuine
non-Euclidean mass, so a forced 3-D MDS embedding discards most of the signal.
This audit reproduces that behavior; it does not change production defaults.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import PI, build_distogram, proximity_to_distance  # noqa: E402
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402
from run_rcsb_holdout import fetch_pdb, sha256_bytes  # noqa: E402

MANIFEST = ROOT / "data" / "rcsb_oriented_backbone_manifest.json"
OUTPUT = ROOT / "data" / "dimensionality_audit.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_oriented_backbone"

# Pinned effective-dimension ladder (FSOTGenetics/ChemLink.lean, full_scalar_law.py).
# These are effective/fractal dimensions, never a count of spatial axes.
CHEM_LINK_D_EFF = {
    "backbone": 8,
    "disulfide": 7,
    "salt_bridge": 9,
    "hydrophobic_pack": 14,
    "hbond_secondary": 8,
    "molecular_sidechain": 9,
    "tertiary_biochem": 13,
}
BASE_LAW_D_EFF = 25
EMBED_DIMENSIONS = (3, 9, 25)


def gram(distances: np.ndarray) -> np.ndarray:
    n = len(distances)
    centering = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * centering @ (distances**2) @ centering


def positive_spectrum(distances: np.ndarray) -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh(gram(distances))
    return np.clip(np.sort(eigenvalues)[::-1], 0.0, None)


def participation_dimension(positive: np.ndarray) -> float:
    total = positive.sum()
    if total <= 0.0:
        return 0.0
    return float(total * total / np.sum(positive * positive))


def top_k_variance_fraction(positive: np.ndarray, k: int) -> float:
    total = positive.sum()
    return float(positive[:k].sum() / total) if total > 0.0 else 0.0


def negative_eigenvalue_mass(distances: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(gram(distances))
    denom = np.abs(eigenvalues).sum()
    return float(-eigenvalues[eigenvalues < 0].sum() / denom) if denom > 0 else 0.0


def embed(distances: np.ndarray, dim: int) -> np.ndarray:
    eigenvalues, vectors = np.linalg.eigh(gram(distances))
    order = np.argsort(eigenvalues)[::-1][:dim]
    return vectors[:, order] @ np.diag(np.sqrt(np.clip(eigenvalues[order], 0.0, None)))


def top_contacts(coordinates: np.ndarray, gate: int, count: int) -> set[tuple[int, int]]:
    n = len(coordinates)
    distances = np.linalg.norm(
        coordinates[:, None, :] - coordinates[None, :, :], axis=2
    )
    pairs = sorted(
        (distances[i, j], i, j)
        for i in range(n)
        for j in range(i + gate, n)
    )
    return {(i, j) for _, i, j in pairs[:count]}


def distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    return np.linalg.norm(
        coordinates[:, None, :] - coordinates[None, :, :], axis=2
    )


def radius_of_gyration(coordinates: np.ndarray) -> float:
    centered = coordinates - coordinates.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum(centered * centered, axis=1))))


def median(values: list[float]) -> float:
    return float(np.median(np.asarray(values, dtype=np.float64)))


def main() -> int:
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    rows = []
    native_part, fsot_part = [], []
    native_top3, fsot_top3 = [], []
    neg_mass = []
    rg_relative_error = []
    recall_by_dim: dict[int, list[float]] = {d: [] for d in EMBED_DIMENSIONS}

    for entry in manifest["entries"]:
        pdb_id, chain = str(entry["pdb_id"]), str(entry["chain"])
        text, _ = fetch_pdb(pdb_id, CACHE)
        sequence, native = parse_pdb_ca(text, chain)
        n = len(sequence)
        native_distances = distance_matrix(native)
        proximity, props, regions, _, interface = build_distogram(sequence)
        target = proximity_to_distance(proximity, props, regions, interface)

        native_positive = positive_spectrum(native_distances)
        fsot_positive = positive_spectrum(target)
        native_pd = participation_dimension(native_positive)
        fsot_pd = participation_dimension(fsot_positive)
        native_v3 = top_k_variance_fraction(native_positive, 3)
        fsot_v3 = top_k_variance_fraction(fsot_positive, 3)
        neg = negative_eigenvalue_mass(target)

        target_rg = PI * (float(n) ** (1.0 / PI))
        native_rg = radius_of_gyration(native)
        rg_rel = abs(target_rg - native_rg) / native_rg

        gate = int(interface["long_range_gate"])
        budget = max(n // 2, 1)
        native_contacts = {
            (i, j)
            for i in range(n)
            for j in range(i + gate, n)
            if native_distances[i, j] < 8.0
        }
        chain_recall = {}
        for dim in EMBED_DIMENSIONS:
            predicted = top_contacts(embed(target, dim), gate, budget)
            recall = (
                len(predicted & native_contacts) / len(native_contacts)
                if native_contacts
                else 0.0
            )
            chain_recall[dim] = recall
            if native_contacts:
                recall_by_dim[dim].append(recall)

        native_part.append(native_pd)
        fsot_part.append(fsot_pd)
        native_top3.append(native_v3)
        fsot_top3.append(fsot_v3)
        neg_mass.append(neg)
        rg_relative_error.append(rg_rel)
        rows.append(
            {
                "pdb_id": pdb_id,
                "chain": chain,
                "length": n,
                "native_participation_dimension": native_pd,
                "fsot_participation_dimension": fsot_pd,
                "native_top3_variance_fraction": native_v3,
                "fsot_top3_variance_fraction": fsot_v3,
                "fsot_negative_eigenvalue_mass": neg,
                "fsot_target_rg_A": target_rg,
                "native_rg_A": native_rg,
                "rg_relative_error": rg_rel,
                "native_contact_recall_by_dim": chain_recall,
            }
        )

    aggregate = {
        "chem_link_d_eff_ladder": CHEM_LINK_D_EFF,
        "base_law_d_eff": BASE_LAW_D_EFF,
        "native_participation_dimension_median": median(native_part),
        "fsot_participation_dimension_median": median(fsot_part),
        "native_top3_variance_fraction_median": median(native_top3),
        "fsot_top3_variance_fraction_median": median(fsot_top3),
        "fsot_negative_eigenvalue_mass_median": median(neg_mass),
        "rg_relative_error_median": median(rg_relative_error),
        "native_contact_recall_median_by_dim": {
            str(dim): median(recall_by_dim[dim]) for dim in EMBED_DIMENSIONS
        },
        "finding": (
            "FSOT proximity participation dimension "
            f"{median(fsot_part):.2f} tracks the D_eff ladder "
            f"[{min(CHEM_LINK_D_EFF.values())}, {max(CHEM_LINK_D_EFF.values())}]; "
            f"native structures sit at {median(native_part):.2f}. Forcing the "
            "FSOT object into 3 axes discards most variance and ignores "
            f"{median(neg_mass):.3f} non-Euclidean mass."
        ),
    }
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "Intrinsic-dimensionality audit; no production formula selection",
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
