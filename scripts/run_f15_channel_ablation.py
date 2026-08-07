#!/usr/bin/env python3
"""Ablate additive F15 channels on the five declared development structures.

This script may select a candidate for a future holdout. It must never read the
frozen RCSB holdout or be used to revise a candidate after holdout evaluation.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import build_distogram  # noqa: E402
from run_fsot_distogram_contact_eval import exp_distance_matrix, pearson  # noqa: E402
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402

OUTPUT = ROOT / "data" / "f15_channel_ablation.json"
PDB_DIR = ROOT / "data" / "pdb_samples"
CHANNELS = ("chemistry", "helix", "sheet", "region")
SAMPLES = (
    ("1UBQ", "A"),
    ("1CRN", "A"),
    ("1VII", "A"),
    ("2GB1", "A"),
    ("1ENH", "A"),
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def conditioned_pearson(score: np.ndarray, distance: np.ndarray, minimum_separation: int = 6) -> float:
    score_residuals: list[float] = []
    truth_residuals: list[float] = []
    for separation in range(minimum_separation, len(score)):
        predicted = np.array([score[i, i + separation] for i in range(len(score) - separation)])
        observed = np.array(
            [1.0 / distance[i, i + separation] for i in range(len(score) - separation)]
        )
        if len(predicted) < 2:
            continue
        score_residuals.extend(predicted - predicted.mean())
        truth_residuals.extend(observed - observed.mean())
    return pearson(np.asarray(score_residuals), np.asarray(truth_residuals))


def long_range_precision(
    score: np.ndarray,
    distance: np.ndarray,
    minimum_separation: int = 24,
    contact_cutoff_A: float = 8.0,
) -> float:
    length = len(score)
    pairs = [
        (float(score[i, j]), float(distance[i, j]))
        for i in range(length)
        for j in range(i + minimum_separation, length)
    ]
    count = min(max(length // 2, 1), len(pairs))
    top = sorted(pairs, key=lambda pair: pair[0], reverse=True)[:count]
    return sum(distance_A < contact_cutoff_A for _, distance_A in top) / max(count, 1)


def load_development_data() -> list[dict[str, Any]]:
    datasets = []
    for pdb_id, chain in SAMPLES:
        path = PDB_DIR / f"{pdb_id}.pdb"
        sequence, coordinates = parse_pdb_ca(
            path.read_text(encoding="utf-8", errors="replace"), chain
        )
        proximity, _, _, _, interface = build_distogram(sequence, collect_channels=True)
        datasets.append(
            {
                "pdb_id": pdb_id,
                "chain": chain,
                "pdb_sha256": sha256_file(path),
                "distance": exp_distance_matrix(coordinates),
                "proximity": proximity,
                "channels": interface["channels"],
            }
        )
    return datasets


def main() -> int:
    datasets = load_development_data()
    candidates = []
    for size in range(1, len(CHANNELS) + 1):
        for subset in itertools.combinations(CHANNELS, size):
            per_structure = []
            for dataset in datasets:
                distance = dataset["distance"]
                score = sum(
                    (dataset["channels"][channel] for channel in subset),
                    start=np.zeros_like(distance),
                )
                per_structure.append(
                    {
                        "pdb_id": dataset["pdb_id"],
                        "conditioned_pearson": conditioned_pearson(score, distance),
                        "long_range_top_L2_precision": long_range_precision(score, distance),
                    }
                )
            correlations = [row["conditioned_pearson"] for row in per_structure]
            precisions = [row["long_range_top_L2_precision"] for row in per_structure]
            candidates.append(
                {
                    "channels": list(subset),
                    "median_conditioned_pearson": float(np.median(correlations)),
                    "worst_conditioned_pearson": float(min(correlations)),
                    "median_long_range_top_L2_precision": float(np.median(precisions)),
                    "per_structure": per_structure,
                }
            )

    candidates.sort(
        key=lambda row: (
            row["median_conditioned_pearson"],
            row["worst_conditioned_pearson"],
            row["median_long_range_top_L2_precision"],
        ),
        reverse=True,
    )
    selected = candidates[0]
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "role": "development_only",
        "selection_rule": [
            "maximize median exact-separation Pearson",
            "then maximize worst-structure exact-separation Pearson",
            "then maximize median separation>=24 Top-L/2 contact precision",
        ],
        "git_commit": git_commit(),
        "engine_sha256": sha256_file(ROOT / "scripts" / "fsot_structure_engine.py"),
        "development_structures": [
            {
                "pdb_id": dataset["pdb_id"],
                "chain": dataset["chain"],
                "pdb_sha256": dataset["pdb_sha256"],
            }
            for dataset in datasets
        ],
        "selected_candidate": selected,
        "candidates": candidates,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(
        f"selected={'+'.join(selected['channels'])} "
        f"median_r={selected['median_conditioned_pearson']:+.4f} "
        f"worst_r={selected['worst_conditioned_pearson']:+.4f} "
        f"LR24={selected['median_long_range_top_L2_precision']:.4f}"
    )
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())