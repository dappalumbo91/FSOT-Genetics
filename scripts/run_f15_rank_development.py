#!/usr/bin/env python3
"""Select a parameter-free F15 rank transform on development structures only."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import build_distogram  # noqa: E402
from run_f15_channel_ablation import conditioned_pearson, long_range_precision  # noqa: E402
from run_fsot_distogram_contact_eval import exp_distance_matrix  # noqa: E402
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402

OUTPUT = ROOT / "data" / "f15_rank_development.json"
PDB_DIR = ROOT / "data" / "pdb_samples"
SAMPLES = (("1UBQ", "A"), ("1CRN", "A"), ("1VII", "A"), ("2GB1", "A"), ("1ENH", "A"))


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> int:
    candidates: dict[str, list[dict[str, float | str]]] = {
        name: []
        for name in (
            "full_f15",
            "locality",
            "sheet_only",
            "locality_plus_sheet",
            "locality_times_one_plus_sheet",
        )
    }
    structure_hashes = []
    for pdb_id, chain in SAMPLES:
        path = PDB_DIR / f"{pdb_id}.pdb"
        sequence, coordinates = parse_pdb_ca(
            path.read_text(encoding="utf-8", errors="replace"), chain
        )
        full_f15, _, _, _, interface = build_distogram(sequence, collect_channels=True)
        channels = interface["channels"]
        distance = exp_distance_matrix(coordinates)
        scores = {
            "full_f15": full_f15,
            "locality": channels["locality"],
            "sheet_only": channels["sheet"],
            "locality_plus_sheet": channels["locality"] + channels["sheet"],
            "locality_times_one_plus_sheet": channels["locality"] * (1.0 + channels["sheet"]),
        }
        structure_hashes.append({"pdb_id": pdb_id, "chain": chain, "pdb_sha256": sha256_file(path)})
        for name, score in scores.items():
            candidates[name].append(
                {
                    "pdb_id": pdb_id,
                    "conditioned_pearson": conditioned_pearson(score, distance),
                    "long_range_top_L2_precision": long_range_precision(score, distance),
                }
            )

    summaries = []
    for name, rows in candidates.items():
        correlations = [float(row["conditioned_pearson"]) for row in rows]
        precisions = [float(row["long_range_top_L2_precision"]) for row in rows]
        summaries.append(
            {
                "name": name,
                "median_conditioned_pearson": float(np.median(correlations)),
                "worst_conditioned_pearson": float(min(correlations)),
                "median_long_range_top_L2_precision": float(np.median(precisions)),
                "per_structure": rows,
            }
        )

    eligible = [row for row in summaries if row["worst_conditioned_pearson"] > 0.0]
    eligible.sort(
        key=lambda row: (
            row["median_long_range_top_L2_precision"],
            row["median_conditioned_pearson"],
        ),
        reverse=True,
    )
    selected = eligible[0]
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "role": "development_only",
        "selection_rule": [
            "require positive exact-separation Pearson on every development structure",
            "maximize median separation>=24 Top-L/2 contact precision",
            "then maximize median exact-separation Pearson",
        ],
        "git_commit": git_commit(),
        "engine_sha256": sha256_file(ROOT / "scripts" / "fsot_structure_engine.py"),
        "development_structures": structure_hashes,
        "selected_candidate": selected,
        "candidates": summaries,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(
        f"selected={selected['name']} "
        f"median_r={selected['median_conditioned_pearson']:+.4f} "
        f"worst_r={selected['worst_conditioned_pearson']:+.4f} "
        f"LR24={selected['median_long_range_top_L2_precision']:.4f}"
    )
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())