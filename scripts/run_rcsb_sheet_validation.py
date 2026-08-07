#!/usr/bin/env python3
"""Validate the preregistered F15 sheet-only candidate exactly once."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import build_distogram  # noqa: E402
from run_f15_channel_ablation import conditioned_pearson, long_range_precision  # noqa: E402
from run_rcsb_holdout import bootstrap_median, fetch_pdb, git_commit, sha256_bytes  # noqa: E402
from run_fsot_distogram_contact_eval import exp_distance_matrix  # noqa: E402
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402

MANIFEST = ROOT / "data" / "rcsb_sheet_validation_manifest.json"
OUTPUT = ROOT / "data" / "rcsb_sheet_validation_eval.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_sheet_validation"
SCORES = ("full_f15", "sheet_only", "locality")


def evaluate_entry(entry: dict[str, Any]) -> dict[str, Any]:
    pdb_id = str(entry["pdb_id"])
    chain = str(entry["chain"])
    text, pdb_sha256 = fetch_pdb(pdb_id, CACHE)
    sequence, coordinates = parse_pdb_ca(text, chain)
    expected_length = int(entry["expected_length"])
    if len(sequence) != expected_length:
        raise ValueError(f"{pdb_id}:{chain} length changed: {len(sequence)} != {expected_length}")

    full_f15, _, _, _, interface = build_distogram(sequence, collect_channels=True)
    distance = exp_distance_matrix(coordinates)
    score_matrices = {
        "full_f15": full_f15,
        "sheet_only": interface["channels"]["sheet"],
        "locality": interface["channels"]["locality"],
    }
    scores = {
        name: {
            "conditioned_pearson": conditioned_pearson(score, distance),
            "long_range_top_L2_precision": long_range_precision(score, distance),
        }
        for name, score in score_matrices.items()
    }
    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "length": len(sequence),
        "release_date": entry["release_date"],
        "resolution_A": entry["resolution_A"],
        "pdb_sha256": pdb_sha256,
        "scores": scores,
    }


def main() -> int:
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    results = []
    for entry in manifest["entries"]:
        result = evaluate_entry(entry)
        results.append(result)
        sheet = result["scores"]["sheet_only"]
        print(
            f"{result['pdb_id']}:{result['chain']} n={result['length']} "
            f"sheet_r={sheet['conditioned_pearson']:+.3f} "
            f"sheet_LR24={sheet['long_range_top_L2_precision']:.3f}"
        )

    summary = {
        score: {
            "conditioned_pearson": bootstrap_median(
                [row["scores"][score]["conditioned_pearson"] for row in results]
            ),
            "long_range_top_L2_precision": bootstrap_median(
                [row["scores"][score]["long_range_top_L2_precision"] for row in results]
            ),
        }
        for score in SCORES
    }
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy": manifest["policy"],
        "candidate": manifest["candidate"],
        "manifest_sha256": sha256_bytes(manifest_bytes),
        "git_commit": git_commit(),
        "benchmark_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "engine_sha256": sha256_bytes((ROOT / "scripts" / "fsot_structure_engine.py").read_bytes()),
        "authority_sha256": sha256_bytes((ROOT / "vendor" / "fsot_compute.py").read_bytes()),
        "summary": summary,
        "results": results,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())