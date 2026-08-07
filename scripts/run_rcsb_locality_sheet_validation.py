#!/usr/bin/env python3
"""Validate the preregistered locality-plus-sheet F15 rank candidate once."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import build_distogram  # noqa: E402
from run_f15_channel_ablation import conditioned_pearson, long_range_precision  # noqa: E402
from run_rcsb_holdout import bootstrap_median, fetch_pdb, git_commit, sha256_bytes  # noqa: E402
from run_fsot_distogram_contact_eval import exp_distance_matrix  # noqa: E402
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402

MANIFEST = ROOT / "data" / "rcsb_locality_sheet_manifest.json"
OUTPUT = ROOT / "data" / "rcsb_locality_sheet_eval.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_locality_sheet"
SCORES = ("full_f15", "locality", "sheet_only", "locality_plus_sheet")


def evaluate_entry(entry: dict) -> dict:
    pdb_id, chain = str(entry["pdb_id"]), str(entry["chain"])
    text, pdb_sha256 = fetch_pdb(pdb_id, CACHE)
    sequence, coordinates = parse_pdb_ca(text, chain)
    if len(sequence) != int(entry["expected_length"]):
        raise ValueError(f"{pdb_id}:{chain} length changed")
    full_f15, _, _, _, interface = build_distogram(sequence, collect_channels=True)
    channels = interface["channels"]
    distance = exp_distance_matrix(coordinates)
    matrices = {
        "full_f15": full_f15,
        "locality": channels["locality"],
        "sheet_only": channels["sheet"],
        "locality_plus_sheet": channels["locality"] + channels["sheet"],
    }
    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "length": len(sequence),
        "release_date": entry["release_date"],
        "resolution_A": entry["resolution_A"],
        "pdb_sha256": pdb_sha256,
        "scores": {
            name: {
                "conditioned_pearson": conditioned_pearson(score, distance),
                "long_range_top_L2_precision": long_range_precision(score, distance),
            }
            for name, score in matrices.items()
        },
    }


def main() -> int:
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    results = []
    for entry in manifest["entries"]:
        row = evaluate_entry(entry)
        results.append(row)
        candidate = row["scores"]["locality_plus_sheet"]
        print(
            f"{row['pdb_id']}:{row['chain']} n={row['length']} "
            f"r={candidate['conditioned_pearson']:+.3f} "
            f"LR24={candidate['long_range_top_L2_precision']:.3f}"
        )
    summary = {
        name: {
            metric: bootstrap_median([row["scores"][name][metric] for row in results])
            for metric in ("conditioned_pearson", "long_range_top_L2_precision")
        }
        for name in SCORES
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