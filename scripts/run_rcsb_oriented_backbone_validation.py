#!/usr/bin/env python3
"""Run the preregistered F19 oriented-backbone coordinate holdout once."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import predict_ca_coords  # noqa: E402
from run_fsot_vs_alphafold_structure import kabsch_rmsd, parse_pdb_ca  # noqa: E402
from run_rcsb_holdout import (  # noqa: E402
    bootstrap_median,
    fetch_pdb,
    git_commit,
    sha256_bytes,
)

MANIFEST = ROOT / "data" / "rcsb_oriented_backbone_manifest.json"
OUTPUT = ROOT / "data" / "rcsb_oriented_backbone_eval.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_oriented_backbone"
MODELS = ("baseline", "f12c", "f12c_f19")


def distance_matrix(coordinates: np.ndarray) -> np.ndarray:
    return np.linalg.norm(
        coordinates[:, None, :] - coordinates[None, :, :], axis=2
    )


def evaluate_entry(entry: dict) -> dict:
    pdb_id, chain = str(entry["pdb_id"]), str(entry["chain"])
    text, pdb_sha256 = fetch_pdb(pdb_id, CACHE)
    sequence, native = parse_pdb_ca(text, chain)
    if len(sequence) != int(entry["expected_length"]):
        raise ValueError(f"{pdb_id}:{chain} length changed")

    baseline = predict_ca_coords(sequence)
    f12c = predict_ca_coords(sequence, cooperative_regions=True)
    f12c_f19 = predict_ca_coords(
        sequence,
        cooperative_regions=True,
        canonicalize_chirality=True,
    )
    predictions = {
        "baseline": baseline,
        "f12c": f12c,
        "f12c_f19": f12c_f19,
    }
    distance_delta = float(
        np.max(
            np.abs(
                distance_matrix(f12c["ca_coords"])
                - distance_matrix(f12c_f19["ca_coords"])
            )
        )
    )
    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "length": len(sequence),
        "release_date": entry["release_date"],
        "resolution_A": entry["resolution_A"],
        "pdb_sha256": pdb_sha256,
        "rmsd_A": {
            name: kabsch_rmsd(result["ca_coords"], native)
            for name, result in predictions.items()
        },
        "f19_reflected": f12c_f19["chirality_reflected"],
        "f19_pair_distance_max_delta_A": distance_delta,
    }


def main() -> int:
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    results = []
    for entry in manifest["entries"]:
        row = evaluate_entry(entry)
        results.append(row)
        print(
            f"{row['pdb_id']}:{row['chain']} n={row['length']} "
            f"F12c={row['rmsd_A']['f12c']:.3f} "
            f"F19={row['rmsd_A']['f12c_f19']:.3f} "
            f"reflected={row['f19_reflected']}"
        )

    summary = {
        "rmsd_A": {
            model: bootstrap_median([row["rmsd_A"][model] for row in results])
            for model in MODELS
        },
        "f19_reflected_count": sum(row["f19_reflected"] for row in results),
        "maximum_pair_distance_delta_A": max(
            row["f19_pair_distance_max_delta_A"] for row in results
        ),
    }
    summary["success_gate_passed"] = bool(
        summary["rmsd_A"]["f12c_f19"]["median"]
        <= summary["rmsd_A"]["f12c"]["median"]
        and summary["maximum_pair_distance_delta_A"] == 0.0
    )
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy": manifest["policy"],
        "candidate": manifest["candidate"],
        "manifest_sha256": sha256_bytes(manifest_bytes),
        "git_commit": git_commit(),
        "benchmark_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "engine_sha256": sha256_bytes(
            (ROOT / "scripts" / "fsot_structure_engine.py").read_bytes()
        ),
        "authority_sha256": sha256_bytes(
            (ROOT / "vendor" / "fsot_compute.py").read_bytes()
        ),
        "summary": summary,
        "results": results,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())