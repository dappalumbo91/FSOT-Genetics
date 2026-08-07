#!/usr/bin/env python3
"""Evaluate frozen FSOT contact predictions against an RCSB PDB holdout.

The holdout is for falsification, not routing or formula selection. Reported
metrics explicitly separate the sequence-separation prior from FSOT residual
signal.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import PI, build_distogram  # noqa: E402
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402

MANIFEST = ROOT / "data" / "rcsb_holdout_manifest.json"
OUTPUT = ROOT / "data" / "rcsb_holdout_eval.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_holdout"
CONTACT_CUTOFF_A = 8.0
SEPARATION_THRESHOLDS = (6, 12, 24)


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x_centered = x - x.mean()
    y_centered = y - y.mean()
    denominator = np.sqrt((x_centered * x_centered).sum() * (y_centered * y_centered).sum())
    if denominator < 1e-15:
        return 0.0
    return float((x_centered * y_centered).sum() / denominator)


def bootstrap_median(values: list[float], samples: int = 10_000) -> dict[str, Any]:
    array = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(0)
    draws = rng.choice(array, size=(samples, len(array)), replace=True)
    medians = np.median(draws, axis=1)
    return {
        "median": float(np.median(array)),
        "bootstrap_95_ci": [
            float(np.quantile(medians, 0.025)),
            float(np.quantile(medians, 0.975)),
        ],
        "values": [float(value) for value in values],
    }


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def fetch_pdb(pdb_id: str, cache: Path) -> tuple[str, str]:
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{pdb_id}.pdb"
    if path.is_file():
        content = path.read_bytes()
    else:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        request = urllib.request.Request(url, headers={"User-Agent": "FSOT-Genetics-holdout/1"})
        with urllib.request.urlopen(request, timeout=90) as response:
            content = response.read()
        path.write_bytes(content)
    return content.decode("utf-8", errors="replace"), sha256_bytes(content)


def exact_separation_pearson(records: np.ndarray) -> float:
    score_residuals: list[float] = []
    truth_residuals: list[float] = []
    for separation in np.unique(records[:, 0]):
        group = records[records[:, 0] == separation]
        if len(group) < 2:
            continue
        score_residuals.extend(group[:, 1] - group[:, 1].mean())
        truth_residuals.extend(group[:, 5] - group[:, 5].mean())
    return pearson(np.asarray(score_residuals), np.asarray(truth_residuals))


def top_precision(records: np.ndarray, score_column: int, length: int, minimum_separation: int) -> float:
    eligible = records[records[:, 0] >= minimum_separation]
    count = min(max(length // 2, 1), len(eligible))
    order = np.argsort(-eligible[:, score_column], kind="stable")[:count]
    return float(np.mean(eligible[order, 4] < CONTACT_CUTOFF_A))


def evaluate_entry(entry: dict[str, Any], cache: Path) -> dict[str, Any]:
    pdb_id = str(entry["pdb_id"])
    chain = str(entry["chain"])
    text, pdb_sha256 = fetch_pdb(pdb_id, cache)
    sequence, coordinates = parse_pdb_ca(text, chain)
    expected_length = int(entry["expected_length"])
    if len(sequence) != expected_length:
        raise ValueError(f"{pdb_id}:{chain} length changed: {len(sequence)} != {expected_length}")

    proximity, _, _, _, interface = build_distogram(sequence)
    records: list[tuple[float, float, float, float, float, float]] = []
    for i in range(len(sequence)):
        for j in range(i + 1, len(sequence)):
            separation = j - i
            distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            locality = separation ** (-1.0 / PI)
            score = float(proximity[i, j])
            records.append(
                (separation, score, score - locality, locality, distance, 1.0 / distance)
            )
    array = np.asarray(records, dtype=np.float64)

    contact_precision: dict[str, Any] = {}
    for threshold in SEPARATION_THRESHOLDS:
        contact_precision[str(threshold)] = {
            "fsot_total": top_precision(array, 1, len(sequence), threshold),
            "fsot_residual": top_precision(array, 2, len(sequence), threshold),
            "separation_only": top_precision(array, 3, len(sequence), threshold),
        }

    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "length": len(sequence),
        "release_date": entry["release_date"],
        "resolution_A": entry["resolution_A"],
        "pdb_sha256": pdb_sha256,
        "routing": interface["routing"],
        "all_pair_pearson": {
            "fsot": pearson(array[:, 1], array[:, 5]),
            "separation_only": pearson(array[:, 3], array[:, 5]),
        },
        "exact_separation_pearson_fsot": exact_separation_pearson(array),
        "top_L2_contact_precision": contact_precision,
    }


def main() -> int:
    manifest_bytes = MANIFEST.read_bytes()
    manifest = json.loads(manifest_bytes)
    results = []
    for entry in manifest["entries"]:
        result = evaluate_entry(entry, CACHE)
        results.append(result)
        print(
            f"{result['pdb_id']}:{result['chain']} n={result['length']} "
            f"r_cond={result['exact_separation_pearson_fsot']:+.3f}"
        )

    summary: dict[str, Any] = {
        "n_structures": len(results),
        "all_pair_pearson_fsot": bootstrap_median(
            [row["all_pair_pearson"]["fsot"] for row in results]
        ),
        "all_pair_pearson_separation_only": bootstrap_median(
            [row["all_pair_pearson"]["separation_only"] for row in results]
        ),
        "exact_separation_pearson_fsot": bootstrap_median(
            [row["exact_separation_pearson_fsot"] for row in results]
        ),
        "top_L2_contact_precision": {},
    }
    for threshold in SEPARATION_THRESHOLDS:
        key = str(threshold)
        summary["top_L2_contact_precision"][key] = {
            score: bootstrap_median(
                [row["top_L2_contact_precision"][key][score] for row in results]
            )
            for score in ("fsot_total", "fsot_residual", "separation_only")
        }

    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy": manifest["policy"],
        "source": manifest["source"],
        "manifest_sha256": sha256_bytes(manifest_bytes),
        "git_commit": git_commit(),
        "benchmark_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "engine_sha256": sha256_bytes((ROOT / "scripts" / "fsot_structure_engine.py").read_bytes()),
        "authority_sha256": sha256_bytes((ROOT / "vendor" / "fsot_compute.py").read_bytes()),
        "contact_cutoff_A": CONTACT_CUTOFF_A,
        "summary": summary,
        "results": results,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())