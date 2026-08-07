#!/usr/bin/env python3
"""Run the frozen F12 candidate once on the pinned Kaggle/PISCES DSSP set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import (  # noqa: E402
    SsPropensity,
    detect_regions_cooperative,
)
from run_f12_secondary_audit import predicted_labels  # noqa: E402

MANIFEST = ROOT / "data" / "kaggle_f12_validation_manifest.json"
OUTPUT = ROOT / "data" / "kaggle_f12_validation_eval.json"
AA20 = set("ARNDCEQGHILKMFPSTWYV")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def candidate_labels(sequence: str) -> list[str]:
    labels = ["C"] * len(sequence)
    props = [SsPropensity.from_expanded_amino_acid(amino_acid) for amino_acid in sequence]
    for region in detect_regions_cooperative(props):
        for position in range(region.start, region.end + 1):
            labels[position] = region.kind
    return labels


def update_metrics(
    observed: str,
    predicted: list[str],
    confusion: Counter,
) -> tuple[float, float]:
    confusion.update(zip(observed, predicted))
    accuracy = sum(truth == guess for truth, guess in zip(observed, predicted)) / len(observed)
    beta_count = observed.count("E")
    beta_recall = (
        sum(truth == "E" and guess == "E" for truth, guess in zip(observed, predicted))
        / beta_count
        if beta_count
        else float("nan")
    )
    return accuracy, beta_recall


def summarize(confusion: Counter) -> dict:
    total = sum(confusion.values())
    recalls = {}
    for kind in ("H", "E", "C"):
        observed = sum(count for (truth, _), count in confusion.items() if truth == kind)
        recalls[kind] = confusion[(kind, kind)] / observed
    return {
        "accuracy": sum(confusion[(kind, kind)] for kind in ("H", "E", "C")) / total,
        "recall": recalls,
        "macro_recall": sum(recalls.values()) / 3.0,
        "confusion": {
            f"{truth}->{guess}": confusion[(truth, guess)]
            for truth in ("H", "E", "C")
            for guess in ("H", "E", "C")
        },
    }


def bootstrap_mean_interval(values: list[float]) -> list[float]:
    array = np.asarray(values, dtype=np.float64)
    generator = np.random.default_rng(314159265)
    means = np.empty(2000, dtype=np.float64)
    for index in range(len(means)):
        means[index] = generator.choice(array, size=len(array), replace=True).mean()
    return [float(value) for value in np.quantile(means, [0.025, 0.975])]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path, help="Path to the pinned Kaggle CSV")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    actual_hash = sha256_file(args.csv)
    expected_hash = manifest["dataset"]["sha256"]
    if actual_hash != expected_hash:
        raise SystemExit(f"dataset SHA-256 mismatch: {actual_hash} != {expected_hash}")

    excluded = set(manifest["excluded_development_pdb_ids"])
    baseline_confusion: Counter = Counter()
    candidate_confusion: Counter = Counter()
    baseline_accuracies = []
    candidate_accuracies = []
    baseline_beta_recalls = []
    candidate_beta_recalls = []
    skipped = Counter()
    residue_count = 0

    with args.csv.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            pdb_id = row["pdb_id"].upper()
            sequence = row["seq"].upper()
            observed = row["sst3"].upper()
            if pdb_id in excluded:
                skipped["development_pdb"] += 1
                continue
            if len(sequence) != len(observed) or not sequence:
                skipped["unaligned"] += 1
                continue
            if set(sequence) - AA20 or set(observed) - set("HEC"):
                skipped["nonstandard"] += 1
                continue

            baseline_accuracy, baseline_beta = update_metrics(
                observed, predicted_labels(sequence), baseline_confusion
            )
            candidate_accuracy, candidate_beta = update_metrics(
                observed, candidate_labels(sequence), candidate_confusion
            )
            baseline_accuracies.append(baseline_accuracy)
            candidate_accuracies.append(candidate_accuracy)
            if not np.isnan(baseline_beta):
                baseline_beta_recalls.append(baseline_beta)
                candidate_beta_recalls.append(candidate_beta)
            residue_count += len(sequence)

    baseline = summarize(baseline_confusion)
    candidate = summarize(candidate_confusion)
    paired_accuracy_delta = [
        candidate_value - baseline_value
        for baseline_value, candidate_value in zip(baseline_accuracies, candidate_accuracies)
    ]
    gates = {
        "aggregate_macro_recall": candidate["macro_recall"] > baseline["macro_recall"],
        "aggregate_beta_recall": candidate["recall"]["E"] > baseline["recall"]["E"],
        "aggregate_accuracy": candidate["accuracy"] >= baseline["accuracy"] - 0.01,
    }
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "role": "frozen_external_validation_result",
        "git_commit": git_commit(),
        "manifest_sha256": sha256_file(MANIFEST),
        "dataset_sha256": actual_hash,
        "evaluated_chains": len(baseline_accuracies),
        "evaluated_residues": residue_count,
        "skipped": dict(skipped),
        "baseline": baseline,
        "candidate": candidate,
        "per_chain": {
            "baseline_mean_accuracy": float(np.mean(baseline_accuracies)),
            "candidate_mean_accuracy": float(np.mean(candidate_accuracies)),
            "mean_paired_accuracy_delta": float(np.mean(paired_accuracy_delta)),
            "paired_accuracy_delta_95pct_bootstrap": bootstrap_mean_interval(
                paired_accuracy_delta
            ),
            "baseline_mean_beta_recall": float(np.mean(baseline_beta_recalls)),
            "candidate_mean_beta_recall": float(np.mean(candidate_beta_recalls)),
        },
        "success_gates": gates,
        "passed": all(gates.values()),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"chains={output['evaluated_chains']} residues={residue_count}")
    print(
        f"macro_recall={baseline['macro_recall']:.4f}->{candidate['macro_recall']:.4f} "
        f"beta_recall={baseline['recall']['E']:.4f}->{candidate['recall']['E']:.4f} "
        f"accuracy={baseline['accuracy']:.4f}->{candidate['accuracy']:.4f}"
    )
    print(f"gates={gates} passed={output['passed']}")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
