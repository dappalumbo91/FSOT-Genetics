#!/usr/bin/env python3
"""Compare baseline F12 with the expanded cooperative candidate on development data."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import (  # noqa: E402
    PHI,
    SsPropensity,
    detect_regions_cooperative,
)
from run_f12_secondary_audit import (  # noqa: E402
    PDB_DIR,
    SAMPLES,
    experimental_labels,
    parse_chain,
    predicted_labels,
)

OUTPUT = ROOT / "data" / "f12_candidate_development.json"


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


def metrics(observed: list[str], predicted: list[str]) -> dict:
    confusion = Counter(zip(observed, predicted))
    recalls = {
        kind: sum(truth == kind and guess == kind for truth, guess in zip(observed, predicted))
        / max(observed.count(kind), 1)
        for kind in ("H", "E", "C")
    }
    return {
        "accuracy": sum(truth == guess for truth, guess in zip(observed, predicted))
        / len(observed),
        "recall": recalls,
        "macro_recall": sum(recalls.values()) / 3.0,
        "confusion": {
            f"{truth}->{guess}": count for (truth, guess), count in sorted(confusion.items())
        },
    }


def main() -> int:
    results = []
    for pdb_id, chain in SAMPLES:
        path = PDB_DIR / f"{pdb_id}.pdb"
        text = path.read_text(encoding="utf-8", errors="replace")
        sequence, residue_numbers = parse_chain(text, chain)
        observed = experimental_labels(text, chain, residue_numbers)
        results.append(
            {
                "pdb_id": pdb_id,
                "chain": chain,
                "pdb_sha256": sha256_file(path),
                "length": len(sequence),
                "observed_counts": dict(Counter(observed)),
                "baseline": metrics(observed, predicted_labels(sequence)),
                "candidate": metrics(observed, candidate_labels(sequence)),
            }
        )

    beta_rows = [row for row in results if row["observed_counts"].get("E", 0)]
    baseline_macro = sum(row["baseline"]["macro_recall"] for row in results) / len(results)
    candidate_macro = sum(row["candidate"]["macro_recall"] for row in results) / len(results)
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "role": "development_only",
        "git_commit_before_candidate_freeze": git_commit(),
        "engine_sha256": sha256_file(ROOT / "scripts" / "fsot_structure_engine.py"),
        "ground_truth": "PDB HELIX and SHEET records",
        "candidate": {
            "beta_topology": "max(branch,0)+abs(aromatic)+abs(hetero)/phi",
            "raw_beta": "exp((volume-polarity+beta_topology)/pi)",
            "same_state_factor": "phi^(1/phi)",
            "selection_family_same_state_exponents": ["1/pi", "1/phi", "1", "phi"],
            "production_enabled": False,
        },
        "gates": {
            "nonzero_beta_recall_each_beta_protein": all(
                row["candidate"]["recall"]["E"] > 0.0 for row in beta_rows
            ),
            "macro_recall_improved": candidate_macro > baseline_macro,
        },
        "summary": {
            "baseline_mean_per_protein_macro_recall": baseline_macro,
            "candidate_mean_per_protein_macro_recall": candidate_macro,
            "baseline_mean_accuracy": sum(row["baseline"]["accuracy"] for row in results)
            / len(results),
            "candidate_mean_accuracy": sum(row["candidate"]["accuracy"] for row in results)
            / len(results),
        },
        "results": results,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(
        f"macro_recall={baseline_macro:.4f}->{candidate_macro:.4f} "
        f"mean_accuracy={output['summary']['baseline_mean_accuracy']:.4f}"
        f"->{output['summary']['candidate_mean_accuracy']:.4f}"
    )
    print(f"gates={output['gates']}")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
