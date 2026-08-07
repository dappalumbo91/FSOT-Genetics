#!/usr/bin/env python3
"""Audit F12 regions against experimental PDB HELIX/SHEET annotations."""

from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import E, SsPropensity, detect_regions  # noqa: E402

OUTPUT = ROOT / "data" / "f12_secondary_audit.json"
PDB_DIR = ROOT / "data" / "pdb_samples"
SAMPLES = (("1UBQ", "A"), ("1CRN", "A"), ("1VII", "A"), ("2GB1", "A"), ("1ENH", "A"))
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_chain(text: str, chain: str) -> tuple[str, list[int]]:
    sequence = []
    residue_numbers = []
    seen = set()
    for line in text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[21].strip() != chain:
            continue
        key = (line[22:26].strip(), line[26].strip())
        residue = AA3_TO_1.get(line[17:20].strip())
        if key in seen or residue is None:
            continue
        seen.add(key)
        sequence.append(residue)
        residue_numbers.append(int(key[0]))
    return "".join(sequence), residue_numbers


def experimental_labels(text: str, chain: str, residue_numbers: list[int]) -> list[str]:
    labels = ["C"] * len(residue_numbers)
    index = {number: position for position, number in enumerate(residue_numbers)}
    for line in text.splitlines():
        try:
            if line.startswith("HELIX ") and line[19].strip() == chain and line[31].strip() == chain:
                kind, start, end = "H", int(line[21:25]), int(line[33:37])
            elif line.startswith("SHEET ") and line[21].strip() == chain and line[32].strip() == chain:
                kind, start, end = "E", int(line[22:26]), int(line[33:37])
            else:
                continue
        except ValueError:
            continue
        for residue_number in range(start, end + 1):
            if residue_number in index:
                labels[index[residue_number]] = kind
    return labels


def predicted_labels(sequence: str) -> list[str]:
    labels = ["C"] * len(sequence)
    regions = detect_regions([SsPropensity.from_amino_acid(amino_acid) for amino_acid in sequence])
    for region in regions:
        for position in range(region.start, region.end + 1):
            labels[position] = region.kind
    return labels


def main() -> int:
    results = []
    aggregate = Counter()
    for pdb_id, chain in SAMPLES:
        text = (PDB_DIR / f"{pdb_id}.pdb").read_text(encoding="utf-8", errors="replace")
        sequence, residue_numbers = parse_chain(text, chain)
        observed = experimental_labels(text, chain, residue_numbers)
        predicted = predicted_labels(sequence)
        confusion = Counter(zip(observed, predicted))
        aggregate.update(confusion)
        recalls = {
            kind: sum(truth == kind and guess == kind for truth, guess in zip(observed, predicted))
            / max(observed.count(kind), 1)
            for kind in ("H", "E", "C")
        }
        results.append(
            {
                "pdb_id": pdb_id,
                "chain": chain,
                "length": len(sequence),
                "observed_counts": dict(Counter(observed)),
                "predicted_counts": dict(Counter(predicted)),
                "recall": recalls,
                "accuracy": sum(truth == guess for truth, guess in zip(observed, predicted)) / len(observed),
                "confusion": {f"{truth}->{guess}": count for (truth, guess), count in confusion.items()},
            }
        )

    beta_gate = 1.0 / E
    propensity_table = {}
    for amino_acid in "ARNDCEQGHILKMFPSTWYV":
        propensity = SsPropensity.from_amino_acid(amino_acid)
        propensity_table[amino_acid] = {
            "p_alpha": propensity.p_alpha,
            "p_beta": propensity.p_beta,
            "p_coil": propensity.p_coil,
            "crosses_beta_gate": propensity.p_beta > beta_gate and propensity.p_beta > propensity.p_alpha,
        }
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "role": "development_diagnostic",
        "ground_truth": "PDB HELIX and SHEET records",
        "beta_gate": beta_gate,
        "amino_acids_crossing_beta_gate": [
            amino_acid for amino_acid, row in propensity_table.items() if row["crosses_beta_gate"]
        ],
        "aggregate_confusion": {
            f"{truth}->{guess}": count for (truth, guess), count in aggregate.items()
        },
        "propensity_table": propensity_table,
        "results": results,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    beta_recalls = [row["recall"]["E"] for row in results if row["observed_counts"].get("E", 0)]
    print(f"beta_gate_residues={output['amino_acids_crossing_beta_gate']}")
    print(f"beta_recall_per_beta_protein={beta_recalls}")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())