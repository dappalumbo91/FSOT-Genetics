#!/usr/bin/env python3
"""Live RCSB API holdout: rule-selected, non-curated structures, one-shot eval.

Selection is programmatic and deterministic (no hand curation): the RCSB Search
API returns X-ray, high-resolution, monomer-sized entries sorted by id; any
structure already present in an existing manifest is excluded so the set is
unseen. Compares the production default fold against the opt-in D_eff bulk
observer. Falsification only; no formula is fit to these results.
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
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

OUTPUT = ROOT / "data" / "rcsb_live_api_holdout_eval.json"
MANIFEST_OUT = ROOT / "data" / "rcsb_live_api_holdout_manifest.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_live_api_holdout"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json="
TARGET_N = 60
MIN_LEN, MAX_LEN = 50, 150
MAX_RESOLUTION = 1.5
BULK_DIM = 25

QUERY = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "exptl.method", "operator": "exact_match",
                "value": "X-RAY DIFFRACTION"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less", "value": MAX_RESOLUTION}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.deposited_polymer_monomer_count",
                "operator": "range", "value": {"from": MIN_LEN, "to": MAX_LEN}}},
        ],
    },
    "return_type": "entry",
    "request_options": {
        "paginate": {"start": 0, "rows": 600},
        "sort": [{"sort_by": "rcsb_id", "direction": "asc"}],
        "results_content_type": ["experimental"],
    },
}


def search_ids() -> list[str]:
    url = SEARCH_URL + urllib.parse.quote(json.dumps(QUERY))
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = json.loads(response.read())
    return [row["identifier"] for row in payload.get("result_set", [])]


def seen_ids() -> set[str]:
    seen: set[str] = set()
    for name in ("rcsb_oriented_backbone_manifest.json", "rcsb_holdout_manifest.json"):
        path = ROOT / "data" / name
        if not path.exists():
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        for entry in data.get("entries", data.get("holdout", [])):
            if isinstance(entry, dict) and entry.get("pdb_id"):
                seen.add(str(entry["pdb_id"]).upper())
    return seen


def first_valid_chain(text: str) -> str | None:
    counts: dict[str, int] = {}
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            counts[line[21]] = counts.get(line[21], 0) + 1
    for chain in sorted(counts):
        if MIN_LEN <= counts[chain] <= MAX_LEN and chain.strip():
            return chain
    return None


def main() -> int:
    excluded = seen_ids()
    candidates = [pid for pid in search_ids() if pid.upper() not in excluded]
    rows, manifest = [], []
    for pdb_id in candidates:
        if len(rows) >= TARGET_N:
            break
        try:
            text, pdb_sha256 = fetch_pdb(pdb_id, CACHE)
        except Exception:
            continue
        chain = first_valid_chain(text)
        if chain is None:
            continue
        try:
            sequence, native = parse_pdb_ca(text, chain)
        except Exception:
            continue
        if not (MIN_LEN <= len(sequence) <= MAX_LEN) or len(native) != len(sequence):
            continue
        baseline = predict_ca_coords(sequence)
        bulk = predict_ca_coords(
            sequence, canonicalize_chirality=True, observer_bulk_dim=BULK_DIM
        )
        rows.append({
            "pdb_id": pdb_id,
            "chain": chain,
            "length": len(sequence),
            "pdb_sha256": pdb_sha256,
            "rmsd_A": {
                "baseline": kabsch_rmsd(baseline["ca_coords"], native),
                "bulk_observer": kabsch_rmsd(bulk["ca_coords"], native),
            },
        })
        manifest.append({"pdb_id": pdb_id, "chain": chain, "length": len(sequence)})

    baseline_rmsd = [r["rmsd_A"]["baseline"] for r in rows]
    bulk_rmsd = [r["rmsd_A"]["bulk_observer"] for r in rows]
    paired_delta = [b - a for a, b in zip(baseline_rmsd, bulk_rmsd)]
    wins = int(sum(1 for d in paired_delta if d < 0))
    aggregate = {
        "n_chains": len(rows),
        "baseline_rmsd_A": bootstrap_median(baseline_rmsd) if rows else None,
        "bulk_observer_rmsd_A": bootstrap_median(bulk_rmsd) if rows else None,
        "paired_delta_median_A": float(np.median(paired_delta)) if rows else None,
        "bulk_wins_fraction": wins / len(rows) if rows else None,
    }
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "Live RCSB API holdout; falsification only, no formula fit",
        "candidate_commit": git_commit(),
        "selection": {
            "source": "RCSB Search API v2",
            "method": "X-RAY DIFFRACTION",
            "max_resolution_A": MAX_RESOLUTION,
            "length_range": [MIN_LEN, MAX_LEN],
            "sort": "rcsb_id asc",
            "excluded_prior_manifest_ids": sorted(excluded),
            "bulk_dim": BULK_DIM,
        },
        "aggregate": aggregate,
        "results": rows,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    MANIFEST_OUT.write_text(
        json.dumps({"query": QUERY, "entries": manifest}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(aggregate, indent=2))
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
