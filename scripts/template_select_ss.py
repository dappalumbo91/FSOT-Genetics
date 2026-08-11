#!/usr/bin/env python3
"""Native-free template re-rank: greedy pool + secondary-structure geometry score.

Uses FSOT F12 regions (sequence-only) to score how well a transferred template
matches helix period distances and sheet-ish contacts — no native coords.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import (  # noqa: E402
    SsPropensity,
    detect_regions,
    CA_CA,
    PI,
    E,
    PHI,
)
from run_rcsb_template_holdout import (  # noqa: E402
    homolog_ids,
    pfam_family_pdbs,
    fetch_template_pdb,
    chains_of,
    nw_align,
    build_from_template,
    model_is_sane,
    MIN_IDENTITY,
    MIN_COVERAGE,
)
from run_fsot_vs_alphafold_structure import parse_pdb_ca  # noqa: E402

# Pauling helix i→i+4 Cα distance ~6.2 A
HELIX_I4 = math.sqrt((4 * 1.5) ** 2 + (2 * 2.3 * math.sin(2 * 100 * PI / 180)) ** 2)


def ss_geometry_score(model: np.ndarray, sequence: str) -> float:
    """Higher = better agreement with sequence-predicted SS geometry (seed scales)."""
    props = [SsPropensity.from_amino_acid(c) for c in sequence]
    regions = detect_regions(props)
    n = len(sequence)
    score = 0.0
    wsum = 0.0
    for r in regions:
        if r.kind == "H" and r.length() >= 5:
            for i in range(r.start, min(r.end - 3, n - 4)):
                d = float(np.linalg.norm(model[i] - model[i + 4]))
                # reward near ideal helix i+4
                err = abs(d - HELIX_I4)
                score += math.exp(-err / PHI)
                wsum += 1.0
        if r.kind == "E" and r.length() >= 3:
            # strand CA step ~3.5
            for i in range(r.start, min(r.end, n - 1)):
                d = float(np.linalg.norm(model[i] - model[i + 1]))
                err = abs(d - CA_CA)
                score += math.exp(-err)
                wsum += 1.0
    if wsum < 1:
        return 0.0
    # Rg fidelity to FSOT globule target
    from fsot_structure_engine import target_rg_fsot, radius_of_gyration  # noqa: WPS433

    rg = radius_of_gyration(model)
    trg = target_rg_fsot(n)
    rg_term = math.exp(-abs(rg - trg) / (PHI * E))
    return (score / wsum) * rg_term


def best_template_ss(
    sequence: str,
    exclude_pdb: str,
    identity_cap: float = 0.95,
    *,
    pool_limit: int = 20,
) -> dict | None:
    """Collect pool, rank by identity×coverage×ss_geometry_score."""
    pool = []
    seen = set()
    for pdb in homolog_ids(sequence) + pfam_family_pdbs(exclude_pdb):
        if pdb == exclude_pdb.upper() or pdb in seen:
            continue
        seen.add(pdb)
        try:
            text = fetch_template_pdb(pdb)
        except Exception:
            continue
        for chain in chains_of(text):
            try:
                tseq, tcoords = parse_pdb_ca(text, chain)
            except Exception:
                continue
            if len(tseq) < 20:
                continue
            pairs = nw_align(sequence, tseq)
            if len(pairs) < 10:
                continue
            identity = sum(1 for a, b in pairs if sequence[a] == tseq[b]) / len(pairs)
            coverage = len(pairs) / len(sequence)
            if identity > identity_cap or identity < MIN_IDENTITY or coverage < MIN_COVERAGE:
                continue
            model = build_from_template(len(sequence), tcoords, pairs)
            if not model_is_sane(model, len(sequence)):
                continue
            ss = ss_geometry_score(model, sequence)
            # combine: identity/coverage primary, SS as native-free quality
            score = coverage * identity * (1.0 + ss / PHI)
            pool.append(
                {
                    "score": score,
                    "ss_score": ss,
                    "pdb_id": pdb,
                    "chain": chain,
                    "model": model,
                    "identity": identity,
                    "coverage": coverage,
                    "selection": "identity_x_ss_geometry",
                }
            )
            if len(pool) >= pool_limit:
                break
        if len(pool) >= pool_limit:
            break
    if not pool:
        return None
    return max(pool, key=lambda r: r["score"])
