#!/usr/bin/env python3
"""Per-structure (structural-consensus) template selection vs global identity.

User's idea: instead of one global identity cap + greedy highest-identity pick,
judge each candidate by whether its SHAPE agrees with the other candidates
(structural medoid), independent of sequence identity. Native is used only to
score the two selectors, never to choose.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import BENCHMARK_SET, fetch_pdb, parse_pdb_ca, kabsch_rmsd  # noqa: E402
from run_rcsb_template_holdout import (  # noqa: E402
    homolog_ids, pfam_family_pdbs, fetch_template_pdb, chains_of, nw_align,
    build_from_template, model_is_sane, MIN_COVERAGE,
)

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
BAND = (0.45, 0.99)  # wide band, NO tight global cap


def candidates(nseq, pdb, cap_lo=BAND[0], cap_hi=BAND[1], limit=16):
    out, seen = [], set()
    for cid in homolog_ids(nseq) + pfam_family_pdbs(pdb):
        if cid == pdb.upper() or cid in seen:
            continue
        seen.add(cid)
        try:
            text = fetch_template_pdb(cid)
        except Exception:
            continue
        for ch in chains_of(text):
            try:
                tseq, tc = parse_pdb_ca(text, ch)
            except Exception:
                continue
            if len(tseq) < 20:
                continue
            pairs = nw_align(nseq, tseq)
            if len(pairs) < 10:
                continue
            ident = sum(1 for a, b in pairs if nseq[a] == tseq[b]) / len(pairs)
            cov = len(pairs) / len(nseq)
            if not (cap_lo <= ident <= cap_hi) or cov < MIN_COVERAGE:
                continue
            model = build_from_template(len(nseq), tc, pairs)
            if not model_is_sane(model, len(nseq)):
                continue
            out.append((ident, cov, cid, model))
        if len(out) >= limit:
            break
    return out


def medoid(models):
    n = len(models)
    if n == 1:
        return 0
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = kabsch_rmsd(models[i], models[j])
    return int(np.argmin(D.sum(1)))


def main() -> int:
    print(f"{'protein':<22}{'#cand':>6}{'greedy-id':>11}{'medoid':>9}  greedy/medoid template")
    print("-" * 78)
    g_list, m_list = [], []
    for acc, pdb, chain, name in BENCHMARK_SET:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            continue
        nseq, nxyz = r
        cand = candidates(nseq, pdb)
        if not cand:
            print(f"{name:<22}{0:>6}  no sane candidate")
            continue
        gi = int(np.argmax([c[0] * c[1] for c in cand]))
        models = [c[3] for c in cand]
        mi = medoid(models)
        g_rmsd = kabsch_rmsd(cand[gi][3], nxyz)
        m_rmsd = kabsch_rmsd(cand[mi][3], nxyz)
        g_list.append(g_rmsd)
        m_list.append(m_rmsd)
        print(f"{name:<22}{len(cand):>6}{g_rmsd:>11.2f}{m_rmsd:>9.2f}  "
              f"{cand[gi][2]}(id{cand[gi][0]:.2f}) / {cand[mi][2]}(id{cand[mi][0]:.2f})")
    print("-" * 78)
    print(f"{'MEDIAN':<22}{'':>6}{np.median(g_list):>11.2f}{np.median(m_list):>9.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
