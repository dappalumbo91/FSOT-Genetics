#!/usr/bin/env python3
"""AlphaFold vs FSOT head-to-head on identical experimental natives.

Three columns per target, all scored as Cα RMSD (Å) to the same PDB chain:
  * AlphaFold DB model (the ML competitor, learned weights)
  * FSOT template   (real homolog transfer, self-excluded, zero trained weights)
  * FSOT bulk       (de-novo single-sequence, zero trained weights)

Shows both regimes side by side so the competitive gap is explicit.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import (  # noqa: E402
    BENCHMARK_SET, parse_pdb_ca, fetch_pdb, fetch_alphafold_pdb, kabsch_rmsd,
)
from run_rcsb_template_holdout import nw_align, best_template  # noqa: E402
from fsot_structure_engine import predict_ca_coords  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CACHE.mkdir(parents=True, exist_ok=True)


def af_rmsd(acc: str, nseq: str, nxyz: np.ndarray):
    r = fetch_alphafold_pdb(acc, CACHE)
    if not r:
        return None
    afseq, afxyz = r
    pairs = nw_align(nseq, afseq)
    if len(pairs) < 10:
        return None
    qi = [a for a, b in pairs]
    ti = [b for a, b in pairs]
    return kabsch_rmsd(afxyz[ti], nxyz[qi])


def med(vals):
    vals = [v for v in vals if v is not None]
    return float(np.median(vals)) if vals else None


def main() -> int:
    rows = []
    print(f"{'protein':<22}{'pdb':<6}{'N':>4}  {'AlphaFold':>10}{'FSOT-tmpl':>11}{'FSOT-bulk':>11}  tmpl")
    print("-" * 78)
    for acc, pdb, chain, name in BENCHMARK_SET:
        r = fetch_pdb(pdb, chain, CACHE)
        if not r:
            continue
        nseq, nxyz = r
        af = af_rmsd(acc, nseq, nxyz)
        bulk = kabsch_rmsd(
            predict_ca_coords(nseq, canonicalize_chirality=True, observer_bulk_dim=25)["ca_coords"],
            nxyz)
        t = best_template(nseq, pdb)
        tr = kabsch_rmsd(t["model"], nxyz) if t else None
        tid = f"{t['pdb_id']} id={t['identity']:.2f}" if t else "-"
        rows.append((name, af, tr, bulk))
        af_s = f"{af:.2f}" if af is not None else "  n/a"
        tr_s = f"{tr:.2f}" if tr is not None else "  n/a"
        print(f"{name:<22}{pdb:<6}{len(nseq):>4}  {af_s:>10}{tr_s:>11}{bulk:>11.2f}  {tid}")
    print("-" * 78)
    af_m, tr_m, bl_m = (med([r[1] for r in rows]), med([r[2] for r in rows]),
                        med([r[3] for r in rows]))
    print(f"{'MEDIAN':<22}{'':<6}{'':>4}  {af_m:>10.2f}{tr_m:>11.2f}{bl_m:>11.2f}")
    within = sum(1 for r in rows if r[2] is not None and r[1] is not None and r[2] - r[1] <= 1.5)
    print(f"\nFSOT-template within 1.5 A of AlphaFold: {within}/{len(rows)}")
    print(f"FSOT-template sub-2 A: {sum(1 for r in rows if r[2] is not None and r[2] < 2.0)}/{len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
