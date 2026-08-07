#!/usr/bin/env python3
"""Probe lawful multi-scale D_eff routings on offline PDB samples.

Does NOT free-fit D. Only named pin-table routings from domain_interface.ROUTINGS.
Reports Cα RMSD vs experimental for residual-at-interface diagnosis.

  python scripts/run_deff_interface_probe.py
"""

from __future__ import annotations

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from domain_interface import ROUTINGS, dump_protein_domains  # noqa: E402
from fsot_structure_engine import predict_ca_coords, clean_sequence  # noqa: E402
from run_fsot_vs_alphafold_structure import parse_pdb_ca, kabsch_rmsd, align_by_sequence  # noqa: E402

OUT = ROOT / "data" / "deff_interface_probe.json"
OUT_MD = ROOT / "predictions" / "reports" / "DEFF_INTERFACE_PROBE.md"
PDB_DIR = ROOT / "data" / "pdb_samples"

# Offline samples: (pdb_file, chain, plain_name, optional_fasta_override)
SAMPLES = [
    ("1UBQ.pdb", "A", "Ubiquitin"),
    ("1CRN.pdb", "A", "Crambin"),
    ("1VII.pdb", "A", "Villin"),
    ("2GB1.pdb", "A", "Protein G B1"),
    ("1ENH.pdb", "A", "Engrailed HD"),
]


def main() -> int:
    print("=" * 64)
    print("D_eff multi-interface probe (named pin domains only)")
    print("  free_parameters = 0  · residual-at-interface diagnosis")
    print("=" * 64)

    # load experimental
    exps = []
    for fname, chain, name in SAMPLES:
        path = PDB_DIR / fname
        if not path.is_file():
            print(f"  skip {fname} (missing)")
            continue
        seq, xyz = parse_pdb_ca(path.read_text(encoding="utf-8", errors="replace"), chain)
        if len(seq) < 20:
            print(f"  skip {fname} short n={len(seq)}")
            continue
        exps.append({"pdb": fname, "name": name, "chain": chain, "seq": seq, "xyz": xyz})
        print(f"  loaded {fname} n={len(seq)} {name}")

    rows = []
    for rname in ROUTINGS:
        print(f"\n--- routing: {rname} ---")
        rmsds = []
        per = []
        t0 = time.perf_counter()
        for exp in exps:
            pred = predict_ca_coords(exp["seq"], rounds=20, routing=rname)
            p_xyz, e_xyz, Ln = align_by_sequence(
                pred["sequence"], pred["ca_coords"], exp["seq"], exp["xyz"]
            )
            if Ln < 20:
                continue
            rmsd = kabsch_rmsd(p_xyz, e_xyz)
            rmsds.append(rmsd)
            per.append(
                {
                    "pdb": exp["pdb"],
                    "name": exp["name"],
                    "n": Ln,
                    "rmsd_A": rmsd,
                    "predict_ms": pred.get("predict_ms"),
                    "gate": pred.get("long_range_gate"),
                    "D_eff_region": pred.get("D_eff_region"),
                }
            )
            print(f"  {exp['pdb']:8s} RMSD={rmsd:7.3f} Å  n={Ln}  gate={pred.get('long_range_gate')}")
        wall = time.perf_counter() - t0
        med = float(np.median(rmsds)) if rmsds else None
        mean = float(np.mean(rmsds)) if rmsds else None
        rows.append(
            {
                "routing": rname,
                "notes": ROUTINGS[rname].notes,
                "median_rmsd_A": med,
                "mean_rmsd_A": mean,
                "n_proteins": len(rmsds),
                "wall_s": wall,
                "per_protein": per,
            }
        )
        print(f"  → median RMSD = {med} Å  wall={wall:.2f}s")

    # rank (diagnosis only — does not change default unless we theory-agree)
    ranked = sorted([r for r in rows if r["median_rmsd_A"] is not None], key=lambda x: x["median_rmsd_A"])
    doc = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "doctrine": "named pin domains only; residual-at-interface; not free D search",
        "free_parameters": 0,
        "default_routing": "multi_scale_v9",
        "domains": dump_protein_domains(),
        "results": rows,
        "ranked_by_median_rmsd": [r["routing"] for r in ranked],
        "literature_note": (
            "Folding free-energy landscapes are often low-dimensional in reaction "
            "coordinates (Onuchic/Wolynes/Dill) while configuration space is high-D. "
            "FSOT D_eff is the domain interface in the 35-domain table — not PCA dim."
        ),
    }
    OUT.write_text(json.dumps(doc, indent=2), encoding="utf-8")

    lines = [
        "# D_eff multi-interface probe",
        "",
        f"*Generated {doc['generated_at']}*",
        "",
        "Named pin-table routings only. **Zero free parameters.** Residual-at-interface diagnosis.",
        "",
        "## Ranked by median Cα RMSD (offline samples)",
        "",
        "| Rank | Routing | Median RMSD Å | Mean | n |",
        "|-----:|---------|--------------:|-----:|--:|",
    ]
    for i, r in enumerate(ranked, 1):
        lines.append(
            f"| {i} | `{r['routing']}` | {r['median_rmsd_A']:.3f} | {r['mean_rmsd_A']:.3f} | {r['n_proteins']} |"
        )
    lines.extend(
        [
            "",
            f"Default theory routing remains **`multi_scale_v9`** (not auto-switched by this table).",
            "",
            "## Notes",
            "",
            doc["literature_note"],
            "",
            "See `docs/DOMAIN_INTERFACE_FOLD.md` and `docs/RUNTIME_STACK.md`.",
            "",
        ]
    )
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print("\n" + "=" * 64)
    print("RANKED (diagnosis)")
    for i, r in enumerate(ranked, 1):
        print(f"  {i}. {r['routing']:20s}  median={r['median_rmsd_A']:.3f} Å")
    print(f"Wrote {OUT}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
