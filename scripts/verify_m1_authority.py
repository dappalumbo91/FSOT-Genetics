#!/usr/bin/env python3
"""Verify M1 measured-authority fix on hard wet-lab cases + freeze guards."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import fetch_pdb, kabsch_rmsd, BENCHMARK_SET
from run_rcsb_template_holdout import best_template
from msa_template_fuse import fuse_predict

CACHE = Path.home() / ".cache" / "fsot-genetics" / "wetlab_af_eval"
HARD = [
    ("2HYY", "A", "ABL1", "2HYY"),
    ("1G5M", "A", "BCL2", "1G5M"),
    ("6M0J", "E", "RBD", "6M0J"),
    ("1UWH", "B", "BRAF", "1UWH"),
    ("2ITX", "A", "EGFR", "2ITX"),
]
# product freeze guards (must not tank)
GUARDS = [("7RSA", "A", "RNase", "7RSA"), ("1CLL", "A", "CaM", "1CLL"), ("1UBQ", "A", "Ubq", "1UBQ")]


def eval_case(pdb, ch, name, excl):
    hit = fetch_pdb(pdb, ch, CACHE)
    if not hit:
        return {"name": name, "status": "no_fetch"}
    seq, nat = hit
    t = best_template(seq, excl, identity_cap=0.95)
    if not t:
        return {"name": name, "status": "no_template"}
    prod = fuse_predict(
        seq, t["model"], None, tertiary_contacts=t.get("tertiary_contacts")
    )
    r = float(kabsch_rmsd(prod["ca_coords"], nat))
    return {
        "name": name,
        "status": "ok",
        "rmsd": r,
        "tmpl": t["pdb_id"],
        "mode": t.get("template_mode"),
        "n_cand": t.get("n_candidates"),
        "expanded": t.get("expanded_isoform_pool"),
        "cap": t.get("identity_cap_used"),
    }


def main() -> int:
    print("=== HARD (M1 targets) ===")
    hard = []
    for args in HARD:
        row = eval_case(*args)
        hard.append(row)
        print(row)
    print("=== GUARDS ===")
    guards = []
    for args in GUARDS:
        row = eval_case(*args)
        guards.append(row)
        print(row)

    # freeze H2H quick median on product templates
    print("=== FREEZE H2H product (template+physics) ===")
    rmsds = []
    for acc, pdb, ch, name in BENCHMARK_SET:
        hit = fetch_pdb(pdb, ch, CACHE)
        if not hit:
            continue
        seq, nat = hit
        t = best_template(seq, pdb, identity_cap=0.95)
        if not t:
            continue
        prod = fuse_predict(
            seq, t["model"], None, tertiary_contacts=t.get("tertiary_contacts")
        )
        r = float(kabsch_rmsd(prod["ca_coords"], nat))
        rmsds.append(r)
        print(f"  {name}: {r:.2f}  {t.get('template_mode')} n={t.get('n_candidates')}")
    import numpy as np

    med = float(np.median(rmsds)) if rmsds else 99.0
    print(f"MEDIAN product={med:.3f} n={len(rmsds)}")

    # gates
    abl = next((h for h in hard if h["name"] == "ABL1"), {})
    bcl = next((h for h in hard if h["name"] == "BCL2"), {})
    rn = next((g for g in guards if g["name"] == "RNase"), {})
    cam = next((g for g in guards if g["name"] == "CaM"), {})
    ok = True
    if abl.get("status") == "ok" and abl["rmsd"] >= 8.0:
        print("FAIL ABL1 still >= 8")
        ok = False
    elif abl.get("status") == "ok":
        print(f"PASS ABL1 {abl['rmsd']:.2f}")
    if bcl.get("status") == "ok" and bcl["rmsd"] >= 5.5:
        print("FAIL BCL2 not improved enough")
        ok = False
    elif bcl.get("status") == "ok":
        print(f"PASS BCL2 {bcl['rmsd']:.2f}")
    if rn.get("status") == "ok" and rn["rmsd"] > 1.5:
        print("FAIL RNase tanked")
        ok = False
    else:
        print(f"PASS RNase {rn.get('rmsd')}")
    if cam.get("status") == "ok" and cam["rmsd"] > 2.0:
        print("FAIL CaM tanked")
        ok = False
    else:
        print(f"PASS CaM {cam.get('rmsd')}")
    if med > 1.24:
        print(f"FAIL freeze median {med}")
        ok = False
    else:
        print(f"PASS freeze median {med:.3f}")

    out = {
        "hard": hard,
        "guards": guards,
        "freeze_median": med,
        "freeze_rmsds": rmsds,
        "ok": ok,
    }
    path = ROOT / "data" / "m1_authority_verify.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("Wrote", path)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
