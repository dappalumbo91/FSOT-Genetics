#!/usr/bin/env python3
"""Diagnose product systems with RMSD > ~1.2 A vs wet-lab experimental.

For each high-error case:
  - residual-at-interface energy of transferred map
  - per-system residual channels (bond/clash/Rg)
  - multi-system ChemLink routing slice
  - primary error mode + hotspot termini fraction
  - template competition (data score vs residual energy)

Writes data/high_rmsd_system_diagnosis.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

from wetlab_benchmark_catalog import STRUCTURE_CASES
from run_fsot_vs_alphafold_structure import fetch_pdb, kabsch_rmsd
from run_rcsb_template_holdout import (
    best_template,
    collect_template_candidates,
    residual_interface_energy,
    CA_CA,
    _R_BOND,
    _R_CLASH,
    _R_FOLD,
    _CLASH,
)
from msa_template_fuse import fuse_predict
from run_error_margin_log import analyze_one
from domain_interface import get_routing, domain_slice
from fsot_structure_engine import target_rg_fsot

CACHE = Path.home() / ".cache" / "fsot-genetics" / "high_rmsd_diag"
CACHE.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data" / "high_rmsd_system_diagnosis.json"
THRESH = 1.20


def channel_breakdown(X: np.ndarray) -> dict:
    n = len(X)
    d = X[1:] - X[:-1]
    L = np.linalg.norm(d, axis=1)
    e_bond = float(_R_BOND * np.sum((L - CA_CA) ** 2))
    e_clash = 0.0
    step = max(1, n // 80)
    for i in range(0, n, step):
        for j in range(i + 2, n, step):
            dist = float(np.linalg.norm(X[i] - X[j]))
            if dist < _CLASH:
                e_clash += (_CLASH - dist) ** 2
    e_clash *= _R_CLASH
    rg = float(np.sqrt(((X - X.mean(0)) ** 2).sum(axis=1).mean()))
    trg = float(target_rg_fsot(n))
    e_fold = float(_R_FOLD * (rg - trg) ** 2)
    return {
        "e_bond_Physical_Chemistry": e_bond,
        "e_clash_Chemistry": e_clash,
        "e_fold_Biochemistry": e_fold,
        "rg_A": rg,
        "target_rg_A": trg,
        "total": e_bond + e_clash + e_fold,
    }


def termini_hotspot_frac(hotspots, n: int) -> dict:
    if not hotspots or n < 4:
        return {"n_term": 0, "c_term": 0, "core": 0, "frac_termini": 0.0}
    band = max(3, n // 10)
    nt = ct = core = 0
    for h in hotspots:
        i = int(h["i"])
        if i < band:
            nt += 1
        elif i >= n - band:
            ct += 1
        else:
            core += 1
    t = nt + ct + core
    return {
        "n_term": nt,
        "c_term": ct,
        "core": core,
        "frac_termini": (nt + ct) / max(t, 1),
    }


def main() -> int:
    routing = get_routing()
    systems = []
    for case in STRUCTURE_CASES:
        hit = fetch_pdb(case["pdb"], case["chain"], CACHE)
        if not hit:
            continue
        seq, nat = hit
        tmpl = best_template(seq, case["pdb"], identity_cap=0.95)
        if not tmpl:
            systems.append(
                {
                    "id": case["id"],
                    "status": "no_template",
                    "category": case["category"],
                }
            )
            continue
        prod = fuse_predict(
            seq, tmpl["model"], None, flip_model=tmpl.get("flip_model")
        )
        X = prod["ca_coords"]
        rmsd = float(kabsch_rmsd(X, nat))
        rmsd_flip = None
        if prod.get("ca_coords_flip") is not None:
            rmsd_flip = float(kabsch_rmsd(prod["ca_coords_flip"], nat))
        rmsd_app = rmsd if rmsd_flip is None else min(rmsd, rmsd_flip)
        if rmsd_app <= THRESH:
            continue  # only high-error systems (apparatus = both collapses)
        margin = analyze_one(
            case["name"],
            f"{case['pdb']}:{case['chain']}",
            seq,
            nat,
            {"sequence": seq, "ca_coords": X, "long_range_gate": 7},
        )
        # competitor templates residual vs data score
        cands = collect_template_candidates(seq, case["pdb"], 0.95)
        if len(cands) < 3:
            cands = collect_template_candidates(seq, case["pdb"], 0.99)
        for c in cands:
            c["residual_energy"] = residual_interface_energy(c["model"])
        by_data = sorted(cands, key=lambda c: -c["score"])[:5]
        by_res = sorted(cands, key=lambda c: c["residual_energy"])[:5]
        # multi-system domain slices active for proteins
        domains = {
            name: {
                "D_eff": domain_slice(name).D_eff,
                "S": domain_slice(name).S,
                "abs_S": domain_slice(name).abs_S,
            }
            for name in (
                "Physical_Chemistry",
                "Chemistry",
                "Molecular_Chemistry",
                "Biochemistry",
                "Condensed_Matter",
                "Atomic_Physics",
            )
        }
        hot = termini_hotspot_frac(margin.get("hotspots") or [], len(seq))
        row = {
            "id": case["id"],
            "name": case["name"],
            "category": case["category"],
            "wetlab": case.get("wetlab"),
            "n": len(seq),
            "rmsd_A": rmsd,
            "rmsd_flip_A": rmsd_flip,
            "rmsd_apparatus_A": rmsd_app,
            "flip_pdb": tmpl.get("flip_pdb"),
            "template_pdb": tmpl.get("pdb_id"),
            "template_mode": tmpl.get("template_mode"),
            "template_identity": tmpl.get("identity"),
            "template_coverage": tmpl.get("coverage"),
            "template_residual_energy": tmpl.get("residual_energy"),
            "n_candidates": tmpl.get("n_candidates"),
            "primary_error_mode": margin.get("primary_error_mode"),
            "ranked_modes": margin.get("ranked_modes"),
            "mode_scores": margin.get("mode_scores"),
            "per_res_p90_A": margin.get("per_res_p90_A"),
            "contact_mae_A": margin.get("contact_mae_A"),
            "rg_err_A": margin.get("rg_err_A"),
            "channel_breakdown_product": channel_breakdown(X),
            "channel_breakdown_template": channel_breakdown(tmpl["model"]),
            "termini_hotspots": hot,
            "top5_data_templates": [
                {
                    "pdb": c["pdb_id"],
                    "score": c["score"],
                    "E_res": c["residual_energy"],
                    "id": c["identity"],
                    "cov": c["coverage"],
                }
                for c in by_data
            ],
            "top5_residual_templates": [
                {
                    "pdb": c["pdb_id"],
                    "score": c["score"],
                    "E_res": c["residual_energy"],
                    "id": c["identity"],
                    "cov": c["coverage"],
                }
                for c in by_res
            ],
            "multi_system_domains": domains,
            "routing": {
                "chem_D": routing["chem"]["D_eff"],
                "ss_D": routing["ss"]["D_eff"],
                "region_D": routing["region"]["D_eff"],
                "pack_D": routing["packing"]["D_eff"],
            },
            "hypothesis": None,
        }
        # auto hypothesis from signals
        hyps = []
        if hot["frac_termini"] >= 0.4:
            hyps.append("termini_disorder")
        if (margin.get("rg_err_A") or 0) > 2.0:
            hyps.append("global_Rg_mismatch_multi_domain_or_state")
        if (margin.get("contact_mae_A") or 0) > 2.0:
            hyps.append("tertiary_contact_map_weak")
        if tmpl.get("identity", 1) < 0.7:
            hyps.append("remote_homolog_map")
        if by_res and by_data and by_res[0]["pdb_id"] != by_data[0]["pdb_id"]:
            hyps.append("data_best_vs_residual_best_disagree")
        ch = row["channel_breakdown_product"]
        dominant_ch = max(
            [
                ("bond", ch["e_bond_Physical_Chemistry"]),
                ("clash", ch["e_clash_Chemistry"]),
                ("fold_Rg", ch["e_fold_Biochemistry"]),
            ],
            key=lambda x: x[1],
        )[0]
        hyps.append(f"dominant_channel_{dominant_ch}")
        row["hypothesis"] = hyps
        systems.append(row)
        print(
            f"{row['id']:16} {rmsd:5.2f}  {row['primary_error_mode']:22}  "
            f"term={hot['frac_termini']:.2f}  hyp={hyps}",
            flush=True,
        )

    # H2H high as well
    h2h = json.loads((ROOT / "data" / "product_vs_alphafold.json").read_text())
    h2h_high = [
        {
            "name": r["name"],
            "rmsd": r["fsot_product_rmsd_A"],
            "template": r.get("template_pdb"),
        }
        for r in h2h.get("results") or []
        if (r.get("fsot_product_rmsd_A") or 0) > THRESH
    ]

    # aggregate fix levers from hypotheses
    from collections import Counter

    all_h = Counter()
    for s in systems:
        for h in s.get("hypothesis") or []:
            all_h[h] += 1

    out = {
        "threshold_A": THRESH,
        "n_high": len(systems),
        "hypothesis_counts": dict(all_h),
        "systems": systems,
        "h2h_above_threshold": h2h_high,
        "fsot_application": "data authority when strong; residual channels bond/clash/Rg; multi-system D_eff from pin",
        "next_levers": [
            "termini_disorder → soft termini rebuild (seed CA_CA walk) only where residual fold stress high at ends",
            "data_best_vs_residual_best_disagree → residual among near-ties only (already gated)",
            "tertiary_contact_map_weak → ChemLink long-range residual at Biochemistry D=13 on measured contacts only",
            "global_Rg_mismatch → domain-split assembly for multi-domain chains",
        ],
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nWrote {OUT}  n_high={len(systems)}")
    print("hypothesis counts:", dict(all_h))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
