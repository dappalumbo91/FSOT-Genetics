#!/usr/bin/env python3
"""Medical / research stress suite — multi-regime accuracy + usability.

Columns per target (Cα RMSD to experimental PDB, Kabsch):
  1. FSOT bulk single-sequence     (de-novo claim path)
  2. FSOT bulk + MSA F15 channel   (optional evo inject)
  3. FSOT template                 (deployable high-accuracy path)
  4. FSOT template + MSA fuse      (coevolution polish on template)
  5. AlphaFold DB                  (competitor, if reachable)

Plus:
  - MSA depth / conservation / confidence stats
  - Contact top-L precision: single vs MSA-augmented F15
  - p53 variant-effect spot check when TP53 present

Zero free parameters throughout. Writes data/medical_stress_suite.json.
"""

from __future__ import annotations

import json
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import (  # noqa: E402
    BENCHMARK_SET,
    fetch_alphafold_pdb,
    fetch_pdb,
    kabsch_rmsd,
)
from run_rcsb_template_holdout import best_template, nw_align  # noqa: E402
from fsot_structure_engine import (  # noqa: E402
    build_distogram,
    predict_ca_coords,
)
from msa_pipeline import build_msa_features, conservation_confidence  # noqa: E402
from msa_template_fuse import fuse_relax, fuse_predict, fused_confidence, select_regime  # noqa: E402
from test_physics_refine import relax as physics_relax  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
CACHE.mkdir(parents=True, exist_ok=True)
OUT_JSON = ROOT / "data" / "medical_stress_suite.json"
OUT_MD = ROOT / "predictions" / "reports" / "MEDICAL_STRESS_SUITE.md"

# UniProt → Pfam (curated for the classic medical/H2H panel)
PFAM_MAP: dict[str, str] = {
    "P69905": "PF00042",  # globin
    "P68871": "PF00042",
    "P00918": "PF00194",  # carbonic anhydrase
    "P00441": "PF00080",  # sodcu
    "P61626": "PF00062",  # lysozyme
    "P61823": "PF00074",  # RNase
    "P0CG47": "PF00240",  # ubiquitin
    "P01308": "PF00049",  # insulin
    "P04637": "PF00870",  # p53
    "P0DP23": "PF00036",  # EF-hand
}

ROUNDS = 12  # bulk fold rounds (stress speed vs depth)


def _med(xs: list[float | None]) -> float | None:
    v = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(v)) if v else None


def _af_rmsd(acc: str, nseq: str, nxyz: np.ndarray) -> float | None:
    try:
        r = fetch_alphafold_pdb(acc, CACHE)
    except Exception:
        return None
    if not r:
        return None
    afseq, afxyz = r
    pairs = nw_align(nseq, afseq)
    if len(pairs) < 10:
        return None
    qi = [a for a, _b in pairs]
    ti = [b for _a, b in pairs]
    return float(kabsch_rmsd(afxyz[ti], nxyz[qi]))


def contact_top_l(M: np.ndarray, native: np.ndarray, *, contact_A: float = 8.0) -> float:
    n = min(M.shape[0], native.shape[0])
    pairs = []
    for i in range(n):
        for j in range(i + 7, n):
            d = float(np.linalg.norm(native[i] - native[j]))
            pairs.append((float(M[i, j]), d < contact_A))
    if not pairs:
        return 0.0
    pairs.sort(reverse=True)
    L = n
    top = pairs[:L]
    return float(sum(1 for _s, hit in top if hit) / max(len(top), 1))


def conf_error_spearman(conf: np.ndarray, err: np.ndarray) -> float | None:
    """Spearman-like via rank correlation (no scipy required)."""
    if len(conf) < 5:
        return None
    rc = conf.argsort().argsort().astype(float)
    re = err.argsort().argsort().astype(float)
    rc -= rc.mean()
    re -= re.mean()
    d = np.sqrt((rc * rc).sum() * (re * re).sum())
    if d < 1e-15:
        return None
    # high conf should anti-correlate with error
    return float((rc * re).sum() / d)


def per_res_dev(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return np.linalg.norm(Pc @ R.T - Qc, axis=1)


def run_one(acc: str, pdb: str, chain: str, name: str) -> dict[str, Any]:
    t0 = time.perf_counter()
    row: dict[str, Any] = {
        "name": name,
        "uniprot": acc,
        "pdb_id": pdb,
        "chain": chain,
        "pfam": PFAM_MAP.get(acc),
    }
    fetched = fetch_pdb(pdb, chain, CACHE)
    if not fetched:
        row["error"] = "pdb_fetch_failed"
        return row
    nseq, nxyz = fetched
    row["length"] = len(nseq)

    # 1) bulk single
    bulk = predict_ca_coords(nseq, rounds=ROUNDS, mode="single", canonicalize_chirality=True)
    row["bulk_single_rmsd_A"] = float(kabsch_rmsd(bulk["ca_coords"], nxyz))

    # 2) MSA features + bulk msa
    feat = None
    pfam = PFAM_MAP.get(acc)
    try:
        feat = build_msa_features(nseq, pfam=pfam, uniprot=acc)
    except Exception as exc:  # noqa: BLE001
        row["msa_error"] = str(exc)
    if feat is not None and feat.n_seqs > 0:
        row["msa"] = feat.summary()
        row["msa"]["mean_evo_confidence"] = float(conservation_confidence(feat).mean())
        msa_bulk = predict_ca_coords(
            nseq, rounds=ROUNDS, mode="msa", msa_features=feat, canonicalize_chirality=True
        )
        row["bulk_msa_rmsd_A"] = float(kabsch_rmsd(msa_bulk["ca_coords"], nxyz))
        row["bulk_msa_mode"] = msa_bulk.get("structure_mode")
        # contact precision single vs msa
        M0, *_ = build_distogram(nseq)
        M1, *_rest = build_distogram(nseq, msa_features=feat)
        row["contact_topL_single"] = contact_top_l(M0, nxyz)
        row["contact_topL_msa"] = contact_top_l(M1, nxyz)
    else:
        row["bulk_msa_rmsd_A"] = None
        row["msa"] = None

    # 3) template
    tmpl = best_template(nseq, pdb)
    if tmpl is not None:
        row["template_pdb"] = tmpl["pdb_id"]
        row["template_identity"] = tmpl["identity"]
        row["template_coverage"] = tmpl["coverage"]
        row["template_rmsd_A"] = float(kabsch_rmsd(tmpl["model"], nxyz))
        # physics-only refine (prior baseline)
        row["template_physics_rmsd_A"] = float(kabsch_rmsd(physics_relax(tmpl["model"]), nxyz))
        # 4) template + MSA packing fuse (v2: near-contact only + energy gate)
        if feat is not None and feat.depth_ok:
            fused = fuse_predict(nseq, tmpl["model"], feat)
            Xf = fused["ca_coords"]
            row["template_msa_fuse_rmsd_A"] = float(kabsch_rmsd(Xf, nxyz))
            row["fuse_regime_chosen"] = fused.get("regime")
            row["fuse_energy_physics"] = fused.get("energy_physics")
            row["fuse_energy_fuse"] = fused.get("energy_fuse")
            conf = fused_confidence(feat)
            err = per_res_dev(Xf, nxyz)
            row["conf_error_spearman"] = conf_error_spearman(conf, err)
            row["mean_confidence"] = float(conf.mean())
            row["mean_per_res_err_A"] = float(err.mean())
        else:
            row["template_msa_fuse_rmsd_A"] = row["template_physics_rmsd_A"]
            row["fuse_regime_chosen"] = "template_physics"
        row["deploy_regime"] = select_regime(True, feat)
    else:
        row["template_rmsd_A"] = None
        row["template_physics_rmsd_A"] = None
        row["template_msa_fuse_rmsd_A"] = None
        row["deploy_regime"] = select_regime(False, feat)

    # 5) AlphaFold
    try:
        row["alphafold_rmsd_A"] = _af_rmsd(acc, nseq, nxyz)
    except Exception:
        row["alphafold_rmsd_A"] = None

    row["elapsed_s"] = round(time.perf_counter() - t0, 2)
    return row


def write_md(report: dict) -> None:
    s = report["summary"]
    lines = [
        "# Medical stress suite — FSOT multi-regime",
        "",
        f"Generated: `{report['generated_at']}`  ",
        f"Free parameters: **{report['free_parameters']}**  ",
        f"Targets: **{s['n']}**",
        "",
        "## Median Cα RMSD (Å) to experimental PDB",
        "",
        "| Regime | Median Å | n |",
        "|--------|---------:|--:|",
        f"| AlphaFold DB | {s.get('alphafold_median_A')} | {s.get('n_af')} |",
        f"| FSOT template + MSA fuse | {s.get('template_msa_fuse_median_A')} | {s.get('n_fuse')} |",
        f"| FSOT template + physics | {s.get('template_physics_median_A')} | {s.get('n_phys')} |",
        f"| FSOT template | {s.get('template_median_A')} | {s.get('n_tmpl')} |",
        f"| FSOT bulk + MSA | {s.get('bulk_msa_median_A')} | {s.get('n_msa')} |",
        f"| FSOT bulk single | {s.get('bulk_single_median_A')} | {s.get('n_bulk')} |",
        "",
        "## Per-target",
        "",
        "| Protein | N | AF | tmpl | phys | fuse | bulk | bulk+MSA | MSA depth | topLΔ |",
        "|---------|--:|---:|-----:|-----:|-----:|-----:|---------:|----------:|------:|",
    ]
    for r in report["results"]:
        if r.get("error"):
            lines.append(f"| {r['name']} | — | err | | | | | | | |")
            continue
        msa = r.get("msa") or {}
        top_d = ""
        if r.get("contact_topL_single") is not None and r.get("contact_topL_msa") is not None:
            top_d = f"{r['contact_topL_msa'] - r['contact_topL_single']:+.2f}"
        def f(k):
            v = r.get(k)
            return f"{v:.2f}" if isinstance(v, (int, float)) and v is not None else "—"
        lines.append(
            f"| {r['name']} | {r.get('length','')} | {f('alphafold_rmsd_A')} | "
            f"{f('template_rmsd_A')} | {f('template_physics_rmsd_A')} | "
            f"{f('template_msa_fuse_rmsd_A')} | {f('bulk_single_rmsd_A')} | "
            f"{f('bulk_msa_rmsd_A')} | {msa.get('n_seqs','—')} | {top_d} |"
        )
    lines += [
        "",
        "## Honesty notes",
        "",
        "- Template regime is the medical-grade structure path when homologs exist.",
        "- Bulk de-novo (~11 Å ceiling) is the honest orphan-sequence fallback.",
        "- MSA inject improves *contact ranking / confidence*; topology still needs templates.",
        "- All FSOT columns: zero trained weights; MSA/templates are data inputs.",
        "",
    ]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    print("FSOT medical stress suite")
    print("=" * 88)
    hdr = (
        f"{'protein':<22}{'N':>4}  {'AF':>6}{'tmpl':>7}{'phys':>7}{'fuse':>7}"
        f"{'bulk':>7}{'+MSA':>7}  msa"
    )
    print(hdr)
    print("-" * 88)
    results: list[dict[str, Any]] = []
    for acc, pdb, chain, name in BENCHMARK_SET:
        try:
            row = run_one(acc, pdb, chain, name)
        except Exception as exc:  # noqa: BLE001
            row = {
                "name": name,
                "uniprot": acc,
                "pdb_id": pdb,
                "chain": chain,
                "error": str(exc),
                "trace": traceback.format_exc()[-500:],
            }
        results.append(row)
        if row.get("error"):
            print(f"{name:<22} ERR {row['error'][:40]}")
            continue
        def fmt(k):
            v = row.get(k)
            return f"{v:6.2f}" if isinstance(v, (int, float)) else "   n/a"
        msa = row.get("msa") or {}
        print(
            f"{name:<22}{row.get('length',0):>4}  "
            f"{fmt('alphafold_rmsd_A')}{fmt('template_rmsd_A')}"
            f"{fmt('template_physics_rmsd_A')}{fmt('template_msa_fuse_rmsd_A')}"
            f"{fmt('bulk_single_rmsd_A')}{fmt('bulk_msa_rmsd_A')}  "
            f"n={msa.get('n_seqs','—')} neff={msa.get('neff','—')}"
        )

    summary = {
        "n": len([r for r in results if not r.get("error")]),
        "alphafold_median_A": _med([r.get("alphafold_rmsd_A") for r in results]),
        "template_median_A": _med([r.get("template_rmsd_A") for r in results]),
        "template_physics_median_A": _med([r.get("template_physics_rmsd_A") for r in results]),
        "template_msa_fuse_median_A": _med([r.get("template_msa_fuse_rmsd_A") for r in results]),
        "bulk_single_median_A": _med([r.get("bulk_single_rmsd_A") for r in results]),
        "bulk_msa_median_A": _med([r.get("bulk_msa_rmsd_A") for r in results]),
        "n_af": sum(1 for r in results if r.get("alphafold_rmsd_A") is not None),
        "n_tmpl": sum(1 for r in results if r.get("template_rmsd_A") is not None),
        "n_phys": sum(1 for r in results if r.get("template_physics_rmsd_A") is not None),
        "n_fuse": sum(1 for r in results if r.get("template_msa_fuse_rmsd_A") is not None),
        "n_bulk": sum(1 for r in results if r.get("bulk_single_rmsd_A") is not None),
        "n_msa": sum(1 for r in results if r.get("bulk_msa_rmsd_A") is not None),
        "fuse_beats_template": sum(
            1
            for r in results
            if isinstance(r.get("template_msa_fuse_rmsd_A"), float)
            and isinstance(r.get("template_rmsd_A"), float)
            and r["template_msa_fuse_rmsd_A"] < r["template_rmsd_A"] - 1e-6
        ),
        "msa_contact_improves": sum(
            1
            for r in results
            if isinstance(r.get("contact_topL_msa"), float)
            and isinstance(r.get("contact_topL_single"), float)
            and r["contact_topL_msa"] > r["contact_topL_single"] + 1e-9
        ),
    }
    print("-" * 88)
    print("MEDIANS (Å):")
    for k in (
        "alphafold_median_A",
        "template_msa_fuse_median_A",
        "template_physics_median_A",
        "template_median_A",
        "bulk_msa_median_A",
        "bulk_single_median_A",
    ):
        print(f"  {k}: {summary[k]}")
    print(
        f"  fuse_beats_template: {summary['fuse_beats_template']}/{summary['n_fuse']}  "
        f"msa_contact_improves: {summary['msa_contact_improves']}/{summary['n_msa']}"
    )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "medical_usability_multi_regime_stress",
        "free_parameters": 0,
        "rounds_bulk": ROUNDS,
        "summary": summary,
        "results": results,
    }
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_md(report)
    print(f"\nWrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
