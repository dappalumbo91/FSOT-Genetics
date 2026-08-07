#!/usr/bin/env python3
"""Error-margin log — localize WHERE the fold fails (not just global RMSD).

Doctrine (same residual discipline as MPCORB / multi-domain campaign):
  1. Predict on experimental chain sequence (not wrong UniProt polyprotein)
  2. Kabsch-align; measure residual by region of structure
  3. Rank error modes by contribution
  4. Attach literature / FSOT fix handle per mode
  5. Emit fix queue — solve modes one by one under full S=K(T1+T2+T3)

  python scripts/run_error_margin_log.py
"""

from __future__ import annotations

import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import predict_ca_coords, clean_sequence  # noqa: E402
from run_fsot_vs_alphafold_structure import (  # noqa: E402
    parse_pdb_ca,
    kabsch_rmsd,
    align_by_sequence,
)

OUT_JSON = ROOT / "data" / "error_margin_log.json"
OUT_MD = ROOT / "predictions" / "reports" / "ERROR_MARGIN_LOG.md"
PDB_DIR = ROOT / "data" / "pdb_samples"
CONTACT_A = math.pi * math.e  # ~8.54 Å F08 contact scale
PHI = (1.0 + math.sqrt(5.0)) / 2.0

# Offline classic folds (experimental chain = ground truth sequence)
SAMPLES = [
    ("1UBQ.pdb", "A", "Ubiquitin", "all-alpha/beta mixed, classic contact test"),
    ("1CRN.pdb", "A", "Crambin", "small SS-rich, disulfide"),
    ("1VII.pdb", "A", "Villin headpiece", "helical bundle"),
    ("2GB1.pdb", "A", "Protein G B1", "β hairpin + helix"),
    ("1ENH.pdb", "A", "Engrailed HD", "helix bundle"),
]

# Literature + FSOT fix handles for ranked modes (evidence-based, not free dials)
ERROR_CATALOG: dict[str, dict[str, str]] = {
    "backbone_sep1_2": {
        "meaning": "Virtual Cα–Cα bond / local geometry wrong",
        "literature": "Standard Cα virtual bond ≈ 3.8 Å; local geometry dominates short-range map (Dill polymer; Flory).",
        "fsot_handle": "F07 bb + CA_CA seed; do NOT residual-scale sep=1,2 (geometry is hard constraint).",
        "fix_class": "local_geometry",
    },
    "helix_period_3_4_7": {
        "meaning": "α-helix i,i+3/4/7 distances wrong",
        "literature": "α-helix rise 1.5 Å/res, 3.6 res/turn; Cα i→i+4 ≈ 6.2 Å (Pauling).",
        "fsot_handle": "F10 helix geometry + Chemistry D=8 interface; enforce ideal helix D when p_alpha high.",
        "fix_class": "secondary_structure",
    },
    "mid_range_5_12": {
        "meaning": "Secondary packing / loops at mid separation",
        "literature": "Secondary structure packing; loop closure; contact order (Plaxco et al.).",
        "fsot_handle": "F11 sheet + F12 regions; SS amp at Chemistry domain; full scalar residual mid-sep.",
        "fix_class": "secondary_packing",
    },
    "long_range_contacts": {
        "meaning": "Tertiary native contacts missing or false",
        "literature": "Contact maps / top-L metrics drive fold quality (CASP; Marks/Sander coevolution; AF distograms).",
        "fsot_handle": "F13–F15 + observer tertiary S at Biochemistry D=13; top-L caps; residual-at-interface on contact set.",
        "fix_class": "tertiary_contacts",
    },
    "global_topology": {
        "meaning": "Chain topology / domain packing globally wrong after Kabsch",
        "literature": "Energy landscape funnel; topology from contact order (Onuchic/Wolynes; Baker).",
        "fsot_handle": "MDS is only as good as D; fix contact D first, then sparse polish; multi-start from SS regions.",
        "fix_class": "embedding",
    },
    "per_residue_hotspots": {
        "meaning": "Local segments high RMSD after global align",
        "literature": "Flexible loops, termini disorder; core vs surface (crystallographic B-factors).",
        "fsot_handle": "Segment-wise residual; coil vs H/E different D_eff; do not over-constrain termini.",
        "fix_class": "local_segments",
    },
    "sequence_alignment_mismatch": {
        "meaning": "Predicted on wrong sequence length vs crystal (e.g. polyubiquitin)",
        "literature": "Always score on experimental construct sequence for structure benchmarks.",
        "fsot_handle": "H2H must use PDB chain sequence for fold input when claiming structure RMSD.",
        "fix_class": "protocol",
    },
}


def _kabsch_apply(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Return p rotated/translated onto q (same Kabsch as RMSD)."""
    p0 = p - p.mean(axis=0)
    q0 = q - q.mean(axis=0)
    H = p0.T @ q0
    U, _S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign = 1.0 if d >= 0 else -1.0
    R = Vt.T @ np.diag([1.0, 1.0, sign]) @ U.T
    return p0 @ R.T + q.mean(axis=0)


def distance_matrix(X: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    D = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(X[i] - X[j]))
            D[i, j] = D[j, i] = d
    return D


def sep_bins(n: int) -> list[tuple[str, int, int]]:
    """(name, lo_sep inclusive, hi_sep exclusive)."""
    return [
        ("backbone_sep1_2", 1, 3),
        ("helix_period_3_4_7", 3, 5),  # plus 7 handled specially
        ("mid_range_5_12", 5, 13),
        ("long_range_contacts", 13, n),
    ]


def analyze_one(
    name: str,
    pdb_file: str,
    exp_seq: str,
    exp_xyz: np.ndarray,
    pred: dict[str, Any],
) -> dict[str, Any]:
    pred_seq = pred["sequence"]
    pred_xyz = pred["ca_coords"]
    p_al, e_al, Ln = align_by_sequence(pred_seq, pred_xyz, exp_seq, exp_xyz)
    if Ln < 20:
        return {"name": name, "error": "align_short", "n": Ln}

    # If we predicted experimental sequence, Ln should ≈ n
    seq_match = pred_seq == exp_seq or exp_seq in pred_seq or pred_seq in exp_seq
    p_fit = _kabsch_apply(p_al, e_al)
    per_res = np.linalg.norm(p_fit - e_al, axis=1)
    rmsd = float(np.sqrt((per_res ** 2).mean()))

    De = distance_matrix(e_al)
    Dp = distance_matrix(p_fit)

    # Distance residual by separation
    bin_stats = []
    mode_scores: dict[str, float] = {}
    n = Ln
    for bname, lo, hi in sep_bins(n):
        errs = []
        for i in range(n):
            for j in range(i + lo, min(n, i + hi)):
                errs.append(abs(Dp[i, j] - De[i, j]))
        # helix period 7
        if bname == "helix_period_3_4_7":
            for i in range(n):
                j = i + 7
                if j < n:
                    errs.append(abs(Dp[i, j] - De[i, j]))
        if not errs:
            continue
        mean_e = float(np.mean(errs))
        med_e = float(np.median(errs))
        # contribution proxy: mean error * count / n²
        contrib = mean_e * len(errs) / max(n * n, 1)
        mode_scores[bname] = contrib
        bin_stats.append(
            {
                "mode": bname,
                "sep_lo": lo,
                "sep_hi": hi,
                "n_pairs": len(errs),
                "mean_abs_dA": mean_e,
                "median_abs_dA": med_e,
                "p90_abs_dA": float(np.percentile(errs, 90)),
                "contribution_proxy": contrib,
            }
        )

    # Contact metrics (native contacts in exp < CONTACT_A, |i-j|>=gate)
    gate = int(pred.get("long_range_gate") or 7)
    native = []
    for i in range(n):
        for j in range(i + gate, n):
            if De[i, j] < CONTACT_A:
                native.append((i, j, De[i, j]))
    evidence_diag = None
    if native:
        contact_mae = float(np.mean([abs(Dp[i, j] - de) for i, j, de in native]))
        # Structure top-L (from folded distances)
        L = n
        scores = []
        for i in range(n):
            for j in range(i + gate, n):
                scores.append((1.0 / max(Dp[i, j], 1e-3), i, j))
        scores.sort(reverse=True)
        top = scores[:L]
        hit = sum(1 for _, i, j in top if De[i, j] < CONTACT_A)
        top_l_prec = hit / max(L, 1)
        # Evidence top-L (pre-fold ranker — the data we need to solidify contacts)
        try:
            from contact_rank import rank_long_range_contacts, top_l_precision_vs_native
            from fsot_structure_engine import build_distogram

            M, _, regions, seq_m, _if = build_distogram(exp_seq)
            ranked = rank_long_range_contacts(seq_m, M, regions, gate)
            # only natives with sep>=gate for De
            De_lr = De.copy()
            evidence_diag = top_l_precision_vs_native(ranked, De_lr, CONTACT_A, L, gate=gate)
            # restrict native set in diagnostic to sep>=gate
            evidence_diag["top_L_precision_structure"] = top_l_prec
        except Exception as ex:
            evidence_diag = {"error": str(ex)}
        mode_scores["long_range_contacts"] = mode_scores.get("long_range_contacts", 0) + contact_mae
        if evidence_diag and evidence_diag.get("top_L_precision_evidence") is not None:
            # invert precision so low top-L increases mode pressure
            mode_scores["long_range_contacts"] += (
                1.0 - float(evidence_diag["top_L_precision_evidence"])
            ) * 5.0
    else:
        contact_mae = None
        top_l_prec = None

    # Per-residue hotspots: top 20% by error
    thr = float(np.percentile(per_res, 80))
    hotspots = [
        {"i": int(i), "aa": exp_seq[i] if i < len(exp_seq) else "?", "err_A": float(per_res[i])}
        for i in np.argsort(-per_res)[: max(5, n // 5)]
        if per_res[i] >= thr
    ]
    hotspot_mean = float(np.mean([h["err_A"] for h in hotspots])) if hotspots else 0.0
    mode_scores["per_residue_hotspots"] = hotspot_mean / max(rmsd, 1e-6)

    # Global topology: radius of gyration error + end-to-end
    def rg(X):
        c = X.mean(axis=0)
        return float(np.sqrt(((X - c) ** 2).sum(axis=1).mean()))

    rg_e, rg_p = rg(e_al), rg(p_fit)
    ee_e = float(np.linalg.norm(e_al[-1] - e_al[0]))
    ee_p = float(np.linalg.norm(p_fit[-1] - p_fit[0]))
    topo_err = abs(rg_p - rg_e) + abs(ee_p - ee_e) / PHI
    mode_scores["global_topology"] = topo_err

    if not seq_match and abs(len(pred_seq) - len(exp_seq)) > 5:
        mode_scores["sequence_alignment_mismatch"] = abs(len(pred_seq) - len(exp_seq))

    # Rank modes
    ranked = sorted(mode_scores.items(), key=lambda x: -x[1])
    primary = ranked[0][0] if ranked else "global_topology"

    return {
        "name": name,
        "pdb": pdb_file,
        "n": Ln,
        "pred_len": len(pred_seq),
        "exp_len": len(exp_seq),
        "seq_match": bool(seq_match),
        "rmsd_A": rmsd,
        "per_res_mean_A": float(per_res.mean()),
        "per_res_median_A": float(np.median(per_res)),
        "per_res_p90_A": float(np.percentile(per_res, 90)),
        "per_res_max_A": float(per_res.max()),
        "distance_bins": bin_stats,
        "native_contacts": len(native) if native else 0,
        "contact_mae_A": contact_mae,
        "top_L_precision": top_l_prec,
        "top_L_precision_structure": top_l_prec,
        "evidence_contact_diag": evidence_diag,
        "rg_exp_A": rg_e,
        "rg_pred_A": rg_p,
        "rg_err_A": abs(rg_p - rg_e),
        "end_to_end_err_A": abs(ee_p - ee_e),
        "hotspots": hotspots[:15],
        "mode_scores": mode_scores,
        "primary_error_mode": primary,
        "ranked_modes": [m for m, _ in ranked],
        "engine": pred.get("engine"),
        "full_law": pred.get("full_law"),
        "S_final": pred.get("S_final_observation"),
        "predict_ms": pred.get("predict_ms"),
    }


def build_fix_queue(protein_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Aggregate primary modes across proteins → prioritized fix queue."""
    counts: dict[str, int] = {}
    score_sum: dict[str, float] = {}
    for r in protein_rows:
        if r.get("error"):
            continue
        for m, sc in (r.get("mode_scores") or {}).items():
            counts[m] = counts.get(m, 0) + 1
            score_sum[m] = score_sum.get(m, 0.0) + float(sc)
        # primary vote
        pm = r.get("primary_error_mode")
        if pm:
            counts[pm] = counts.get(pm, 0) + 2  # weight primary

    queue = []
    for mode, cnt in sorted(counts.items(), key=lambda x: -x[1]):
        cat = ERROR_CATALOG.get(mode, {})
        queue.append(
            {
                "priority": len(queue) + 1,
                "mode": mode,
                "votes": cnt,
                "score_sum": score_sum.get(mode, 0.0),
                "meaning": cat.get("meaning", ""),
                "literature": cat.get("literature", ""),
                "fsot_handle": cat.get("fsot_handle", ""),
                "fix_class": cat.get("fix_class", "unknown"),
                "status": "open",
            }
        )
    # re-number priority by score_sum then votes
    queue.sort(key=lambda x: (-x["score_sum"], -x["votes"]))
    for i, q in enumerate(queue, 1):
        q["priority"] = i
    return queue


def main() -> int:
    print("=" * 64)
    print("ERROR MARGIN LOG — localize failure modes (FSOT residual discipline)")
    print("=" * 64)

    rows = []
    t0 = time.perf_counter()
    for fname, chain, name, note in SAMPLES:
        path = PDB_DIR / fname
        if not path.is_file():
            print(f"  SKIP {fname} missing")
            continue
        exp_seq, exp_xyz = parse_pdb_ca(path.read_text(encoding="utf-8", errors="replace"), chain)
        if len(exp_seq) < 20:
            print(f"  SKIP {fname} short")
            continue
        # CRITICAL: fold the experimental construct sequence (not UniProt polyprotein)
        pred = predict_ca_coords(exp_seq, rounds=24)
        rec = analyze_one(name, fname, exp_seq, exp_xyz, pred)
        rec["note"] = note
        rows.append(rec)
        ev = rec.get("evidence_contact_diag") or {}
        ev_p = ev.get("top_L_precision_evidence")
        print(
            f"  {fname:8s} n={rec.get('n')}  RMSD={rec.get('rmsd_A', float('nan')):6.2f} Å  "
            f"primary={rec.get('primary_error_mode')}  "
            f"contact_mae={rec.get('contact_mae_A')}  "
            f"topL_struct={rec.get('top_L_precision')}  topL_evid={ev_p}"
        )

    queue = build_fix_queue(rows)
    med_rmsd = float(np.median([r["rmsd_A"] for r in rows if "rmsd_A" in r])) if rows else None

    doc = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "doctrine": (
            "Error margin log: localize residual by structural mode, map literature, "
            "fix one mode at a time under full S=K(T1+T2+T3). Same residual discipline "
            "as multi-domain / MPCORB campaigns. Zero free parameters."
        ),
        "protocol": {
            "sequence_source": "experimental PDB chain (not UniProt polyprotein)",
            "metric": "Kabsch Cα RMSD + distance residual bins + contacts",
            "contact_cutoff_A": CONTACT_A,
        },
        "summary": {
            "n_proteins": len(rows),
            "median_rmsd_A": med_rmsd,
            "wall_s": time.perf_counter() - t0,
            "primary_modes_seen": list({r.get("primary_error_mode") for r in rows if r.get("primary_error_mode")}),
        },
        "proteins": rows,
        "fix_queue": queue,
        "error_catalog": ERROR_CATALOG,
        "next_action": queue[0] if queue else None,
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(doc, indent=2), encoding="utf-8")

    lines = [
        "# Error margin log",
        "",
        f"*Generated {doc['generated_at']}*",
        "",
        "## Protocol",
        "",
        "- Fold **experimental PDB sequence** (not UniProt polyprotein).",
        "- Kabsch Cα RMSD + **distance residual by sep bin** + native contact MAE + top-L precision.",
        "- Rank modes → **fix queue** with literature + FSOT handle.",
        "- Full law \(S=K(T_1+T_2+T_3)\) only; **0 free parameters**.",
        "",
        f"**Median RMSD (this set):** {med_rmsd} Å",
        "",
        "## Per protein",
        "",
        "| Protein | n | RMSD Å | Primary mode | Contact MAE Å | Top-L prec |",
        "|---------|--:|-------:|:-------------|--------------:|-----------:|",
    ]
    for r in rows:
        if r.get("error"):
            lines.append(f"| {r.get('name')} | — | err | {r.get('error')} | — | — |")
            continue
        lines.append(
            f"| {r['name']} | {r['n']} | {r['rmsd_A']:.2f} | `{r['primary_error_mode']}` | "
            f"{r['contact_mae_A'] if r['contact_mae_A'] is not None else '—'} | "
            f"{r['top_L_precision'] if r['top_L_precision'] is not None else '—'} |"
        )

    lines.extend(["", "## Fix queue (priority)", ""])
    for q in queue:
        lines.append(f"### {q['priority']}. `{q['mode']}` (votes={q['votes']})")
        lines.append("")
        lines.append(f"- **Meaning:** {q['meaning']}")
        lines.append(f"- **Literature:** {q['literature']}")
        lines.append(f"- **FSOT handle:** {q['fsot_handle']}")
        lines.append(f"- **Status:** {q['status']}")
        lines.append("")

    if queue:
        top = queue[0]
        lines.extend(
            [
                "## Next solve (do not skip)",
                "",
                f"**Mode:** `{top['mode']}`",
                "",
                top["fsot_handle"],
                "",
            ]
        )

    lines.append(f"Full JSON: `{OUT_JSON.relative_to(ROOT)}`")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")

    print("\n" + "=" * 64)
    print(f"Median RMSD: {med_rmsd} Å")
    print("FIX QUEUE:")
    for q in queue[:6]:
        print(f"  {q['priority']}. {q['mode']:28s} votes={q['votes']}  {q['fsot_handle'][:60]}...")
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
