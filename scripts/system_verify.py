#!/usr/bin/env python3
"""Full-system claim verify for FSOT-Genetics.

Checks every numbered claim we publish (product Å, hop splits, homologs,
plants, analog pointer, PGx, pin) against live JSON + the scalar engine.
Does not re-fold the 10-protein freeze (that JSON is the freeze).
Does not invent connectomes.

  python scripts/system_verify.py
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402
from plant_signal import CLASSES  # noqa: E402
from trinary_syntax import aa_pair_weight, uniqueness_report  # noqa: E402

DATA = ROOT / "data"
PHI = float(fc.PHI)
CLOSE = 1.0 / PHI
LEFTOVER = 1.0 / (PHI * PHI)
PIN = "D1D38A"

rows: list[dict[str, Any]] = []
failed: list[str] = []


def load(name: str) -> dict[str, Any]:
    p = DATA / name
    if not p.is_file():
        raise FileNotFoundError(p)
    return json.loads(p.read_text(encoding="utf-8"))


def rec(name: str, ok: bool, detail: str) -> None:
    rows.append({"name": name, "ok": ok, "detail": detail})
    tag = "OK  " if ok else "FAIL"
    print(f"{tag}  {name}  {detail}", flush=True)
    if not ok:
        failed.append(name)


def near(a: float, b: float, eps: float) -> bool:
    return abs(float(a) - float(b)) <= eps


def hop(doc: dict[str, Any], n: int, program: str, field: str) -> float:
    recd = ((doc.get("compare") or {}).get(f"hop_{n}") or {}).get(program) or {}
    return float(recd.get(field) or 0.0)


def main() -> int:
    print("=" * 64)
    print("FSOT-Genetics system-wide claim verify")
    print(f"  root = {ROOT}")
    print("=" * 64)

    raw = (ROOT / "vendor" / "fsot_compute.py").read_bytes()
    sha = hashlib.sha256(raw).hexdigest().upper()
    rec("pin", sha.startswith(PIN), f"sha={sha[:16]}… bytes={len(raw)}")
    cert = json.loads((ROOT / "vendor" / "fsot_compute_AUTHORITY_PIN.json").read_text(encoding="utf-8"))
    rec(
        "pin_cert",
        str(cert.get("certificate_authority") or "").upper() == sha,
        "certificate matches file",
    )

    r = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))
    rec("residual_Biochemistry", near(r, 1.091958932349586, 1e-9), f"r={r:.9f}")
    rec("seeds", float(fc.PHI) > 1.6 and float(fc.PI) > 3.1, f"φ={float(fc.PHI):.6f} π={float(fc.PI):.6f}")

    u = uniqueness_report()
    rec("trinary_20_unique", bool(u.get("all_unique")), str(u.get("n_aa")))
    w = aa_pair_weight("F", "W", 8)
    rec("pair_FW_d8", near(w, 1.6700524010190179, 1e-9), f"w={w:.9f}")

    prod = load("product_vs_alphafold.json")
    s = prod.get("summary") or {}
    rec("product_n", int(s.get("n") or 0) == 10, f"n={s.get('n')}")
    rec("product_median", near(float(s["fsot_product_median_A"]), 0.13, 0.01), f"{s['fsot_product_median_A']:.4f} Å")
    rec("af_median", near(float(s["alphafold_median_A"]), 0.47, 0.01), f"{s['alphafold_median_A']:.4f} Å")
    rec(
        "product_beats_af",
        float(s["fsot_product_median_A"]) < float(s["alphafold_median_A"]),
        "0.13 < 0.47",
    )
    rec("product_sub2A", int(s.get("product_sub2A") or 0) == 10, str(s.get("product_sub2A")))
    rec("product_free_params", int(s.get("free_parameters", -1)) == 0, str(s.get("free_parameters")))
    rec(
        "bulk_not_product",
        float(s["fsot_bulk_median_A"]) > 10.0,
        f"bulk={s['fsot_bulk_median_A']:.2f} Å",
    )
    n_win = sum(
        1
        for r0 in prod.get("results") or []
        if float(r0.get("fsot_product_rmsd_A") or 99) < float(r0.get("alphafold_rmsd_A") or 0)
    )
    rec("product_beats_af_each", n_win == 10, f"{n_win}/10")

    male = load("male_cns_boot.json")
    banc = load("banc_connectome_boot.json")
    larva = load("larva_connectome_boot.json")
    rec("male_n", int(male["n_neurons"]) == 165122, str(male["n_neurons"]))
    rec("male_gaba", int(male["n_gaba"]) == 22055, str(male["n_gaba"]))
    rec("banc_n", int(banc["n_neurons"]) == 175401, str(banc["n_neurons"]))
    male_vnc = hop(male, 2, "vnc_sensory", "vnc_motor")
    male_jo = hop(male, 2, "JO", "vnc_motor")
    male_olf = hop(male, 2, "olfactory", "vnc_motor")
    banc_vnc = hop(banc, 2, "vnc_sensory", "vnc_motor")
    banc_jo = hop(banc, 2, "JO", "vnc_motor")
    banc_olf = hop(banc, 2, "olfactory", "vnc_motor")
    rec("male_vnc", near(male_vnc, 10.74, 0.02), f"{male_vnc:.4f}")
    rec("male_jo", near(male_jo, 7.77, 0.02), f"{male_jo:.4f}")
    rec("male_olf_dark", male_olf < 0.01, f"{male_olf:.6f}")
    rec("male_split", male_vnc > male_jo > 1 > male_olf, "VNC > JO >> olfactory")
    rec("banc_vnc", near(banc_vnc, 10.16, 0.02), f"{banc_vnc:.4f}")
    rec("banc_jo", near(banc_jo, 6.64, 0.05), f"{banc_jo:.4f}")
    rec("banc_olf_dark", banc_olf < 0.01, f"{banc_olf:.6f}")
    rec("banc_split", banc_vnc > banc_jo > 1 > banc_olf, "VNC > JO >> olfactory")
    larva_m = hop(larva, 2, "mechano", "DN_VNC")
    larva_o = hop(larva, 2, "olfactory", "DN_VNC")
    rec("larva_split", larva_m > 10 and larva_o < 1, f"mech={larva_m:.2f} olf={larva_o:.2f}")

    hemi_p = DATA / "hemibrain_connectome_boot.json"
    if hemi_p.is_file():
        hemi = load("hemibrain_connectome_boot.json")
        hj = hop(hemi, 2, "JO", "descending")
        ho = hop(hemi, 2, "olfactory", "descending")
        rec("hemibrain_split", hj > 1 and ho < 0.5, f"JO desc={hj:.2f} olf={ho:.3f}")
        rec("hemibrain_no_dng29", int(hemi.get("n_dng29") or 0) == 0, "DNg29 absent")
        rec("hemibrain_unsigned", int(hemi.get("n_gaba") or 0) == 0, "no invented GABA")
    else:
        rec("hemibrain_json", False, "missing data/hemibrain_connectome_boot.json")

    worm = load("worm_connectome_boot.json")
    ciona = load("ciona_connectome_boot.json")
    platy = load("platynereis_connectome_boot.json")
    rec("worm_nonempty", int(worm.get("n_cells") or 0) > 0, str(worm.get("n_cells")))
    rec("ciona_nonempty", int(ciona.get("n_cells") or 0) > 0, str(ciona.get("n_cells")))
    rec("platy_nonempty", int(platy.get("n_cells") or 0) > 0, str(platy.get("n_cells")))
    rec("ciona_unsigned", ciona.get("gaba_inhibitory") is False, str(ciona.get("gaba_inhibitory")))
    rec("platy_unsigned", platy.get("gaba_inhibitory") is False, str(platy.get("gaba_inhibitory")))

    homo = load("homolog_correspondence.json")
    rec("homolog_measured", int(homo.get("n_measured_homologs") or 0) == 26, str(homo.get("n_measured_homologs")))
    rec("homolog_folded", int(homo.get("n_folded") or 0) == 26, str(homo.get("n_folded")))
    misses = (homo.get("miss_audit") or {}).get("true_1to1_misses") or []
    rec("homolog_true_miss", len(misses) == 1, str(misses))
    analog = load("analog_pointer.json")
    rec("analog_jobs", len(analog.get("jobs") or []) == 2, str(len(analog.get("jobs") or [])))

    plant = load("plant_product.json")
    rec("plant_arabidopsis", int(plant.get("n_folded") or 0) == 6, str(plant.get("n_folded")))
    ph = load("plant_homolog.json")
    rec("plant_crop_folded", int(ph.get("n_folded") or 0) == 24, str(ph.get("n_folded")))
    rec("plant_crop_miss", int(ph.get("n_miss") or 0) == 0, str(ph.get("n_miss")))

    sig = load("plant_signal_boot.json")
    rec("plant_signal_nodes", int(sig.get("n_nodes") or 0) == 8329, str(sig.get("n_nodes")))
    rec("plant_signal_edges", int(sig.get("n_undirected_edges") or 0) == 39664, str(sig.get("n_undirected_edges")))
    light1 = ((sig.get("compare") or {}).get("hop_1") or {}).get("photoreceptor") or {}
    aba2 = ((sig.get("compare") or {}).get("hop_2") or {}).get("aba") or {}
    rec(
        "plant_light_PIF3",
        (light1.get("top") or {}).get("gene") == "PIF3" and float(light1.get("light_tf") or 0) > 0.99,
        str((light1.get("top") or {}).get("gene")),
    )
    rec(
        "plant_slac_unlit",
        float(aba2.get("stomata_channel") or 0) < 0.01,
        f"aba hop-2 stomata={aba2.get('stomata_channel')}",
    )
    rec("plant_classes_n", len(CLASSES) == 22, str(len(CLASSES)))

    sp = load("plant_signal_product.json")
    rec("plant_signal_product_n", int(sp.get("n_folded") or 0) >= 8, str(sp.get("n_folded")))
    pif3 = next((g for g in sp.get("genes") or [] if g.get("symbol") == "PIF3"), {})
    rec(
        "pif3_not_human",
        pif3.get("uniprot") == "O80536",
        f"PIF3 acc={pif3.get('uniprot')} mode={pif3.get('structure_mode')}",
    )
    phot1 = next((g for g in sp.get("genes") or [] if g.get("symbol") == "PHOT1"), {})
    rec("phot1_O48963", phot1.get("uniprot") == "O48963", str(phot1.get("uniprot")))

    pgx = load("experimental_pgx.json")
    rec("pgx_n", int(pgx.get("n") or 0) == 10, str(pgx.get("n")))
    rec("pgx_concordant", int(pgx.get("n_concordant") or 0) == 10, str(pgx.get("n_concordant")))
    rec("pgx_disclosure", "EXPERIMENTAL" in str(pgx.get("disclosure") or ""), "disclosure present")

    out = {
        "pin": PIN,
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3)",
        "n": len(rows),
        "n_fail": len(failed),
        "overall_ok": len(failed) == 0,
        "failed": failed,
        "checks": rows,
    }
    dest = DATA / "system_verify.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("=" * 64)
    print(f"overall_ok={out['overall_ok']}  n={out['n']}  fail={out['n_fail']}")
    print(f"wrote {dest}")
    return 0 if out["overall_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
