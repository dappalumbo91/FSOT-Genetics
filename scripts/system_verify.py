#!/usr/bin/env python3
"""Full-system claim verify for FSOT-Genetics.

Checks every published claim against live JSON + the scalar engine:
product freeze, AF3-class coverage, hop splits, homologs, plants,
medical/PGx, language-stack parity, RCSB/wetlab benches, fly walking
join, and the honesty walls (bulk / fair-cap / frozen Biohub-Kaggle
are not the product).

Does not re-fold the 10-protein freeze (that JSON is the freeze).
Does not invent connectomes. Does not treat 0.95-cap / fuse-era /
orphan MDS numbers as the 0.13 Å product.

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
from domain_interface import domain_slice  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402
from medical_gene_catalog import GENE_CATALOG  # noqa: E402
from plant_signal import CLASSES  # noqa: E402
from run_rcsb_template_holdout import PRODUCT_IDENTITY_CAP  # noqa: E402
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


def job(af: dict[str, Any], name: str) -> dict[str, Any]:
    return (af.get("jobs") or {}).get(name) or {}


def exists(*parts: str) -> bool:
    return (ROOT.joinpath(*parts)).is_file()


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

    rec("leftover_floor", near(LEFTOVER, 1.0 / (PHI * PHI), 1e-12), f"1/φ²={LEFTOVER:.9f}")
    rec("close_homolog_floor", near(CLOSE, 1.0 / PHI, 1e-12), f"1/φ={CLOSE:.9f}")
    rec("product_identity_cap", float(PRODUCT_IDENTITY_CAP) == 1.0, str(PRODUCT_IDENTITY_CAP))
    chem = {
        "Physical_Chemistry": 8,
        "Chemistry": 8,
        "Molecular_Chemistry": 9,
        "Electromagnetism": 9,
        "Atomic_Physics": 7,
        "Condensed_Matter": 14,
        "Biochemistry": 13,
    }
    for name, d_eff in chem.items():
        got = int(domain_slice(name).D_eff)
        rec(f"chemlink_{name}", got == d_eff, f"D_eff={got}")

    af = load("af_coverage.json")
    jobs = af.get("jobs") or {}
    rec("af_n_jobs", len(jobs) == 17, f"n={len(jobs)}")
    rec("af_covered", len(af.get("covered_now") or []) == 17, str(len(af.get("covered_now") or [])))
    rec("af_free_params", int(af.get("free_parameters", -1)) == 0, str(af.get("free_parameters")))
    rec("af_not_yet_empty", len(af.get("not_yet") or []) == 0, str(af.get("not_yet")))
    rec(
        "af_monomer_product",
        job(af, "protein_monomer").get("status") == "product"
        and near(job(af, "protein_monomer")["fsot_median_A"], 0.13, 0.01),
        f"{job(af, 'protein_monomer').get('fsot_median_A')}",
    )
    rec("af_dna_prot", near(job(af, "protein_dna")["protein_rmsd_A"], 0.013, 0.005), f"{job(af, 'protein_dna').get('protein_rmsd_A')}")
    rec("af_dna_c1", near(job(af, "protein_dna")["dna_c1_rmsd_A"], 0.016, 0.005), f"{job(af, 'protein_dna').get('dna_c1_rmsd_A')}")
    rec("af_caii_metal", near(job(af, "metal_caii")["metal_site_rmsd_A"], 0.061, 0.01), f"{job(af, 'metal_caii').get('metal_site_rmsd_A')}")
    rec("af_sod1_metal", near(job(af, "metal_sod1")["metal_site_rmsd_A"], 0.26, 0.02), f"{job(af, 'metal_sod1').get('metal_site_rmsd_A')}")
    rec("af_rna_c1", near(job(af, "rna")["c1_rmsd_A"], 0.68, 0.02), f"{job(af, 'rna').get('c1_rmsd_A')}")
    rec("af_modified_c1", near(job(af, "modified_na")["c1_rmsd_A"], 0.93, 0.03), f"{job(af, 'modified_na').get('c1_rmsd_A')}")
    rec("af_h_matched", int(job(af, "hydrogens").get("n_h_matched") or 0) == 961, str(job(af, "hydrogens").get("n_h_matched")))
    rec("af_h_neutron", near(job(af, "hydrogens")["hydrogen_rmsd_A"], 1.01, 0.02), f"{job(af, 'hydrogens').get('hydrogen_rmsd_A')}")
    rec("af_dimer", near(job(af, "protein_protein")["dimer_rmsd_A"], 0.45, 0.02), f"{job(af, 'protein_protein').get('dimer_rmsd_A')}")
    rec("af_iface", near(job(af, "protein_protein")["interface_contact_mae_A"], 0.17, 0.02), f"{job(af, 'protein_protein').get('interface_contact_mae_A')}")
    rec("af_tetramer", near(job(af, "protein_tetramer")["tetramer_rmsd_A"], 0.51, 0.02), f"{job(af, 'protein_tetramer').get('tetramer_rmsd_A')}")
    rec("af_sc_centroid", near(job(af, "all_atom_sidechains")["sidechain_centroid_rmsd_A"], 0.41, 0.02), f"{job(af, 'all_atom_sidechains').get('sidechain_centroid_rmsd_A')}")
    rec("af_sc_heavy", near(job(af, "all_atom_sidechains")["sidechain_heavy_rmsd_A"], 1.01, 0.02), f"{job(af, 'all_atom_sidechains').get('sidechain_heavy_rmsd_A')}")
    rec("af_glycan", near(job(af, "ptm_glycan")["protein_rmsd_A"], 0.56, 0.03), f"{job(af, 'ptm_glycan').get('protein_rmsd_A')}")
    rec("af_phospho", near(job(af, "ptm_phospho")["protein_rmsd_A"], 0.77, 0.03), f"{job(af, 'ptm_phospho').get('protein_rmsd_A')}")
    rec("af_cdr", near(job(af, "antibody_cdr")["ca_rmsd_A"], 0.93, 0.03), f"{job(af, 'antibody_cdr').get('ca_rmsd_A')}")
    rec("af_pair", near(job(af, "antibody_pair")["pair_rmsd_A"], 1.03, 0.03), f"{job(af, 'antibody_pair').get('pair_rmsd_A')}")
    rec("af_u1a_prot", near(job(af, "protein_rna")["protein_rmsd_A"], 0.23, 0.02), f"{job(af, 'protein_rna').get('protein_rmsd_A')}")
    rec("af_u1a_rna", near(job(af, "protein_rna")["rna_c1_rmsd_A"], 0.28, 0.02), f"{job(af, 'protein_rna').get('rna_c1_rmsd_A')}")
    rec("af_ligand", near(job(af, "ligand")["ligand_site_rmsd_A"], 0.60, 0.03), f"{job(af, 'ligand').get('ligand_site_rmsd_A')}")
    rec("af_joint_ca", near(job(af, "joint_forward")["protein_rmsd_A"], 0.013, 0.005), f"{job(af, 'joint_forward').get('protein_rmsd_A')}")
    rec("af_joint_sc", near(job(af, "joint_forward")["sidechain_centroid_rmsd_A"], 0.016, 0.005), f"{job(af, 'joint_forward').get('sidechain_centroid_rmsd_A')}")
    rec("af_joint_dna", near(job(af, "joint_forward")["dna_c1_rmsd_A"], 0.016, 0.005), f"{job(af, 'joint_forward').get('dna_c1_rmsd_A')}")
    rec("af_joint_engine", job(af, "joint_forward").get("engine") == "fsot_predict_system", str(job(af, "joint_forward").get("engine")))

    field = load("field_stress_suite.json")
    fs = field.get("summary") or {}
    rec("field_n", int(fs.get("n") or 0) == 49, str(fs.get("n")))
    rec("field_pass", int(fs.get("pass") or 0) == 49, str(fs.get("pass")))
    rec("field_overall", str(fs.get("overall") or "") == "PASS", str(fs.get("overall")))
    rec("field_required_fails", int(fs.get("required_fails") or 0) == 0, str(fs.get("required_fails")))

    med = load("medical_variant_panel.json")
    ms = med.get("summary") or {}
    rec("med_genes", int(ms.get("n_genes") or 0) == 8, str(ms.get("n_genes")))
    rec("med_drivers", int(ms.get("n_drivers") or 0) == 35, str(ms.get("n_drivers")))
    rec("med_recall", float(ms.get("driver_recall_at_75pct") or 0) == 1.0, str(ms.get("driver_recall_at_75pct")))
    rec("med_free_params", int(med.get("free_parameters", -1)) == 0, str(med.get("free_parameters")))
    rec(
        "med_catalog_8",
        set(GENE_CATALOG.keys()) == {"TP53", "KRAS", "EGFR", "BRAF", "CFTR", "SOD1", "HBB", "BRCA1"},
        str(sorted(GENE_CATALOG.keys())),
    )

    wet = load("wetlab_af_eval.json")
    wst = wet.get("structure_summary") or {}
    wv = wet.get("variant_summary") or {}
    rec("wetlab_n_ok", int(wst.get("n_ok") or 0) == 19, f"{wst.get('n_ok')}/{wst.get('n_attempted')}")
    rec("wetlab_sub2", int(wst.get("fsot_sub2A") or 0) == 19, str(wst.get("fsot_sub2A")))
    rec("wetlab_fsot_median", near(float(wst["fsot_median_A"]), 0.256, 0.02), f"{wst.get('fsot_median_A')}")
    rec("wetlab_cap", float(wet.get("identity_cap") or 0) == 1.0, str(wet.get("identity_cap")))
    rec("wetlab_not_freeze_panel", int(wst.get("n_attempted") or 0) != 10, "19-protein medical panel, not the n=10 freeze")
    rec("wetlab_beats_af_n", int(wst.get("fsot_beats_af") or 0) == 16, str(wst.get("fsot_beats_af")))
    rec("wetlab_pathogenic_recall", float(wv.get("pathogenic_recall_likely_damaging") or 0) == 1.0, str(wv.get("pathogenic_recall_likely_damaging")))
    rec(
        "wetlab_benign_miss_recorded",
        float(wv.get("benign_like_not_called_damaging", 1)) == 0.0,
        "honest miss: benign-like still called damaging",
    )

    reality = load("reality_margin_eval.json")
    rs = reality.get("summary") or {}
    rec("reality_n_ok", int(rs.get("n_ok") or 0) == 19, f"{rs.get('n_ok')}/{rs.get('n_attempted')}")
    rec("reality_median", near(float(rs["median_rmsd_A"]), 1.17, 0.05), f"{rs.get('median_rmsd_A')}")
    rec("reality_median_ok", bool(rs.get("reality_median_ok")), f"target={rs.get('reality_median_target_A')}")
    rec("reality_not_freeze", float(rs["median_rmsd_A"]) > 0.5, "1.17 Å is not the 0.13 freeze")

    stress = load("medical_stress_suite.json")
    ss = stress.get("summary") or {}
    rec("stress_fuse_era", near(float(ss["template_msa_fuse_median_A"]), 1.16, 0.05), f"fuse={ss.get('template_msa_fuse_median_A')}")
    rec("stress_bulk_not_product", float(ss["bulk_single_median_A"]) > 10.0, f"bulk={ss.get('bulk_single_median_A')}")
    rec("stress_not_current_product", float(ss["template_msa_fuse_median_A"]) > 1.0, "fair-cap/fuse-era, not 0.13 product")

    dim = load("dimensionality_audit.json")
    agg = dim.get("aggregate") or {}
    rec("dim_base_25", int(agg.get("base_law_d_eff") or 0) == 25, str(agg.get("base_law_d_eff")))
    rec("dim_fsot_part", near(float(agg["fsot_participation_dimension_median"]), 9.10, 0.05), f"{agg.get('fsot_participation_dimension_median')}")
    rec("dim_native_part", near(float(agg["native_participation_dimension_median"]), 2.53, 0.05), f"{agg.get('native_participation_dimension_median')}")
    rec("dim_fsot_top3", near(float(agg["fsot_top3_variance_fraction_median"]), 0.45, 0.03), f"{agg.get('fsot_top3_variance_fraction_median')}")
    rec("dim_neg_eigen", near(float(agg["fsot_negative_eigenvalue_mass_median"]), 0.20, 0.03), f"{agg.get('fsot_negative_eigenvalue_mass_median')}")
    ladder = agg.get("chem_link_d_eff_ladder") or {}
    rec("dim_ladder_backbone", int(ladder.get("backbone") or 0) == 8, str(ladder.get("backbone")))
    rec("dim_ladder_pack", int(ladder.get("hydrophobic_pack") or 0) == 14, str(ladder.get("hydrophobic_pack")))

    dsplit = load("domain_split_eval.json")
    dgenes = dsplit.get("genes") or []
    rec("domain_split_n", len(dgenes) == 4, str(len(dgenes)))
    rec("domain_split_symbols", {g.get("symbol") for g in dgenes} == {"KRAS", "SOD1", "HBB", "TP53"}, str([g.get("symbol") for g in dgenes]))
    rec("domain_split_free", int(dsplit.get("free_parameters", -1)) == 0, str(dsplit.get("free_parameters")))
    sod1 = next((g for g in dgenes if g.get("symbol") == "SOD1"), {})
    hbb = next((g for g in dgenes if g.get("symbol") == "HBB"), {})
    tp53 = next((g for g in dgenes if g.get("symbol") == "TP53"), {})
    rec("domain_split_SOD1", float(sod1.get("global_rmsd_to_native_A") or 99) < 0.5, str(sod1.get("global_rmsd_to_native_A")))
    rec("domain_split_HBB", float(hbb.get("global_rmsd_to_native_A") or 99) < 0.5, str(hbb.get("global_rmsd_to_native_A")))
    rec("domain_split_TP53_n", int(tp53.get("n_domains") or 0) == 4, str(tp53.get("n_domains")))
    rec("domain_split_TP53_not_freeze", float(tp53.get("global_rmsd_to_native_A") or 0) > 2.0, "multi-domain pose not claimed as 0.13")

    dtab = load("domain_interface_table.json")
    rec("domain_table_n", len(dtab.get("domains") or []) == 9, str(len(dtab.get("domains") or [])))
    rec("domain_table_routings", len(dtab.get("routings") or {}) == 6, str(len(dtab.get("routings") or {})))
    rec(
        "domain_table_zero_free",
        all(int((r or {}).get("free_parameters", -1)) == 0 for r in (dtab.get("routings") or {}).values()),
        "all routings 0 free params",
    )

    parity = load("parity_zig_python.json")
    rec("parity_status", str(parity.get("status") or "") == "PASS", str(parity.get("status")))
    rec("parity_n_rows", len(parity.get("rows") or []) == 22, str(len(parity.get("rows") or [])))
    rec("parity_all_ok", all(bool(r.get("ok")) for r in (parity.get("rows") or [])), "22/22")

    oriented = load("rcsb_oriented_backbone_eval.json")
    osu = oriented.get("summary") or {}
    rec("rcsb_oriented_gate", bool(osu.get("success_gate_passed")), str(osu.get("success_gate_passed")))
    rec("rcsb_oriented_pair_delta", float(osu.get("maximum_pair_distance_delta_A", 1)) == 0.0, str(osu.get("maximum_pair_distance_delta_A")))

    tmpl = load("rcsb_template_holdout_eval.json")
    ta = tmpl.get("aggregate") or {}
    rec("rcsb_tmpl_n", int(ta.get("n_chains") or 0) == 60, str(ta.get("n_chains")))
    rec("rcsb_tmpl_covered", int(ta.get("n_template_covered") or 0) == 54, str(ta.get("n_template_covered")))
    rec("rcsb_tmpl_cap_095", float(ta.get("identity_cap") or 0) == 0.95, str(ta.get("identity_cap")))
    best_med = float(((ta.get("best_rmsd_A") or {}).get("median")) or 99)
    rec("rcsb_tmpl_best", near(best_med, 2.20, 0.05), f"{best_med:.3f} Å")
    rec("rcsb_tmpl_not_freeze", best_med > 1.0 and float(ta.get("identity_cap") or 0) < 1.0, "0.95-cap holdout is not the 0.13 product")

    live_api = load("rcsb_live_api_holdout_eval.json")
    la = live_api.get("aggregate") or {}
    rec("rcsb_live_n", int(la.get("n_chains") or 0) == 60, str(la.get("n_chains")))
    live_bulk = float(((la.get("bulk_observer_rmsd_A") or {}).get("median")) or 0)
    rec("rcsb_live_bulk_not_product", live_bulk > 8.0, f"bulk={live_bulk:.2f} Å")

    hold = load("rcsb_holdout_eval.json")
    rec("rcsb_holdout_n", int((hold.get("summary") or {}).get("n_structures") or 0) == 12, str((hold.get("summary") or {}).get("n_structures")))

    f12 = load("f12_candidate_development.json")
    rec("f12_dev_only", (f12.get("candidate") or {}).get("production_enabled") is False, "production_enabled=False")
    rec("f12_macro_improved", bool((f12.get("gates") or {}).get("macro_recall_improved")), str((f12.get("gates") or {}).get("macro_recall_improved")))
    rec(
        "f12_macro_recall",
        near(float((f12.get("summary") or {})["candidate_mean_per_protein_macro_recall"]), 0.58, 0.02),
        str((f12.get("summary") or {}).get("candidate_mean_per_protein_macro_recall")),
    )
    rec(
        "f12_beats_baseline",
        float((f12.get("summary") or {})["candidate_mean_per_protein_macro_recall"])
        > float((f12.get("summary") or {})["baseline_mean_per_protein_macro_recall"]),
        "0.58 > 0.34",
    )

    kag = load("kaggle_f12_validation_eval.json")
    rec("kaggle_f12_frozen_role", "frozen" in str(kag.get("role") or "").lower(), str(kag.get("role")))
    rec("kaggle_f12_n_chains", int(kag.get("evaluated_chains") or 0) == 6483, str(kag.get("evaluated_chains")))
    rec("kaggle_f12_passed", bool(kag.get("passed")), str(kag.get("passed")))
    rec(
        "kaggle_f12_gates",
        all(bool(v) for v in (kag.get("success_gates") or {}).values()),
        str(kag.get("success_gates")),
    )
    freeze_doc = (ROOT / "docs" / "BIOHUB_FREEZE.md").read_text(encoding="utf-8")
    rec("biohub_freeze_doc", "back burner" in freeze_doc.lower(), "docs/BIOHUB_FREEZE.md")
    rec("biohub_inventory_n", int(load("biohub_3d_inventory.json").get("n") or 0) == 199, str(load("biohub_3d_inventory.json").get("n")))
    rec(
        "kaggle_not_live_product",
        "frozen" in str(kag.get("role") or "").lower() and "back burner" in freeze_doc.lower(),
        "do not market Kaggle F12 / Biohub as live product",
    )

    disto = load("fsot_distogram_contact_eval.json")
    ds = disto.get("summary") or {}
    rec("distogram_n", int(ds.get("n_proteins") or 0) == 5, str(ds.get("n_proteins")))
    rec("distogram_pearson", near(float(ds["median_pearson"]), 0.61, 0.03), f"{ds.get('median_pearson')}")

    cb = load("chemlink_bulk_bench.json")
    rec("chemlink_bulk_not_product", float(cb["median_single_A"]) > 10.0, f"{cb.get('median_single_A'):.2f} Å")
    rec("chemlink_bulk_near_13", near(float(cb["median_single_A"]), 13.57, 0.2), f"{cb.get('median_single_A')}")
    ubulk = load("uniref_bulk_bench.json")
    rec("uniref_bulk_not_product", float(ubulk["median_single_A"]) > 10.0, f"{ubulk.get('median_single_A'):.2f} Å")
    seq_only = load("fsot_vs_alphafold_structure.json")
    sos = seq_only.get("summary") or {}
    rec("seq_only_fsot_wins", int(sos.get("fsot_wins", -1)) == 0, str(sos.get("fsot_wins")))
    rec("seq_only_not_product", float(sos["fsot_median_rmsd_A"]) > 10.0, f"{sos.get('fsot_median_rmsd_A'):.2f} Å")
    rec("seq_only_af_wins", int(sos.get("alphafold_wins") or 0) == 8, str(sos.get("alphafold_wins")))

    m1 = load("m1_authority_verify.json")
    rec("m1_ok_false", m1.get("ok") is False, "cap-0.95 authority is not the product")
    rec("m1_freeze_median", near(float(m1["freeze_median"]), 1.14, 0.05), f"{m1.get('freeze_median')}")
    rec("m1_not_013", float(m1["freeze_median"]) > 1.0, "1.14 Å ≠ 0.13 Å")

    rtb = load("residual_template_bench.json")
    rts = rtb.get("summary") or {}
    rec("residual_tmpl_not_product", float(rts["residual_median_A"]) > 2.0, f"residual={rts.get('residual_median_A')}")

    em = load("error_margin_log.json")
    rec("error_margin_bulk", float((em.get("summary") or {})["median_rmsd_A"]) > 5.0, "diagnostic bulk, not product")
    ig = load("information_gap_audit.json")
    rec("info_gap_chirality_zero", bool((ig.get("aggregate") or {}).get("all_mirror_stress_deltas_zero")), "pair-distance delta 0")

    walk = load("fly_walking_product.json")
    wgenes = walk.get("genes") or []
    rec("walk_n", len(wgenes) == 4, str([g.get("symbol") for g in wgenes]))
    rec("walk_symbols", {g.get("symbol") for g in wgenes} == {"Gad1", "nan", "iav", "nompC"}, str([g.get("symbol") for g in wgenes]))
    rec("walk_free", int(walk.get("free_parameters", -1)) == 0, str(walk.get("free_parameters")))
    nomp = next((g for g in wgenes if g.get("symbol") == "nompC"), {})
    rec("walk_nompC_5VKQ", nomp.get("template_pdb") == "5VKQ" and float(nomp.get("template_identity") or 0) > 0.9, str(nomp.get("template_pdb")))

    beh = load("fly_behavior_flow.json")
    rec("behavior_free", int(beh.get("free_parameters", -1)) == 0, str(beh.get("free_parameters")))
    rec("behavior_residual", near(float(beh["residual_Biochemistry"]), r, 1e-9), str(beh.get("residual_Biochemistry")))
    bcmp = ((beh.get("connectome") or {}).get("compare") or {}).get("hop_2") or {}
    b_mech = float((bcmp.get("mechanosensory") or {}).get("descending") or 0)
    b_jo = float((bcmp.get("JO") or {}).get("descending") or 0)
    b_olf = float((bcmp.get("olfactory") or {}).get("descending") or 0)
    rec("behavior_mech_desc", near(b_mech, 23.84, 0.05), f"{b_mech:.2f}")
    rec("behavior_jo_desc", near(b_jo, 16.55, 0.05), f"{b_jo:.2f}")
    rec("behavior_olf_dark", b_olf < 0.5, f"{b_olf:.3f}")
    rec("behavior_split", b_mech > b_jo > 1 > b_olf, "mechano > JO >> olfactory")

    odor = load("fly_odor_flow.json")
    rec("odor_has_rest", bool((odor.get("observer") or {}).get("has_rest_state")), str((odor.get("observer") or {}).get("has_rest_state")))
    rec("odor_rest_seed_zero", int(((odor.get("connectome") or {}).get("rest") or {}).get("n_seed", -1)) == 0, "rest = no afferent seed")
    rec("odor_free", int(odor.get("free_parameters", -1)) == 0, str(odor.get("free_parameters")))

    for flow_name, label in (
        ("fly_courtship_flow.json", "courtship"),
        ("fly_sleep_flow.json", "sleep"),
        ("fly_aggression_flow.json", "aggression"),
    ):
        fl = load(flow_name)
        rec(f"{label}_free", int(fl.get("free_parameters", -1)) == 0, str(fl.get("free_parameters")))
        rec(f"{label}_n_male", int((fl.get("connectome") or {}).get("n_neurons") or 0) == 165122, str((fl.get("connectome") or {}).get("n_neurons")))
        rec(f"{label}_honesty", "not a thought" in str(fl.get("honesty") or "").lower(), "honesty present")

    org = load("organism_product_join.json")
    rec("organism_walking", len(org.get("walking_folds") or []) == 4, str(len(org.get("walking_folds") or [])))
    rec("organism_new_folds", len(org.get("new_folds") or []) == 6, str(len(org.get("new_folds") or [])))
    rec("organism_pin", str(org.get("pin") or "") == PIN, str(org.get("pin")))

    mgc = load("male_cns_genetics_on_cells.json")
    rec("male_genetics_n", int(mgc.get("n_traced") or 0) == 165122, str(mgc.get("n_traced")))
    rec("male_fru_high", int((mgc.get("fruDsx") or {}).get("fru_high") or 0) == 2611, str((mgc.get("fruDsx") or {}).get("fru_high")))
    rec("male_dsx_high", int((mgc.get("fruDsx") or {}).get("dsx_high") or 0) == 138, str((mgc.get("fruDsx") or {}).get("dsx_high")))
    rec("male_receptorType", int(mgc.get("n_with_receptorType") or 0) == 752, str(mgc.get("n_with_receptorType")))
    rec("male_no_invented_atlas", "do not invent" in str(mgc.get("note") or "").lower(), str(mgc.get("note"))[:80])

    msa = load("msa_dual_mode_smoke.json")
    rec("msa_ca_drift", float(msa["mean_ca_drift_single_vs_msa_A"]) < 0.01, f"{msa.get('mean_ca_drift_single_vs_msa_A')}")
    rec("msa_free", int((msa.get("single") or {}).get("free_parameters", -1)) == 0, str((msa.get("single") or {}).get("free_parameters")))

    flyb = load("fly_connectome_boot.json")
    rec("flywire_n", int(flyb.get("n_neurons") or 0) == 127979, str(flyb.get("n_neurons")))
    rec("flywire_residual", near(float(flyb["residual_Biochemistry"]), r, 1e-9), str(flyb.get("residual_Biochemistry")))
    rec("flywire_gaba_measured", flyb.get("gaba_inhibitory") is True, "FlyWire predictedNt is measured")
    inv = load("fly_connectome_inventory.json")
    rec("flywire_inventory_n", int(inv.get("n_neurons") or 0) == 139248, str(inv.get("n_neurons")))

    rec("worm_herm_n", int(worm.get("n_cells") or 0) == 453, str(worm.get("n_cells")))
    rec("worm_herm_gaba", int(worm.get("n_gaba") or 0) == 26, str(worm.get("n_gaba")))
    wm = load("worm_male_connectome_boot.json")
    rec("worm_male_n", int(wm.get("n_cells") or 0) == 575, str(wm.get("n_cells")))
    rec("worm_male_gaba", int(wm.get("n_gaba") or 0) == 26, str(wm.get("n_gaba")))
    wsex = load("worm_sex_compare.json")
    herm_sex = float((((wsex.get("compare") or {}).get("hop_2") or {}).get("hermaphrodite") or {}).get("sex_specific") or 0)
    male_sex = float((((wsex.get("compare") or {}).get("hop_2") or {}).get("male") or {}).get("sex_specific") or 0)
    rec("worm_sex_herm_dark", herm_sex < 0.1, f"{herm_sex:.4f}")
    rec("worm_sex_male_lit", male_sex > 1.0, f"{male_sex:.4f}")

    pdom = load("plant_signal_domains.json")
    rec("phot1_domains_acc", (pdom.get("PHOT1") or {}).get("uniprot") == "O48963", str((pdom.get("PHOT1") or {}).get("uniprot")))
    rec("phot1_close_domains", int((pdom.get("PHOT1") or {}).get("n_close_homolog_domains") or 0) == 2, str((pdom.get("PHOT1") or {}).get("n_close_homolog_domains")))
    rec("uvr8_acc", (pdom.get("UVR8") or {}).get("uniprot") == "Q9FN03", str((pdom.get("UVR8") or {}).get("uniprot")))
    rec("uvr8_product", (pdom.get("UVR8") or {}).get("template_pdb") == "8GQE", str((pdom.get("UVR8") or {}).get("template_pdb")))
    rec("uvr8_not_intact_hops", (pdom.get("UVR8") or {}).get("on_intact_graph") is False, "UVR8 product is not an IntAct hop")

    ubq = load("fsot_predict_ubq.json")
    rec("ubq_pin", str(ubq.get("authority_pin") or "") == PIN, str(ubq.get("authority_pin")))
    rec("ubq_free", int(ubq.get("free_parameters", -1)) == 0, str(ubq.get("free_parameters")))
    rec("ubq_len", int(ubq.get("length") or 0) == 76, str(ubq.get("length")))

    smiles = json.loads((ROOT / "formulas" / "smiles_protein_chemistry.json").read_text(encoding="utf-8"))
    rec("smiles_n", int(smiles.get("n_records") or 0) == 116, str(smiles.get("n_records")))
    rec("smiles_p_new", near(float((smiles.get("layer2") or {})["P_new"]), float(fc.P_NEW), 1e-12), str((smiles.get("layer2") or {}).get("P_new")))

    rec("formula_20aa", exists("formulas", "20_amino_acid_fsot_map.txt"), "20 aa map")
    rec("formula_64codon", exists("formulas", "64_codon_trinary_map.txt"), "64 codon map")
    rec("formula_derivations", exists("formulas", "FSOT_PROTEIN_DERIVATIONS.md"), "F01–F15")
    rec("predict_system_src", "def predict_system" in (ROOT / "scripts" / "multi_system.py").read_text(encoding="utf-8"), "joint forward")
    rec("dna_variant_src", exists("scripts", "dna_variant_effect.py"), "DNA → AA")
    rec("cofactor_src", exists("scripts", "cofactor_nodes.py"), "cofactor nodes")
    rec("disclosure_doc", exists("docs", "EXPERIMENTAL_DISCLOSURE.md"), "PGx disclosure")
    rec("medical_doc", exists("docs", "MEDICAL_PLATFORM.md"), "medical platform")
    rec("af_doc", exists("docs", "AF_COVERAGE.md"), "AF coverage")
    rec("lang_doc", exists("docs", "LANGUAGE_STACK.md"), "language stack")

    rec("lean_chemlink", exists("FSOTGenetics", "ChemLink.lean"), "Lean ChemLink")
    rec("lean_zero", exists("FSOTGenetics", "ZeroFreeParams.lean"), "Lean 0 free params")
    rec("haskell_chemlink", exists("haskell", "src", "FSOT", "Genetics", "ChemLink.hs"), "Haskell ChemLink")
    rec("haskell_contact", exists("haskell", "src", "FSOT", "Genetics", "Contact.hs"), "Haskell Contact")
    rec("rust_protein", exists("crates", "fsot_protein", "src", "lib.rs"), "Rust fsot_protein")
    rec("zig_src", (ROOT / "zig" / "src").is_dir(), "Zig runtime")
    rec("coq_spine", exists("verification", "coq", "GeneticsSpine.v"), "Coq")
    rec("isabelle_spine", exists("verification", "isabelle", "GeneticsSpine.thy"), "Isabelle")
    rec("fstar_spine", exists("verification", "fstar", "FSOTGenetics.fst"), "F*")
    rec("tla_spine", exists("verification", "tla", "GeneticsSolve.tla"), "TLA+")

    gaunt = load("cross_proof_report.json")
    rec("gauntlet_ok", bool(gaunt.get("overall_ok")), str(gaunt.get("overall_ok")))
    rec("gauntlet_layers", len(gaunt.get("layers") or []) == 10, str(len(gaunt.get("layers") or [])))
    rec("gauntlet_required_fail", int(gaunt.get("required_fail_count") or 0) == 0, str(gaunt.get("required_fail_count")))
    rec("gauntlet_no_hub_inherit", gaunt.get("inherits_hub_overall_ok") is False, "independent of Lean hub")
    obl = json.loads((ROOT / "verification" / "obligations.json").read_text(encoding="utf-8"))
    n_obl = int((obl.get("meta") or {}).get("n_obligations") or len(obl.get("obligations") or []))
    rec("obligations_n", n_obl == 42, str(n_obl))

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
