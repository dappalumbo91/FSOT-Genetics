#!/usr/bin/env python3
"""Wet-lab / clinical / AlphaFold-style benchmark catalog (data only).

Curated cases with **experimental** anchors:
  - structure: UniProt + experimental PDB chain (wet lab crystal/NMR/cryo-EM)
  - AlphaFold DB available for UniProt accession (fair structure competitor)
  - variants: literature pathogenic drivers with experimental/clinical support

Zero free parameters — this is an evaluation *set*, not fitted weights.
Residue numbers are UniProt canonical unless noted.
"""

from __future__ import annotations

from typing import Any

# ---------------------------------------------------------------------------
# Structure panel: FSOT product vs AlphaFold DB vs experimental Cα
# identity_cap 0.95 fair H2H (exclude near-identical redeposits of same PDB)
# ---------------------------------------------------------------------------
STRUCTURE_CASES: list[dict[str, Any]] = [
    # --- Cancer / oncogenes (experimental structures exist) ---
    {
        "id": "p53_dbd",
        "category": "cancer",
        "name": "p53 DNA-binding domain",
        "uniprot": "P04637",
        "pdb": "1TUP",
        "chain": "A",
        "wetlab": "X-ray p53–DNA complex (Cho et al.)",
        "note": "hotspot R175/R248/R273 in structure",
    },
    {
        "id": "kras",
        "category": "cancer",
        "name": "KRAS G-domain",
        "uniprot": "P01116",
        "pdb": "4OBE",
        "chain": "A",
        "wetlab": "X-ray KRAS (GDP)",
        "note": "G12/G13/Q61 oncogenic sites",
    },
    {
        "id": "egfr_kinase",
        "category": "cancer",
        "name": "EGFR kinase domain",
        "uniprot": "P00533",
        "pdb": "2ITX",
        "chain": "A",
        "wetlab": "X-ray EGFR kinase",
        "note": "L858R / T790M drug sites in kinase",
    },
    {
        "id": "braf_kinase",
        "category": "cancer",
        "name": "BRAF kinase",
        "uniprot": "P15056",
        "pdb": "1UWH",
        "chain": "B",
        "wetlab": "X-ray BRAF kinase",
        "note": "V600E melanoma hotspot",
    },
    {
        "id": "abl1_kinase",
        "category": "cancer",
        "name": "ABL1 kinase (imatinib target)",
        "uniprot": "P00519",
        "pdb": "2HYY",
        "chain": "A",
        "wetlab": "X-ray ABL–imatinib",
        "note": "CML drug target",
    },
    {
        "id": "bcl2",
        "category": "cancer",
        "name": "BCL-2",
        "uniprot": "P10415",
        "pdb": "1G5M",
        "chain": "A",
        "wetlab": "NMR BCL-2",
        "note": "apoptosis / venetoclax target family",
    },
    # --- Vaccine / pathogen antigens ---
    {
        "id": "sars2_rbd",
        "category": "vaccine",
        "name": "SARS-CoV-2 spike RBD",
        "uniprot": "P0DTC2",
        "pdb": "6M0J",
        "chain": "E",
        "wetlab": "X-ray RBD–ACE2 (Lan et al. Nature 2020)",
        "note": "vaccine antigen / neutralizing Ab target",
    },
    {
        "id": "ha_h3",
        "category": "vaccine",
        "name": "Influenza HA (H3)",
        "uniprot": "P03437",
        "pdb": "4WE8",
        "chain": "A",
        "wetlab": "X-ray hemagglutinin",
        "note": "flu vaccine antigen",
    },
    {
        "id": "hiv_pr",
        "category": "drug",
        "name": "HIV-1 protease",
        "uniprot": "P03367",
        "pdb": "1HVR",
        "chain": "A",
        "wetlab": "X-ray HIV protease",
        "note": "classic antiviral drug target",
    },
    # --- Drug targets / medical enzymes ---
    {
        "id": "ace2",
        "category": "drug",
        "name": "Human ACE2 (peptidase domain)",
        "uniprot": "Q9BYF1",
        "pdb": "1R42",
        "chain": "A",
        "wetlab": "X-ray ACE2",
        "note": "SARS-CoV-2 host receptor / drug target",
    },
    {
        "id": "dhfr",
        "category": "drug",
        "name": "Human DHFR",
        "uniprot": "P00374",
        "pdb": "1U72",
        "chain": "A",
        "wetlab": "X-ray DHFR–methotrexate",
        "note": "classic antifolate drug target",
    },
    {
        "id": "cox2",
        "category": "drug",
        "name": "COX-2 (PTGS2)",
        "uniprot": "P35354",
        "pdb": "5KIR",
        "chain": "A",
        "wetlab": "X-ray COX-2",
        "note": "NSAID / celecoxib target",
    },
    # --- Controls (already strong product path) ---
    {
        "id": "ubiquitin",
        "category": "control",
        "name": "Ubiquitin",
        "uniprot": "P0CG47",
        "pdb": "1UBQ",
        "chain": "A",
        "wetlab": "X-ray ubiquitin (Vijay-Kumar)",
        "note": "classic fold control",
    },
    {
        "id": "lysozyme",
        "category": "control",
        "name": "Human lysozyme",
        "uniprot": "P61626",
        "pdb": "1LZ1",
        "chain": "A",
        "wetlab": "X-ray lysozyme",
        "note": "classic structural biology control",
    },
    {
        "id": "sod1",
        "category": "cancer",  # ALS — neurodegenerative medical
        "name": "SOD1",
        "uniprot": "P00441",
        "pdb": "2C9V",
        "chain": "A",
        "wetlab": "X-ray SOD1",
        "note": "ALS drivers; medical structure",
    },
]

# ---------------------------------------------------------------------------
# Variant panel with wet-lab / clinical labels (positive controls)
# label: pathogenic | benign_like | drug_resistance
# ---------------------------------------------------------------------------
VARIANT_CASES: list[dict[str, Any]] = [
    # TP53 — Cho / IARC hotspot literature
    {"gene": "TP53", "uniprot": "P04637", "pos": 175, "wt": "R", "mut": "H", "label": "pathogenic", "evidence": "IARC hotspot; structural Zn", "category": "cancer"},
    {"gene": "TP53", "uniprot": "P04637", "pos": 248, "wt": "R", "mut": "Q", "label": "pathogenic", "evidence": "DNA contact; IARC", "category": "cancer"},
    {"gene": "TP53", "uniprot": "P04637", "pos": 273, "wt": "R", "mut": "H", "label": "pathogenic", "evidence": "DNA contact; IARC", "category": "cancer"},
    {"gene": "TP53", "uniprot": "P04637", "pos": 72, "wt": "P", "mut": "R", "label": "benign_like", "evidence": "common polymorphism", "category": "cancer"},
    # KRAS — COSMIC classic
    {"gene": "KRAS", "uniprot": "P01116", "pos": 12, "wt": "G", "mut": "D", "label": "pathogenic", "evidence": "COSMIC codon 12", "category": "cancer"},
    {"gene": "KRAS", "uniprot": "P01116", "pos": 12, "wt": "G", "mut": "C", "label": "pathogenic", "evidence": "sotorasib-sensitive G12C", "category": "cancer"},
    {"gene": "KRAS", "uniprot": "P01116", "pos": 61, "wt": "Q", "mut": "H", "label": "pathogenic", "evidence": "switch II", "category": "cancer"},
    # EGFR — NSCLC wet-lab / trials
    {"gene": "EGFR", "uniprot": "P00533", "pos": 858, "wt": "L", "mut": "R", "label": "pathogenic", "evidence": "TKI-sensitive; crystal/clinic", "category": "cancer"},
    {"gene": "EGFR", "uniprot": "P00533", "pos": 790, "wt": "T", "mut": "M", "label": "drug_resistance", "evidence": "gatekeeper; TKI resistance", "category": "cancer"},
    # BRAF
    {"gene": "BRAF", "uniprot": "P15056", "pos": 600, "wt": "V", "mut": "E", "label": "pathogenic", "evidence": "melanoma; vemurafenib", "category": "cancer"},
    # CFTR — classic Mendelian wet-lab genetics
    {"gene": "CFTR", "uniprot": "P13569", "pos": 551, "wt": "G", "mut": "D", "label": "pathogenic", "evidence": "gating; ivacaftor responsive", "category": "drug"},
    {"gene": "CFTR", "uniprot": "P13569", "pos": 508, "wt": "F", "mut": "*", "label": "pathogenic", "evidence": "ΔF508 most common CF (del); scored as severe if present", "category": "drug", "skip_if_not_missense": True},
    # HBB
    {"gene": "HBB", "uniprot": "P68871", "pos": 7, "wt": "E", "mut": "V", "label": "pathogenic", "evidence": "sickle cell (HbS; UniProt pos 7 = Hb 6)", "category": "drug"},
    # SOD1 ALS
    {"gene": "SOD1", "uniprot": "P00441", "pos": 94, "wt": "G", "mut": "A", "label": "pathogenic", "evidence": "ALS familial", "category": "cancer"},
]


def by_category(cases: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list] = {}
    for c in cases:
        out.setdefault(c.get("category", "other"), []).append(c)
    return out
