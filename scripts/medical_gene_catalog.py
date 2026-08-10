#!/usr/bin/env python3
"""Curated multi-gene medical catalog for FSOT variant + domain work.

Residue numbers are **UniProt canonical** (1-based, includes initiator Met).
Drivers are well-known pathogenic missense (or classic Mendelian) used as
positive controls — not a clinical database dump. Synonymous / mild controls
are included where available for calibration.

Zero free parameters: this file is *data*, not fitted weights.
"""

from __future__ import annotations

from typing import Any

# Each gene:
#   uniprot, symbol, name, pfam (primary MSA for conservation),
#   structure_pdb/chain (optional experimental chain for sanity),
#   drivers: list of {hgvs_p, pos, wt, mut, note}
#   dna_examples: optional list of DNA-routed variants (pos, wt_codon, cpos, alt, hgvs_c)
#   domains: optional static Pfam ranges (start,end inclusive UniProt) — runtime
#            InterPro fetch preferred when network available

GENE_CATALOG: dict[str, dict[str, Any]] = {
    "TP53": {
        "symbol": "TP53",
        "name": "Tumor protein p53",
        "uniprot": "P04637",
        "pfam": "PF00870",  # DNA-binding domain
        "structure_pdb": "1TUP",
        "structure_chain": "A",
        "indication": "cancer hotspot / Li-Fraumeni",
        "drivers": [
            {"hgvs_p": "p.R175H", "pos": 175, "wt": "R", "mut": "H", "note": "structural Zn node"},
            {"hgvs_p": "p.G245S", "pos": 245, "wt": "G", "mut": "S", "note": "structural"},
            {"hgvs_p": "p.R248Q", "pos": 248, "wt": "R", "mut": "Q", "note": "DNA contact"},
            {"hgvs_p": "p.R248W", "pos": 248, "wt": "R", "mut": "W", "note": "DNA contact"},
            {"hgvs_p": "p.R249S", "pos": 249, "wt": "R", "mut": "S", "note": "structural"},
            {"hgvs_p": "p.R273H", "pos": 273, "wt": "R", "mut": "H", "note": "DNA contact"},
            {"hgvs_p": "p.R282W", "pos": 282, "wt": "R", "mut": "W", "note": "structural"},
        ],
        "controls": [
            {"hgvs_p": "p.P72R", "pos": 72, "wt": "P", "mut": "R", "note": "common polymorphism (mild)"},
        ],
        "dna_examples": [
            (175, "CGC", 1, "A", "c.524G>A"),
            (245, "GGC", 0, "A", "c.733G>A"),
            (248, "CGG", 0, "T", "c.742C>T"),
            (273, "CGT", 1, "A", "c.818G>A"),
            (282, "CGG", 0, "T", "c.844C>T"),
            (248, "CGG", 2, "A", "c.744G>A"),  # synonymous control
        ],
        "domains_static": [
            {"pfam": "PF08563", "name": "TAD1", "start": 6, "end": 30},
            {"pfam": "PF18521", "name": "TAD2", "start": 35, "end": 59},
            {"pfam": "PF00870", "name": "DBD", "start": 100, "end": 288},
            {"pfam": "PF07710", "name": "TET", "start": 319, "end": 357},
        ],
    },
    "KRAS": {
        "symbol": "KRAS",
        "name": "GTPase KRas",
        "uniprot": "P01116",
        "pfam": "PF00071",  # Ras family
        "structure_pdb": "4OBE",
        "structure_chain": "A",
        "indication": "oncogene / MAPK",
        "drivers": [
            {"hgvs_p": "p.G12D", "pos": 12, "wt": "G", "mut": "D", "note": "codon 12 classic"},
            {"hgvs_p": "p.G12V", "pos": 12, "wt": "G", "mut": "V", "note": "codon 12 classic"},
            {"hgvs_p": "p.G12C", "pos": 12, "wt": "G", "mut": "C", "note": "codon 12 (sotorasib)"},
            {"hgvs_p": "p.G13D", "pos": 13, "wt": "G", "mut": "D", "note": "codon 13"},
            {"hgvs_p": "p.Q61H", "pos": 61, "wt": "Q", "mut": "H", "note": "switch II"},
            {"hgvs_p": "p.Q61L", "pos": 61, "wt": "Q", "mut": "L", "note": "switch II"},
        ],
        "controls": [],
        "dna_examples": [],
        "domains_static": [
            {"pfam": "PF00071", "name": "Ras", "start": 5, "end": 164},
        ],
    },
    "EGFR": {
        "symbol": "EGFR",
        "name": "Epidermal growth factor receptor",
        "uniprot": "P00533",
        "pfam": "PF07714",  # protein kinase domain (L858R lives here)
        "structure_pdb": "2ITX",
        "structure_chain": "A",
        "indication": "NSCLC / kinase inhibitors",
        "drivers": [
            {"hgvs_p": "p.G719S", "pos": 719, "wt": "G", "mut": "S", "note": "exon 18"},
            {"hgvs_p": "p.T790M", "pos": 790, "wt": "T", "mut": "M", "note": "gatekeeper resistance"},
            {"hgvs_p": "p.L858R", "pos": 858, "wt": "L", "mut": "R", "note": "exon 21 classic"},
            {"hgvs_p": "p.L861Q", "pos": 861, "wt": "L", "mut": "Q", "note": "exon 21"},
        ],
        "controls": [],
        "dna_examples": [],
        "domains_static": [
            {"pfam": "PF01030", "name": "Recep_L_domain", "start": 57, "end": 168},
            {"pfam": "PF00757", "name": "Furin-like", "start": 177, "end": 338},
            {"pfam": "PF14843", "name": "GF_recep_IV", "start": 505, "end": 637},
            {"pfam": "PF07714", "name": "Pkinase_Tyr", "start": 712, "end": 968},
        ],
    },
    "BRAF": {
        "symbol": "BRAF",
        "name": "Serine/threonine-protein kinase B-raf",
        "uniprot": "P15056",
        "pfam": "PF07714",
        "structure_pdb": "1UWH",
        "structure_chain": "B",
        "indication": "melanoma / MAPK",
        "drivers": [
            {"hgvs_p": "p.V600E", "pos": 600, "wt": "V", "mut": "E", "note": "activation loop classic"},
            {"hgvs_p": "p.V600K", "pos": 600, "wt": "V", "mut": "K", "note": "activation loop"},
            {"hgvs_p": "p.G469A", "pos": 469, "wt": "G", "mut": "A", "note": "P-loop"},
        ],
        "controls": [],
        "dna_examples": [],
        "domains_static": [
            {"pfam": "PF00069", "name": "Pkinase", "start": 457, "end": 717},
        ],
    },
    "CFTR": {
        "symbol": "CFTR",
        "name": "Cystic fibrosis transmembrane conductance regulator",
        "uniprot": "P13569",
        "pfam": "PF00664",  # ABC transporter
        "structure_pdb": "5UAK",
        "structure_chain": "A",
        "indication": "cystic fibrosis",
        "drivers": [
            {"hgvs_p": "p.G551D", "pos": 551, "wt": "G", "mut": "D", "note": "gating (ivacaftor)"},
            {"hgvs_p": "p.R117H", "pos": 117, "wt": "R", "mut": "H", "note": "mild/variable"},
            {"hgvs_p": "p.N1303K", "pos": 1303, "wt": "N", "mut": "K", "note": "NBD2 classic"},
            {"hgvs_p": "p.G542X", "pos": 542, "wt": "G", "mut": "*", "note": "nonsense (special)"},
        ],
        "controls": [],
        "dna_examples": [],
        "domains_static": [
            {"pfam": "PF00664", "name": "ABC_membrane", "start": 81, "end": 350},
            {"pfam": "PF00005", "name": "ABC_tran NBD1", "start": 389, "end": 670},
            {"pfam": "PF00664", "name": "ABC_membrane2", "start": 859, "end": 1155},
            {"pfam": "PF00005", "name": "ABC_tran NBD2", "start": 1207, "end": 1436},
        ],
    },
    "SOD1": {
        "symbol": "SOD1",
        "name": "Superoxide dismutase [Cu-Zn]",
        "uniprot": "P00441",
        "pfam": "PF00080",
        "structure_pdb": "2C9V",
        "structure_chain": "A",
        "indication": "ALS",
        "drivers": [
            # UniProt numbering includes Met1; classic "A4V" = p.A5V, "G93A" = p.G94A
            {"hgvs_p": "p.A5V", "pos": 5, "wt": "A", "mut": "V", "note": "A4V historic / aggressive ALS"},
            {"hgvs_p": "p.G94A", "pos": 94, "wt": "G", "mut": "A", "note": "G93A historic mouse model"},
            {"hgvs_p": "p.H47R", "pos": 47, "wt": "H", "mut": "R", "note": "Cu ligand"},
            {"hgvs_p": "p.G38R", "pos": 38, "wt": "G", "mut": "R", "note": "ALS"},
        ],
        "controls": [],
        "dna_examples": [],
        "domains_static": [
            {"pfam": "PF00080", "name": "Sod_Cu", "start": 2, "end": 150},
        ],
    },
    "HBB": {
        "symbol": "HBB",
        "name": "Hemoglobin subunit beta",
        "uniprot": "P68871",
        "pfam": "PF00042",
        "structure_pdb": "1A3N",
        "structure_chain": "B",
        "indication": "sickle cell / hemoglobinopathy",
        "drivers": [
            # UniProt: sickle is E7V (historic E6V after Met removal)
            {"hgvs_p": "p.E7V", "pos": 7, "wt": "E", "mut": "V", "note": "sickle cell HbS"},
            {"hgvs_p": "p.E7K", "pos": 7, "wt": "E", "mut": "K", "note": "HbC"},
            {"hgvs_p": "p.E122Q", "pos": 122, "wt": "E", "mut": "Q", "note": "HbD Punjab"},
        ],
        "controls": [],
        "dna_examples": [],
        "domains_static": [
            {"pfam": "PF00042", "name": "Globin", "start": 3, "end": 146},
        ],
    },
    "BRCA1": {
        "symbol": "BRCA1",
        "name": "Breast cancer type 1 susceptibility protein",
        "uniprot": "P38398",
        "pfam": "PF00533",  # BRCA1 C-terminus (BRCT) — RING is PF00097
        "structure_pdb": "1T15",
        "structure_chain": "A",
        "indication": "hereditary breast/ovarian cancer",
        "drivers": [
            {"hgvs_p": "p.C61G", "pos": 61, "wt": "C", "mut": "G", "note": "RING Zn finger pathogenic"},
            {"hgvs_p": "p.C64Y", "pos": 64, "wt": "C", "mut": "Y", "note": "RING Zn finger"},
            {"hgvs_p": "p.R1699W", "pos": 1699, "wt": "R", "mut": "W", "note": "BRCT pathogenic"},
            {"hgvs_p": "p.A1708E", "pos": 1708, "wt": "A", "mut": "E", "note": "BRCT pathogenic"},
            {"hgvs_p": "p.M1775R", "pos": 1775, "wt": "M", "mut": "R", "note": "BRCT pathogenic"},
        ],
        "controls": [],
        "dna_examples": [],
        "domains_static": [
            {"pfam": "PF00097", "name": "zf-C3HC4 RING", "start": 24, "end": 64},
            {"pfam": "PF00533", "name": "BRCT", "start": 1642, "end": 1736},
            {"pfam": "PF00533", "name": "BRCT2", "start": 1756, "end": 1855},
        ],
    },
}


def list_genes() -> list[str]:
    return sorted(GENE_CATALOG.keys())


def get_gene(symbol: str) -> dict[str, Any]:
    key = symbol.upper()
    if key not in GENE_CATALOG:
        raise KeyError(f"unknown gene {symbol}; known={list_genes()}")
    return GENE_CATALOG[key]
