#!/usr/bin/env python3
"""Variant-effect pilot: does FSOT + structural annotation flag pathogenic mutations?

Medical direction step 1. For p53 (1TUP: protein + its DNA + structural Zn), build a
zero-parameter variant-impact score for every possible missense mutation:
  impact = criticality(position) * chemistry_change(wt -> mut)
where criticality is derived ONLY from the structure (non-circular): DNA-contact,
Zn-coordination, and burial; and chemistry_change is the FSOT/trinary propensity
change (charge, hydrophobicity, volume). Validation: do the six well-known Cho'94
cancer hotspot mutations rank near the top of all ~3700 possible missense variants?
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import chemical_propensity  # noqa: E402
from trinary_syntax import aa_opcode  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}
DNA = {"DA", "DT", "DG", "DC", "A", "T", "G", "C"}
AA20 = "ARNDCQEGHILKMFPSTWYV"
# Cho'94 p53 core cancer hotspots (resnum, wt, mut)
HOTSPOTS = [(175, "R", "H"), (245, "G", "S"), (248, "R", "Q"),
            (249, "R", "S"), (273, "R", "H"), (282, "R", "W")]


def chem_change(wt, mut):
    a, b = chemical_propensity(wt), chemical_propensity(mut)
    return (abs(a.q - b.q) + abs(a.h - b.h) / 4.5 + abs(a.vol - b.vol))


def trinary_delta(wt, mut):
    """Mutation expressed in the trinary substrate: distance over the 7-trit opcode word."""
    a = np.array(aa_opcode(wt).word(), dtype=float)
    b = np.array(aa_opcode(mut).word(), dtype=float)
    return float(np.linalg.norm(a - b))


def parse_1tup():
    text = (CACHE / "1TUP.pdb").read_text(encoding="utf-8", errors="replace")
    ca, nums, prot_atoms, dna_atoms, zn = [], [], [], [], []
    seen = set()
    for line in text.splitlines():
        tag = line[:6].strip()
        if tag not in ("ATOM", "HETATM"):
            continue
        res = line[17:20].strip()
        try:
            p = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        if tag == "ATOM" and line[21] == "A" and res in AA3:
            prot_atoms.append((line[22:26].strip(), p))
            if line[12:16].strip() == "CA" and (line[21], line[22:26].strip()) not in seen:
                seen.add((line[21], line[22:26].strip()))
                ca.append(p)
                nums.append(int(line[22:26].strip()))
        elif res in DNA:
            dna_atoms.append(p)
        elif tag == "HETATM" and res == "ZN":
            zn.append(p)
    return np.array(ca), nums, prot_atoms, np.array(dna_atoms), np.array(zn)


def main() -> int:
    ca, nums, prot_atoms, dna, zn = parse_1tup()
    n = len(nums)
    idx = {num: i for i, num in enumerate(nums)}
    # sequence from prot_atoms (one letter per resnum, via CA order)
    seqmap = {}
    for line_num, p in prot_atoms:
        pass
    # build wt seq by resnum from CA residues
    text = (CACHE / "1TUP.pdb").read_text(encoding="utf-8", errors="replace")
    wt = {}
    for line in text.splitlines():
        if line.startswith("ATOM") and line[21] == "A" and line[12:16].strip() == "CA":
            r = line[17:20].strip()
            if r in AA3:
                wt[int(line[22:26].strip())] = AA3[r]

    # structural criticality per residue (non-circular): DNA-contact, Zn, burial
    w_dna = np.zeros(n)
    w_zn = np.zeros(n)
    burial = np.zeros(n)
    # group protein atoms by resnum
    by_res: dict = {}
    for num, p in prot_atoms:
        by_res.setdefault(int(num), []).append(p)
    for num, i in idx.items():
        atoms = by_res.get(num, [])
        if dna.size and atoms:
            dmin = min(float(np.linalg.norm(a - d)) for a in atoms for d in dna)
            w_dna[i] = 1.0 if dmin <= 4.5 else 0.0
        if zn.size and atoms:
            zmin = min(float(np.linalg.norm(a - z)) for a in atoms for z in zn)
            w_zn[i] = 1.0 if zmin <= 3.0 else 0.0
    d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=2)
    burial = ((d < 10.0).sum(1) - 1).astype(float)
    burial /= burial.max()
    crit = 1.0 + w_dna + w_zn + burial

    # score every possible missense variant, chemistry vs trinary-opcode terms
    chem_scores, tri_scores = [], []
    for num, i in idx.items():
        w = wt[num]
        for mut in AA20:
            if mut == w:
                continue
            chem_scores.append(crit[i] * chem_change(w, mut))
            tri_scores.append(crit[i] * trinary_delta(w, mut))
    chem_scores = np.array(chem_scores)
    tri_scores = np.array(tri_scores)

    print(f"p53 1TUP/A  n={n}  DNA-contact={int(w_dna.sum())} Zn={int(w_zn.sum())}")
    print(f"scored {len(chem_scores)} possible missense variants\n")
    print(f"{'hotspot':<9}{'wt-opcode':>22}{'mut-opcode':>22}{'triDelta':>9}{'chem%':>7}{'tri%':>7}")
    print("-" * 76)
    cpct, tpct = [], []
    for num, w, mut in HOTSPOTS:
        if num not in idx:
            continue
        i = idx[num]
        cimp = crit[i] * chem_change(w, mut)
        timp = crit[i] * trinary_delta(w, mut)
        cp = float((chem_scores < cimp).mean()) * 100
        tp = float((tri_scores < timp).mean()) * 100
        cpct.append(cp)
        tpct.append(tp)
        print(f"{w}{num}{mut:<5}{aa_opcode(w).as_string():>22}{aa_opcode(mut).as_string():>22}"
              f"{trinary_delta(w, mut):>9.2f}{cp:>6.0f}%{tp:>6.0f}%")
    print("-" * 76)
    print(f"mean hotspot percentile:  chemistry {np.mean(cpct):.1f}%   "
          f"trinary {np.mean(tpct):.1f}%   (random 50%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
