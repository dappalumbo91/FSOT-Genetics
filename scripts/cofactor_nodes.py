#!/usr/bin/env python3
"""Extract every metal/cofactor node from real structures and their coordinating residues.

Step 1 of the cofactor program: assign & associate the mineral/chemistry composites
(Fe-heme, Zn, Cu, Ca, Mg, Mn, ...) that each protein carries, and find which residues
coordinate them (any protein atom within a coordination cutoff of the cofactor). This
is the real-data basis for adding cofactor-coordination D_eff rungs to the FSOT ladder.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))
from fsot_compute import domain_scalar  # noqa: E402
import json as _json  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
FSOT_MAP = _json.loads((ROOT / "formulas" / "cofactor_fsot_map.json").read_text(encoding="utf-8"))["cofactors"]

# metals + biological cofactors of interest (exclude water and common cryo additives)
METALS = {"ZN", "FE", "CU", "CA", "MG", "MN", "CO", "NI", "MO", "FE2", "K", "NA"}
COFACTORS = {"HEM", "HEC", "NAD", "NAP", "FAD", "FMN", "PLP", "SAM", "BTN", "COA"}
IGNORE = {"HOH", "GOL", "EDO", "SO4", "PO4", "ACT", "CL", "BR", "IOD", "MPD", "PEG", "DMS", "TRS"}
AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
       "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"}
BENCHMARK = [
    ("1A3N", "Hemoglobin"), ("1CA2", "Carbonic anhydrase II"), ("2C9V", "SOD1"),
    ("1LZ1", "Lysozyme"), ("7RSA", "RNase A"), ("1UBQ", "Ubiquitin"),
    ("4INS", "Insulin"), ("1TUP", "p53"), ("1CLL", "Calmodulin"),
]


def atoms(text):
    prot, het = [], []
    for line in text.splitlines():
        tag = line[:6].strip()
        if tag not in ("ATOM", "HETATM"):
            continue
        res = line[17:20].strip()
        if res == "HOH":
            continue
        try:
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
        rec = (res, line[21], line[22:26].strip(), line[12:16].strip(), xyz)
        if tag == "ATOM" and res in AA3:
            prot.append(rec)
        elif tag == "HETATM":
            het.append(rec)
    return prot, het


def coord_cutoff(res):
    return 3.0 if res in METALS else 4.0  # metals coordinate tighter than organic cofactors


def main() -> int:
    for pdb, name in BENCHMARK:
        fp = CACHE / f"{pdb}.pdb"
        if not fp.exists():
            print(f"{name} ({pdb}): not cached")
            continue
        prot, het = atoms(fp.read_text(encoding="utf-8", errors="replace"))
        # group cofactor atoms by (resname, chain, resnum)
        groups: dict = {}
        for res, ch, num, atom, xyz in het:
            if res in IGNORE:
                continue
            if res in METALS or res in COFACTORS:
                groups.setdefault((res, ch, num), []).append(xyz)
        if not groups:
            print(f"{name:<22}({pdb}): no metal/cofactor")
            continue
        parts = []
        for (res, ch, num), xyzs in sorted(groups.items()):
            cut = coord_cutoff(res)
            coord_res = {}
            for pres, pch, pnum, patom, pxyz in prot:
                d = min(float(np.linalg.norm(pxyz - c)) for c in xyzs)
                if d <= cut:
                    key = (pnum, pres)
                    coord_res[key] = min(coord_res.get(key, 9.9), d)
            ligs = ",".join(f"{pres}{pnum}" for (pnum, pres) in sorted(coord_res, key=lambda k: int(k[0])))
            fm = FSOT_MAP.get(res)
            if fm:
                s = abs(float(domain_scalar(fm["domain"])))
                lig_types = {pres for (pnum, pres) in coord_res}
                consistent = bool(lig_types & set(fm["typical_ligands"])) or not coord_res
                tag = (f"[FSOT {fm['domain']} D_eff={fm['D_eff']} |S|={s:.3f} "
                       f"CN_obs={len(coord_res)}/exp{fm['coordination']} "
                       f"{'chem-OK' if consistent else 'chem-MISMATCH'}]")
            else:
                tag = "[no FSOT route]"
            parts.append(f"{res}[{ch}{num}] {tag}\n        <- {ligs or 'none'}")
        print(f"{name:<22}({pdb}):")
        for p in parts:
            print(f"    {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
