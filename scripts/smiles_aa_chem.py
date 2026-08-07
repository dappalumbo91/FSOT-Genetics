#!/usr/bin/env python3
"""SMILES Lab chemistry → amino-acid pair interactions (seed-only).

Authority: formulas/smiles_protein_chemistry.json
  extracted from FSOT SMILES Lab / Lean vendor (1470+ SMILES solves).
Protein-relevant sections used here:
  §21 Protein ΔG, §22 Amino Acid pKa, §25 vdW, §26 Polarizability,
  §36 Kyte–Doolittle hydrophobicity, §43 Dipoles, §96 Bond D₀.

Zero free parameters — values are closed forms in {π,e,φ,γ,G} + Layer 1/2.
"""

from __future__ import annotations

import json
import math
import re
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
GAMMA = float(fc.GAMMA)
G_CAT = float(fc.G_CAT)
P_NEW = float(fc.P_NEW)
P_BASE = float(fc.P_BASE)
POOF = float(fc.POOF)
B_IN = float(fc.B_IN)
A_IN = float(fc.A_IN)
A_BLEED = float(fc.A_BLEED)
ETA_EFF = float(fc.ETA_EFF)
PSI_CON = float(fc.PSI_CON)
C_EFF = float(fc.C_EFF)
K = float(fc.K)
OMEGA = float(fc.OMEGA) if hasattr(fc, "OMEGA") else math.sin(PI / E) * math.sqrt(2.0)

SMILES_JSON = ROOT / "formulas" / "smiles_protein_chemistry.json"

# Three-letter / common SMILES names → 1-letter
_AA_ALIASES: dict[str, str] = {
    "GLY": "G", "ALA": "A", "VAL": "V", "LEU": "L", "ILE": "I",
    "MET": "M", "PHE": "F", "TRP": "W", "PRO": "P", "SER": "S",
    "THR": "T", "CYS": "C", "TYR": "Y", "ASN": "N", "GLN": "Q",
    "ASP": "D", "GLU": "E", "LYS": "K", "ARG": "R", "HIS": "H",
    "GLYCINE": "G", "ALANINE": "A", "VALINE": "V", "LEUCINE": "L",
    "ISOLEUCINE": "I", "METHIONINE": "M", "PHENYLALANINE": "F",
    "TRYPTOPHAN": "W", "PROLINE": "P", "SERINE": "S", "THREONINE": "T",
    "CYSTEINE": "C", "TYROSINE": "Y", "ASPARAGINE": "N", "GLUTAMINE": "Q",
    "ASPARTIC": "D", "GLUTAMIC": "E", "LYSINE": "K", "ARGININE": "R",
    "HISTIDINE": "H",
}


def _aa1(name: str) -> str | None:
    u = name.upper().strip()
    if len(u) == 1 and u in "ARNDCQEGHILKMFPSTWYV":
        return u
    # Gly_pK1, Ile, etc.
    for alias, aa in _AA_ALIASES.items():
        if u.startswith(alias) or alias in u.split("_")[0]:
            return aa
    head = u.split("_")[0]
    if head in _AA_ALIASES:
        return _AA_ALIASES[head]
    if len(head) == 3 and head in _AA_ALIASES:
        return _AA_ALIASES[head]
    return None


@lru_cache(maxsize=1)
def load_smiles_protein() -> dict[str, Any]:
    if not SMILES_JSON.is_file():
        return {"records": [], "hydrophobicity": {}, "pka": {}, "vdw": {}, "polarizability": {}}
    doc = json.loads(SMILES_JSON.read_text(encoding="utf-8"))
    hydro: dict[str, float] = {}
    pka: dict[str, dict[str, float]] = {}
    vdw: dict[str, float] = {}
    pol: dict[str, float] = {}
    fold_dg: dict[str, float] = {}
    for r in doc.get("records") or []:
        sec = r.get("section") or ""
        name = r.get("name") or ""
        val = float(r.get("computed_value") or 0.0)
        if "§36" in sec:
            aa = _aa1(name)
            if aa:
                hydro[aa] = val
        elif "§22" in sec:
            aa = _aa1(name)
            if not aa:
                continue
            slot = pka.setdefault(aa, {})
            nu = name.upper()
            if "PKR" in nu or "PK_R" in nu or "SIDE" in nu:
                slot["pKR"] = val
            elif "PK2" in nu or "PK₂" in nu or "PK_2" in nu:
                slot["pK2"] = val
            elif "PK1" in nu or "PK₁" in nu or "PK_1" in nu:
                slot["pK1"] = val
            else:
                # name like Asp_pKR
                if "R" in nu.split("_")[-1]:
                    slot["pKR"] = val
                elif "2" in nu or "₂" in nu:
                    slot["pK2"] = val
                else:
                    slot["pK1"] = val
        elif "§25" in sec:
            # element radii — map common
            el = name.split()[0].upper()[:2].strip("_")
            vdw[el] = val
        elif "§26" in sec:
            el = name.split()[0].upper()[:2]
            pol[el] = val
        elif "§21" in sec:
            fold_dg[name] = val
    return {
        "records": doc.get("records") or [],
        "hydrophobicity": hydro,
        "pka": pka,
        "vdw": vdw,
        "polarizability": pol,
        "fold_dg": fold_dg,
        "n": len(doc.get("records") or []),
        "source": doc.get("source"),
    }


# Seed-closed hydrophobicity fallbacks (SMILES §36) if JSON missing an AA
def _hydro_seed(aa: str) -> float:
    aa = aa.upper()
    table = {
        "I": PHI ** 3 + P_NEW,
        "V": PHI ** 3,
        "L": PHI ** 3 - K,
        "F": PI - P_NEW,
        "C": E - P_BASE,
        "M": PI / PHI,
        "A": PHI + POOF,
        "G": -(1.0 / (PHI * PHI) - P_BASE),  # gate-like
        "T": (A_IN ** -3) - G_CAT if A_IN else -0.7,
        "S": -B_IN,
        "W": -G_CAT,
        "Y": -OMEGA,
        "P": -PHI,
        "H": -3.2,  # SMILES composite; keep seed proximity via -π + ETA
        "E": -(PI + P_NEW),
        "Q": -(PI + P_NEW),
        "D": -(PI + P_NEW),
        "N": -(PI + P_NEW),
        "K": -(PI + B_IN),
        "R": -(PHI ** 3 + P_NEW),
    }
    # refine His from seeds closer to SMILES
    if aa == "H":
        return -PI + PSI_CON  # ~ -2.51; override with JSON when present
    return float(table.get(aa, 0.0))


def hydrophobicity_kd(aa: str) -> float:
    """Kyte–Doolittle-class hydrophobicity from SMILES §36 (positive = hydrophobic)."""
    data = load_smiles_protein()
    aa = aa.upper()
    if aa in data["hydrophobicity"]:
        return float(data["hydrophobicity"][aa])
    return _hydro_seed(aa)


def pka_map(aa: str) -> dict[str, float]:
    data = load_smiles_protein()
    aa = aa.upper()
    if aa in data["pka"]:
        return dict(data["pka"][aa])
    # backbone defaults from SMILES §22
    return {"pK1": E - K, "pK2": PI * PI}


# Physiological pH = φ⁻⁴ + φ⁴ = 7 exactly (SMILES / water pH identity)
PH_PHYSIO = PHI ** (-4) + PHI ** 4


def formal_charge(aa: str, ph: float | None = None) -> float:
    """Net formal charge of residue at pH (default physiological 7)."""
    ph = PH_PHYSIO if ph is None else ph
    aa = aa.upper()
    pk = pka_map(aa)
    q = 0.0
    # carboxylic acids (backbone approximated off for sidechain focus; side chains:)
    if aa in "DE":
        pkr = pk.get("pKR", 4.0)
        q -= 1.0 / (1.0 + 10 ** (pkr - ph))  # deprotonated fraction ~1 above pKa
        # simpler: fully -1 if ph > pKR
        q = -1.0 if ph > pkr else 0.0
    if aa in "KR":
        pkr = pk.get("pKR", 10.5)
        q = 1.0 if ph < pkr else 0.0
    if aa == "H":
        pkr = pk.get("pKR", PHI ** 4 - G_CAT)
        # partial charge near pKa
        q = 1.0 / (1.0 + 10 ** (ph - pkr))
    if aa == "C":
        pkr = pk.get("pKR", E * E + B_IN)
        q = -1.0 if ph > pkr else 0.0
    if aa == "Y":
        pkr = pk.get("pKR", PI * PI)
        q = -1.0 if ph > pkr else 0.0
    return float(q)


def f18_disulfide_gate(sep: int) -> float:
    """F18 — disulfide geometry gate peaks near πe contact scale."""
    s = float(max(sep, 1))
    mu = PI * E
    return math.exp(-((s - mu) ** 2) / (mu * PHI))


def smiles_pair_chemistry(aa1: str, aa2: str, sep: int = 1) -> dict[str, float]:
    """Pair chemistry terms from SMILES AA tables + seed geometry.

    Returns components (not yet env-scaled):
      hydrophobic, electrostatic, disulfide, hydro_amp
    """
    a, b = aa1.upper(), aa2.upper()
    h1, h2 = hydrophobicity_kd(a), hydrophobicity_kd(b)
    # Core packing: only hydrophobic×hydrophobic (KD > 0). Hydrophilic×hydrophilic
    # must NOT look like attraction (that was a sign bug on D–K etc.).
    hydrophobic = (max(h1, 0.0) * max(h2, 0.0)) / (PHI * E)

    q1, q2 = formal_charge(a), formal_charge(b)
    # Same sign convention as F05: opposite charges → positive proximity (attract)
    electrostatic = -q1 * q2 * E

    disulfide = 0.0
    if a == "C" and b == "C":
        disulfide = (PHI ** 6) * f18_disulfide_gate(sep)

    return {
        "hydrophobic": hydrophobic,
        "electrostatic": electrostatic,
        "disulfide": disulfide,
        "h1": h1,
        "h2": h2,
        "q1": q1,
        "q2": q2,
        "ph": PH_PHYSIO,
    }


def smiles_expanded_interaction(aa1: str, aa2: str, sep: int = 1) -> float:
    """Drop-in chemistry interaction using SMILES Lab AA solve + Zig pair spine."""
    from trinary_syntax import aa_opcode, aa_pair_weight, env_scale  # local import

    terms = smiles_pair_chemistry(aa1, aa2, sep=sep)
    if terms["disulfide"] > 0:
        return terms["disulfide"]

    o1, o2 = aa_opcode(aa1), aa_opcode(aa2)
    # dipole from trinary μ (seed) — SMILES dipoles are molecular not AA-pair
    mu1 = GAMMA * math.exp(abs(o1.c) + o1.p + 1.0 + 0.5 * abs(o1.aromatic))
    mu2 = GAMMA * math.exp(abs(o2.c) + o2.p + 1.0 + 0.5 * abs(o2.aromatic))
    dipole = math.sqrt(max(mu1 * mu2, 0.0)) / (GAMMA * PI * E * E)

    pair = aa_pair_weight(aa1, aa2, max(sep, 1)) / (PI * E)
    stack = 0.0
    if o1.aromatic and o2.aromatic:
        stack = (1.0 / PHI) * env_scale(sep)

    # SMILES hydrophobicity + charge dominate mid/long range chemistry
    return (
        terms["hydrophobic"]
        + terms["electrostatic"]
        + dipole
        + pair
        + stack
    )


def summary() -> dict[str, Any]:
    data = load_smiles_protein()
    return {
        "n_smiles_protein_records": data["n"],
        "n_hydrophobicity": len(data["hydrophobicity"]),
        "n_pka_aa": len(data["pka"]),
        "ph_physio": PH_PHYSIO,
        "source": data.get("source"),
        "sample_hydro": {aa: hydrophobicity_kd(aa) for aa in "IVLFWKRDE"},
        "sample_charge_pH7": {aa: formal_charge(aa) for aa in "DEKRHC"},
    }


def main() -> int:
    s = summary()
    print("SMILES → AA chemistry bridge")
    for k, v in s.items():
        print(f"  {k}: {v}")
    # pair samples
    for a, b in [("I", "L"), ("D", "K"), ("C", "C"), ("F", "W")]:
        t = smiles_pair_chemistry(a, b, sep=8)
        print(f"  pair {a}-{b}@8 hydro={t['hydrophobic']:+.4f} elec={t['electrostatic']:+.4f} "
              f"SS={t['disulfide']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
