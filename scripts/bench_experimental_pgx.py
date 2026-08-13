#!/usr/bin/env python3
"""Experimental PGx cards vs public known outcomes.

Not a medical device. See docs/EXPERIMENTAL_DISCLOSURE.md.

Each case is a measured co-crystal (or DNA/metal complex) plus a published
mechanism class. FSOT only asks: does this residue sit on the ChemLink
observer (ligand / metal / DNA)? Concordance confirms the interface.
Discordance is a diagnostic for refinement — same loop as Cα RMSD.
0 free parameters. Pin D1D38A.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from multi_system import (  # noqa: E402
    CONTACT,
    DNA3,
    E,
    PHI,
    _get_pdb,
    ca_with_nums,
    na_chains,
    parse_ligands,
    parse_metals,
    parse_na_c1,
    parse_pdb_ca,
)

OUT = ROOT / "data" / "experimental_pgx.json"
DISCLOSURE = (
    "EXPERIMENTAL RESEARCH — not a diagnosis, prescription, or medical device. "
    "See docs/EXPERIMENTAL_DISCLOSURE.md."
)

# Public, textbook mechanism classes. Labels are evaluation data (like a native PDB).
CASES: list[dict[str, Any]] = [
    {
        "id": "kras_g12c_sotorasib",
        "gene": "KRAS",
        "pdb": "6OIM",
        "chain": "A",
        "pdb_resnum": "12",
        "drug": "sotorasib (AMG 510)",
        "known_class": "on_drug_site",
        "known_outcome": "covalent Cys12 engagement; G12C-selective inhibitor",
        "source": "PDB 6OIM (canon co-crystal)",
    },
    {
        "id": "egfr_l858r_gefitinib",
        "gene": "EGFR",
        "pdb": "2ITY",
        "chain": "A",
        "pdb_resnum": "858",
        "drug": "gefitinib",
        "known_class": "allosteric_state",
        "known_outcome": "sensitizing L858R is activation-loop state, not the TKI first shell",
        "source": "PDB 2ITY (gefitinib binds cleft; L858 is ~16 Å away)",
    },
    {
        "id": "egfr_t790m_gefitinib",
        "gene": "EGFR",
        "pdb": "2ITY",
        "chain": "A",
        "pdb_resnum": "790",
        "drug": "gefitinib",
        "known_class": "on_drug_site",
        "known_outcome": "gatekeeper; steric resistance to first-gen TKI",
        "source": "PDB 2ITY residue 790 (ATP/TKI cleft)",
    },
    {
        "id": "braf_v600e_vemurafenib",
        "gene": "BRAF",
        "pdb": "3OG7",
        "chain": "B",
        "pdb_resnum": "600",
        "drug": "vemurafenib",
        "known_class": "allosteric_state",
        "known_outcome": "V600E activation-loop state; vemurafenib binds the cleft, not V600 first-shell",
        "source": "PDB 3OG7 chain B",
    },
    {
        "id": "sod1_h47r_cu",
        "gene": "SOD1",
        "pdb": "2C9V",
        "chain": "A",
        "pdb_resnum": "46",
        "drug": "Cu (metal observer)",
        "known_class": "on_metal_site",
        "known_outcome": "Cu ligand His; ALS-associated",
        "source": "PDB 2C9V (historic H46 / UniProt H47)",
    },
    {
        "id": "tp53_r248w_dna",
        "gene": "TP53",
        "pdb": "1TUP",
        "chain": "A",
        "pdb_resnum": "248",
        "drug": "DNA (observer)",
        "known_class": "on_dna_site",
        "known_outcome": "DNA-contact hotspot; damaging",
        "source": "PDB 1TUP / Cho et al.",
    },
    {
        "id": "tp53_p72r_control",
        "gene": "TP53",
        "pdb": "1TUP",
        "chain": "A",
        "pdb_resnum": "72",
        "drug": None,
        "known_class": "off_site_polymorphism",
        "known_outcome": "common P72R polymorphism; not a DNA-contact driver",
        "source": "population / catalog control",
    },
    {
        "id": "cyp2c9_i359l_warfarin",
        "gene": "CYP2C9",
        "pdb": "1OG5",
        "chain": "A",
        "pdb_resnum": "359",
        "drug": "warfarin",
        "known_class": "allosteric_state",
        "known_outcome": "CYP2C9*3 I359L reduces warfarin clearance; SRS5 ~15 Å from the warfarin first shell (CPIC-class metabolizer, not the binder)",
        "source": "PDB 1OG5 (CYP2C9–warfarin)",
    },
    {
        "id": "caii_zn_site",
        "gene": "CA2",
        "pdb": "1CA2",
        "chain": "A",
        "pdb_resnum": "94",
        "drug": "Zn (metal observer)",
        "known_class": "on_metal_site",
        "known_outcome": "catalytic His94–Zn",
        "source": "PDB 1CA2",
    },
]


def _res_index(nums: list[str], pdb_resnum: str) -> int | None:
    key = str(pdb_resnum).strip()
    for i, n in enumerate(nums):
        if n.strip() == key:
            return i
    return None


def _atom_hits(
    text: str, chain: str, nums: list[str], cutoff: float
) -> tuple[set[int], set[int], set[int]]:
    """Residue indices contacting ligand / metal / DNA C1'."""
    idx = {n.strip(): i for i, n in enumerate(nums)}
    atoms: list[tuple[str, np.ndarray]] = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        num = line[22:26].strip()
        if num not in idx:
            continue
        atom = line[12:16].strip()
        if atom.startswith("H") or atom.startswith("D"):
            continue
        atoms.append(
            (num, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
        )
    lig, met, dna = set(), set(), set()
    cut = cutoff
    for lg in parse_ligands(text):
        for num, p in atoms:
            if any(float(np.linalg.norm(p - a)) <= cut for a in lg["atoms"]):
                lig.add(idx[num])
    for m in parse_metals(text):
        for num, p in atoms:
            if float(np.linalg.norm(p - m["xyz"])) <= 3.0:
                met.add(idx[num])
    for ch in na_chains(text, DNA3):
        _s, nx = parse_na_c1(text, ch, DNA3)
        if len(nx) < 2:
            continue
        for num, p in atoms:
            if float(np.linalg.norm(nx - p, axis=1).min()) <= CONTACT * PHI:
                dna.add(idx[num])
    return lig, met, dna


def classify(i: int | None, lig: set[int], met: set[int], dna: set[int]) -> str:
    if i is None:
        return "not_on_chain"
    if i in met:
        return "on_metal_site"
    if i in dna:
        return "on_dna_site"
    if i in lig:
        return "on_drug_site"
    return "off_site"


def concordant(known: str, got: str) -> bool:
    if known == "off_site_polymorphism":
        return got in ("off_site", "not_on_chain")
    if known == "allosteric_state":
        # Mutation changes collapse; it is not the ligand first shell.
        return got in ("off_site", "activation_loop_near_site", "allosteric_state")
    if known == "activation_loop_near_site":
        return got in ("on_drug_site", "activation_loop_near_site")
    return got == known


def run_case(case: dict[str, Any]) -> dict[str, Any]:
    text = _get_pdb(case["pdb"])
    seq, xyz, nums = ca_with_nums(text, case["chain"])
    if len(seq) < 8:
        _s, xyz = parse_pdb_ca(text, case["chain"])
        nums = [str(i + 1) for i in range(len(_s))]
        seq = _s
    i = _res_index(nums, case["pdb_resnum"])
    cutoff = float(E + PHI)
    lig, met, dna = _atom_hits(text, case["chain"], nums, cutoff)
    got = classify(i, lig, met, dna)
    # Activation-loop class: on-site or CA within φ of any ligand-contact CA.
    if case["known_class"] == "activation_loop_near_site" and i is not None and got == "off_site":
        if lig and xyz is not None and len(xyz) == len(nums):
            site = xyz[sorted(lig)]
            dmin = float(np.linalg.norm(site - xyz[i], axis=1).min())
            if dmin <= PHI * PHI:
                got = "activation_loop_near_site"
    ok = concordant(case["known_class"], got)
    aa = seq[i] if i is not None and i < len(seq) else None
    card = {
        "id": case["id"],
        "disclosure": DISCLOSURE,
        "gene": case["gene"],
        "pdb": case["pdb"],
        "residue": case["pdb_resnum"],
        "aa_on_chain": aa,
        "drug": case["drug"],
        "fsot_class": got,
        "known_class": case["known_class"],
        "known_outcome": case["known_outcome"],
        "source": case["source"],
        "concordant": ok,
        "n_ligand_site": len(lig),
        "n_metal_site": len(met),
        "n_dna_site": len(dna),
        "domain": (
            "Molecular_Chemistry"
            if "drug" in (got + case["known_class"])
            else "Electromagnetism"
            if "dna" in got or "dna" in case["known_class"]
            else "Atomic_Physics"
            if "metal" in got or "metal" in case["known_class"]
            else "Biochemistry"
        ),
        "diagnostic": None if ok else (
            "interface_miss — residue not first-shell on this crystal; "
            "refine numbering, pose, or ChemLink cutoff"
        ),
        "free_parameters": 0,
    }
    return card


def main() -> int:
    print(DISCLOSURE, flush=True)
    print("Experimental PGx bench (FSOT ChemLink vs public mechanism)", flush=True)
    cards = []
    for case in CASES:
        print(f"  {case['id']}…", flush=True)
        try:
            card = run_case(case)
        except Exception as exc:
            card = {
                "id": case["id"],
                "disclosure": DISCLOSURE,
                "status": "error",
                "error": str(exc),
                "concordant": False,
                "diagnostic": "fetch_or_parse_failed",
            }
        cards.append(card)
        print(
            f"    {card.get('fsot_class')} vs {case['known_class']} "
            f"{'OK' if card.get('concordant') else 'MISS'} "
            f"({case['pdb']})",
            flush=True,
        )
    n = len(cards)
    hits = sum(1 for c in cards if c.get("concordant"))
    misses = [c["id"] for c in cards if not c.get("concordant")]
    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "disclosure": DISCLOSURE,
        "doctrine": "Public dose/tox/mechanism labels are measured authority. "
        "FSOT reports whether the variant occupies that ChemLink. "
        "Hits confirm the interface; misses diagnose it. Not a prescription.",
        "n": n,
        "n_concordant": hits,
        "concordance": hits / n if n else None,
        "misses": misses,
        "cards": cards,
        "free_parameters": 0,
        "pin": "D1D38A",
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Wrote {OUT}", flush=True)
    print(f"concordant={hits}/{n} misses={misses}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
