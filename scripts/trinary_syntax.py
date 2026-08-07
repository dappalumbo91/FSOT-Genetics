#!/usr/bin/env python3
"""FSOT trinary syntax expansion — genetics as executable code.

Doctrine (from zig/docs/GENETICS_AS_TRINARY_CODE.md):
  Codons → trinary trips → AA opcodes → expanded syntax word → pair laws.

F01 alone collapses five large-nonpolar AAs to the same (c,p,v).
Expansion adds *lawful structural trits* (not free lookup tables):

  word[7] = (c, p, v, aromatic, branch, hetero, detail)
    aromatic ∈ {-1,0,+1}   ring chemistry
    branch   ∈ {-1,0,+1}   side-chain topology
    hetero   ∈ {-1,0,+1}   S / OH / special N
    detail   ∈ {-1,0,+1}   isomer / functional-group detail (I≠L, R≠K, …)

Continuous spin/charge for Zig pair weight are seed-composites of the word.
Zero free parameters — only {π, e, φ, γ} + discrete chemistry facts as trits.

Also mirrors neuron-zig:
  geometricScaleDist = φ · dist^(-1/π)
  fsotPairWeight = geom · (base + 0.15·elec) · (0.35 + 0.65·env)
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
GAMMA = float(fc.GAMMA)

# ── F01 base triple (c, p, v) — authority ─────────────────────────────────
F01: dict[str, tuple[int, int, int]] = {
    "A": (0, -1, -1),
    "R": (1, 1, 1),
    "N": (0, 1, 0),
    "D": (-1, 1, 0),
    "C": (0, 0, -1),
    "Q": (0, 1, 1),
    "E": (-1, 1, 1),
    "G": (0, -1, -1),
    "H": (1, 1, 1),
    "I": (0, -1, 1),
    "L": (0, -1, 1),
    "K": (1, 1, 1),
    "M": (0, -1, 1),
    "F": (0, -1, 1),
    "P": (0, -1, 0),
    "S": (0, 1, -1),
    "T": (0, 1, 0),
    "W": (0, -1, 1),
    "Y": (0, 1, 1),
    "V": (0, -1, 0),
}

AA20 = "ARNDCQEGHILKMFPSTWYV"


def _trit(x: int | float) -> int:
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


# Structural expansion trits — discrete chemistry, not fitted weights.
# aromatic: F,Y,W = +1; H has imidazole (partial) → 0 with charge already in c
# branch: I,V,L,T = +1; P = −1 (ring constraint); else 0
# hetero: C,M = +1 (S); S,T,Y = −1 (OH); N,Q = 0 already polar; W ring N → +1 with aromatic
def aromatic_trit(aa: str) -> int:
    return 1 if aa in "FYW" else 0


def branch_trit(aa: str) -> int:
    if aa in "IVLT":
        return 1
    if aa == "P":
        return -1
    return 0


def hetero_trit(aa: str) -> int:
    """S / OH / special heteroatom marker (lawful structural trit)."""
    if aa in "CM":
        return 1  # sulfur
    if aa in "STY":
        return -1  # hydroxyl
    if aa == "W":
        return 1  # indole N (with aromatic)
    if aa == "H":
        return 1  # imidazole N
    if aa == "G":
        return -1  # null side chain (syntax NOP)
    return 0


def detail_trit(aa: str) -> int:
    """Isomer / functional-group detail — splits remaining F01 twins.

    Structural facts only:
      I vs L: β-branch (sec-butyl) vs γ-branch (isobutyl)
      R vs K: guanidino vs primary amine
      F vs others already split by aromatic; F detail=+1 pure phenyl
      N vs Q length already in v; keep 0
    """
    if aa == "I":
        return 1  # β-branched
    if aa == "L":
        return -1  # γ-branched
    if aa == "R":
        return 1  # guanidino
    if aa == "K":
        return -1  # primary amine alkyl
    if aa == "F":
        return 1  # pure phenyl (vs W indole already hetero+)
    if aa == "V":
        return 1  # gem-dimethyl β
    if aa == "T":
        return -1  # β-OH already hetero-; mark β-methyl
    return 0


@dataclass(frozen=True)
class AaOpcode:
    """Executable trinary syntax word for one amino acid."""

    aa: str
    c: int
    p: int
    v: int
    aromatic: int
    branch: int
    hetero: int
    detail: int

    def word(self) -> tuple[int, int, int, int, int, int, int]:
        return (
            self.c,
            self.p,
            self.v,
            self.aromatic,
            self.branch,
            self.hetero,
            self.detail,
        )

    def as_string(self) -> str:
        return "[" + ", ".join(f"{t:+d}" for t in self.word()) + "]"

    def spin(self) -> float:
        """Composite spin ∈ [-1,1] — mean of structural trits."""
        parts = (self.p, self.v, self.branch, self.aromatic, self.detail)
        return sum(parts) / float(len(parts))

    def charge(self) -> float:
        """Composite charge — F01 c + hetero/φ + detail/(φ²) (seed scales)."""
        return float(self.c) + float(self.hetero) / PHI + float(self.detail) / (PHI * PHI)

    def side_volume(self) -> float:
        """F02 volume refined by expansion trits (seed-only)."""
        base = PI * E * (PHI ** self.v)
        # aromatic rings expand volume by φ; branch by √φ; S by e^(1/π)
        base *= PHI ** (0.5 * abs(self.aromatic))
        base *= math.sqrt(PHI) ** abs(self.branch) if self.branch > 0 else 1.0
        if self.branch < 0:  # proline ring constraint
            base /= PHI
        if self.hetero > 0:
            base *= math.exp(1.0 / PI)
        if self.hetero < 0 and self.aa != "G":
            base *= math.exp(-1.0 / (PI * PHI))
        if self.aa == "G":
            base /= PHI
        # isomer detail: slight φ-scale volume shift (I vs L, R vs K)
        if self.detail != 0:
            base *= PHI ** (self.detail / PI)
        return base

    def hydrophobicity(self) -> float:
        """F02 h refined: center + aromatic / branch lift."""
        h = (PHI ** (-self.p)) * math.exp(self.v / PI)
        if self.aromatic:
            h *= PHI  # π-stacking / ring hydrophobicity scale
        if self.branch > 0:
            h *= math.sqrt(PHI)
        if self.hetero < 0:  # OH
            h /= PHI
        if self.hetero > 0 and self.aa in "CM":
            h *= math.exp(GAMMA / PI)
        if self.detail != 0:
            h *= math.exp(self.detail / (PI * PHI))
        return h


def aa_opcode(aa: str) -> AaOpcode:
    aa = aa.upper()
    if aa not in F01:
        return AaOpcode(aa, 0, 0, 0, 0, 0, 0, 0)
    c, p, v = F01[aa]
    return AaOpcode(
        aa,
        c,
        p,
        v,
        aromatic_trit(aa),
        branch_trit(aa),
        hetero_trit(aa),
        detail_trit(aa),
    )


def all_opcodes() -> dict[str, AaOpcode]:
    return {aa: aa_opcode(aa) for aa in AA20}


def uniqueness_report() -> dict:
    ops = all_opcodes()
    words = {}
    collisions = []
    for aa, op in ops.items():
        w = op.word()
        if w in words:
            collisions.append((words[w], aa, w))
        else:
            words[w] = aa
    f01_collisions = []
    f01_map: dict[tuple[int, int, int], list[str]] = {}
    for aa, op in ops.items():
        key = (op.c, op.p, op.v)
        f01_map.setdefault(key, []).append(aa)
    for k, aas in f01_map.items():
        if len(aas) > 1:
            f01_collisions.append({"phase": k, "aas": aas})
    return {
        "n_aa": len(ops),
        "unique_expanded_words": len(words),
        "expanded_collisions": collisions,
        "f01_collisions": f01_collisions,
        "all_unique": len(collisions) == 0 and len(words) == 20,
    }


# ── Codon syntax (PRIMARY / SECONDARY) — twin of codon.zig ────────────────
def base_primary(b: str) -> int:
    b = b.upper()
    if b in "AG":
        return 1
    if b in "CTU":
        return -1
    return 0


def base_secondary(b: str) -> int:
    b = b.upper()
    if b == "A":
        return 1
    if b in "TU":
        return -1
    if b in "GC":
        return 0
    return 0


def codon_primary(codon: str) -> tuple[int, int, int]:
    c = codon.upper().replace("U", "T")
    assert len(c) == 3
    return base_primary(c[0]), base_primary(c[1]), base_primary(c[2])


def codon_secondary(codon: str) -> tuple[int, int, int]:
    c = codon.upper().replace("U", "T")
    return base_secondary(c[0]), base_secondary(c[1]), base_secondary(c[2])


# Standard genetic code DNA
_GENETIC_CODE: dict[str, str] = {}


def _init_code() -> None:
    # Built to match codon.zig dnaToAa
    bases = "TCAG"
    # Explicit table from codon.zig structure
    table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    _GENETIC_CODE.update(table)


_init_code()


def dna_to_aa(codon: str) -> str:
    c = codon.upper().replace("U", "T")
    return _GENETIC_CODE.get(c, "?")


def all_codons() -> list[str]:
    return [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"]


# ── Zig pair geometry (genetic_pair.zig / genetic.zig) ────────────────────
def trinary_pair_interaction(tau_i: float, tau_j: float) -> float:
    ti = max(-1.0, min(1.0, float(tau_i)))
    tj = max(-1.0, min(1.0, float(tau_j)))
    prod = ti * tj
    return prod * E + (1.0 - abs(prod)) * PI


def geometric_scale_dist(dist: int) -> float:
    d = float(max(int(dist), 1))
    return PHI * (d ** (-1.0 / PI))


def electrostatic_term(q_i: float, q_j: float) -> float:
    return -q_i * q_j * E


def env_scale(dist: int) -> float:
    d = float(max(int(dist), 1))
    return d / (d + PI * E)


def fsot_pair_weight(
    spin_i: float,
    spin_j: float,
    charge_i: float,
    charge_j: float,
    dist: int,
) -> float:
    base = trinary_pair_interaction(spin_i, spin_j)
    geom = geometric_scale_dist(dist)
    elec = electrostatic_term(charge_i, charge_j)
    env = env_scale(dist)
    return geom * (base + 0.15 * elec) * (0.35 + 0.65 * env)


def aa_pair_weight(aa1: str, aa2: str, dist: int) -> float:
    o1, o2 = aa_opcode(aa1), aa_opcode(aa2)
    return fsot_pair_weight(o1.spin(), o2.spin(), o1.charge(), o2.charge(), dist)


def expanded_chemical_interaction(aa1: str, aa2: str, sep: int = 1) -> float:
    """F03–F06 chemistry + Zig pair geometry (higher precision path).

    Combines classic F03–F06 terms with fsotPairWeight scaled by 1/(πe)
    so magnitude stays formula-native (not a free fit weight).
    """
    c1, c2 = aa1.upper(), aa2.upper()
    if c1 == "C" and c2 == "C":
        return PHI ** 6  # F03 disulfide still dominates
    o1, o2 = aa_opcode(c1), aa_opcode(c2)
    # F04 hydrophobic (expanded h)
    h1 = (o1.hydrophobicity() - 1.0) / PHI
    h2 = (o2.hydrophobicity() - 1.0) / PHI
    hydrophobic = h1 * h2
    # F05 electrostatic on composite charge
    electrostatic = -o1.charge() * o2.charge() * E
    # F06 dipole from expanded volume + polarity
    mu1 = GAMMA * math.exp(abs(o1.c) + o1.p + 1.0 + 0.5 * abs(o1.aromatic))
    mu2 = GAMMA * math.exp(abs(o2.c) + o2.p + 1.0 + 0.5 * abs(o2.aromatic))
    dipole = math.sqrt(max(mu1 * mu2, 0.0)) / (GAMMA * PI * E * E)
    # Zig pair geometry contribution (syntax-level)
    pair = aa_pair_weight(c1, c2, max(sep, 1)) / (PI * E)
    # aromatic stacking bonus when both aromatic (φ-law)
    stack = 0.0
    if o1.aromatic and o2.aromatic:
        stack = (1.0 / PHI) * env_scale(sep)
    return hydrophobic + electrostatic + dipole + pair + stack


def write_expanded_maps(out_dir: Path | None = None) -> dict:
    out_dir = out_dir or (ROOT / "formulas")
    out_dir.mkdir(parents=True, exist_ok=True)
    ops = all_opcodes()
    uniq = uniqueness_report()

    # 20-AA expanded map
    lines = [
        "FSOT 20-AMINO-ACID EXPANDED TRINARY SYNTAX (6-trit opcode)",
        "==========================================================================================================",
        " AA | F01 [C,P,V]     | EXPANDED [C,P,V,Aro,Br,Het,Det] | SPIN     | CHARGE   | H_EXP    | VOL_EXP",
        "----------------------------------------------------------------------------------------------------------",
    ]
    rows = []
    for aa in AA20:
        op = ops[aa]
        f01 = f"[{op.c:+d}, {op.p:+d}, {op.v:+d}]"
        lines.append(
            f" {aa}  | {f01:15s} | {op.as_string():28s} | "
            f"{op.spin():+8.4f} | {op.charge():+8.4f} | "
            f"{op.hydrophobicity():8.4f} | {op.side_volume():8.4f}"
        )
        rows.append(
            {
                "aa": aa,
                "f01": [op.c, op.p, op.v],
                "expanded": list(op.word()),
                "spin": op.spin(),
                "charge": op.charge(),
                "h": op.hydrophobicity(),
                "vol": op.side_volume(),
            }
        )
    lines.append("==========================================================================================================")
    lines.append(f"unique_expanded_words={uniq['unique_expanded_words']}/20  all_unique={uniq['all_unique']}")
    lines.append("F01 collisions (pre-expansion): " + json.dumps(uniq["f01_collisions"]))
    lines.append("Law: aromatic/branch/hetero are structural trits; continuous fields from {π,e,φ,γ} only.")
    (out_dir / "20_amino_acid_expanded_trinary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 64-codon with AA opcode
    clines = [
        "FSOT 64-CODON TRINARY + AA EXPANDED OPCODE",
        "================================================================================",
        "CODON | PRIMARY      | SECONDARY    | AA | EXPANDED WORD",
        "--------------------------------------------------------------------------------",
    ]
    codon_rows = []
    for codon in all_codons():
        prim = codon_primary(codon)
        sec = codon_secondary(codon)
        aa = dna_to_aa(codon)
        if aa == "*":
            exp = "STOP"
            word = None
        else:
            op = aa_opcode(aa)
            exp = op.as_string()
            word = list(op.word())
        clines.append(
            f"{codon}  | {list(prim)} | {list(sec)} | {aa:1s}  | {exp}"
        )
        codon_rows.append(
            {
                "codon": codon,
                "primary": list(prim),
                "secondary": list(sec),
                "aa": aa,
                "expanded": word,
            }
        )
    clines.append("================================================================================")
    (out_dir / "64_codon_expanded_trinary.txt").write_text("\n".join(clines) + "\n", encoding="utf-8")

    payload = {
        "authority": "F01 + structural expansion trits + neuron-zig pair geometry",
        "free_parameters": 0,
        "pin": "D1D38A",
        "uniqueness": uniq,
        "amino_acids": rows,
        "codons": codon_rows,
        "pair_law": "fsotPairWeight = φ·d^{-1/π} · (τ_iτ_j e + (1-|prod|)π + 0.15·elec) · (0.35+0.65·env)",
    }
    (out_dir / "trinary_syntax_expanded.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    return payload


def main() -> int:
    payload = write_expanded_maps()
    u = payload["uniqueness"]
    print("FSOT trinary syntax expansion")
    print(f"  unique expanded AA words: {u['unique_expanded_words']}/20  all_unique={u['all_unique']}")
    print(f"  F01 collision groups: {len(u['f01_collisions'])}")
    for g in u["f01_collisions"]:
        print(f"    F01 {g['phase']} → {g['aas']}  (resolved by expansion)")
    if u["expanded_collisions"]:
        print("  FAIL expanded collisions:", u["expanded_collisions"])
        return 1
    # sample pair weights
    for d in (1, 4, 8, 20):
        w = aa_pair_weight("F", "W", d)
        print(f"  pair F–W dist={d:2d}  w={w:.6f}")
    print("  wrote formulas/20_amino_acid_expanded_trinary.txt")
    print("  wrote formulas/64_codon_expanded_trinary.txt")
    print("  wrote formulas/trinary_syntax_expanded.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
