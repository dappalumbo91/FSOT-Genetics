#!/usr/bin/env python3
"""Does FSOT's sequence-only scalar importance flag p53's functional mechanism?

Literature mechanism (1TUP): Zn structural node C176/H179/C238/C242 shapes the
DNA-binding surface; DNA-contact hotspots R248 (minor groove) and R273 (phosphate
backbone) grip DNA; R175 is the structural hotspot. We compute a per-residue FSOT
importance = mean |S| over long-range pairs (center-line deviation = functional
cost) and test whether the known functional residues are enriched at the top,
against a permutation null. Sequence in, structure-defined function tested.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import fetch_pdb  # noqa: E402
from fsot_structure_engine import SsPropensity, LONG_RANGE_GATE  # noqa: E402
from trinary_syntax import aa_opcode  # noqa: E402
from full_scalar_law import pair_full_scalar  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "C"}

ZN = [176, 179, 238, 242]
DNA_CONTACT = [248, 273, 280, 120, 241, 277]
HOTSPOTS = [175, 245, 248, 249, 273, 282]


def parse_ca_resnums(text, chain):
    seq, nums, seen = [], [], set()
    for line in text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[21] != chain:
            continue
        num = line[22:26].strip()
        if (chain, num) in seen:
            continue
        seen.add((chain, num))
        aa = AA3.get(line[17:20].strip())
        if aa:
            seq.append(aa)
            nums.append(int(num))
    return "".join(seq), nums


def importance(seq):
    chars = [c for c in seq.upper() if c in "ARNDCEQGHILKMFPSTWYV"]
    n = len(chars)
    ops = [aa_opcode(c) for c in chars]
    spins = [op.spin() for op in ops]
    charges = [op.charge() for op in ops]
    branches = [op.branch for op in ops]
    aros = [op.aromatic for op in ops]
    props = [SsPropensity.from_amino_acid(c) for c in chars]
    gate = int(LONG_RANGE_GATE)
    acc = np.zeros(n)
    cnt = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if abs(i - j) < gate:
                continue
            fs = pair_full_scalar(
                abs(i - j), spins[i], spins[j], charges[i], charges[j],
                branch_i=branches[i], branch_j=branches[j], aro_i=aros[i], aro_j=aros[j],
                long_range_gate=gate, chain_len=n, recent_hits=0.0, aa1=chars[i], aa2=chars[j],
                p_alpha_i=props[i].p_alpha, p_alpha_j=props[j].p_alpha,
                p_beta_i=props[i].p_beta, p_beta_j=props[j].p_beta)
            acc[i] += abs(float(fs["S"]))
            cnt[i] += 1
    return acc / np.maximum(cnt, 1)


def enrich(imp, idxs, label, rng):
    pct = np.array([(imp < imp[i]).mean() for i in range(len(imp))])  # 0..1, high=important
    sel = [i for i in idxs if 0 <= i < len(imp)]
    obs = float(np.mean([pct[i] for i in sel]))
    null = np.array([np.mean(rng.choice(pct, size=len(sel), replace=False)) for _ in range(20000)])
    p = float((null >= obs).mean())
    print(f"{label:<26} n={len(sel):2d}  mean-percentile={obs*100:5.1f}%  perm-p={p:.4f}"
          f"  {'ENRICHED' if p < 0.05 else ''}")


def enrich_typed(imp, seq, idxs, label, rng):
    """AA-type-controlled null: each functional residue vs same-AA residues only."""
    pct = np.array([(imp < imp[i]).mean() for i in range(len(imp))])
    sel = [i for i in idxs if 0 <= i < len(imp)]
    pool = {}
    for i, c in enumerate(seq):
        pool.setdefault(c, []).append(i)
    obs = float(np.mean([pct[i] for i in sel]))
    null = np.empty(20000)
    for t in range(20000):
        picks = [rng.choice(pool[seq[i]]) for i in sel]
        null[t] = np.mean([pct[i] for i in picks])
    p = float((null >= obs).mean())
    print(f"{label:<26} n={len(sel):2d}  obs={obs*100:5.1f}%  same-AA-null p={p:.4f}"
          f"  {'POSITIONAL' if p < 0.05 else 'type-driven'}")


def main() -> int:
    text, _ = None, None
    r = fetch_pdb("1TUP", "A", CACHE)
    # fetch_pdb returns (seq, xyz); re-read raw for residue numbers
    raw = (CACHE / "1TUP.pdb").read_text(encoding="utf-8", errors="replace")
    seq, nums = parse_ca_resnums(raw, "A")
    num_to_idx = {num: i for i, num in enumerate(nums)}
    print(f"p53 core 1TUP/A  len={len(seq)}  resnum {nums[0]}..{nums[-1]}")
    imp = importance(seq)
    pct = np.array([(imp < imp[i]).mean() for i in range(len(imp))])
    order = np.argsort(-imp)
    funcset = set(ZN + DNA_CONTACT + HOTSPOTS)
    print("\nTop-15 FSOT-important residues (resnum aa pctile  *=known functional):")
    for i in order[:15]:
        num = nums[i]
        mark = "*" if num in funcset else " "
        print(f"  {num:>4} {seq[i]}  {pct[i]*100:5.1f}%  {mark}")
    print()
    rng = np.random.default_rng(0)
    enrich(imp, [num_to_idx[x] for x in ZN if x in num_to_idx], "Zn node C176/H179/C238/C242", rng)
    enrich(imp, [num_to_idx[x] for x in DNA_CONTACT if x in num_to_idx], "DNA-contact residues", rng)
    enrich(imp, [num_to_idx[x] for x in HOTSPOTS if x in num_to_idx], "Cancer hotspots (Cho94)", rng)
    enrich(imp, [num_to_idx[x] for x in funcset if x in num_to_idx], "All functional (union)", rng)
    print("\nAA-type-controlled (does position beat same-AA residues?):")
    enrich_typed(imp, seq, [num_to_idx[x] for x in ZN if x in num_to_idx], "Zn node (vs other C/H)", rng)
    enrich_typed(imp, seq, [num_to_idx[x] for x in HOTSPOTS if x in num_to_idx], "Cancer hotspots (vs same-AA)", rng)
    enrich_typed(imp, seq, [num_to_idx[x] for x in funcset if x in num_to_idx], "All functional (vs same-AA)", rng)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
