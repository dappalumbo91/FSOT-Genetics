#!/usr/bin/env python3
"""Genetics-native front door: DNA coding variant -> codon -> AA -> effect.

A real clinical variant is a nucleotide change. This routes it through the trinary
codon layer (dna_to_aa / codon_primary, twin of codon.zig), classifies it
(synonymous / missense / nonsense / start-loss), and for missense runs the working
conservation variant-effect predictor. Output is explainable end to end, from the
DNA change to the predicted impact, expressed in trinary codon opcodes.
"""

from __future__ import annotations

import sys
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from trinary_syntax import dna_to_aa, codon_primary  # noqa: E402
from run_fsot_vs_alphafold_structure import fetch_pdb  # noqa: E402
from variant_conservation import conservation_profile, resnums  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
AA20 = "ARNDCQEGHILKMFPSTWYV"

# p53 cancer variants as real DNA changes (protein_pos, wt_codon, codon_pos0, alt, hgvs)
# plus a synonymous control that must read out benign.
VARIANTS = [
    (175, "CGC", 1, "A", "c.524G>A"),   # R175H structural
    (245, "GGC", 0, "A", "c.733G>A"),   # G245S structural
    (248, "CGG", 0, "T", "c.742C>T"),   # R248W DNA-contact
    (273, "CGT", 1, "A", "c.818G>A"),   # R273H DNA-contact
    (282, "CGG", 0, "T", "c.844C>T"),   # R282W structural
    (248, "CGG", 2, "A", "c.744G>A"),   # R248R synonymous control (CGG->CGA, still Arg)
]


def classify(wt_codon, codon_pos, alt):
    mut = list(wt_codon)
    mut[codon_pos] = alt
    mut_codon = "".join(mut)
    wt_aa, mut_aa = dna_to_aa(wt_codon), dna_to_aa(mut_codon)
    if mut_aa == "*" or mut_aa == "Stop":
        kind = "nonsense (truncating)"
    elif wt_aa == mut_aa:
        kind = "synonymous (silent)"
    else:
        kind = "missense"
    return mut_codon, wt_aa, mut_aa, kind


def main() -> int:
    fetch_pdb("1TUP", "A", CACHE)
    raw = (CACHE / "1TUP.pdb").read_text(encoding="utf-8", errors="replace")
    seq, nums = resnums(raw, "A")
    idx = {num: i for i, num in enumerate(nums)}
    cons, freq, nrows, pfam = conservation_profile(seq, "1TUP", pfam="PF00870")

    var_scores = []
    for i in range(len(seq)):
        for m in AA20:
            if m == seq[i]:
                continue
            fm = freq[i].get(m, 0.0) if i < len(freq) else 0.0
            var_scores.append(cons[i] * (1.0 - fm))
    var_scores = np.array(var_scores)

    print(f"TP53 variant interpreter  (Pfam {pfam}, MSA {nrows} seqs)\n")
    for pos, wt_codon, cpos, alt, hgvs in VARIANTS:
        mut_codon, wt_aa, mut_aa, kind = classify(wt_codon, cpos, alt)
        tag = f"{wt_aa}{pos}{mut_aa if mut_aa != wt_aa else '='}"
        line = (f"{hgvs:<10} {wt_codon}->{mut_codon}  trit{list(codon_primary(wt_codon))}->"
                f"{list(codon_primary(mut_codon))}  {tag:<7} {kind}")
        if kind == "missense" and pos in idx:
            i = idx[pos]
            fm = freq[i].get(mut_aa, 0.0) if i < len(freq) else 0.0
            imp = cons[i] * (1.0 - fm)
            pct = float((var_scores < imp).mean()) * 100
            call = "LIKELY DAMAGING" if pct >= 75 else ("uncertain" if pct >= 40 else "likely tolerated")
            line += f"  | conservation={cons[i]:.2f} impact-pctile={pct:.0f}% -> {call}"
        elif kind == "synonymous (silent)":
            line += "  | no AA change -> likely benign"
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
