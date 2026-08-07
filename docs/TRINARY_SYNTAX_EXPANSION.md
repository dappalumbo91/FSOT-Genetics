# Trinary syntax expansion — genetics as code

## Idea

Genetics is **executable code**. Codons are instructions. Amino acids are
**opcodes**. Structure is what happens when that code runs under FSOT law.

F01 alone maps each AA to a 3-trit word `(charge, polarity, volume)`.  
That is necessary but not sufficient: five large non-polars share one word.

| F01 phase | AAs (collision) |
|-----------|-----------------|
| `[0,-1,+1]` | I, L, M, F, W |
| `[0,-1,-1]` | A, G |
| `[1,+1,+1]` | R, H, K |
| … | … |

## Expansion (7-trit opcode)

```text
word = (c, p, v, aromatic, branch, hetero, detail)
```

| Trit | Meaning | Lawful basis |
|------|---------|--------------|
| c | charge | F01 / chemistry |
| p | polarity | F01 |
| v | volume class | F01 |
| aromatic | ring side chain | F,Y,W = +1 |
| branch | side-chain topology | I,V,L,T = +1; P = −1 |
| hetero | S / OH / special | C,M,W,H = +1; S,T,Y,G = −1 |
| detail | isomer / FG detail | I≠L (β/γ), R≠K (guanidino/amine), … |

**Not free parameters** — structural chemistry facts encoded as trits, same
spirit as encoding charge as a trit.

Continuous fields used in pair geometry:

- `spin` = mean(p, v, branch, aromatic) ∈ [−1, 1]
- `charge` = c + hetero/φ
- `h`, `vol` refined from F02 with aromatic/branch/hetero multipliers from {π,e,φ,γ}

## Zig twin (neuron)

```text
geometricScaleDist(d) = φ · d^(-1/π)
fsotPairWeight = geom · (trinaryPair + 0.15·elec) · (0.35 + 0.65·env)
env(d) = d / (d + πe)
```

Sources under `zig/src/` (pulled from `fsot-neuron-zig`).

## Pipeline

```text
DNA codon
  → PRIMARY / SECONDARY trit trips   (codon.zig)
  → AA
  → 6-trit opcode                     (trinary_syntax.py)
  → spin, charge, h, vol
  → F15 proximity M_ij (+ Zig pair)
  → distance D → MDS → Cα
```

## Commands

```powershell
python scripts/trinary_syntax.py     # write expanded maps
python scripts/verify_cross.py       # uniqueness + pin + fold
python scripts/run_fsot_vs_alphafold_structure.py --max-proteins 8 --rounds 24
```

## Files

| Path | Role |
|------|------|
| `scripts/trinary_syntax.py` | Expansion + pair law + map writers |
| `formulas/20_amino_acid_expanded_trinary.txt` | Human map |
| `formulas/64_codon_expanded_trinary.txt` | Codon → opcode |
| `formulas/trinary_syntax_expanded.json` | Machine map |
| `zig/src/genetic_pair.zig` | Standalone Zig pair law |
| `zig/docs/GENETICS_AS_TRINARY_CODE.md` | Doctrine |
