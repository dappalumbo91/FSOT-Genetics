# FSOT-Genetics

**Mathematical genetics & protein structure under Fluid Space-Time Omni Theory (FSOT).**  
Zero free parameters. No neural-network claim path.

**Repository:** [github.com/dappalumbo91/FSOT-Genetics](https://github.com/dappalumbo91/FSOT-Genetics)  
**Authority pin:** `D1D38A` (`vendor/fsot_compute.py`)  
**Law:** \(S = K(T_1+T_2+T_3)\)

---

## Why this repo exists

[FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean) is the multi-domain formal hub.  
**This repo** is the **genetics / structure formula branch**: sequence → F01–F15 → Cα coordinates, scored against experimental PDB and AlphaFold DB — still **0 free parameters**.

| | AlphaFold | FSOT-Genetics |
|--|-----------|---------------|
| Free parameters | ~tens of millions of weights | **0** |
| Method | Trained network | Closed-form seeds + domain scalars |
| Home PC fold time (this path) | N/A (precomputed DB) | **~0.1–0.4 s / chain** formula path |
| Accuracy (current H2H median Cα RMSD) | **~0.4 Å** (wins most) | **~15 Å** (open — next work) |

We are **not** training nets. We are **solving** with math and publishing honest scoreboards.

---

## Quick start

```powershell
git clone https://github.com/dappalumbo91/FSOT-Genetics.git
cd FSOT-Genetics

# Cross-verification (must pass — Lean-style gate)
python scripts/verify_cross.py

# Head-to-head vs AlphaFold (needs network for UniProt/PDB/AF)
python scripts/run_fsot_vs_alphafold_structure.py --max-proteins 8 --rounds 24 --sleep 0.2

# Distogram contact metrics (F15 design metrics)
python scripts/run_fsot_distogram_contact_eval.py
```

Optional Rust formula crates:

```powershell
cargo test --workspace
```

---

## Layout

```
vendor/fsot_compute.py          # D1D38A scalar engine (byte-pinned)
scripts/trinary_syntax.py        # 7-trit AA opcodes + codon syntax + Zig pair law
scripts/fsot_structure_engine.py # F01–F15 + expanded trinary + MDS + sparse polish
scripts/run_fsot_vs_alphafold_*  # dual scoreboard
scripts/verify_cross.py          # hard green gate (includes 20/20 unique opcodes)
formulas/                        # F01–F15 + expanded trinary maps
zig/                             # neuron-zig genetics twin (seeds, trit, codon, pair)
crates/                          # codon_core, fsot_core, fsot_protein (Rust)
docs/DESIGN.md                   # architecture
docs/TRINARY_SYNTAX_EXPANSION.md # genetics-as-code syntax
docs/BEAT_ALPHAFOLD_PLAN.md      # campaign bar
docs/CROSS_VERIFICATION.md       # gate policy
data/                            # scoreboard snapshots + sample PDBs
```

### Genetics as code (trinary)

```text
codon → PRIMARY/SECONDARY trits → AA → 7-trit opcode → spin/charge → F15 pair → fold
```

F01 alone collides 6 AA groups. Expanded opcodes are **20/20 unique** (still zero free params).

---

## Cross-verification (required)

Same spirit as Lean green gates in FSOT-2.1-Lean:

1. SHA-256 of `vendor/fsot_compute.py` starts with **D1D38A**
2. Engine `free_parameters == 0`
3. Formula fold finishes under hard time gate
4. Finite Cα coordinates
5. Derivations document present

```text
python scripts/verify_cross.py   → exit 0 only if all pass
```

CI: `.github/workflows/ci.yml`

---

## Related embodiments

| Project | Role |
|---------|------|
| [FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean) | Law + multi-domain verification |
| [fsot-neuron-zig](https://github.com/dappalumbo91/fsot-neuron-zig) | Neural mind; genetic pair geometry |
| [FSOT-2.1-Neural](https://github.com/dappalumbo91/FSOT-2.1-Neural) | Wet-lab / neural monorepo |
| **FSOT-Genetics** | Genetics structure formula branch (this repo) |

---

## License

MIT (unless otherwise noted in subcrates).

## Author

Damian Arthur Palumbo — FSOT (Fluid Space-Time Omni Theory).
