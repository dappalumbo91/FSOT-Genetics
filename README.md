# FSOT-Genetics

**Mathematical genetics & protein structure under Fluid Space-Time Omni Theory (FSOT).**  
Zero free parameters. No neural-network claim path.

**Repository:** [github.com/dappalumbo91/FSOT-Genetics](https://github.com/dappalumbo91/FSOT-Genetics)  
**Authority pin:** `D1D38A` (`vendor/fsot_compute.py`)  
**Law:** \(S = K(T_1+T_2+T_3)\)

**Mathematical authority:** [FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean). This repository derives its genetics and protein formulas from that formal hub and carries a byte-identical `D1D38A` scalar engine pin.

---

## Why this repo exists

[FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean) is the multi-domain formal hub.  
**This repo** is the **genetics / structure formula branch**: sequence → F01–F15 → Cα coordinates, scored against experimental PDB and AlphaFold DB — still **0 free parameters**.

| | AlphaFold | FSOT-Genetics (this repo) |
|--|-----------|---------------|
| Free / trained parameters | ~tens of millions of weights | **0** |
| Method | Trained network + MSA | Closed-form scalar law; real data used as **input**, never training |
| Cα RMSD — **product** (multi-template + residual physics) | ~0.47 Å median | **1.16 Å median; 10/10 within 1.5 Å of AF; 9/10 sub-2 Å** |
| Cα RMSD — **de-novo, single sequence** | — | **~11–14 Å** (proven single-sequence information ceiling) |

**Product freeze (2026-08-11):** further AF RMSD grinding paused.  
Runtime moves to **Zig bare metal** (`zig/` host + QEMU Multiboot).  
See `docs/PRODUCT_FREEZE.md` and `docs/BARE_METAL_GENETICS_ROADMAP.md`.

We are **not** training nets. We use the zero-parameter scalar law as the map and
real observed data (homolog structures, evolutionary alignments) as **input** —
fully deterministic, auditable, and reproducible.

---

## Current results (validated, zero trained weights)

**Protein structure (Cα RMSD to experimental PDB):**
- **AlphaFold head-to-head** (10 classic proteins, `data/alphafold_headtohead.json`): FSOT-template **1.2 Å** median vs AlphaFold 0.47 Å — **9/10 within 1.5 Å**, 8/10 sub-2 Å. FSOT-template *beats* AlphaFold on flexible calmodulin (0.77 Å vs 6.45 Å).
- **60-chain live holdout** (`data/rcsb_template_holdout_eval.json`): median **2.20 Å**, 54/60 template-covered.
- **De-novo single sequence**: ~11 Å — the proven information ceiling (native full-distance reconstruction → 0.0 Å; perfect contacts → 11 Å).

**RNA** (`scripts/run_rna_template_probe.py`): homolog C1′ transfer → **0.6–0.9 Å** near-native. Flexible multi-domain RNA interdomain angle solved from FSOT coaxial-stacking (`scripts/solve_hinge_angle.py`): 38.7 → 10.4 Å native-free.

**Per-residue confidence** (pLDDT analog, `scripts/test_confidence.py`): identical 0.70 Å / mutated 0.93 Å / gap 2.91 Å — reliably flags low-confidence regions.

**Cofactor chemistry** (`scripts/cofactor_nodes.py`, `formulas/cofactor_fsot_map.json`): every metal/cofactor routed to its FSOT physics domain (Zn/Cu/Fe → Atomic_Physics D7; Ca/Mg → Electromagnetism D9), validated against textbook active sites.

**Medical — variant-effect prediction** (`scripts/variant_conservation.py`, `scripts/dna_variant_effect.py`): evolutionary conservation (diverse Pfam MSA) + SIFT-style substitution-specificity flags all six p53 cancer-driver mutations at mean **84th percentile**, zero trained weights, fully explainable. A DNA coding variant is routed through the trinary codon layer (`c.742C>T → CGG→TGG → R248W → LIKELY DAMAGING`).

---

## Quick start

### Shipping runtime (Zig — preferred)

```powershell
cd zig
zig build host          # product residual + codon + scalar gate
zig build kernel        # freestanding Multiboot image
.\run_qemu.ps1          # QEMU serial gate (if QEMU installed)
```

Host residual must match pin: `r_bond≈1.100 r_clash≈1.122 r_anchor≈1.092`.

### Research oracle (Python — metrics only)

```powershell
git clone https://github.com/dappalumbo91/FSOT-Genetics.git
cd FSOT-Genetics
python -m pip install -r requirements.txt

# Cross-verification (must pass — Lean-style gate)
python scripts/verify_cross.py

# Head-to-head vs AlphaFold (needs network for UniProt/PDB/AF)
python scripts/run_alphafold_headtohead.py

# Medical: variant-effect prediction (conservation + substitution-specificity)
python scripts/variant_conservation.py

# Medical: DNA coding variant -> codon -> amino acid -> effect (trinary front door)
python scripts/dna_variant_effect.py

# Distogram contact metrics (F15 design metrics)
python scripts/run_fsot_distogram_contact_eval.py

# Dual-mode fold: pure single-sequence vs optional MSA-augmented F15
python scripts/run_msa_augmented_fold.py --id 1UBQ

# Frozen real-data reproduction audit (downloads 12 experimental RCSB structures)
python scripts/run_rcsb_holdout.py
```

**MSA (optional usability layer):** default fold is pure single-sequence (published
claims). Pass `mode="msa"` or use `scripts/run_msa_augmented_fold.py` to inject
Pfam / local JackHMMER–HHblits / file MSA features into long-range F15. See
`docs/MSA_AUGMENTATION.md`. MSA is **data input**, not training — still 0 free parameters.

**Medical front door + stress suite:**

```powershell
python scripts/fsot_predict.py --id 1UBQ --pdb-out model.pdb
python scripts/run_medical_stress_suite.py   # AF vs template vs fuse vs bulk±MSA
```

Latest multi-regime stress (`data/medical_stress_suite.json`): FSOT **template+MSA
packing fuse median ~1.16 Å** (7/9 beats raw template); bulk remains the orphan
fallback. Roadmap: `docs/CAPABILITY_ROADMAP.md`.

**Multi-gene medical variant panel** (TP53, KRAS, EGFR, BRAF, CFTR, SOD1, HBB, BRCA1):

```powershell
python scripts/run_medical_variant_panel.py
```

**UniRef50 protein-specific MSAs** (real UniProt homolog clusters) + Pfam gap-fill +
absolute evolutionary gates (seed-closed). Latest scoreboard: **34/35 drivers
LIKELY DAMAGING (97%)**. See `data/medical_variant_panel.json`.

**Domain-split + joint multi-domain templates:**

```powershell
python scripts/domain_split_assemble.py TP53 KRAS SOD1 EGFR
```

When a single homolog covers multiple domains, experimental inter-domain pose is
kept; otherwise per-domain templates + FSOT interface pad.


The RCSB manifest is a no-tuning holdout: do not select formulas, routings, or
thresholds from its outcomes. After a formula change is frozen in Git, validate
it on a newly preregistered holdout rather than re-optimizing against this set.

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

### Runtime honesty

| Layer | Status |
|-------|--------|
| Trinary **syntax / law** | Yes — opcodes + pair geometry |
| Python H2H runner | **Floats emulating** the law (lab path) |
| Zig / trit bare metal | Twin sources under `zig/` — **not** yet the H2H claim runtime |

See `docs/RUNTIME_STACK.md`. Bare metal is the direction; Python is the current formula lab.

### D_eff interfaces

Protein fold is multi-scale. We route **named pin domains** only (no free D):

```powershell
python scripts/domain_interface.py
python scripts/run_deff_interface_probe.py
```

Docs: `docs/DOMAIN_INTERFACE_FOLD.md`

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

Apache License 2.0. See [LICENSE](LICENSE).

## Author

Damian Arthur Palumbo — FSOT (Fluid Space-Time Omni Theory).
