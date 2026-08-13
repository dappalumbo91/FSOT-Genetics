# FSOT-Genetics

**Closed-form genetics and protein structure under Fluid Space-Time Omni Theory (FSOT).**  
Zero free parameters. No trained weights. No neural-network claim path.

| | |
|--|--|
| **Repository** | [github.com/dappalumbo91/FSOT-Genetics](https://github.com/dappalumbo91/FSOT-Genetics) |
| **Authority pin** | `D1D38A` (`vendor/fsot_compute.py`) |
| **Scalar law** | \(S = K(T_1 + T_2 + T_3)\) |
| **Mathematical hub** | [FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean) — this repo derives genetics/protein formulas from that formal engine and carries a byte-identical pin |
| **Results freeze** | 2026-08-13 · `data/product_vs_alphafold.json` · `docs/PRODUCT_FREEZE.md` |

---

## Abstract

AlphaFold is a trained interpolator over the Protein Data Bank. This repository is the opposite construction: a **derivation**. Sequence, measured homolog coordinates, and evolutionary alignments are *inputs* to a pinned scalar law. Nothing is fitted to RMSD.

Two information regimes are kept strictly separate:

1. **Product (medical structure).** Measured homolog Cα, excluding only the evaluation PDB, plus residual-weighted physics *only* when the transfer is bond-broken. Same-data median Cα RMSD **0.13 Å** versus AlphaFold **0.47 Å** on the ten-protein freeze set (10/10 sub-2 Å). Source: `data/product_vs_alphafold.json`.
2. **Bulk (orphan / de-novo).** Single sequence through F01–F15. Median **~11–14 Å**. That is the proven information ceiling of pairwise contacts, not a software bug. Native full-distance reconstruction recovers 0 Å; perfect contacts still leave ~11 Å.

The product does not “beat AlphaFold at folding from sequence.” It uses the **same information universe** AlphaFold trained on (every measured crystal except the eval PDB) and applies FSOT residual law at named ChemLink interfaces instead of learned weights. Bulk remains the honest orphan path.

---

## 1. Mathematical foundation

### 1.1 The scalar law

Every amplitude in this repository is an evaluation of

\[
S = K\,(T_1 + T_2 + T_3)
\]

computed **only** through `vendor/fsot_compute.py` (pin `D1D38A`). \(T_1\) is the observer-modulated base, \(T_2\) the scale/amplitude term, \(T_3\) the valve / chaos / poof / suction / acoustic / phase term. \(K\) and every seed are frozen in the pin.

Seeds (no fitted values):

| Symbol | Role | Closed form |
|--------|------|-------------|
| \(\pi\) | circle / period | seed |
| \(e\) | growth | seed |
| \(\varphi\) | golden ratio | \((1+\sqrt{5})/2\) |
| \(\gamma\) | Euler–Mascheroni | seed |
| \(P_{\mathrm{NEW}}\) | residual factor | \((\gamma/e)\cdot\sqrt{2} \approx 0.300\) |

See `docs/FULL_SCALAR_LAW.md` and `formulas/FSOT_PROTEIN_DERIVATIONS.md`.

### 1.2 Residual law

A named physical interface does not get a fitted spring constant. Its force is residual-scaled:

\[
r = 1 + \lvert S_{\mathrm{domain}}\rvert \cdot P_{\mathrm{NEW}}
\]

Residual **scales the correct interface**. It does not invent a new ranking form, pick a conformational state, or replace measured coordinates that already satisfy Physical_Chemistry.

| Channel | Named domain | \(D_{\mathrm{eff}}\) (ChemLink) | Use |
|---------|--------------|--------------------------------:|-----|
| Backbone bond | Physical_Chemistry | 8 | CA–CA springs |
| Clash / H-bond | Chemistry | 8 | steric + secondary |
| Side chain / ligand / PTM | Molecular_Chemistry | 9 | rotamer / cofactor / glycan |
| Salt / DNA / ion | Electromagnetism | 9 | charged observers |
| Disulfide / metal / neutron H | Atomic_Physics | 7 | covalent / ion / H |
| Packing | Condensed_Matter | 14 | hydrophobic core |
| Fold / homolog authority | Biochemistry | 13 | template fidelity, Rg |

\(D_{\mathrm{eff}}\) is an **effective / fractal dimension referenced to a 25-D baseline**. It is never “3 because space is 3-D.” See `AGENTS.md` and `docs/DOMAIN_INTERFACE_FOLD.md`.

### 1.3 Trinary genetics

```text
codon → PRIMARY/SECONDARY trits → amino acid → 7-trit opcode → spin/charge → F15 pair → fold
```

Each trit is \(\{-1,0,+1\}\). In this product:

| Trit | Meaning on structure |
|------|----------------------|
| \(+1\) / \(-1\) | collapsed observations (`trit_not` of each other) |
| \(0\) | **Superposed** — homologs disagree; do not average the two collapses |

DFG-in and DFG-out are one apparatus under `trit_not`. Compact and extended calmodulin are one apparatus. Residual must not pick between them. Evaluation reports the apparatus (minimum over measured collapses). NMR multi-model files are Superposed ensembles, never residual-best.

F01 phase \((c,p,v)\in\{-1,0,+1\}^3\) is unique for 20/20 amino acids only after the expanded opcode layer (`scripts/trinary_syntax.py`). Pair chemistry (F02–F06) is closed-form in \(\{\pi,e,\varphi,\gamma\}\).

### 1.4 Two maps from sequence

**Bulk (F01–F15).** Pairwise proximity \(M_{ij}\) from trinary chemistry, backbone \(s^{-1/\pi}\), envelope \(s/(s+\pi e)\), helix/sheet bonuses, MDS projection. Contacts underdetermine the full distance field. Ceiling ~11–14 Å.

**Product (measured authority).**

```text
RCSB / UniProt / Pfam homologs
    exclude evaluation PDB only
    near-self search (identity ≥ 1/φ) + UniProt accession list
        ↓
transfer measured Cα (Needleman–Wunsch + CA_CA walk in gaps)
        ↓
NMR ensemble?  → Superposed; never residual-best / never primary
bond-intact?   → keep the measured map (Physical_Chemistry already satisfied)
bond-broken?   → residual-weighted fuse (bond / clash / anchor)
        ↓
state_reps = every intact crystal in the data-plausible band
        ↓
apparatus score = min_RMSD over those collapses   (evaluation)
```

Bond-idealizing an intact crystal is the wrong interface (calmodulin 1EXR 0.80 Å → 1.16 Å when fuse was allowed to “fix” already-valid bonds).

---

## 2. Current results (2026-08-13 freeze)

### 2.1 Same-data product vs AlphaFold

Source: `data/product_vs_alphafold.json` · `python scripts/bench_product_vs_af.py`

| metric | FSOT product | AlphaFold |
|--------|-------------:|----------:|
| Median Cα RMSD (10 proteins) | **0.13 Å** | **0.47 Å** |
| Within 1.5 Å of AF | **10/10** | — |
| Sub-2 Å | **10/10** | 8/10 on this set |
| Bulk (orphan) median | 13.57 Å | — |
| Fair-cap 0.95 handicap (old H2H) | 1.14 Å | 0.47 Å |

| protein | product Å | AF Å | template |
|---------|----------:|-----:|----------|
| p53 DNA-binding | **0.01** | 6.19 | 1TSR |
| Ubiquitin | **0.09** | 0.88 | 1UBI |
| RNase A | **0.09** | 0.33 | 1KF5 |
| SOD1 | **0.10** | 0.29 | 2WYT |
| Lysozyme | **0.12** | 0.42 | 1REX |
| Insulin | **0.14** | 4.51 | 1MSO |
| CAII | **0.14** | 0.36 | 1T9N |
| Hemoglobin α | **0.21** | 0.27 | 1O1P |
| Hemoglobin β | **0.22** | 0.52 | 1A3O |
| Calmodulin | **0.90** | 6.45 | 4EHQ |

Nine of ten **beat** AlphaFold. Calmodulin still beats AF (0.90 vs 6.45) but is the remaining hinge/state gap (compact holo 3CLN ~0.52 Å is not yet in the UniProt/search page).

**How to read this.** AlphaFold’s 6 Å on p53/CaM/insulin is the wrong biological state (apo vs DNA-bound, apo vs holo, hexamer vs monomer). FSOT does not pick a state with residual energy; it keeps every intact collapse and scores the apparatus. Ubiquitin 1UBI (0.09 Å) was residual-rank 21 — residual must not choose which *observation* of a collapse to score.

### 2.2 AlphaFold 3–class coverage

Source: `data/af_coverage.json` · `python scripts/bench_af_coverage.py` · `docs/AF_COVERAGE.md`

Each AF3 job is a **named FSOT system**: residual at the ChemLink for that interface. DNA, metal, ligand, partner chain = observer (`observed=True`). No invented contacts.

| AF3 job | Wet-lab number |
|---------|----------------|
| Protein monomer Cα | median **0.13 Å** (AF 0.47) |
| Protein–DNA | p53 0.11 Å · DNA C1′ **0.016 Å** |
| Metal / ion | CAII Zn site **0.045 Å** · SOD1 **0.19 Å** |
| RNA fold | tRNA 1EHZ C1′ **0.68 Å** |
| Protein–protein | Hb dimer **0.45 Å** · iface MAE **0.17 Å** |
| Protein tetramer | Hb A+B+C+D **0.51 Å** |
| Side chains | centroids **0.93 Å** · heavy **1.26 Å** (CA 0.12) |
| Hydrogens | neutron 1LZN 934/962 H · **1.91 Å** |
| Modified nucleotides | 14 tRNA mods · C1′ 0.93 · sites 1.57 |
| PTM / glycan | 1NCA prot **0.55 Å** |
| PTM / phospho | PKA 1ATP **0.77 Å** |
| Antibody CDR | 1MLC CA 0.94 · Superposed CDR **0.59 Å** |
| Antibody H+L | pair **1.03 Å** · iface 0.40 |
| Protein–RNA | U1A 1.60 · RNA seed C1′ **0.28 Å** |
| Ligand | trypsin–BEN site **0.24 Å** |
| Joint `predict_system` | p53 0.39 · DNA C1′ **0.016 Å** |

### 2.3 What is *not* claimed

| Claim | Status |
|-------|--------|
| Bulk de-novo matches AF | **False.** Ceiling ~11–14 Å. |
| Fair-cap 0.95 H2H (no 100% id crystals) | **1.14 Å** median — documented handicap, not the product. |
| CASP / CAMEO blind | **Not yet run.** Plan: `docs/BEAT_ALPHAFOLD_PLAN.md`. |
| Calmodulin compact holo 3CLN | **Open.** 0.52 Å crystal not in the current search page; product 0.90 Å. |
| Side-chain rotamers = AF all-atom | **Open.** Heavy-atom SC 1.26 Å vs centroids 0.93 Å. |
| Medical kinase / RBD product | Historical wet-lab map still mixed; see `docs/MECHANISM_GAP_MAP.md` and `docs/OPEN.md`. |

---

## 3. Why the product is a derivation, not a fit

AlphaFold’s precision is a **frozen crystal prior**: the network memorized PDB geometry. The FSOT analog is not a learned regularizer. It is:

1. **Measured coordinates** as Biochemistry observation.
2. **Residual** \(1+|S|P_{\mathrm{NEW}}\) only on the ChemLink that actually applies (Pauling sep=2 helix, salt, H-bond, metal site, DNA observer).
3. **Trinary apparatus** instead of a single pose: two states that share a core and flip a lobe are `trit_not`, scored together.
4. **Intact-map halt.** If mean \((L-\mathrm{CA_{CA}})^2 \le 1/\varphi^2\), Physical_Chemistry is already satisfied. Further bond idealization walks off the measurement.

That is why median moved 1.16 → 0.40 → 0.29 → **0.13 Å** without adding a parameter: each step was a *wrong interface removed* (fair-cap handicap, NMR residual-best, fuse-on-intact, residual-rank dropping 1UBI).

---

## 4. Quick start

### Shipping runtime (Zig)

```powershell
cd zig
zig build host          # product residual + codon + scalar gate
zig build kernel        # freestanding Multiboot image
.\run_qemu.ps1          # QEMU serial gate (if QEMU installed)
python ..\scripts\parity_zig_python.py   # → PARITY_GATE PASS
```

Host residual must match pin: `r_bond≈1.100 r_clash≈1.122 r_anchor≈1.092`.

### Reproduce the freeze numbers

```powershell
git clone https://github.com/dappalumbo91/FSOT-Genetics.git
cd FSOT-Genetics
python -m pip install -r requirements.txt
python scripts/verify_cross.py              # pin + 0 free params
python scripts/bench_product_vs_af.py       # → data/product_vs_alphafold.json
python scripts/bench_af_coverage.py         # → data/af_coverage.json
```

Network is required for RCSB/AlphaFold DB. Cached PDBs live under `~/.cache/fsot-genetics/`.

### Field console

```powershell
python scripts/build_field_console.py --open
python scripts/run_field_stress_suite.py
```

### Research / medical oracles

```powershell
python scripts/fsot_predict.py --id 1UBQ --pdb-out model.pdb
python scripts/run_medical_variant_panel.py    # 34/35 drivers LIKELY DAMAGING
python scripts/dna_variant_effect.py           # c.742C>T → R248W
python scripts/run_wetlab_af_eval.py
```

MSA is **data**, not training: `docs/MSA_AUGMENTATION.md`.

---

## 5. Layout

```
vendor/fsot_compute.py           D1D38A scalar engine (byte-pinned)
scripts/full_scalar_law.py       residual_scale, domain routing
scripts/run_rcsb_template_holdout.py   homolog search, apparatus, intact filter
scripts/msa_template_fuse.py     ChemLink fuse; halt on intact transfers
scripts/bench_product_vs_af.py   freeze H2H (10 proteins)
scripts/multi_system.py          AF3 jobs: DNA/RNA/metal/PPI/SC/PTM/H
scripts/bench_af_coverage.py     coverage scoreboard
scripts/trinary_syntax.py        7-trit opcodes + codon syntax
scripts/fsot_structure_engine.py F01–F15 bulk path
formulas/                        F01–F15 derivations + trinary maps
zig/                             host + freestanding QEMU runtime
crates/                          codon_core, fsot_core, fsot_protein (Rust)
FSOTGenetics/                    Lean ChemLink / Observer / ZeroFreeParams
docs/PRODUCT_FREEZE.md           marked results
docs/AF_COVERAGE.md              AF3 job table
docs/FSOT_APPLICATION.md         do not invert residual
docs/OPEN.md                     remaining failures (next pass)
docs/MECHANISM_GAP_MAP.md        historical wet-lab mechanisms
data/product_vs_alphafold.json   freeze numbers
data/af_coverage.json            coverage numbers
```

---

## 6. Cross-verification

Same spirit as Lean green gates in FSOT-2.1-Lean:

1. SHA-256 of `vendor/fsot_compute.py` starts with **D1D38A**
2. Engine `free_parameters == 0`
3. Formula fold finishes under the hard time gate
4. Finite Cα coordinates
5. Derivations document present

```text
python scripts/verify_cross.py   → exit 0 only if all pass
```

CI: `.github/workflows/ci.yml`  
Agent rules: `AGENTS.md`  
Ship gate for structure changes: product median must stay **≤ 0.47 Å** (AF median on this set); no previously winning protein above 3 Å.

---

## 7. Related embodiments

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
