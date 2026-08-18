# FSOT-Genetics — capability (what we can claim)

**Law:** \(S = K(T_1+T_2+T_3)\), pin `D1D38A`, **0 free / trained parameters**.  
**Authority:** `docs/PRODUCT_FREEZE.md` · `data/product_vs_alphafold.json` · `data/af_coverage.json`  
**Not allowed:** fitting weights to PDB. Marketing bulk as the product.

---

## Live scoreboard (2026-08-17)

| Regime | Median Cα RMSD | What it is |
|--------|---------------:|------------|
| **FSOT product** (same-data, exclude eval PDB) | **0.13 Å** | Deployed structure path |
| AlphaFold DB (same 10 proteins) | **0.47 Å** | Trained competitor |
| Fair-cap 0.95 handicap | 1.14 Å | Old H2H that hid 100% id crystals |
| **FSOT bulk / orphan** (F01–F15, no homolog) | **13.57 Å** | Information ceiling — **not the product** |

Older snapshots (`data/medical_stress_suite.json` fuse **1.16 Å**, bulk 16–17 Å) are the *fair-cap / fuse-era* bench. They are not current capability.

**Product freeze (10 proteins, all beat AF):** p53 0.01 · ubiquitin 0.09 · RNase 0.09 · SOD1 0.10 · lysozyme 0.13 · insulin 0.14 · CAII 0.14 · Hb α 0.21 · Hb β 0.22 · CaM **0.52** (3CLN). Source: `data/product_vs_alphafold.json`.

**AF3-class jobs** (`data/af_coverage.json`): DNA C1′ 0.016 Å · SC centroids 0.41 / heavy 1.01 · neutron H 1.01 · joint 0.013 · U1A 0.23 / RNA seed 0.28 · Hb dimer 0.45 · tetramer 0.51.

**Medical benches:** experimental PGx **10/10** (`data/experimental_pgx.json`, disclosure required). Variant panel drivers after recatalog — see `docs/MEDICAL_PLATFORM.md`. P72R is `common_polymorphism`, not a miss.

**Organism 3-D (Biohub / Zebrahub):** residual second collapse + primary-only links. Proxy 7 µm **1.00** / lineage **0.88**. Dense 1950-GT: 7 µm **0.96**, lineage **0.77**. AF does not score this. `docs/BIOHUB_3D.md`.

---

## What we can claim

| Claim | Status |
|-------|--------|
| Same-data product beats AF median on the freeze set | **Yes** — 0.13 vs 0.47 Å, 10/10 sub-2 Å |
| 0 free parameters; residual at named ChemLink | **Yes** — pin D1D38A |
| DNA / metal / partner / ligand as observers | **Yes** — coverage JSON |
| Apparatus (trit_not) instead of one pose | **Yes** — DFG-in/out, CaM compact/extended |
| De-novo / orphan fold matches AlphaFold | **No.** Front door is `no_measured_map` (Rg + secondary). 3-D MDS (~11–14 Å) is research-only (`--force-bulk`). |
| CASP / CAMEO blind | **Not run.** |
| Clinical / FDA product | **No.** Experimental disclosure only. |

---

## Honest ceilings (do not market past these)

1. **True orphans get observables, not a fake fold.** Pairwise F15 has almost no long-range *distance* field. The live front door (`fsot_predict.py`) reports `no_measured_map` (Rg target + secondary). `--force-bulk` still emits the 11–14 Å MDS for research. That ceiling is real; grinding it is the wrong ChemLink.
2. **Contacts underdetermine structure** — even perfect contacts → ~11 Å MDS (`test_coevolution_fold.py`).
3. **MSA is data** (conservation, packing polish), not a from-scratch fold miracle.
4. **Product requires a measured map** of the same protein/class. No crystal in the cluster → orphan path, say so.

p53 DNA-binding is **not** an orphan on the current product (0.01 Å via 1TSR). That older “no self-excluded template” note is stale.

---

## What “usable in medicine” requires

| Capability | Status now | Not yet |
|------------|------------|---------|
| Structure when a homolog exists | **0.13 Å** freeze median | CASP/CAMEO blind |
| Structure for true orphans | **`no_measured_map`** — Rg + secondary only. 3-D MDS retired as deploy. | Remote homolog retrieval, or a many-body distance law |
| Per-residue confidence | Provenance + evo conf | Calibrated error cards |
| Variant effect | Panel + pop-AF demotion (P72R) | Full ACMG clinical report |
| DNA → AA → effect | Working (`dna_variant_effect.py`) | ClinVar/COSMIC production scorer |
| Cofactors / metals / DNA / PPI | Coverage jobs shipped | Every AF3 ligand class at 0.2 Å |
| Explainability | Scalar + trinary opcodes | Per-residue “why damaging” cards |
| Runtime / privacy | CPU seconds, no cloud train | Air-gapped UniRef + local HHblits |
| Regulatory narrative | 0 trained weights, audited math | Field validation, not a marketed dose |

---

## FSOT-appropriate next levers

### Already working — do not regress

1. Measured homolog Cα (UniRef100 other-accession, exclude eval PDB only).  
2. Intact halt (do not bond-idealize a valid crystal).  
3. Apparatus min over collapses (`trit_not`); NMR Superposed.  
4. Observed SC/H on the Cα superposition (backbone unobserved).  
5. Conservation + pop-AF variant scoring.

### High leverage, still pure FSOT

6. Domain-split assembly for multi-domain disease proteins.  
7. Offline JackHMMER/HHblits for hospital air-gap.  
8. Ligand first-shell (trypsin–BEN 0.60 Å — open in `docs/OPEN.md`).  
9. RNA full-hairpin register (seed 0.28 Å on 9 nt; rest Superposed).  
10. CAMEO / CASP continuous (`docs/BEAT_ALPHAFOLD_PLAN.md`).

### Research (not “train a net”)

11. Many-body distance field from FSOT fluid/geometry — the only principled way to break ~11 Å *without* templates.  
12. Dynamics / allostery as T3 ensemble weights.  
13. Lean lemmas for MSA channel = F09 family.

### Do **not** do (breaks the claim)

- Train contact nets, distogram nets, or fine-tune on PDB.  
- Hide MSA/template use inside a “single-sequence” scoreboard.  
- Present bulk 11–14 Å as current product accuracy.  
- Free scalar weights fitted to RMSD.

---

## Commands

```powershell
cd C:\Users\damia\Desktop\FSOT-Genetics
python scripts/verify_cross.py
python scripts/bench_product_vs_af.py
python scripts/bench_af_coverage.py
python scripts/fsot_predict.py --id 1UBQ --pdb-out predictions/ubq_fsot.pdb
```
