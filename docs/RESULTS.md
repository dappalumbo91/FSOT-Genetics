# What FSOT-Genetics actually does (live results)

Pin `D1D38A`. Law \(S = K(T_1+T_2+T_3)\). **0 free parameters.**  
Data stays on `D:\FlyWire_Connectome` (not git). Do not mix product Å with hop mass.

## 1. Protein product (accuracy that is Å)

Same-data homolog Cα vs eval PDB. Source: `docs/PRODUCT_FREEZE.md`, `data/product_vs_alphafold.json`.

| Regime | Median Cα RMSD |
|--------|---------------:|
| **FSOT product** | **0.13 Å** (n=10) |
| AlphaFold on the same 10 | 0.47 Å |
| Bulk / no measured map | ~13.6 Å — **not the product** |

Product requires a measured homolog. No map → `no_measured_map` (Rg + secondary only). We do not invent coordinates.

## 2. Cross-species homologs of residual-mass genes

`python scripts/homolog_correspondence.py` → `data/homolog_correspondence.json`

Bee / mosquito / beetle transfers of proteins that already sit on the live fly and worm graphs.

| | n |
|--|--:|
| Measured homologs folded | **26** |
| Same-OrthoDB covers (not a new protein) | 3 |
| True 1:1 miss | **1** (`mec-4` → *Anopheles*) |
| Template identity (24 with a map) | median **0.605**; **10/24** at close-homolog ≥ 1/φ |
| Template coverage | median **0.74** |
| `no_measured_map` | 2 (bee/beetle *mec-4* sequence hits, no structure) |

Notable transfers (same 5VKQ template as fly nompC):

| Species | nompC | identity |
|---------|-------|----------:|
| *A. mellifera* | A0A7M7IF52 | 0.81 |
| *A. gambiae* | AGAP008559 / A0A1S4GZD0 | **0.85** |
| *T. castaneum* | NCBI 662890 / XP_015838654.2 | **0.81** |

A 1:1 miss is not a license to invent. The law is the same; the gene family may not be.

## 3. Measured graphs (accuracy that is a hop split)

Not RMSD. Seed a sensory class → residual hops → does motor/descending light?

Independent sexes, same split:

| Seed | Male CNS hop-2 `vnc_motor` | BANC (female) hop-2 `vnc_motor` |
|------|---------------------------:|--------------------------------:|
| VNC sensory | **10.74** | **10.16** |
| JO | **7.77** | **6.64** |
| olfactory | **0.001** | **0.0008** |

JO hop-1 peaks on **DNg29** (male and BANC). Olfactory stays at antennal-lobe LNs (`il3LN6` / `v2LN30`). Larva: mechanosensory hop-2 **DN-VNC 16.6** vs olfactory **0.43**.

Whole-animal: worm (both sexes) sensory → command interneuron → muscle NMJ; Ciona MGIN; Platynereis prototroch.

## 4. When the 1:1 is blank, math points at the mapped job

`python scripts/analog_pointer.py` → `data/analog_pointer.json`

| Blank 1:1 | Why | Mapped analog (already folded) |
|-----------|-----|--------------------------------|
| *mec-4* in *An. gambiae* | UniRef50/90 is **Nematoda only**; fly *ppk* ~18% | **nompC** AGAP008559 (0.85 on 5VKQ). Worm ALM uses a degenerin; fly walking/JO uses TRPN on the measured graph. |
| *unc-25* UniRef50 in *An. gambiae* | worm GAD cluster | **Gad1** Q7PNL7 (same OrthoDB) |

Do not fold a random mosquito *ppk* as *mec-4*.

## 5. Plants

No fly-class EM connectome. Genome + crystals only.

`python scripts/plant_product.py` — *Arabidopsis thaliana* (proteome `UP000006548`). First panel, all product (not bulk):

| Gene | UniProt | Template | id | cov | System |
|------|---------|----------|---:|----:|--------|
| *rbcL* | O03042 | 5IU0 | **1.00** | 0.96 | Calvin cycle |
| *psbA* | P83755 | 9LK5 | **1.00** | 0.95 | PSII D1 |
| *LHCB1.3* | P04778 | 8J6Z | **1.00** | 0.82 | LHCII antenna |
| *GAPA1* | P25856 | 6KEZ | **1.00** | 0.85 | chloroplast GAPDH |
| *ACT2* | Q96292 | 6IUG | **0.94** | 0.98 | actin |
| *CESA3* | Q941L0 | 8VI0 | **0.80** | 0.67 | cellulose synthase |

Source: `data/plant_product.json`. No plant synapses invented.

`python scripts/plant_homolog.py` — same six systems into rice (Japonica `39947`), maize (`4577`), soybean (`3847` / `UP000008827`), wheat (`4565` / `UP000019116`). **24/24** measured homologs folded, all `template_raw`. True 1:1 miss: **0**.

| Gene | Rice | Maize | Soybean | Wheat |
|------|------|-------|---------|-------|
| *rbcL* | P0C512 **1.00** | P00874 **0.95** | P27066 **0.96** | P11383 **0.98** |
| *psbA* | P0C434 **1.00** | P48183 **1.00** | P02957 **1.00** | P12463 **1.00** |
| *LHCB1.3* | P12330 *CAB1R* **1.00** | P12329 *CAB1* **1.00** | Q43437 *LHCB1-7* **0.97** | A0A3B6AWZ1 **0.95** |
| *GAPA1* | Q7X8A1 **0.92** | P09315 **0.91** | I1N843 **0.94** | A0A3B6B1C5 **0.91** |
| *ACT2* | A3C6D7 **0.95** | P02582 *ACT1* **0.89** | P02581 *SAC1* **0.90** | A0A3B6HWJ7 **0.94** |
| *CESA3* | Q69V23 **0.83** | A0A1D6P4I8 **0.74** | I1KTE1 **0.81** | A0A3B6HZ43 **0.81** |

Named UniProt gene is the measured recover when UniRef50 of chloroplast *rbcL*/*psbA* is a split cluster. Soybean actin is reviewed *SAC1* (not the 336 aa Soy115 fragment). Wheat *acT2* 442 aa has no protein name — dropped; UniRef50 actin used instead. Still not a plant connectome.

Source: `data/plant_homolog.json`.

## What this is capable of

| Can | Cannot |
|-----|--------|
| Product Cα when a homolog exists (0.13 Å freeze) | Orphan 3-D at AlphaFold grade |
| Residual hops on a **measured** synapse graph | Invent a bee / mosquito / plant connectome |
| Point a clade-restricted 1:1 at the mapped residual job | Treat a genome as synapses |
| Same law on animals and plants (proteins) | Skip VNC and call a fly brain a whole animal |

## 6. Multi-prover stamp

`python verification/run_cross_proof.py` → `data/cross_proof_report.json`. **42** obligations (13 engine, 29 measured). **overall_ok = true.** Does not inherit the hub report.

| Layer | Status |
|-------|--------|
| Python D1D38A + formula path | PASS |
| SMT python + Z3 | PASS |
| Lean 4 + Mathlib (`Catalog.lean`) | PASS |
| Coq / Rocq | PASS |
| Isabelle/HOL | PASS |
| F* | PASS |
| Rust f64 kernel | PASS |
| TLA+ TLC routing | PASS |

Labeled archive: `docs/VERIFIED_SOLVES.md`. Biohub/Kaggle stay frozen (`docs/BIOHUB_FREEZE.md`).
