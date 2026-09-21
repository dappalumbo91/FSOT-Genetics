# What FSOT-Genetics actually does (live results)

Pin `D1D38A`. Law \(S = K(T_1+T_2+T_3)\). **0 free parameters.**  
Data stays on `D:\FlyWire_Connectome` (not git). Do not mix product Å with hop mass.

Full claim inventory + re-run: `docs/SYSTEM_VERIFY.md`. Last stamp: `python scripts/system_verify.py` **247/247**; gauntlet **overall_ok** (42 obligations). The 52-check stamp was freeze+hops only.

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

Independent EM volume, no VNC — neuPrint hemibrain:v1.2.1 typed neurons including cropped Leaves (compact GCS dump drops JO). Same residual, **unsigned** (no predictedNt). Score **descending**, not vnc_motor. DNg29 is not in this cut.

`python scripts/hemibrain_connectome.py` → `data/hemibrain_connectome_boot.json`

| Seed | n | hop-2 `descending` | hop-1 top |
|------|--:|-------------------:|-----------|
| JO (`JO-A/B/C`) | 78 | **2.68** | **Giant Fiber** |
| olfactory (`ORN_*`) | 2577 | **0.072** | antennal-lobe LN `lLN2T_c` |

Same split: JO lights a descending escape neuron; olfactory stays in the AL LN leftover. Do not mix 2.68 descending with 7.77 vnc_motor — different effector pools.

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

### Plant information graph (not a connectome)

Plants do not use synapses. The measured equivalent is experimental protein–protein interactions.

`python scripts/plant_signal.py` — IntAct PSICQUIC *Arabidopsis* (taxid 3702), physical edges only. **8,329** proteins, **39,664** undirected edges. Same residual \(r = 1+|S|\cdot P_{\mathrm{NEW}}\) as the animal graphs.

Hop-1 / hop-2 (which class lights):

| Seed | hop-1 top | hop-1 `light_tf` | hop-2 top |
|------|-----------|-----------------:|-----------|
| photoreceptor | **PIF3** | **1.00** | PHYB |
| ABA | **ABI1** | 0 | PYL9 |
| auxin | ABCB19 | 0 | **TIR1** |
| calvin / PSII | CML9 | 0 | **psbA** |

Light seed lights the light transcription job (PIF3). Calvin seed stays on PSII. Auxin seed goes to TIR1/IAA. Stomatal channels (SLAC1) are almost unlit on this PPI graph (hop-2 ~0.002 from ABA) — kinase→channel is often phosphorylation, not a binary IntAct edge. We do not invent that edge.

Live UniProt gene-name audit of the class table (taxid 3702 reviewed). Dropped wrong accessions (same class as wheat *acT2* / PHOT1 Q2V2M9):

| Was labeled | Wrong acc | Live gene | Correct acc |
|-------------|-----------|-----------|-------------|
| PHOT1 | Q2V2M9 | human FHOD3 | **O48963** |
| PHYA | P42497 | PHYD | **P14712** |
| PIF3 | Q495N3 | human SLC36A3 | **O80536** |
| PYR1 | Q8VZS9 | FZR1 | **O49686** |
| PYL4 | Q8S8E3 | PYL6 | **O80920** |
| ABI2 | P25042 | yeast ROX1 | **O04719** |
| SLAC1 | Q9FLV9 | SLAH3 | **Q9LD83** |
| KAT1 | Q39153 | MYB13 | **Q39128** |
| GORK | Q94KI8 | TPC1 | **Q94A76** |

Hops re-run after the fix: same split (PIF3 / PHYB / ABI1 / TIR1 / psbA). OST1 is UniProt gene *SRK2E* (synonym OST1) — kept.

Source: `data/plant_signal_boot.json`.

Proteins on that class table (`python scripts/plant_signal_product.py`). Close-homolog floor 1/φ. Not a connectome.

| Gene | UniProt | Template | id | cov | Mode |
|------|---------|----------|---:|----:|------|
| *PHOT1* | O48963 | 5HZI | 0.61 | 0.46 | leftover (id &lt; 1/φ). LOV1 **2Z6C** / LOV2 **4HHD** id 1.00 |
| *PHOT2* | P93025 | 5HZK | 0.60 | 0.49 | leftover |
| *PHYA* | P14712 | 8IFF | **1.00** | 0.77 | product |
| *PHYB* | P14713 | 7RZW | **1.00** | 0.74 | product |
| *CRY1* | Q43125 | 1U3C | **1.00** | 0.71 | product |
| *CRY2* | Q96524 | 6K8K | **0.99** | 0.80 | product |
| *UVR8* | Q9FN03 | 8GQE | **1.00** | 0.88 | product; **not on IntAct hops** |
| *PIF3* | O80536 | — | — | — | **no_measured_map** (old 9V3V was human SLC36A3) |
| *HY5* | O24646 | — | — | — | **no_measured_map** |
| *OST1* | Q940H6 | 3UC4 | **0.99** | 0.81 | product (gene SRK2E) |
| *PYR1* | O49686 | 3ZVU | **0.99** | 0.95 | product |
| *PYL4* | O80920 | 8AY6 | 0.55 | 0.89 | not close-homolog |
| *ABI1* | P49597 | 3NMN | **1.00** | 0.64 | product |
| *ABI2* | O04719 | 3UJK | **1.00** | 0.70 | product |
| *SLAC1* | Q9LD83 | 8J0J | **1.00** | 0.68 | product; hop-2 still unlit |
| *KAT1* | Q39128 | 6V1X | **1.00** | 0.66 | product |
| *GORK* | Q94A76 | 9J0X | **1.00** | 0.84 | product |
| *PIN1* | Q9C6B8 | 7Y9T | **1.00** | 0.61 | product |
| *TIR1* | Q570C0 | 2P1N | **1.00** | 0.96 | product |
| *ARF5* | P93024 | — | — | — | **no_measured_map** |

A crystal on SLAC1 is not an IntAct edge. UVR8 product is not a hop seed. Inter-domain PHOT1 pose is not claimed.

Source: `data/plant_signal_product.json`, `data/plant_signal_domains.json`.

## Where this stands against AlphaFold

Same law on every organism. Competitive when a measured homolog exists. Not competitive as an orphan fold.

| Physiology | What is scored | Against AlphaFold |
|------------|----------------|-------------------|
| Human protein / variant | Medical panel median **0.26 Å** vs AF **3.98 Å** (19/19 sub-2 Å). Drivers recalled. P72R demoted by measured allele frequency. | Wins on this panel when the crystal exists. AF still covers sequences with no homolog; that path here is `no_measured_map` (~13.6 Å bulk). |
| Animal | Fly, worm, Ciona, Platynereis, larva, hemibrain: residual hops on measured synapses. Named proteins (nompC, Gad1) get product Cα. | Not a second AlphaFold. The graph is the measurement. |
| Plant | Arabidopsis + rice, maize, soybean, wheat product Cα. IntAct physical PPI is the information graph (8,329 nodes). | Same protein law. Not an invented plant connectome. |

## What this is capable of

| Can | Cannot |
|-----|--------|
| Product Cα when a homolog exists (0.13 Å freeze) | Orphan 3-D at AlphaFold grade |
| AF3-class jobs at the named ChemLink (DNA, metal, SC, H, ligand, PPI) | Invent contacts or a from-sequence AF clone |
| Residual hops on a **measured** synapse graph, or on **measured** plant PPIs | Invent a bee / mosquito / plant connectome |
| Point a clade-restricted 1:1 at the mapped residual job | Treat a genome as synapses |
| Same law on animals and plants (proteins) | Skip VNC and call a fly brain a whole animal |
| Variant / PGx concordance on public labels (research) | A diagnosis, prescription, or FDA device |
| Zig/Lean/Haskell/Rust/Coq/Isabelle/F*/TLA+ on the same pin | Treating 1.14 Å fair-cap or 13.6 Å bulk as the product |

## 6. AF3-class coverage

`python scripts/bench_af_coverage.py` → `data/af_coverage.json`. **17/17** jobs. Same pin, 0 free parameters, measured homolog except the eval PDB. Details: `docs/AF_COVERAGE.md`.

| Job | Number |
|-----|-------:|
| Monomer (freeze) | **0.13 Å** vs AF 0.47 |
| p53–DNA protein / C1′ | **0.013** / **0.016** |
| CAII Zn site / SOD1 metal | **0.061** / **0.26** |
| tRNA C1′ / modified C1′ | **0.68** / **0.93** |
| Neutron H (961/962) | **1.01** |
| Hb dimer / iface / tetramer | **0.45** / **0.17** / **0.51** |
| SC centroids / heavy | **0.41** / **1.01** |
| Glycan / phospho | **0.56** / **0.77** |
| Antibody CDR / H+L pair | **0.93** / **1.03** |
| U1A prot / RNA seed | **0.23** / **0.28** |
| Trypsin–BEN ligand site | **0.60** |
| Joint `predict_system` | CA **0.013** · SC **0.016** · DNA **0.016** |

## 7. Medical / variant (disclosure required)

Not a device. `docs/EXPERIMENTAL_DISCLOSURE.md`, `docs/MEDICAL_PLATFORM.md`.

| Bench | Live |
|-------|------|
| Experimental PGx | **10/10** |
| Variant panel | 8 genes, 35 drivers, recall@75% **1.0** |
| Wet-lab structure | **19/19** sub-2 Å · median **0.26 Å** vs AF **3.98** · beats AF 16/16. Product cap 1.0 + apparatus min. Not the n=10 freeze |
| Wet-lab variants | pathogenic recall **1.0**; P72R **common_polymorphism** (pop AF 0.46 ≥ 1/φ³) |
| Reality margin | 19/19 · median **1.17 Å** vs target 2.5 |
| Domain-split | SOD1 **0.29** / HBB **0.30**; TP53 4 domains, inter-domain pose not freeze |
| Catalog | TP53, KRAS, EGFR, BRAF, CFTR, SOD1, HBB, BRCA1 |

Same law, measured homolog except the eval PDB, apparatus minimum over `trit_not` collapses. The search no longer spends its PDB budget on same-protein redeposits and then keeps the partner chain.

| Case | Was | Now |
|------|----:|----:|
| p53 DBD | 2.60 | **0.01** |
| EGFR kinase | no template | **1.00** |
| ABL1 | 12.7 | **0.81** |
| BRAF | 2.76 | **0.94** |
| BCL-2 | 5.80 | **1.90** |
| SARS-CoV-2 RBD | 5.68 | **0.52** |
| ACE2 peptidase | 5.20 | **0.48** |
| HIV-1 RT | 5.46 | **0.85** |

BCL-2 stops at **1.90 Å** because every close measured map is another NMR model. No crystal in that pool is under 1.5 Å. That is the ensemble width, not a missing homolog.

Fair-cap / fuse-era (`medical_stress_suite.json` fuse **1.16 Å**, bulk 16–17 Å) and the 0.95 handicap (`m1_authority_verify.json`, **1.14 Å**) stay on the honesty wall. They are not this panel.

## 8. Language stack and engine identities

| Layer | Stamp |
|-------|-------|
| Zig ↔ Python | **PASS** 22/22 (`data/parity_zig_python.json`) |
| ChemLink \(D_{\mathrm{eff}}\) | 8 / 8 / 9 / 9 / 7 / 14 / 13 |
| Dimensionality | FSOT participation **9.10** vs native **2.53**; base law **25**; neg-eigen mass **0.20** |
| SMILES | 116 records; \(P_{\mathrm{NEW}}\) matches pin |
| Leftover / close-homolog | \(1/\varphi^2\), \(1/\varphi\) |
| Product identity cap | **1.0** (0.95 is a handicap bench) |

Haskell ChemLink/Contact, Lean ChemLink + ZeroFreeParams, Rust `fsot_protein`, Zig host/kernel, Coq/Isabelle/F*/TLA+ spines: present and gauntlet-stamped. `docs/LANGUAGE_STACK.md`.

## 9. Walking proteins and organism join

Named cell on a measured graph → UniProt → product Cα. No MDS fly-protein brain.

`data/fly_walking_product.json`: **Gad1, nan, iav, nompC**. nompC template **5VKQ** identity **0.97**.

Walking observer (`fly_behavior_flow.json`) hop-2 descending: mechanosensory **23.84** / JO **16.55** / olfactory **0.245**. Odor rest is **no afferent seed**. Courtship / sleep / aggression seed Male CNS 165,122 with measured fru/dsx / ER-FB-LNv / pC1 — not invented fight neurons.

`organism_product_join.json`: 4 walking folds + 6 new folds. `male_cns_genetics_on_cells.json`: fru_high **2611**, dsx_high **138**, receptorType **752**. Fly Cell Atlas bodyId join does not exist — do not invent one.

Worm whole-animal: hermaphrodite **453** cells, male **575**; hop-2 sex-specific mass lights in the male (**6.85**) and stays dark in the hermaphrodite (**0.07**).

## 10. Field, F12, RCSB — labeled correctly

| Stamp | Number | Label |
|-------|--------|-------|
| Field stress | **49/49 PASS** | deploy gate |
| F12 SS candidate | macro **0.58** vs **0.34** | development; `production_enabled=False` |
| Kaggle F12 | 6,483 chains, gates passed | **frozen** (`docs/BIOHUB_FREEZE.md`) |
| Distogram Pearson | **0.61** (n=5) | contact ranking |
| RCSB oriented F19 | pair-distance delta **0**, gate passed | chirality |
| RCSB template holdout | 60 chains / 54 covered / cap **0.95** / best **2.20 Å** | not the freeze |
| RCSB live API bulk | ~11 Å | orphan path |
| MSA dual-mode | CA drift **0.0012 Å** | MSA is data |

## 11. Honesty walls (verified as not the product)

| Stamp | Number |
|-------|--------|
| ChemLink / UniRef bulk | **~13.57 Å** |
| Sequence-only vs AF | FSOT **13.94 Å**, **0/8** wins |
| M1 authority cap 0.95 | `ok=False`, **1.14 Å** |
| Residual-template | **6.77 Å** |
| Error-margin diagnostic | **8.59 Å** |
| Biohub GEFF inventory | 199 videos, **back burner** |

Do not cross-cite product freeze **0.13 Å** · medical panel **0.26 Å** · AF on that panel **3.98 Å** · fair-cap handicap **1.14 Å** · bulk **13.6 Å**. The older wet-lab **2.06 Å** was the 0.95-cap snapshot.

## 12. System verify + multi-prover stamp

`python scripts/system_verify.py` → `data/system_verify.json`. **247/247** live claims vs JSON + engine (was 52 when only freeze+hops were wired).

`python scripts/verify_cross.py` → pin D1D38A, 0 free parameters, Lean chem-link. **PASS.**

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

Labeled archive: `docs/VERIFIED_SOLVES.md`. Biohub/Kaggle stay frozen (`docs/BIOHUB_FREEZE.md`). Fly pack (`C:\Users\damia\Desktop\fsot fly nuron net`, pin AEB2AD) is a separate system — this stamp is Genetics only.
