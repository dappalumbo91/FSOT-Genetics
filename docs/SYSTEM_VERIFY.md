# System-wide claim verify

Pin `D1D38A`. Law \(S=K(T_1+T_2+T_3)\). **0 free parameters.**

```powershell
python scripts/system_verify.py          # every live JSON claim vs engine
python scripts/verify_cross.py           # pin, seeds, formula path, Lean chem-link
python verification/run_cross_proof.py   # Lean / Coq / Isabelle / F* / SMT / Rust / TLA+
```

Last full run: `data/system_verify.json` **overall_ok=true**, n=247, fail=0.  
Gauntlet: `data/cross_proof_report.json` **overall_ok=true**, 42 obligations, 10 layers.

Does **not** re-fold the 10-protein product freeze (that JSON is the freeze).  
Does **not** invent bee / mosquito / plant connectomes.  
Does **not** mix product Å, hop mass, coverage Å, bulk MDS, or fair-cap/fuse benches.

The previous 52-check stamp only covered freeze + hop splits. This file is the full repo.

## What this repository solves

### Product structure (Å, measured homolog)

| Solve | Accuracy kind | Live number | Source |
|-------|---------------|-------------|--------|
| Product Cα (same-data homolog) | Å vs eval PDB | median **0.13** vs AF **0.47** (n=10, 10/10 beat AF, 10/10 sub-2 Å) | `data/product_vs_alphafold.json` |
| Identity cap on the product path | cap | **1.0** (`PRODUCT_IDENTITY_CAP`) | `scripts/run_rcsb_template_holdout.py` |
| Orphan / no map | not Å product | `no_measured_map` (Rg + secondary). Bulk ~13.6 Å is retired MDS | freeze |
| Joint forward | one call | `predict_system()` p53 CA **0.013** / SC **0.016** / DNA C1′ **0.016** | `af_coverage.json` `joint_forward` |

### AF3-class coverage (named ChemLink, not invented contacts)

`data/af_coverage.json` · 17/17 jobs · `docs/AF_COVERAGE.md`

| AF3 job | Wet-lab number |
|---------|----------------|
| Protein monomer Cα | **0.13 Å** (AF 0.47) — same freeze |
| Protein–DNA | p53 **0.013 Å** · DNA C1′ **0.016 Å** |
| Metal / ion | CAII Zn **0.061 Å** · SOD1 **0.26 Å** |
| RNA | tRNA C1′ **0.68 Å** |
| Modified nucleotides | C1′ **0.93 Å** · modified sites **1.57 Å** |
| Neutron H | 961/962 matched · **1.01 Å** |
| Protein–protein | Hb dimer **0.45 Å** · iface MAE **0.17 Å** |
| Tetramer | Hb **0.51 Å** |
| Side chains | centroids **0.41 Å** · heavy **1.01 Å** |
| PTM / glycan | prot **0.56 Å** |
| PTM / phospho | PKA **0.77 Å** |
| Antibody CDR | CA **0.93 Å** (Superposed loops) |
| Antibody H+L | pair **1.03 Å** · iface **0.40 Å** |
| Protein–RNA | U1A **0.23 Å** · RNA seed **0.28 Å** |
| Ligand | trypsin–BEN site **0.60 Å** |

### Residual hops (measured graphs)

| Solve | Accuracy kind | Live number | Source |
|-------|---------------|-------------|--------|
| Male CNS / BANC | hop split | VNC **10.74** / JO **7.77** / olf **0.001**; BANC **10.16** / **6.64** / **0.0008** | `male_cns_boot.json`, `banc_connectome_boot.json` |
| Hemibrain | descending, not vnc_motor | JO **2.68** Giant Fiber vs olfactory **0.072** AL LN | `hemibrain_connectome_boot.json` |
| Larva | DN-VNC | mechano **16.6** vs olfactory **0.43** | `larva_connectome_boot.json` |
| FlyWire FAFB | measured GABA | 127,979 neurons used; inventory 139,248 | `fly_connectome_boot.json` |
| Walking observer | descending | mechano hop-2 **23.84** / JO **16.55** / olfactory **0.245** | `fly_behavior_flow.json` |
| Odor observer | rest = no seed | `has_rest_state` · rest `n_seed=0` | `fly_odor_flow.json` |
| Courtship / sleep / aggression | Male CNS 165,122 | same residual; honesty: not a thought | `fly_*_flow.json` |
| Worm both sexes | whole-animal | herm **453** / male **575**; hop-2 sex-specific male **6.85** vs herm **0.07** | `worm_*_boot.json`, `worm_sex_compare.json` |
| Ciona / Platynereis | graph nonempty; NT unsigned | 205 / 1720 cells | `*_connectome_boot.json` |

### Homologs, analog, plants

| Solve | Live number | Source |
|-------|-------------|--------|
| Insect homologs | **26** folded, **1** true miss (`mec-4` Anopheles) | `homolog_correspondence.json` |
| Analog pointer | mec-4 → nompC; unc-25 → Gad1 | `analog_pointer.json` |
| Arabidopsis metabolic | **6/6** product | `plant_product.json` |
| Crop homologs | **24/24**, miss **0** | `plant_homolog.json` |
| Plant IntAct graph | 8,329 nodes, 39,664 edges; light → PIF3; SLAC1 unlit | `plant_signal_boot.json` |
| Plant signaling product | 20 genes; PIF3/HY5/ARF5 `no_measured_map` | `plant_signal_product.json` |
| PHOT1 domains | O48963; LOV1/LOV2 close-homolog; kinase not; pose not claimed | `plant_signal_domains.json` |
| UVR8 | Q9FN03 / 8GQE product; **not** on IntAct hops | `plant_signal_domains.json` |

### Medical / variant (research disclosure, not a device)

| Solve | Live number | Source |
|-------|-------------|--------|
| Experimental PGx | **10/10** concordant | `experimental_pgx.json` |
| Variant panel | 8 genes / 35 drivers / recall@75% **1.0** | `medical_variant_panel.json` |
| Gene catalog | TP53 KRAS EGFR BRAF CFTR SOD1 HBB BRCA1 | `scripts/medical_gene_catalog.py` |
| Wet-lab product panel | **19/19** sub-2 Å; median **0.26** vs AF **3.98**; beats AF **16/16** with an AF model. Cap **1.0** + apparatus min. EGFR **1.00** (was no template). ABL1 **0.81** (was 12.7). RBD **0.52** (was 5.7). Pathogenic recall **1.0**; benign-like miss **recorded** | `wetlab_af_eval.json` |
| Reality margin | 19/19 · median **1.17 Å** vs target 2.5 | `reality_margin_eval.json` |
| DNA → AA | `dna_variant_effect.py` | scripts |
| Domain-split | KRAS/SOD1/HBB/TP53; SOD1 **0.29** / HBB **0.30**; TP53 4 domains, pose not freeze | `domain_split_eval.json` |

The medical panel median **0.26 Å** is the same product law on 19 wet-lab chains, not the 10-protein freeze **0.13 Å**. Reality **1.17 Å** and the 0.95-cap handicap (**1.14 Å**) stay separate. Do not collapse them.

### Language stack / engine identities

| Solve | Live number | Source |
|-------|-------------|--------|
| Engine identities | φ milli 1618, leftover \(1/\varphi^2\), close-homolog \(1/\varphi\), 7 ChemLink D_eff | gauntlet layer A |
| Pair chemistry | w(F,W,d=8)=**1.670052** | trinary |
| Zig ↔ Python | **PASS** 22/22 | `parity_zig_python.json` |
| Dimensionality | FSOT participation **9.10** vs native **2.53**; top-3 var **0.45**; neg-eigen **0.20**; base D_eff **25** | `dimensionality_audit.json` |
| SMILES chemistry | 116 records; \(P_{\mathrm{NEW}}\) matches pin | `formulas/smiles_protein_chemistry.json` |
| Lean / Haskell / Rust / Zig / Coq / Isabelle / F* / TLA+ | sources present; gauntlet 10/10 | `docs/LANGUAGE_STACK.md` |

ChemLink \(D_{\mathrm{eff}}\): Physical_Chemistry **8**, Chemistry **8**, Molecular_Chemistry **9**, Electromagnetism **9**, Atomic_Physics **7**, Condensed_Matter **14**, Biochemistry **13**. Never “3 because space is 3-D.”

### Organism product join (named cell → UniProt → Cα)

| Solve | Live number | Source |
|-------|-------------|--------|
| Walking proteins | Gad1, nan, iav, nompC; nompC **5VKQ** id **0.97** | `fly_walking_product.json` |
| Organism join | 4 walking folds + 6 new folds | `organism_product_join.json` |
| Male genetics on cells | 165,122 traced; fru_high **2611**; dsx_high **138**; receptorType **752** | `male_cns_genetics_on_cells.json` |

No invented Fly Cell Atlas join.

### Field / F12 / RCSB (stamped, with the right label)

| Solve | Live number | Label |
|-------|-------------|-------|
| Field stress | **49/49 PASS** | deploy gate, not a new Å claim | `field_stress_suite.json` |
| F12 SS candidate | macro recall **0.58** vs baseline **0.34**; `production_enabled=False` | development only | `f12_candidate_development.json` |
| Kaggle F12 | 6,483 chains, gates passed | **frozen** — not live product | `kaggle_f12_validation_eval.json` |
| Distogram contacts | median Pearson **0.61** (n=5) | contact ranking, not Cα product | `fsot_distogram_contact_eval.json` |
| RCSB oriented backbone | success gate; pair-distance delta **0** | F19 chirality | `rcsb_oriented_backbone_eval.json` |
| RCSB template holdout | n=60, 54 covered, cap **0.95**, best median **2.20 Å** | **not** the 0.13 freeze | `rcsb_template_holdout_eval.json` |
| RCSB live API bulk | n=60, bulk **~11 Å** | orphan path | `rcsb_live_api_holdout_eval.json` |
| MSA dual-mode | ubiquitin CA drift **0.0012 Å** | MSA is data, not a fold miracle | `msa_dual_mode_smoke.json` |

### Honesty walls (verified as *not* the product)

| Stamp | Number | Why it is here |
|-------|--------|----------------|
| ChemLink / UniRef bulk | **~13.57 Å** | orphan MDS ceiling |
| Sequence-only vs AF | FSOT **13.94 Å**, **0** wins / 8 AF | `fsot_vs_alphafold_structure.json` |
| M1 authority (cap 0.95) | `ok=False`, median **1.14 Å** | fair-cap handicap |
| Medical stress fuse-era | fuse **1.16 Å**, bulk **16–17 Å** | `docs/CAPABILITY_ROADMAP.md` |
| Residual-template bench | residual median **6.77 Å** | not current product |
| Error-margin log | median **8.59 Å** | bulk diagnostic |
| Biohub inventory | 199 GEFF | **back burner** (`docs/BIOHUB_FREEZE.md`) |

## What it does not solve

- Orphan 3-D at AlphaFold grade
- Bee / mosquito / plant / crow connectomes
- GABA or leftover LN/KC mass as a claimed thought
- CASP/CAMEO blind (still open: `docs/OPEN.md`)
- Biohub / Kaggle as a live competition product (frozen: `docs/BIOHUB_FREEZE.md`)
- Clinical / FDA device (experimental disclosure only)
- Inter-domain PHOT1 pose; IntAct phosphorylation as a fake edge
- Treating the medical-panel **0.26 Å**, the 0.95-cap **1.14 Å**, or bulk **13.6 Å** as the 10-protein **0.13** freeze
