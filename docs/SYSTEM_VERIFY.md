# System-wide claim verify

Pin `D1D38A`. Law \(S=K(T_1+T_2+T_3)\). **0 free parameters.**

```powershell
python scripts/system_verify.py          # 52 published claims vs live JSON + engine
python scripts/verify_cross.py           # pin, seeds, formula path, Lean chem-link
python verification/run_cross_proof.py   # Lean / Coq / Isabelle / F* / SMT / Rust / TLA+
```

Last full run: `data/system_verify.json` **overall_ok=true**, n=52, fail=0.  
Gauntlet: `data/cross_proof_report.json` **overall_ok=true**, 42 obligations, 10 layers.

Does **not** re-fold the 10-protein product freeze (that JSON is the freeze).  
Does **not** invent bee / mosquito / plant connectomes.

## What this repository solves

| Solve | Accuracy kind | Live number | Source |
|-------|---------------|-------------|--------|
| Product Cα (same-data homolog) | Å vs eval PDB | median **0.13** vs AF **0.47** (n=10, 10/10 beat AF, 10/10 sub-2 Å) | `data/product_vs_alphafold.json` |
| Orphan / no map | not Å product | `no_measured_map` (Rg + secondary). Bulk ~13.6 Å is retired MDS | freeze |
| Residual hop on measured synapses | hop split | Male VNC **10.74** / JO **7.77** / olfactory **0.001**; BANC **10.16** / **6.64** / **0.0008** | `male_cns_boot.json`, `banc_connectome_boot.json` |
| Independent fly EM (hemibrain) | descending, not vnc_motor | JO **2.68** Giant Fiber vs olfactory **0.072** AL LN | `hemibrain_connectome_boot.json` |
| Larva | DN-VNC | mechano **16.6** vs olfactory **0.43** | `larva_connectome_boot.json` |
| Worm / Ciona / Platynereis | graph nonempty; NT unsigned where unannotated | 453 / 205 / 1720 cells | `*_connectome_boot.json` |
| Insect homologs | count | **26** folded, **1** true miss (mec-4 Anopheles) | `homolog_correspondence.json` |
| Analog pointer | mapped job | mec-4 → nompC; unc-25 → Gad1 | `analog_pointer.json` |
| Arabidopsis metabolic panel | product Cα | **6/6** | `plant_product.json` |
| Crop homologs | product Cα | **24/24**, miss **0** | `plant_homolog.json` |
| Plant information graph | hop split, not a connectome | 8,329 nodes, 39,664 edges; light → PIF3; SLAC1 unlit | `plant_signal_boot.json` |
| Plant signaling product | homolog Cα | audited accessions; PIF3/HY5/ARF5 `no_measured_map` | `plant_signal_product.json` |
| Experimental PGx | concordance | **10/10** (research disclosure) | `experimental_pgx.json` |
| Engine identities | milliscale Nat | φ milli 1618, leftover 1/φ², close-homolog 1/φ, 7 ChemLink D_eff | gauntlet layer A |
| Pair chemistry | scalar | w(F,W,d=8)=**1.670052** | trinary |

## What it does not solve

- Orphan 3-D at AlphaFold grade
- Bee / mosquito / plant / crow connectomes
- GABA or leftover LN/KC mass as a claimed thought
- CASP/CAMEO blind (still open: `docs/OPEN.md`)
- Biohub / Kaggle (frozen: `docs/BIOHUB_FREEZE.md`)
