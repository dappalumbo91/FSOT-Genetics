# Open failures (marked 2026-08-13)

Results in `docs/PRODUCT_FREEZE.md` and `data/product_vs_alphafold.json` are **frozen**.  
This list is what to clear next. One mechanism at a time. Pin `D1D38A`. 0 free parameters.

## Product H2H

| Item | Now | Target / handle |
|------|-----|-----------------|
| Calmodulin 1CLL | **0.90 Å** (4EHQ) | Compact holo **3CLN ~0.52 Å** is not in the UniProt/search page. Need a lawful same-sequence retrieval that surfaces 1988-era crystals without paginating the 0.25 search (that pulled p53 4MZI). |
| Side-chain heavy atoms | **1.26 Å** vs centroids 0.93 | χ1 flips in the residue frame. Not a Cα problem (CA 0.12). |
| Hydrogens | **1.91 Å** (934/962 on 1LZN) | First coverage number; depth is local H geometry, not missing data. |

## Coverage / AF3 depth

| Item | Note |
|------|------|
| Joint `predict_system` | 0.39 Å vs DNA-job 0.11 — still not taking the full apparatus min in one forward without native. |
| Protein–RNA | Seed C1′ 0.28 Å on 9 nt; full hairpin register still Superposed. |
| 3CLN / 4CLN | Same CaM sequence class; search ranking, not residual. |

## Historical wet-lab (not the 10-protein freeze)

See `docs/MECHANISM_GAP_MAP.md`. Still open there: EGFR coverage, variant P72R specificity, vaccine-antigen domain scope, CASP/CAMEO blind.

Medical expansion (not Ångström work): `docs/MEDICAL_PLATFORM.md`. First gate there is **P72R must not be LIKELY DAMAGING**.

## Anti-goals (do not “fix” these)

- Bulk 11–14 Å (information ceiling).
- Residual picking DFG-in vs DFG-out.
- Bond-idealizing intact crystals.
- Geometric shotgun (medoid-all, invented contacts).
