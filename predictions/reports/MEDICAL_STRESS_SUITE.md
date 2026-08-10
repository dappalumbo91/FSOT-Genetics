# Medical stress suite — FSOT multi-regime

Generated: `2026-08-10T21:49:23.756993+00:00`  
Free parameters: **0**  
Targets: **10**

## Median Cα RMSD (Å) to experimental PDB

| Regime | Median Å | n |
|--------|---------:|--:|
| AlphaFold DB | 0.4706172512401575 | 10 |
| FSOT template + MSA fuse | 1.1560920992862558 | 9 |
| FSOT template + physics | 1.1856870735796325 | 9 |
| FSOT template | 1.2218645296216386 | 9 |
| FSOT bulk + MSA | 16.806315999127698 | 9 |
| FSOT bulk single | 17.368488270101338 | 10 |

## Per-target

| Protein | N | AF | tmpl | phys | fuse | bulk | bulk+MSA | MSA depth | topLΔ |
|---------|--:|---:|-----:|-----:|-----:|-----:|---------:|----------:|------:|
| Hemoglobin alpha | 141 | 0.27 | 1.22 | 1.19 | 1.16 | 14.07 | 14.07 | 12000 | +0.00 |
| Hemoglobin beta | 145 | 0.52 | 1.32 | 1.32 | 1.32 | 14.92 | 14.92 | 12000 | +0.00 |
| Carbonic anhydrase II | 256 | 0.36 | 1.37 | 1.31 | 1.31 | 18.13 | 18.13 | 12000 | +0.00 |
| SOD1 | 153 | 0.29 | 1.11 | 1.18 | 1.13 | 16.92 | 16.81 | 12000 | +0.00 |
| Lysozyme human | 130 | 0.42 | 1.38 | 1.34 | 1.33 | 22.69 | 22.69 | 3573 | +0.00 |
| RNase A | 124 | 0.33 | 0.44 | 0.44 | 0.44 | 19.13 | 19.14 | 2733 | +0.00 |
| Ubiquitin | 76 | 0.88 | 2.12 | 2.12 | 2.12 | 11.10 | 11.10 | 12000 | +0.00 |
| Insulin | 21 | 4.51 | 1.13 | 1.13 | 1.13 | 5.69 | 5.22 | 6194 | +0.00 |
| p53 DNA-binding | 196 | 6.19 | — | — | — | 18.74 | 18.72 | 4272 | +0.00 |
| Calmodulin | 144 | 6.45 | 0.77 | 0.76 | 0.76 | 17.82 | — | — |  |

## Honesty notes

- Template regime is the medical-grade structure path when homologs exist.
- Bulk de-novo (~11 Å ceiling) is the honest orphan-sequence fallback.
- MSA inject improves *contact ranking / confidence*; topology still needs templates.
- All FSOT columns: zero trained weights; MSA/templates are data inputs.
