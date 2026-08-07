# RCSB Frozen Holdout Evaluation

## Protocol

- Source: 12 experimental X-ray structures from RCSB PDB, resolution at most 2.0 A.
- Length range: 87-219 residues.
- Release range: 2024-03-27 through 2026-08-05.
- Selection: one usable chain from each preregistered stratum of the date-sorted RCSB query.
- Frozen protocol commit: `c5affb18bba37ddbc453817521d7ad7dfbc62fcd`.
- Policy: these outcomes must not select formulas, routings, or thresholds.

Every downloaded PDB, the manifest, benchmark, structure engine, and D1D38A
authority implementation are SHA-256 identified in `data/rcsb_holdout_eval.json`.

## Result

| Metric | FSOT | Separation-only baseline |
| --- | ---: | ---: |
| All-pair proximity Pearson, median | 0.0447 | 0.6696 |
| Top-L/2 contact precision, separation >= 6 | 0.0172 | 0.0642 |
| Top-L/2 contact precision, separation >= 12 | 0.0172 | 0.0285 |
| Top-L/2 contact precision, separation >= 24 | 0.0048 | 0.0340 |

After subtracting the mean score and experimental inverse distance separately
within every exact sequence separation, FSOT has median Pearson correlation
`0.0381` with bootstrap 95% interval `[0.0248, 0.0445]`. The correlation is
positive on all 12 structures.

## Interpretation

The exact-separation result reproduces a small signal that cannot be attributed
to sequence separation alone. The current F15 score does not convert that signal
into useful tertiary-contact ordering: the locality baseline wins at every
reported contact threshold.

This localizes the next engineering problem to the mapping from the FSOT
residual channels into contact rank, not to the existence of any measurable
residual signal. Develop that mapping only on a separately declared development
set, freeze it in Git, and evaluate it once on a newly preregistered holdout.
