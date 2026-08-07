# RCSB Sheet-Channel Validation

## Preregistration

- Candidate: F15 `sheet` channel alone.
- Candidate freeze: `dd899a0`.
- Validation protocol freeze: `0516d6d`.
- Data: 12 previously unused RCSB X-ray structures, 1.15-1.97 A resolution.
- Policy: evaluate once and do not revise the candidate from these outcomes.

## Result

| Score | Conditioned Pearson, median | 95% bootstrap interval | LR Top-L/2 precision |
| --- | ---: | ---: | ---: |
| Full F15 | 0.0517 | [0.0385, 0.0608] | 0.0197 |
| Sheet only | **0.0705** | **[0.0386, 0.0828]** | **0.0303** |
| Locality | 0.0000 | approximately zero by construction | 0.0509 |

The sheet channel had positive conditioned correlation on all 12 structures.
It improved both preregistered FSOT metrics over full F15, validating that this
channel carries reproducible information beyond exact sequence separation.

The candidate did not beat locality on long-range contact precision. It is a
validated residual feature, not yet a complete tertiary-contact ranker.

## Next Iteration

The validation failure now becomes development information. A parameter-free
`locality + sheet` score improves median LR precision from 0.0789 to 0.1053 on
the original five-structure development set while preserving sheet's
conditioned correlation. That new candidate requires a fresh holdout and must
not be evaluated on this validation set for selection evidence.
