# Oriented Backbone Development Result

## Scope

F19 resolves the reflection ambiguity proven in the information-gap audit. It
uses only the physical handedness of alpha-helices formed from L-amino acids.
It does not use native coordinates, fitted constants, or learned weights.

## Frozen Development Set

| Model | Median C-alpha RMSD |
| --- | ---: |
| Baseline production fold | 8.590 A |
| Cooperative F12c | 10.172 A |
| Cooperative F12c plus F19 | 8.432 A |

F19 reflected three of five F12c predictions. Engrailed improved from
10.751 A to 8.365 A. Protein G worsened slightly from 10.172 A to 10.261 A.
The maximum pair-distance change from F19 was exactly 0.0 A, and baseline
output with the option disabled had exactly zero coordinate delta.

This is a development result, not evidence of 1-2 A prediction accuracy. The
candidate must be committed before evaluation on a fresh nonredundant PDB
holdout. Failed mutual-partner and coordinate-nearest F13 sparsification probes
were not added to the model.

Machine-readable measurements are in `data/information_gap_audit.json`.
