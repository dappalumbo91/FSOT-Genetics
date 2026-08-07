# F15 Ranking and F12 Failure Localization

## Measurement Scope

The 1.06-1.95 A values in this experiment are experimental X-ray resolutions
of the RCSB reference structures. They are not FSOT coordinate RMSDs.

## Locality Plus Sheet Validation

- Candidate freeze: `957017b`.
- Third-holdout protocol freeze: `1da885c`.
- Holdout: 12 new, nonredundant RCSB chains, 48-219 residues.

| Score | Conditioned Pearson, median | 95% bootstrap interval | LR Top-L/2 precision |
| --- | ---: | ---: | ---: |
| Full F15 | 0.0376 | [0.0274, 0.0465] | 0.0000 |
| Locality | approximately zero | approximately zero | **0.0281** |
| Sheet only | **0.0538** | **[0.0352, 0.0908]** | 0.0149 |
| Locality + sheet | **0.0538** | **[0.0352, 0.0908]** | 0.0149 |

The candidate preserved the validated sheet residual and beat full F15, but it
did not beat locality for contact ranking. Adding locality did not materially
reorder sheet candidates on this holdout.

## Root Failure

F11 sheet propensity carries reproducible continuous signal. F17 register is
implemented only after F12 has produced distinct beta regions. Direct comparison
with PDB `HELIX` and `SHEET` annotations shows:

- Only `I`, `L`, `M`, `F`, and `W` can cross F12's beta gate.
- Beta-strand recall is 0% on all three beta-containing development proteins.
- Ubiquitin has five annotated strands; F12 predicts no beta region.
- Protein G has four annotated strands; F12 predicts no beta region.

Therefore F17 is usually starved of beta regions and cannot supply strand
orientation or register. The next target is F12 secondary-structure inference,
not another scalar weight or contact-score combination.

## Next Engineering Gate

Any F12 replacement must be developed only on declared development structures
and must satisfy all of these before another structure holdout:

1. Nonzero beta recall on every beta-containing development protein.
2. Improved macro H/E/C recall over current F12.
3. No fitted amino-acid lookup table or learned weight.
4. Preserve the validated F11 conditioned signal.
5. Freeze the formula and code before evaluating a fresh RCSB holdout.
