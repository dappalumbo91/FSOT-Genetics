# Structure Behavior Audit

## Scope

This audit reuses the frozen twelve-chain oriented-backbone holdout to inspect
observed behavior. It does not select a formula, tune a value, or change the
production fold. The audit separates ranked-contact constraints, target-distance
error, coordinate-realization error, and structural context.

## Ranked Constraints

| Measurement | Result |
| --- | ---: |
| Tight Top-L/2 native-contact precision | 2.749% |
| Loose next-L/2 native-contact precision | 2.604% |
| Strict separation-24 Top-L recall | 2.060% |
| False capped-pair signed distance error, median | -17.450 A |

The universal Top-L cap is not identifying native contacts on this set. It
forces almost all selected pairs much closer than their observed separation.
These false constraints can still serve as accidental global scale anchors,
which explains why simply deleting them destabilizes the current globally
rescaled MDS path.

## Distance Behavior

| Sequence separation | Target signed median | Realized signed median |
| --- | ---: | ---: |
| 1-2 | -0.088 A | -0.028 A |
| 3-6 | -4.494 A | -4.883 A |
| 7-11 | -5.565 A | -8.380 A |
| 12-23 | -9.422 A | -8.237 A |
| 24+ | -6.979 A | +6.224 A |

The target matrix compresses native distances outside the immediate covalent
neighborhood. Coordinate reconstruction does not merely inherit that error: it
over-compresses separation 7-11 while reversing the long-range bias. The target
matrix and its three-dimensional realization are therefore separate bottlenecks.

The production path obtains a mean adjacent spacing from raw classical MDS and
then scales every coordinate to make that local statistic 3.8 A. When raw MDS
collapses adjacent spacing, the operation magnifies the entire embedding. A
diagnostic path that removes ranked caps, skips global rescaling, locally
rebonds the chain, and canonicalizes chirality reduced median RMSD from
15.066 A to 13.230 A. This is a diagnosis of reconstruction behavior, not a
holdout-valid production candidate.

## Structural Context

Six of twelve references contain another protein chain, ten contain a non-water
hetero component, and all twelve contain at least one of those contexts. The
set includes metal ions, cofactors, bound metabolites, inhibitors, and
multichain interfaces. The current engine receives only the target-chain
sequence, so it cannot represent those causal inputs.

This does not make the references invalid. It means sequence-only autonomous
folding and context-dependent crystallographic conformation must be reported as
different benchmark strata.

## Consequences

1. Do not add another scalar or ranking lever before isolating these failure modes.
2. Replace global bond-statistic scaling with local covalent constraints in a
   separately frozen development candidate.
3. Stratify future structures into autonomous monomer, oligomer/partner-bound,
   ligand/cofactor/metal-bound, and elongated or conditionally folded classes.
4. Evaluate removal or conditioning of universal Top-L caps only on a new
   preregistered holdout.

Machine-readable per-chain measurements are in
`data/structure_behavior_audit.json`.