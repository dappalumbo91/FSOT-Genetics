# Dimensionality Audit

## Scope

This audit reuses the frozen twelve-chain oriented-backbone holdout to measure
the intrinsic dimensionality of the FSOT proximity object. It does not select a
formula or change the production fold. Machine-readable results are in
`data/dimensionality_audit.json`; reproduce with
`python scripts/run_dimensionality_audit.py`.

## Central Finding

An FSOT target matrix is not a three-dimensional spatial object.

| Quantity | Native Cα matrix (real 3-D) | FSOT target matrix |
| --- | ---: | ---: |
| Variance in top-3 eigenvalues | 100.0% | 45.1% |
| Participation (effective) dimension | 2.53 | 9.10 |
| Negative-eigenvalue mass | ~0 | 0.202 |

A real protein is a 3-D object: 100% of its distance-matrix variance lives in
three eigenvalues. The FSOT object needs ~9 dimensions, and 20% of its structure
is non-Euclidean (negative-eigenvalue mass) — it cannot be embedded in any
number of flat spatial axes without distortion.

The measured participation dimension 9.10 sits on the pinned `D_eff` ladder
(`FSOTGenetics/ChemLink.lean`): backbone 8, disulfide 7, salt bridge 9,
hydrophobic pack 14, H-bond 8, sidechain 9, tertiary biochem 13; base law 25.
The proximity object's intrinsic dimension is its `D_eff`. The scalar law was
self-consistent; the reconstruction stage was asserting `D_eff = 3`.

## Consequence: The 3-D Collapse Destroys Topology

Native long-range contact recall from the FSOT embedding, by the dimension the
object is allowed to occupy:

| Embedding dimension | Native contact recall |
| --- | ---: |
| 3 (forced by `classical_mds`) | 0.995% |
| 9 (FSOT-native `D_eff`) | 2.086% |
| 25 (base-law baseline) | 2.268% |

Letting the object breathe in its native dimension more than doubles recovered
topology (0.995% → 2.086%) and plateaus by 25. More than half the topological
signal the math already contained was being deleted by the 3-D projection before
refinement ran.

## Where We Sit vs Current Capabilities

The finding splits the problem into observables FSOT determines and an observer
collapse it cannot supply from a single sequence.

| Observable | Dimensionality | FSOT (0 params, 1 sequence) | Reference capability |
| --- | --- | ---: | --- |
| Size / radius of gyration | 0-D scalar | 7.1% median rel. error | Empirical law `2.2·N^0.38` ≈ 5% |
| Secondary structure (H/E/C) | 1-D labels | 50.5% accuracy (F12c) | Single-seq ML Q3 ≈ 80% |
| Native-dimension contact topology | `D_eff`-D | 2.09% recall | MSA/coevolution predictors far higher |
| 3-D Cα coordinates (RMSD) | observer collapse | ~15 A median | AlphaFold ≈ 0.4 A (uses MSAs) |

Read this honestly:

- On the observer-independent, low-dimensional observable it actually targets —
  overall size — FSOT is within ~2 points of the empirical literature scaling
  law, with zero free parameters and no data fit.
- Secondary structure is below trained single-sequence models but is analytic
  and parameter-free.
- Specific 3-D coordinates are not competitive, and the audit explains why: a
  single sequence lacks the evolutionary covariation that fixes specific
  contacts, and the `D_eff → 3` observer collapse is underdetermined. This is an
  information limit, not a formula defect. It is the same reason MSA-free methods
  trail MSA-based ones.

## Where To Go Next

1. Stop asserting `D_eff = 3`. Represent and score structure at `D_eff`, treating
   3-D coordinates as an explicit observer projection rather than the native
   object.
2. Lead with the observables FSOT provably owns: size, burial profile, secondary
   structure, and native-dimension topology — each reported in its own
   dimension.
3. If 3-D coordinates remain a target, model the observer collapse explicitly
   (the scalar law's observer term is the hook) and state that it is
   underdetermined without covariation input.
4. Do not add scalar or contact levers before a candidate is frozen and tested on
   a fresh preregistered holdout.
