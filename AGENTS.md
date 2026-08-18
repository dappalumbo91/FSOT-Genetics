# FSOT-Genetics — Agent Adherence Charter

This repository derives protein and genetics predictions from the **FSOT scalar
law** with **zero free parameters**. Any agent — human or AI — working here must
follow these rules. They are not style preferences; violating them produces
wrong results that look plausible. When something breaks, the cause is **almost
always `D_eff`** (see Rule 3).

## 1. The scalar engine is frozen

- The law is `S = K(T1 + T2 + T3)`, computed **only** through
  `vendor/fsot_compute.py` at authority pin `D1D38A`
  (SHA-256 `d1d38a185487b452e470ac68ece2eb45aeb1ca9ce25fc9bf9564c19633ffbe70`).
- Never reimplement, approximate, inline, or "simplify" the scalar engine.
- Import and call it. Do not reproduce its arithmetic by hand.

## 2. Zero free parameters — never fit

- Every constant must trace to the seeds (`π`, `e`, `φ`) or the pinned domain
  table. No fitted, tuned, or regressed values, ever.
- Scoring against experimental datasets (PDB, AlphaFold DB) is allowed. Using
  dataset labels to *set* any value is forbidden — that would make it ML and
  break the zero-parameter claim.
- New formulas (`Fxx`) must be derived and recorded in
  `formulas/FSOT_PROTEIN_DERIVATIONS.md` with their seed/domain provenance
  **before** they are used in code.

## 3. `D_eff` is NOT a spatial dimension (the rule that is always broken)

`D_eff` is an **effective / fractal dimension referenced to a 25-D baseline**. It
is never a count of Cartesian axes. Concretely:

- **Never set `D_eff = 3`.** Never treat `D_eff` as "the number of x/y/z axes."
- Use only the pinned ladder (`FSOTGenetics/ChemLink.lean`,
  `scripts/full_scalar_law.py`): backbone 8, disulfide 7, salt bridge 9,
  hydrophobic pack 14, H-bond 8, sidechain 9, tertiary biochem 13; base law 25.
  Do not invent new `D_eff` values.
- The FSOT proximity object is intrinsically **~9-dimensional and non-Euclidean**
  (measured participation dimension 9.10, matching the ladder; negative-
  eigenvalue mass 0.20). See `predictions/reports/DIMENSIONALITY_AUDIT.md`.
- **Do not collapse the proximity object into a 3-D Euclidean MDS embedding and
  call it "the structure."** That deletes >50% of the topological signal. A 3-D
  coordinate set is an *observer projection*, not the native FSOT object.
  Reconstruct and evaluate at `D_eff`.
- **When a result degrades, check `D_eff` first.** It is the most common failure.

## 4. Separate the observables from the observer collapse

- FSOT determines **observer-independent, low-dimensional observables**: radius
  of gyration (0-D), hydrophobic burial profile (1-D radial), secondary structure
  (1-D labels), and native-dimension contact topology (`D_eff`-D).
- A 3-D Cα RMSD requires an **observer collapse** `D_eff → 3` that is
  underdetermined from a single sequence. Do not claim *orphan / bulk*
  AlphaFold-grade coordinate accuracy (~11–14 Å ceiling). The **product**
  (measured homologs, exclude eval PDB) is a different information regime
  and is currently **0.13 Å** median — say which regime you mean.
  Report each observable in its own dimension.

## 5. Validation is falsifiable and frozen

- Freeze the candidate commit and the protocol commit **before** running any
  holdout. Record outcomes unchanged — wins and losses alike.
- Experimental X-ray resolution is not prediction RMSD. Do not conflate them.
- Do not tune defaults on a holdout you have already observed; start a new
  freeze cycle.

## 6. Repository boundary

- Edit **only this repository** (`FSOT-Genetics`).
- `FSOT-2.1-Lean` is the **read-only** mathematical authority. Never modify it,
  its parameter audits, or any sibling repository.

## 7. Gates and hygiene

- `python scripts/verify_cross.py` must pass (authority pin, seeds, zero free
  parameters, Lean + Mathlib). Do not merge or commit if it fails.
- Never bypass hooks (`--no-verify`) or discard unfamiliar in-progress files.
- Preserve user-edited files; never fold them into unrelated commits.

## Quick `D_eff` checklist before any structure/embedding change

1. Am I using a pinned ladder value, not an invented one?
2. Did I avoid setting any dimension to 3 "because space is 3-D"?
3. Am I embedding/scoring at `D_eff`, treating 3-D coordinates as a projection?
4. Does `scripts/verify_cross.py` still pass?

If any answer is "no," stop and fix it before proceeding.
