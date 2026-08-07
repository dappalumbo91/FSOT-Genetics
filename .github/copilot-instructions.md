# FSOT-Genetics — Copilot Instructions

This repo derives protein/genetics predictions from the **FSOT scalar law** with
**zero free parameters**. The full rules are in [`AGENTS.md`](../AGENTS.md). Read
it. The non-negotiables below are repeated here because Copilot loads this file
automatically.

## Non-negotiables

1. **Frozen engine.** `S = K(T1 + T2 + T3)`, computed only via
   `vendor/fsot_compute.py` at pin `D1D38A`. Never reimplement or approximate it.
2. **Zero free parameters.** Every constant traces to seeds (`π`, `e`, `φ`) or the
   pinned domain table. Never fit, tune, or regress against benchmark data.
   Scoring on datasets is fine; setting values from labels is not.
3. **`D_eff` is NOT a spatial dimension — this is the rule that always gets
   broken.** It is an effective/fractal dimension referenced to a 25-D baseline.
   - Never set `D_eff = 3`; never treat it as a count of x/y/z axes.
   - Use only the pinned ladder (`FSOTGenetics/ChemLink.lean`): backbone 8,
     disulfide 7, salt bridge 9, hydrophobic pack 14, H-bond 8, sidechain 9,
     tertiary biochem 13; base law 25.
   - The FSOT proximity object is intrinsically ~9-D and non-Euclidean
     (participation dim 9.10, negative-eigenvalue mass 0.20). Do not collapse it
     into 3-D MDS and call it "the structure." Reconstruct/evaluate at `D_eff`;
     3-D coordinates are an observer projection.
   - **When results degrade, check `D_eff` first.**
4. **Observer split.** FSOT determines low-D observables (radius of gyration,
   burial, secondary structure, native-dimension topology). A 3-D Cα RMSD needs
   an observer collapse that is underdetermined from a single sequence — do not
   claim single-sequence AlphaFold-grade coordinates.
5. **Frozen validation.** Freeze candidate + protocol commits before any holdout;
   record wins and losses unchanged. Experimental resolution ≠ prediction RMSD.
6. **Boundary.** Edit only this repo. `FSOT-2.1-Lean` is read-only authority.
7. **Gates.** `python scripts/verify_cross.py` must pass. Never use `--no-verify`.
   Preserve user-edited files.

See [`AGENTS.md`](../AGENTS.md) for the full charter and the `D_eff` checklist.
