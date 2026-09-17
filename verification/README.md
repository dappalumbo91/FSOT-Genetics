# FSOT-Genetics cross-proof gauntlet

Independent re-proof of Genetics engine identities **and measured catalog gates** across Lean, Coq/Rocq, Isabelle, F*, SMT (Z3), Rust f64, TLA+, and the Python D1D38A oracle.

This tree does **not** inherit `FSOT-2.1-Lean` `overall_ok`. Every layer re-runs here.

## Run

```powershell
cd FSOT-Genetics
python scripts/system_verify.py
python verification/export_obligations.py
python verification/run_cross_proof.py
```

Expect `data/cross_proof_report.json` → `overall_ok: true`.

## What is stamped

| Layer | What | Source |
|-------|------|--------|
| **A** | Pin D1D38A, 0 free params, φ, leftover `1/φ²`, close-homolog `1/φ`, 7 chem-link D_eff, observer (backbone unobserved), F13 gate = 7 | Lean Mathlib + Nat catalog |
| **B** | Product median 0.13 Å < AF 0.47 Å; hop splits (JO/VNC light motor, olfactory does not); insect homologs 26/1 miss; plants 6+12; analog pointer; measured graphs nonempty | live `data/*.json` |

Layer B triangulates **exported numeric gates**, not a re-derivation of FlyWire pixels from type theory.

## Tools (used if present)

| Framework | How |
|-----------|-----|
| Python oracle | `scripts/verify_cross.py` |
| SMT python | `verification/smt_replay.py` (always) |
| Z3 | `verification/smt/genetics_bounds.smt2` |
| Lean 4 + Mathlib | `lake build` (`FSOTGenetics/Catalog.lean`) |
| Coq / Rocq | `verification/coq/GeneticsSpine.v` |
| Isabelle/HOL | `verification/isabelle` session `FSOTGenetics` |
| F* | `verification/fstar/FSOTGenetics.fst` |
| Rust | `verification/rust_genetics_kernel` |
| TLA+ TLC | `verification/tla/GeneticsSolve.tla` |

Missing optional tools SKIP. A present tool that fails FAIL.

Labeled archive: `docs/VERIFIED_SOLVES.md`.
