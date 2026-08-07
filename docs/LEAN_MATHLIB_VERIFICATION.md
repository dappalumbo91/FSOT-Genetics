# Lean + Mathlib verification — FSOT-Genetics

## Purpose

Formal lock on the genetics fold **interface law**:

- chem-link → pin-table `D_eff` (no free continuous D)
- observer on/off by chemical system
- free parameters = 0
- seed / residual inequalities on Mathlib `ℝ`

Runtime multiprecision `S` remains **Python pin D1D38A** / sibling
[FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean) `FSOT.Scalar`.

## Toolchain (aligned with FSOT-2.1-Lean)

| Component | Version |
|-----------|---------|
| Lean | `leanprover/lean4:v4.31.0` (`lean-toolchain`) |
| Mathlib | `mathlib4 @ v4.31.0` via `lakefile.lean` |

**Yes — Mathlib is required.** This is not a std-only stub.

## Build

```powershell
cd FSOT-Genetics
lake update          # first time: fetches Mathlib (large)
lake build           # formal gate (typecheck all modules)
```

On Windows, a Mathlib-linked `lean_exe` can fail the PE **65k export limit**.  
The gate is **`lake build`** (typecheck), not a fat CLI binary.

## Modules

| File | Role |
|------|------|
| `FSOTGenetics/Seeds.lean` | Pin string, freeParameterCount=0, F13 gate Nat |
| `FSOTGenetics/SeedsReal.lean` | `phiR`/`eR`/`piR` on `ℝ`; residual ≥ 1; chaos at D=25 |
| `FSOTGenetics/ChemLink.lean` | 7 chem links → D_eff; Finset pin range |
| `FSOTGenetics/Observer.lean` | Observer policy theorems |
| `FSOTGenetics/ZeroFreeParams.lean` | **Lock** claim path freeParameters=0 (does *not* introduce free params) |

## Python dual

`scripts/verify_cross.py` runs:

1. Python full-law / SMILES / fold smoke  
2. `lake build` + `lake exe fsot_genetics_check` when Lake is on PATH  
3. Parity of chem-link D_eff table vs Lean constants  

## CI

`.github/workflows/ci.yml` — Python job + Lean/Mathlib job.
