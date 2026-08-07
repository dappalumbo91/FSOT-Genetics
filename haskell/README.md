# FSOT-Genetics — Haskell layer

**Role:** compile-time lock for chem-link → `D_eff`, claim path (zero free params),
and contact-ranking types. **Not** a full rewrite of the fold (yet).

| Language | Job on this project |
|----------|---------------------|
| **Python** | Fast Omen lab: F15, MDS, error log, SMILES bridge |
| **Lean + Mathlib** | Prove interface Nat/`ℝ` facts |
| **Haskell** | Types that **will not compile** if free-D / free weights sneak in |
| **Zig / Rust** | Neuron twin + formula crates |
| **Roc** | Interesting for pure functional apps; weaker scientific ecosystem here |

## Build

```powershell
cd haskell
cabal build
cabal run fsot-genetics-check
```

## Doctrine

Same as Lean: named pin domains only, observer by chemical system, free parameters = 0.
Python must stay parity with `ChemLink` / `Contact` tables exported here.
