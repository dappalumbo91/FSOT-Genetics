# Language stack — what each language is for

Python was never chosen because it is “the best language for genetics.”
It was chosen as a **fast lab** on the Omen. That has limits. Here is the honest split.

| Language | Strength for *this* project | Weakness | Role now |
|----------|----------------------------|----------|----------|
| **Python** | Speed of iteration, NumPy MDS, PDB I/O, error logs | Silent runtime nonsense; easy to drift | Fold lab + scoreboards |
| **Haskell** | Types that refuse illegal programs at compile time; pure laws | Rewrite cost; fewer PDB/NumPy-class libs | **Chem-link + contact law locks** (`haskell/`) |
| **Lean + Mathlib** | Proofs of interface / residual facts | Not a fold runtime | Formal D_eff / zero free params |
| **Zig** | Bare-metal / fixed-point twin (neuron) | Less science ecosystem | Pair geometry twin |
| **Rust** | Formula crates (`fsot_protein`) | Rebuild time | F01–F15 authority port |
| **Roc** | Friendly pure functional, good for apps | Young ecosystem for this science stack | Optional later; not blocking |

## Why not “all Haskell” tomorrow

1. The **error modes** (topology, long-range contacts) are scientific, not only language bugs.  
2. A full rewrite pauses the residual campaign.  
3. We already pay Lean for *proof* and Haskell can pay for *typed law* without throwing away the lab.

## Intended growth path

```text
Haskell ChemLink / Contact types  →  compile fails if free-D appears
        ↓ parity
Python implements same tables     →  fold + error_margin_log
        ↓ proofs
Lean ChemLink D_eff               →  Mathlib residual lemmas
```

When contact ranking stabilizes, **port the hot path** (consensus score + maybe MDS) to Haskell or Zig — not the whole monorepo at once.

## Commands

```powershell
# Haskell
cd haskell
cabal build
cabal run fsot-genetics-check

# Lean
cd ..
lake build

# Python residual campaign
python scripts/run_error_margin_log.py
```
