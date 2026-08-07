# Cross-verification — genetics formula branch

Modeled on **FSOT-2.1-Lean** multi-prover / green-gate discipline:

- Same **D1D38A** scalar pin as Lean and neuron-zig.
- Hard CI exit on pin drift, free parameters, or grind-time regression.
- Scoreboards are evidence, not marketing.

## Local

```powershell
cd FSOT-Genetics
python scripts/verify_cross.py
# optional (network):
python scripts/run_fsot_vs_alphafold_structure.py --max-proteins 8 --rounds 24 --sleep 0.2
python scripts/run_fsot_distogram_contact_eval.py
```

## What is *not* a green gate

| Item | Status |
|------|--------|
| Median Cα RMSD beating AlphaFold | **Open research goal** — tracked, not required for CI green |
| Full proteome dump | Out of scope (storage) |
| Neural weight training | Forbidden on claim path |

CI green means: **law pin + zero free params + formula path runs**.  
Campaign green (AF) means: **median RMSD / lDDT competitive** — see `BEAT_ALPHAFOLD_PLAN.md`.

## Sibling verification map

```
FSOT-2.1-Lean     → Lean / multi-prover / domain margins
fsot-neuron-zig   → seed + genetic pair geometry parity
FSOT-Genetics     → F01–F15 + fold + AF H2H (this repo)
```
