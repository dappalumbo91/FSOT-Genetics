# FSOT-Genetics — Design

**Pin:** D1D38A (`vendor/fsot_compute.py`)  
**Law:** \(S = K(T_1 + T_2 + T_3)\)  
**Free parameters:** **0**

## Mission

Mathematically solve sequence → structure (and related genetics observables)  
**without** neural-network weights, MSA-trained free dials, or O(n²) grind that apes deep learning cost.

Sibling of:

| Repo | Role |
|------|------|
| [FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean) | Formal law hub, multi-domain verification |
| [fsot-neuron-zig](https://github.com/dappalumbo91/fsot-neuron-zig) | Neural mind + genetic pair geometry |
| **This repo** | Genetics / protein formula branch + AF head-to-head |

## Architecture

```
sequence
   │
   ▼
F01 trinary AA phases
   │
   ├─ F02–F06 chemistry (h, V, q, μ, disulfide, hydrophobic, elec, dipole)
   ├─ F07 backbone bb(s) = s^(-1/π)
   ├─ F08 env(s) = s/(s + πe)
   ├─ F10–F11 helix / sheet bonuses
   ├─ F12–F14 regions + register
   └─ F15 proximity matrix M_ij
           │
           ▼
    proximity → distance D (seed geometry + contact caps)
           │
           ▼
    classical MDS embed + sparse polish (O(n·k), not O(n²)×N)
           │
           ▼
    Cα coordinates  ·  timed (predict_ms)
           │
           ▼
    Kabsch RMSD vs experimental PDB  ·  dual scoreboard vs AlphaFold DB
```

That diagram is the **bulk / orphan** path only (information ceiling ~11–14 Å). It is **not** the product and must not be quoted as “FSOT is 15 Å off.”

**Product path (medical structure, freeze 2026-08-17):**

```
sequence + exclude-eval-PDB
   │
   ▼
RCSB near-self (id ≥ 1/φ) + query UniProt + UniRef100 other-accession + Pfam
   │
   ▼
measured Cα transfer  →  intact? keep  :  residual fuse
   │
   ▼
state_reps = intact maps of each collapse (NMR Superposed, not residual-best)
   │
   ▼
apparatus = min over collapses   ·   0.13 Å median vs AF 0.47 Å
```

See `docs/PRODUCT_FREEZE.md`, `docs/FSOT_APPLICATION.md`.

## Non-negotiables

1. **Zero free parameters** — only {π, e, φ, γ} + domain scalars from the pin.
2. **Formula speed** — fold path stays seconds-scale on home hardware (Omen-class).
3. **Honest metrics** — Cα RMSD to experiment; AF win/loss table published as-is.
4. **Cross-verification** — `scripts/verify_cross.py` is a hard green gate (like Lean margin gates).
5. **No neural nets** for the claim path. Optional lab tools stay offline / out of the repo claim.

## AlphaFold targeting order

Priority is **accuracy on proteins AF currently wins**, without re-introducing free params:

1. Hemoglobin α/β, carbonic anhydrase, SOD1, lysozyme, RNase, ubiquitin (H2H set).
2. Improve F15 → D map and long-range contact ranking first.
3. Region packing / chirality only if stress-selected (no dials).
4. Expand benchmark only after median RMSD moves.

See `docs/BEAT_ALPHAFOLD_PLAN.md` and live scoreboard under `data/`.

## Cross-verification (Lean-style)

| Gate | What |
|------|------|
| Authority pin | SHA-256 of `vendor/fsot_compute.py` starts with `D1D38A` |
| Free params | Engine returns `free_parameters: 0` |
| Speed | Short-chain predict hard-capped (ms-scale) |
| Finite fold | Cα coords finite, correct shape |
| Derivations | F01–F15 doc present |

```powershell
python scripts/verify_cross.py
```

CI runs the same gate on every push.

## Rust crates

| Crate | Role |
|-------|------|
| `codon_core` | Trinary codon core (no_std) |
| `fsot_core` | Scalar engine + domain table (Rust twin) |
| `fsot_protein` | F01–F15 authority port (secondary, chemical, distogram, regions) |

Python `scripts/fsot_structure_engine.py` is the **fast fold runtime** used for AF H2H.  
Rust remains the **formula authority** for chemistry / SS / distogram parity work.

## Storage / hardware

Designed for HP Omen-class desktop. No multi-GB optical catalogs.  
PDB/AF downloads cache locally; sample PDBs under `data/pdb_samples/` for offline smoke.
