# Product freeze — structure path (2026-08-11)

**Status:** Python research ceiling for this cycle. Further RMSD chasing vs AlphaFold is **paused**.  
**Authority pin:** `D1D38A` · law \(S = K(T_1+T_2+T_3)\) · **0 free parameters**.

## Frozen product numbers (H2H, 10 proteins)

Source: `data/product_vs_alphafold.json` (commit lineage through multi-template coverage).

| metric | value |
|--------|------:|
| AlphaFold median Cα RMSD | **0.47 Å** |
| FSOT **product** median | **1.15 Å** (reality-first M1 + length-sim) |
| FSOT template median | ~1.20 Å |
| FSOT bulk (orphan) median | 13.57 Å |
| Product within 1.5 Å of AF | **10/10** |
| Product sub-2 Å | **9/10** |
| **Wet-lab reality panel median** | **1.17 Å** (`reality_margin_eval`) |

Notable product Cα RMSDs:

| protein | product Å | AF Å |
|---------|----------:|-----:|
| Hemoglobin α | 1.02 | 0.27 |
| RNase A | 0.77 | 0.33 |
| Calmodulin | 0.82 | 6.45 (product wins) |
| Insulin | 1.14 | 4.51 (product wins) |
| Ubiquitin | 1.78 | 0.88 |
| p53 DNA-binding | 2.60 | 6.19 (product wins) |

## What the product *is*

```
measured homolog Cα (multi-template ensemble)
  → residual-weighted physics
      bond   ← Physical_Chemistry · residual (1+|S|·P_NEW)
      clash  ← Chemistry
      anchor ← Biochemistry (template authority)
  → Cα model
```

- Templates = **measured authority** (RCSB/Pfam homologs; identity cap 0.95 fair H2H).
- MSA = **data**, not trained weights.
- Bulk de-novo ≈ 11–14 Å is the single-sequence information ceiling — not a bug.

## What we are *not* doing next

- Geometric shotgun (medoid-all, false contacts, residual springs that tank median).
- More Python RMSD loops as the primary product.
- Claiming AlphaFold-beating bulk fold without more measured coverage / hardware path.

## What we *are* doing next

**Polish + bare metal.** See `docs/BARE_METAL_GENETICS_ROADMAP.md`.

Python under `scripts/` remains the **research bench and oracle** for metrics.  
**Runtime product** moves to Zig (host + freestanding QEMU) / optional Rust `no_std`, same pin, same residual law, trinary codon front door.

Gate for any structure change still: **do not ship if product median > 1.16 Å** on this freeze set (or a newly preregistered holdout).
