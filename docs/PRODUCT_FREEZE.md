# Product freeze — structure path (2026-08-12)

**Status:** Residual-fit transfer filter on measured authority. Gate still holds.  
**Authority pin:** `D1D38A` · law \(S = K(T_1+T_2+T_3)\) · **0 free parameters**.

## Frozen product numbers (H2H, 10 proteins)

Source: `data/product_vs_alphafold.json` (residual-fit data authority).

| metric | value |
|--------|------:|
| AlphaFold median Cα RMSD | **0.47 Å** |
| FSOT **product** median | **1.14 Å** |
| FSOT template median | **1.17 Å** |
| FSOT bulk (orphan) median | 13.57 Å |
| Product within 1.5 Å of AF | **10/10** |
| Product sub-2 Å | **10/10** |
| **Wet-lab reality panel median** | **1.17 Å** (`reality_margin_eval`) |

Notable product Cα RMSDs:

| protein | product Å | AF Å |
|---------|----------:|-----:|
| Insulin | 1.10 | 4.51 (product wins) |
| RNase A | 0.47 | 0.33 |
| Calmodulin | 0.75 | 6.45 (product wins) |
| Lysozyme | 1.16 | 0.42 |
| Hemoglobin α | 1.12 | 0.27 |
| Ubiquitin | 1.48 | 0.88 |
| p53 DNA-binding | 1.57 | 6.19 (product wins) |
| Carbonic anhydrase II | 1.30 | 0.36 (best fair homolog; id cap) |

## What the product *is*

```
measured homolog Cα
  (id×cov authority; residual-fit filter if data-best transfer is unfit)
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

Gate for any structure change still: **do not ship if product median > 1.15 Å** on this freeze set (or a newly preregistered holdout).
