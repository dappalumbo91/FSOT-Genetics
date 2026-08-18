# FSOT vs AlphaFold — structure head-to-head

> **Historical (2026-08-07): sequence-only / bulk F01–F15.** The **13.9 Å** median below is the orphan path. It is **not** the product. Current same-data product: **0.13 Å** vs AF **0.47 Å** (`docs/PRODUCT_FREEZE.md`, `data/product_vs_alphafold.json`). Do not cite this report as present-day capability.

*Generated 2026-08-07T01:48:48.639942+00:00*

## Mission

FSOT **sequence-only (bulk)** structure prediction vs AlphaFold, scored on experimental PDB Cα RMSD. This is the orphan ceiling, not the measured-authority product.

- Engine: `fsot_protein_FULL_SCALAR_v10` · **free parameters: 0**
- Metric: Cα RMSD (Å) after Kabsch alignment to experimental PDB (lower is better)
- Hardware: HP Omen-class desktop; formula path (MDS+sparse), storage-capped cache
- FSOT median fold time: **299.16949998005293 ms**/chain · wall **17.673401099978946 s** (incl. downloads)

## Scoreboard

| Side | Median Cα RMSD (Å) | Wins |
|------|-------------------:|-----:|
| **FSOT** | **13.942122413158094** | **0** |
| AlphaFold | 0.43535089145183603 | 8 |
| Ties | — | 0 |

Paired proteins: **8** · FSOT win rate: **0.0**

## Per protein

| UniProt | Name | PDB | FSOT RMSD Å | AF RMSD Å | predict_ms | Winner |
|---------|------|-----|------------:|----------:|-----------:|:------:|
| P69905 | Hemoglobin alpha | 1A3N | 14.281 | 0.2697832925375286 | 272 | AlphaFold |
| P68871 | Hemoglobin beta | 1A3N | 14.532 | 0.5196952732909565 | 290 | AlphaFold |
| P00918 | Carbonic anhydrase II | 1CA2 | 16.826 | 0.3617948214686656 | 823 | AlphaFold |
| P00441 | SOD1 | 2C9V | 13.641 | 0.28634210236027413 | 314 | AlphaFold |
| P61626 | Lysozyme human | 1LZ1 | 13.660 | 0.43535089145183603 | 291 | AlphaFold |
| P61823 | RNase A | 7RSA | 13.942 | 0.3314052297214222 | 299 | AlphaFold |
| P0CG47 | Ubiquitin | 1UBQ | 11.326 | 1.6959974037473304 | 646 | AlphaFold |
| P01308 | Insulin | 4INS | 7.505 | 6.619892197937176 | 176 | AlphaFold |

## How to run

```powershell
python scripts/run_fsot_vs_alphafold_structure.py --max-proteins 8 --rounds 24 --sleep 0.2
```

## Next

- Accuracy: improve F15→D map / contact caps until RMSD competitive with AF
- Keep formula path seconds-scale (no O(n²) free-param grind)
- Optional Zig CLI for distogram parity + speed
- Expand benchmark set only after median RMSD moves
