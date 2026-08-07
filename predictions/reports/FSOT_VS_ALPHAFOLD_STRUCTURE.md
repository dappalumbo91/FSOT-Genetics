# FSOT vs AlphaFold — structure head-to-head

*Generated 2026-08-07T01:43:59.003982+00:00*

## Mission

FSOT sequence-only structure prediction vs AlphaFold, scored on experimental PDB Cα RMSD

- Engine: `fsot_protein_F01_F15_deff_v9` · **free parameters: 0**
- Metric: Cα RMSD (Å) after Kabsch alignment to experimental PDB (lower is better)
- Hardware: HP Omen-class desktop; formula path (MDS+sparse), storage-capped cache
- FSOT median fold time: **258.9278000232298 ms**/chain · wall **19.080174100003205 s** (incl. downloads)

## Scoreboard

| Side | Median Cα RMSD (Å) | Wins |
|------|-------------------:|-----:|
| **FSOT** | **13.858567478931747** | **0** |
| AlphaFold | 0.43535089145183603 | 8 |
| Ties | — | 0 |

Paired proteins: **8** · FSOT win rate: **0.0**

## Per protein

| UniProt | Name | PDB | FSOT RMSD Å | AF RMSD Å | predict_ms | Winner |
|---------|------|-----|------------:|----------:|-----------:|:------:|
| P69905 | Hemoglobin alpha | 1A3N | 13.651 | 0.2697832925375286 | 234 | AlphaFold |
| P68871 | Hemoglobin beta | 1A3N | 14.010 | 0.5196952732909565 | 245 | AlphaFold |
| P00918 | Carbonic anhydrase II | 1CA2 | 16.653 | 0.3617948214686656 | 711 | AlphaFold |
| P00441 | SOD1 | 2C9V | 13.859 | 0.28634210236027413 | 269 | AlphaFold |
| P61626 | Lysozyme human | 1LZ1 | 13.146 | 0.43535089145183603 | 251 | AlphaFold |
| P61823 | RNase A | 7RSA | 14.148 | 0.3314052297214222 | 259 | AlphaFold |
| P0CG47 | Ubiquitin | 1UBQ | 11.423 | 1.6959974037473304 | 563 | AlphaFold |
| P01308 | Insulin | 4INS | 6.724 | 6.619892197937176 | 151 | AlphaFold |

## How to run

```powershell
python scripts/run_fsot_vs_alphafold_structure.py --max-proteins 8 --rounds 24 --sleep 0.2
```

## Next

- Accuracy: improve F15→D map / contact caps until RMSD competitive with AF
- Keep formula path seconds-scale (no O(n²) free-param grind)
- Optional Zig CLI for distogram parity + speed
- Expand benchmark set only after median RMSD moves
