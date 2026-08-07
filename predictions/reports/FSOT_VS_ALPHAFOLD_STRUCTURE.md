# FSOT vs AlphaFold — structure head-to-head

*Generated 2026-08-07T01:37:03.938192+00:00*

## Mission

FSOT sequence-only structure prediction vs AlphaFold, scored on experimental PDB Cα RMSD

- Engine: `fsot_protein_F01_F15_trinary_v8` · **free parameters: 0**
- Metric: Cα RMSD (Å) after Kabsch alignment to experimental PDB (lower is better)
- Hardware: HP Omen-class desktop; formula path (MDS+sparse), storage-capped cache
- FSOT median fold time: **267.58060001884587 ms**/chain · wall **16.770663200004492 s** (incl. downloads)

## Scoreboard

| Side | Median Cα RMSD (Å) | Wins |
|------|-------------------:|-----:|
| **FSOT** | **13.901767727888158** | **0** |
| AlphaFold | 0.43535089145183603 | 8 |
| Ties | — | 0 |

Paired proteins: **8** · FSOT win rate: **0.0**

## Per protein

| UniProt | Name | PDB | FSOT RMSD Å | AF RMSD Å | predict_ms | Winner |
|---------|------|-----|------------:|----------:|-----------:|:------:|
| P69905 | Hemoglobin alpha | 1A3N | 13.689 | 0.2697832925375286 | 232 | AlphaFold |
| P68871 | Hemoglobin beta | 1A3N | 14.006 | 0.5196952732909565 | 272 | AlphaFold |
| P00918 | Carbonic anhydrase II | 1CA2 | 16.626 | 0.3617948214686656 | 701 | AlphaFold |
| P00441 | SOD1 | 2C9V | 13.902 | 0.28634210236027413 | 268 | AlphaFold |
| P61626 | Lysozyme human | 1LZ1 | 13.168 | 0.43535089145183603 | 249 | AlphaFold |
| P61823 | RNase A | 7RSA | 14.150 | 0.3314052297214222 | 257 | AlphaFold |
| P0CG47 | Ubiquitin | 1UBQ | 11.406 | 1.6959974037473304 | 551 | AlphaFold |
| P01308 | Insulin | 4INS | 6.705 | 6.619892197937176 | 148 | AlphaFold |

## How to run

```powershell
python scripts/run_fsot_vs_alphafold_structure.py --max-proteins 8 --rounds 24 --sleep 0.2
```

## Next

- Accuracy: improve F15→D map / contact caps until RMSD competitive with AF
- Keep formula path seconds-scale (no O(n²) free-param grind)
- Optional Zig CLI for distogram parity + speed
- Expand benchmark set only after median RMSD moves
