# SMILES Lab → protein pair chemistry

## Source

**FSOT SMILES Lab** (~1470 seed-only chemistry solves) + Lean vendor twin:

- Desktop: `FSOT SMILES Lab/FSOT_SMILES_Lab_Dataset.json`
- Lean: `vendor/smiles/FSOT_SMILES_Lab_Dataset.json`
- Genetics extract: `formulas/smiles_protein_chemistry.json` (protein-relevant sections)

## Sections wired into fold chemistry

| Section | Use in structure path |
|---------|------------------------|
| **§36 Hydrophobicity** (20 AA) | Core packing: hydrophobic×hydrophobic only |
| **§22 Amino acid pKa** | Formal charge at pH = φ⁻⁴+φ⁴ = 7; electrostatic F05 |
| **§21 Protein ΔG** | Reference scale (stability context) |
| **F18 disulfide gate** | C–C modulated by sep envelope at πe |
| **§25/§26/§43/§96** | Stored for next layers (radii, polarizability, bonds) |

## Law

Still **zero free parameters**: SMILES formulas are closed forms in {π,e,φ,γ,G} + Layer 1/2.  
We load **computed** pin values from the lab dataset (same authority as Lean verification).

## Code

- `scripts/smiles_aa_chem.py` — tables + `smiles_expanded_interaction`
- `scripts/trinary_syntax.py` — `expanded_chemical_interaction` prefers SMILES bridge
- `scripts/fsot_structure_engine.py` — F15 uses that chemistry

## Rebuild extract

```powershell
# re-run extract from Lean vendor SMILES dataset if lab updates
python -c "..."  # see git history / smiles_aa_chem docs
```
