# High-RMSD systems (>1.2 Å vs wet lab) — multi-system diagnosis

Source: `data/high_rmsd_system_diagnosis.json` + product path after multi-system residual levers.

## Systems above 1.2 Å (reality panel)

| System | RMSD | Primary mode | Termini frac | Dominant residual channel | Hypotheses |
|--------|-----:|--------------|-------------:|---------------------------|------------|
| sars2_rbd | 5.68 | global_topology | 0.53 | bond | termini + tertiary + data≠residual |
| hiv_rt | 5.46 | long_range | 0.00 | fold_Rg | multi-domain / Rg; data≠residual |
| abl1 | 4.93 | long_range | 0.00 | bond | data≠residual; kinase state |
| ace2 | 4.74 | long_range | 0.40 | bond | termini + remote + data≠residual |
| bcl2 | 3.75 | long_range | 0.27 | bond | data≠residual |
| p53_dbd | 2.60 | long_range | 0.67 | bond | **termini_disorder** dominant |
| egfr | 2.51 | long_range | 0.87 | bond | **termini_disorder** |
| ha_h3 | 2.33 | long_range | 0.07 | fold_Rg | multi-domain packing |
| ubiquitin | 2.03 | long_range | 0.40 | fold_Rg | C-term flexibility |
| kras | 1.59 | long_range | 0.20 | fold_Rg | switch-region / state |

## Cross-system variables (shared pin)

| Domain | Role for high-error set |
|--------|-------------------------|
| Physical_Chemistry | Bond residual — elevated when gap interpolation / termini wrong |
| Chemistry | Clash residual — secondary here |
| Biochemistry | Fold/Rg residual — multi-domain & conformation state |
| Molecular_Chemistry / Condensed_Matter | ChemLink packing in structure engine (bulk path) |

## Behaviors → FSOT handle (applied)

1. **termini_disorder** (p53, EGFR, Ubq, RBD, ACE2)  
   Soft CA_CA termini rebuild **only if** termini bond stress > core × `r_bond` (Physical_Chemistry residual-at-interface). Candidate competes on residual energy in `fuse_predict`.

2. **data_best residual-unfit** (ABL1, BCL2, …)  
   Residual override **only if** `E_data/E_res > φ` and residual-best still data-plausible (`score ≥ best/φ`). Does not free-rank CaM/RNase.

3. **fold_Rg / multi-domain** (HA, HIV-RT, KRAS)  
   Next: domain-split product path (existing `domain_split_assemble`) for multi-domain wet-lab chains.

4. **tertiary long_range** (almost all)  
   Next: ChemLink residual on measured long-range contacts only (Biochemistry D=13) — not bulk invent.

## Gate after multi-system product v6

| Guard | RMSD |
|-------|-----:|
| Product H2H median | **1.15 Å** |
| RNase | 0.77 |
| CaM | 0.89 |
| Insulin (override) | **0.87** |

Re-run: `python scripts/diagnose_high_rmsd_systems.py`
