# FSOT-Genetics — capability push roadmap (medical + research)

**Law:** \(S = K(T_1+T_2+T_3)\), pin `D1D38A`, **0 free / trained parameters**.  
**Data inputs allowed:** experimental homolog structures, MSAs (Pfam / JackHMMER / HHblits).  
**Not allowed:** fitting weights to PDB.

---

## Live stress scoreboard (this campaign)

Source: `data/medical_stress_suite.json` (10 classic medical/H2H proteins).

| Regime | Median Cα RMSD (Å) | Role |
|--------|-------------------:|------|
| AlphaFold DB | **0.47** | Competitor (trained + MSA) |
| **FSOT template + MSA packing fuse** | **1.16** | **Best FSOT deploy path** |
| FSOT template + physics | 1.19 | Prior best polish |
| FSOT template raw | 1.22 | Homolog transfer |
| FSOT bulk + MSA | 16.8 | Orphan fallback + evo channel |
| FSOT bulk single | 17.4 | Pure single-sequence claim path |

**Fuse beats raw template on 7/9** templated targets (packing-only coevolution clamps).  
**p53 variant panel:** 6/6 drivers flagged, mean **83.8th** percentile (conservation).  
**Cross-verify:** ALL GATES PASSED.

### Honest ceilings (do not market past these)

1. **De-novo bulk ~11–17 Å** — FSOT pairwise scalar has almost no long-range *distance* field (`diagnose_distance_wall.py`). Topology needs templates or a full many-body distance law.
2. **Contacts underdetermine structure** — even perfect contacts → ~11 Å MDS ceiling (`test_coevolution_fold.py`).
3. **MSA helps most as confidence + packing polish + medical conservation**, not as a from-scratch fold miracle.
4. **p53 DNA-binding** often lacks a clean self-excluded template in the current search → bulk/MSA/variant path matters medically.

---

## What “usable in medicine” requires

| Capability | Status now | Next push |
|------------|------------|-----------|
| Structure when homolog exists | **~1.2 Å, fuse ~1.16** | Multi-template medoid + domain assembly |
| Structure for orphans | ~11–17 Å bulk | Domain split + stronger D_eff observer; accept ceiling |
| Per-residue confidence | Provenance + evo conf | Fuse provenance×conservation; calibrate vs error |
| Variant effect | **p53 drivers 84th pctile** | Multi-gene panel (BRCA1/2, CFTR, …) + DNA front door |
| DNA → AA → effect | Working (`dna_variant_effect.py`) | ClinVar/COSMIC batch scorer |
| Cofactors / metals | Mapped to FSOT domains | Constrain template around validated sites |
| Explainability | Full scalar + trinary opcodes | Per-residue “why damaging” cards |
| Runtime / privacy | CPU seconds, no cloud train | Offline UniRef + local HHblits |
| Regulatory narrative | 0 trained weights, audited math | Repro kit + frozen holdouts |

---

## FSOT-appropriate accuracy levers (ranked)

### A. Already working — double down

1. **Template transfer (real homolog Cα)** — primary medical structure product.  
2. **Physics relax** (bond + clash, template-anchored).  
3. **MSA packing fuse v2** — only near-contact coevolution springs + intrinsic energy gate.  
4. **Conservation variant scoring** — medical win path independent of de-novo RMSD.  
5. **Regime auto-select** (`select_regime` / `fsot_predict.py`).

### B. High leverage, still pure FSOT

6. **Domain-aware fold** — split multi-domain chains by Pfam/InterPro ranges; fold/template each domain; assemble with FSOT coaxial / interface D_eff (p53, receptors).  
7. **Template ranking by coevolution agreement** — `score = coverage·identity·(1 + agreement/φ)` (data agreement, not fit).  
8. **Inter-template structural agreement confidence** — residual AF gap is partly domain orientation; multi-template variance → confidence (already flagged in prior commits).  
9. **DNA/RNA + protein joint** — regulatory interfaces (p53–DNA) as second system in chem-link D_eff (archive doctrine: multi-system observer).  
10. **Full-law residual as ranking energy** — use \(S=K(T_1+T_2+T_3)\) observer stress to pick among template candidates (already used in refine).  
11. **Local JackHMMER/HHblits + UniRef90** — orphan MSA depth without Pfam membership.  
12. **Codon-aware somatic panels** — batch `c.XXX` → trinary delta → conservation impact for tumor boards.

### C. Research depth (not “train a net”)

13. **Many-body distance field from FSOT fluid/geometry** — the missing non-contact distance law (wall diagnosis). This is the only principled way to break ~11 Å *without* templates.  
14. **Complexes / PPI** — interface D_eff routing between two chains.  
15. **Dynamics / allostery** — T3 chaos + observer hits as conformational ensemble weights.  
16. **CAMEO / CASP continuous** — external public scoreboard (BEAT_ALPHAFOLD_PLAN).  
17. **Formal Lean lemmas** for MSA channel amplitude = F09 family (Mathlib parity with chem-link).

### D. Do **not** do (breaks the claim)

- Train contact nets, distogram nets, or fine-tune on PDB.  
- Hide MSA/template use inside a “single-sequence” scoreboard.  
- Free scalar weights fitted to RMSD.

---

## Recommended execution order (next sessions)

```text
1. Expand medical gene panel (variant + structure) beyond p53
2. Domain-split template assembly for multi-domain disease proteins
3. Offline MSA tools (JackHMMER/HHblits) for hospital air-gapped use
4. Confidence card JSON for every residue (provenance × evo × clamp support)
5. Batch ClinVar driver test (precision/recall at fixed FSOT threshold)
6. Research track: many-body distance field from FSOT seeds only
```

---

## Commands

```powershell
cd C:\Users\damia\Desktop\FSOT-Genetics
python scripts/verify_cross.py
python scripts/run_medical_stress_suite.py
python scripts/fsot_predict.py --id 1UBQ --pdb-out predictions/ubq_fsot.pdb
python scripts/variant_conservation.py
python scripts/dna_variant_effect.py
```
