# MSA augmentation (optional evolutionary channel)

**Status:** optional evolutionary *data* channel on top of the zero-free-parameter engine.  
**Claim path:** the **product** (measured homologs, exclude eval PDB) — median **0.13 Å**.  
Single-sequence F01–F15 is the **orphan / bulk** fallback (~11–14 Å), not the product.

## Why this exists

The published product is measured authority + residual law (`docs/PRODUCT_FREEZE.md`).
The older “template ~1.2 Å / de-novo ~11 Å” sentence was the fair-cap + bulk era.
Those numbers are still true *as those regimes*; they are not current product accuracy.

For **practical** use (medical, design, shallow-MSA orphans vs deep families),
evolutionary signal is still on the table:

| Need | Single-sequence | MSA-augmented |
|------|-----------------|---------------|
| Orphan / novel fold | Primary path | Falls back cleanly |
| Deep family | Works, leaves coevolution unused | Injects MI+APC × conservation |
| Confidence | Geometry / template provenance | + per-residue evo confidence |
| Free parameters | 0 | 0 (MSA is **data input**) |

Prior engineering in this repo already measured that **contact-only coevolution MDS**
does not break the ~11 Å de-novo ceiling (`scripts/test_coevolution_fold.py`).
That is expected: contacts underdetermine full distance structure. The dual-mode
path here is therefore framed as **usability + robustness + confidence**, not as
a promise to beat the information ceiling with MI alone.

## Architecture

```
sequence
   │
   ├─ mode=single ──────────────────────► F01–F15 ► MDS/refine ► Cα
   │
   └─ mode=msa
         │
         ├─ MSA source (first available)
         │    1. --msa FILE  (Stockholm / A3M / FASTA)
         │    2. jackhmmer / hhblits  (if on PATH + FSOT_*_DB)
         │    3. Pfam full via InterPro  (--pfam / --uniprot)
         │
         ├─ features (closed-form stats)
         │    conservation, gap_frac, entropy, aa_freq, MI−APC, Neff
         │
         └─ inject long-range F15 term
              M ← M + evo_amp · Ĉ_ij · √(c_i c_j) · gap_penalty
              evo_amp = |S_biochem| · P_NEW · C_EFF / φ   (F09 family)
```

No trained weights. Amplitude reuses the same domain-scalar family as F09 region amplitude.

## API

```python
from fsot_structure_engine import predict_ca_coords
from msa_pipeline import build_msa_features

# Pure single-sequence (default — published claims)
r0 = predict_ca_coords(seq, mode="single")

# MSA-augmented
feat = build_msa_features(seq, pfam="PF00240")  # or msa_path=..., uniprot=...
r1 = predict_ca_coords(seq, mode="msa", msa_features=feat)
# r1["structure_mode"] in {"single_sequence", "msa_augmented"}
# r1["msa"]  → backend, n_seqs, neff, mean_conservation, ...
```

CLI:

```powershell
python scripts/run_msa_augmented_fold.py --id 1UBQ
python scripts/run_msa_augmented_fold.py --seq MQIFV... --pfam PF00240
python scripts/run_msa_augmented_fold.py --msa path\to\family.sto --seq ...
```

Local databases (optional):

```powershell
$env:FSOT_JACKHMMER_DB = "D:\dbs\uniref90"
$env:FSOT_HHBLITS_DB   = "D:\dbs\uniclust30"
```

## Relation to existing scripts

| Script | Role |
|--------|------|
| `msa_pipeline.py` | Shared MSA obtain + features + evo boost |
| `run_msa_augmented_fold.py` | Dual-mode smoke / report |
| `variant_conservation.py` | Medical conservation scoring (Pfam MSA) |
| `test_coevolution_fold.py` | Ceiling measurement (contact MDS only) |
| `test_confidence.py` | Template-provenance confidence |

## Honesty bar

- Do **not** mix MSA-augmented medians into the published single-sequence de-novo scoreboard without labelling the mode.
- Template regime remains the high-accuracy product path when homolog structures exist.
- MSA channel strengthens **when evolutionary support exists**; depth failure returns single-sequence automatically.
