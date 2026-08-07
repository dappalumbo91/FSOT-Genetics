# Contributing

## Doctrine

1. **Zero free parameters** on the claim path.
2. Prefer closed forms in {π, e, φ, γ} + domain scalars from `vendor/fsot_compute.py`.
3. Do not add trained weights, MSA-fitted dials, or Chou-Fasman lookup tables.
4. Keep the fold path **fast** (seconds-scale). O(n²)×huge multi-start is a regression.
5. Publish honest RMSD / contact metrics — no cherry-picking.

## Before push

```powershell
python scripts/verify_cross.py
```

## PR checklist

- [ ] `verify_cross.py` green
- [ ] No large binaries / caches
- [ ] Scoreboard updated if structure math changed
- [ ] Design/docs updated if F-layer math changed
