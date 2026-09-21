# Genetics system verify (fly pack is a separate system)

Pin **D1D38A**. Do not copy this pin onto the fly pack.

Full claim inventory: `docs/SYSTEM_VERIFY.md`.

- `python scripts/system_verify.py` → **259/259** (Genetics claims; fly pack is separate)
- `python scripts/verify_cross.py` → PASS
- `python verification/run_cross_proof.py` → **overall_ok**, 42 obligations, 10 layers

Plant accession audit remains in `docs/RESULTS.md`. CASP/CAMEO still open (`docs/OPEN.md`). Biohub/Kaggle frozen.
