# Biohub / Kaggle — freeze (2026-09-07)

**Status: back burner.** Do not grind this for the competition clock. Pick up from here when organism tracking is the job again.

Competition: [Biohub – Cell Tracking During Development](https://github.com/dappalumbo91/FSOT-Genetics). Data stays on `D:\Kaggle_Biohub_Data` (199 train GEFF+Zarr, 4 test Zarr). Not in git.

## What is shipped (photon claim, 0 free params)

Pin `D1D38A`. Native observer: MAD+φ gate, φ³ NMS, residual second collapse, product graph (primary + isolated residual, in-shell fold), leftover-yield linker (2-opt at φ⁵, leftover-length dests yield).

| Video | Find 7 | Jaccard (adj.) | Notes |
|-------|-------:|----------------|-------|
| `44b6_0113de3b` proxy | 1.00 | **1.00** | closed |
| `6bba_09961292` dense | 0.96 | **0.82 (0.83)** | leftover-yield `4dbacd6` |

Git: `4dbacd6` last photon linker ship. Docs: `docs/BIOHUB_3D.md`.

## Eye / relay (lab apparatus — not the claim)

ML paints a **field**; FSOT still reports the center. Net must not emit annotator xyz.

| Path | Dense Jac | Proxy | Note |
|------|----------:|------:|------|
| Steal-shell correspondence fill | **0.848** | 1.00 | holdout `6bba_09961292`; weights `D:\_fsot_eye_cache\correspondence\edge_predictor_correspondence.pth` |
| Low-contrast curriculum | 0.839 (holdout drop) | 1.00 | **0.67** on `0b24845f` (native find 0.10); **0.61** on `05db0fb1` (steal-shell was 0.64) |

Do **not** replace steal-shell with low-contrast. Extreme gate-miss only (`0b248`: GT 1337 vs MAD-gate 2432).

Variety panel: `data/biohub_eye_variety.json`. Scripts: `scripts/_eye_relay.py`, `_eye_correspondence_train.py`, `_eye_panel.py`.

## Why we stopped

Held-out dense 0.848 is the Kaggle **Final number on one train video**, not public LB and not live top 0.962. Native find is not 0.96 on low-contrast volumes (105/197 train videos paint-frac < 1/φ³). A working T4 kernel was not proven (v64 leftover `import tracksdata` died). Do not Final until public > 0.848.

## Anti-goals (still)

No U-Net wrap as the product. No leftover-all-tier mix. No greedy steal. No 7 µm fold / third NMS collapse (all measured regressions).

## Resume command

```text
python -u scripts/biohub_3d.py --voxels --dataset 6bba_09961292
python -u scripts/_eye_relay.py 6bba_09961292 corr
```

Photon claim first. Eye fill only with steal-shell weights unless the volume is the `0b248` gate-miss class.
