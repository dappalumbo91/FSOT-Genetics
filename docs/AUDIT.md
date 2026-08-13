# Repository audit — 2026-08-13

Compared **code + `data/product_vs_alphafold.json` + `data/af_coverage.json`** to every public claim surface. This file is the punch list that the same-day README / freeze / application updates close.

## Authority (current)

| Surface | Status |
|---------|--------|
| Pin `D1D38A` | Current. `vendor/fsot_compute.py`, `AGENTS.md`, Lean `FSOTGenetics/`. |
| Product H2H median **0.13 Å** vs AF **0.47 Å** | Current. `data/product_vs_alphafold.json` (`generated_at` 2026-08-13). |
| Coverage jobs including H, modified NA, tetramer | Current. `data/af_coverage.json`. |
| `docs/PRODUCT_FREEZE.md` | Brought to 0.13 Å table. |
| `docs/AF_COVERAGE.md` | Brought to coverage JSON. |
| `docs/OPEN.md` | New. Remaining failures only. |

## Was stale (fixed in this pass)

| Surface | Old claim | Truth |
|---------|-----------|-------|
| `README.md` | Product median **1.16 Å**, 9/10 sub-2 Å; freeze “paused” 2026-08-11 | **0.13 Å**, 10/10; work continued |
| `README.md` | H2H template 1.2 Å; CaM 0.77 | Fair-cap handicap era; product CaM is **0.90**, ubq **0.09** |
| `docs/FSOT_APPLICATION.md` | “PRODUCT physics always”; ship gate **≤ 1.16 Å** | Intact maps stay raw; gate is AF median **0.47 Å** |
| `docs/DESIGN.md` | Architecture ends at F15 MDS | Product path is measured authority + apparatus |
| `docs/FIELD_READY.md` | Freeze gate ≤ 1.16 Å | ≤ 0.47 Å (AF); freeze number 0.13 |
| `docs/BEAT_ALPHAFOLD_PLAN.md` | Implied no H2H win yet | Same-data product **beats AF median** on the freeze set (not CASP) |
| `docs/MECHANISM_GAP_MAP.md` | Controls 1.47 Å; p53 2.60 | Historical wet-lab snapshot — labeled as such; H2H is 0.13 |
| `field/console_data.json` | Product 1.16 Å (git 59c889f, 2026-08-11) | Rebuild from current JSON |

## Intentionally historical (do not overwrite)

| Artifact | Why keep |
|----------|----------|
| `data/alphafold_headtohead.json` | Fair-cap 0.95 H2H (handicap). |
| `data/wetlab_af_eval.json` | Broader medical set, older protocol. |
| `data/m1_authority_verify.json` | ABL1/BCL2/EGFR M1 experiment. |
| `docs/MECHANISM_GAP_MAP.md` body | Failure-mode ledger; still the kinase/RBD map. |
| `predictions/reports/*` | Dated reports. |

## Code orientation (current entry points)

| Job | Script |
|-----|--------|
| Freeze H2H | `scripts/bench_product_vs_af.py` |
| AF3 coverage | `scripts/bench_af_coverage.py` + `scripts/multi_system.py` |
| Authority / apparatus | `scripts/run_rcsb_template_holdout.py` |
| Fuse / intact halt | `scripts/msa_template_fuse.py` |
| Bulk F01–F15 | `scripts/fsot_structure_engine.py` |
| Residual / domains | `scripts/full_scalar_law.py` |
| Pin gate | `scripts/verify_cross.py` |

Probe scripts `scripts/_*.py` are **not** product. Do not ship them.

## Claim language

Say: *same-data product, exclude eval PDB, 0 free parameters, residual at named ChemLink, apparatus min over intact crystals.*

Do not say: *de-novo fold beats AlphaFold* or *we trained nothing and still recover 0.13 Å from sequence alone.*
