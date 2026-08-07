# D_eff multi-interface probe

*Generated 2026-08-07T01:43:18.934852+00:00*

Named pin-table routings only. **Zero free parameters.** Residual-at-interface diagnosis.

## Ranked by median Cα RMSD (offline samples)

| Rank | Routing | Median RMSD Å | Mean | n |
|-----:|---------|--------------:|-----:|--:|
| 1 | `legacy_v7` | 8.287 | 8.429 | 5 |
| 2 | `atomic_to_biochem` | 8.844 | 8.773 | 5 |
| 3 | `polymer_physchem` | 9.147 | 8.725 | 5 |
| 4 | `packing_dense` | 9.377 | 8.887 | 5 |
| 5 | `multi_scale_v9` | 9.408 | 8.939 | 5 |
| 6 | `bio_context` | 9.758 | 8.869 | 5 |

Default theory routing remains **`multi_scale_v9`** (not auto-switched by this table).

## Notes

Folding free-energy landscapes are often low-dimensional in reaction coordinates (Onuchic/Wolynes/Dill) while configuration space is high-D. FSOT D_eff is the domain interface in the 35-domain table — not PCA dim.

See `docs/DOMAIN_INTERFACE_FOLD.md` and `docs/RUNTIME_STACK.md`.
