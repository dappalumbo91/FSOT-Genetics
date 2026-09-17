# Genetics plant panel (fly pack is a separate system)

Fly experiment folder stays on its own pin. Do not copy Genetics `D1D38A` over it.

## Plant findings (this repo)

IntAct hops unchanged (gene-symbol class match already used real PHOT1 O48963 on the graph).

Product accession bug: Q2V2M9 is **human FHOD3**, not Arabidopsis PHOT1. Fixed to **O48963**.

| Item | Result |
|------|--------|
| PHOT1 O48963 full chain | leftover 5HZI id 0.61 cov 0.46 (below close-homolog 1/φ) |
| PHOT1 LOV1 205–301 | **2Z6C** id 1.00 |
| PHOT1 LOV2 485–577 | **4HHD** id 1.00 |
| PHOT1 kinase 665–952 | 4L3J id 0.55 — not close-homolog |
| UVR8 Q9FN03 | **8GQE** id 1.00 / 0.88; still missing from IntAct hops |
| SLAC1 | still ~0.001 — do not invent a phospho edge |

Not a plant connectome. Inter-domain PHOT1 pose is not a claimed fold.
