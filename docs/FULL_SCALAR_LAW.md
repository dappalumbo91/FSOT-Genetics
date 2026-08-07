# Full scalar law on the fold path

## Why this exists

Using only frozen `|domain_scalar|` amplitudes throws away most of FSOT:

- observer quirk on T1  
- chaos / poof / suction on T3  
- live \(\delta\psi\), hits, \(D_{\mathrm{eff}}\)  

That was wrong. The fold now runs the **whole** formula.

## Formula

\[
S = K(T_1 + T_2 + T_3)
\]

| Term | Role on the fold |
|------|------------------|
| **T1** | Observer-modulated base; **`observed=True`** on tertiary pairs and every refine round |
| **T2** | scale · amplitude (pin defaults) |
| **T3** | Valve + **chaos** \((D-25)/25\) + poof/suction + acoustic \(\delta\theta\) |

**Residual scale** (same spirit as residual packs):

\[
\text{term} \leftarrow \text{term}\cdot(1 + |S|\cdot P_{\mathrm{NEW}})
\]

## Where it runs

| Stage | Full law usage |
|-------|----------------|
| Distogram pair \((i,j)\) | `pair_full_scalar(sep, spin, charge, …)` → S, residual, chaos, domain interface |
| Backbone sep≤2 | S with `observed=False` (geometry) |
| Local / SS / tertiary | Domain table D_eff + trinary-modulated \(\delta\psi\) |
| Refine round \(r\) | `refine_observation_scalar` — **observer ON**, hits\(=r\) |
| Output | Reports T1, T2, T3, observer_mod, chaos_factor |

## Code

- `scripts/full_scalar_law.py` — float twin of vendor `compute_scalar` + pair/refine helpers  
- `scripts/fsot_structure_engine.py` — engine `fsot_protein_FULL_SCALAR_v10`  
- Parity gate in `verify_cross.py` vs `vendor/fsot_compute.domain_scalar`

## Zero free parameters

No fitted weights. Inputs are pin seeds, DomainConfig rows, and trinary opcodes only.
