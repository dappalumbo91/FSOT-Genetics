# Fly pack + plant panel (2026-09-17)

Fly experiment folder: `C:\Users\damia\Desktop\fsot fly nuron net`

Hemibrain hops recorded: `data/hemibrain_connectome_boot.json` (JO hop-2 descending 2.68 Giant Fiber vs olfactory 0.072 AL LN). Script `scripts/hemibrain_connectome.py`. Cache on `D:\FlyWire_Connectome\hemibrain` (not git).

## Plant findings to pick up (HEAD started at 61872b4)

| Seed | hop-1 top | hop-2 top | note |
|------|-----------|-----------|------|
| photoreceptor | PIF3 (`light_tf=1.00`) | PHYB | UVR8 missing from IntAct graph |
| ABA | ABI1 | PYL9 | SLAC1 ~0.001 |
| auxin | ABCB19 | TIR1 | |
| calvin / PSII | CML9 | psbA | stays on photosystem |

Signaling product: 7 folded, **PHOT1 Q2V2M9 `no_measured_map`**. PIF3 identity 0.69 ≥ 1/φ.

Next here: re-stamp gauntlet (`python verification/run_cross_proof.py`). Do not invent plant synapses or plasmodesmata.
