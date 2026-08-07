# All-atom MD: lab tool, not cognitive runtime

## Status

**Implemented** as a host-side classical MD lab (`src/allatom_md.zig` in this product; monorepo path `embodiment/zig/src/allatom_md.zig`):

| Piece | What it is |
|-------|------------|
| Integrator | Velocity-Verlet, Δt ~ 0.5–1 fs |
| Bonded | Harmonic bonds + angles |
| Nonbonded | LJ 12-6 + Coulomb, cutoff, cubic PBC (minimum image) |
| Thermostat | Berendsen weak coupling |
| Demo systems | TIP3P-like water box (~27 H₂O); K⁺ + carbonyl “filter” |
| Entry | `fsot_mind md` / `allatom` / `all-atom` |
| Numeric path | **Host f64** — never on the Fixed cognitive lattice |

Gate markers: `FSOT_ALLATOM_MD PASS` · `FSOT_MD_LAB_OK`.

## What “all-atom molecular dynamics” is

**All-atom MD** simulates **every atom** (protein, lipid, water, ion) with classical force fields, timesteps of **~1–2 femtoseconds** (10⁻¹⁵ s), for **nanoseconds to microseconds** of wall-clock physics.

Typical use:

- How one NMDA pore or CaMKII domain moves at atomic scale  
- Drug docking, membrane thickness, hydrogen bonds  

A **single synapse** with all atoms of receptors + spine + water is already **huge**. A **whole brain / whole mind** at all-atom MD is **not** how any biological AI or neural sim of learning works at scale. Supercomputers run MD for *one protein* for *nanoseconds*, not “think about plants and math.”

## What the mind uses (correct biological *process* scale)

| Scale | Biology | This Zig mind |
|-------|---------|----------------|
| Atoms (fs) | MD of one protein / water / ion filter | **`allatom_md`** — **lab only**, not mind step |
| Single channel (µs–ms) | Markov open/close | **`channel_stoch_fixed`** 50 µs, 48 AMPA / 16 NMDA |
| Spine chemistry (ms) | Ca, kinases, traffic | **`molecular_fixed`** ODEs + consolidate W |
| Glia (ms–s) | Uptake, prune, myelo | **`glia_fixed`** |
| Plasticity | STDP / LTP / LTD | **`stdp_fixed`** × eligibility |
| Genetics | Expression → synapse law | **codon → genotype → fsotPairWeight** |
| Learning / curriculum | Memory, concepts | bank + depth + pathway thought |

**Credibility for a genetic FSOT mind** comes from **lawful process at channel → synapse → network → concept**. All-atom MD is the offline structural biophysics tool when we need atomic motion; it does **not** replace the wet stack or curriculum.

## When all-atom MD matters here

- Offline lab: energy drift / force sanity, ion–carbonyl geometry probes  
- Future: calibrate a rate constant or pore geometry against a short trajectory  
- Publish / audit: show we can run real all-atom classical MD in-repo  

**Not** as the online “thinking loop.” Using MD as the mind would be a different project (and still wouldn’t give grade-school understanding).

## Doctrine

We **implement** all-atom MD because structural biophysics may require it and we refuse the “we don’t need it so we never build it” shortcut.  
We **do not** route cognition through femtosecond atoms.  
We **did** refuse the shortcut *up* into fake 4-field tags for wet chemistry.  
**Channel-scale stochastic wet biophysics + FSOT Fixed + curriculum** remains the correct stack for learning and depth; **`allatom_md` sits beside it as lab.**
