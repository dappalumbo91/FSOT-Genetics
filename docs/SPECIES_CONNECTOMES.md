# Which animals are actually mapped

Pin `D1D38A`. Same product rule as protein Cα and the fly boot: **measured cells and measured edges are authority**. Residual does not invent synapses. A genome is not a connectome. An fMRI “connectome” is not this object.

## Short answer

**Three species have a whole-body synaptic wiring diagram** (every reconstructed cell, including effectors):

| Species | What | Cells / neurons | Year |
|---------|------|----------------:|------|
| *Caenorhabditis elegans* | Whole adult (hermaphrodite + male), neurons **and muscles** | ~1,000 somatic cells; 302 / 385 neurons | 1986, Cook 2019 |
| *Ciona intestinalis* tadpole larva | Whole larva | 301 cells, 177 neurons | Ryan 2016 |
| *Platynereis dumerilii* 3-day larva | Whole segmented larva | 9,162 cells, ~966 neurons, 294 types | Verasztó 2025 |

**One species has a complete complex CNS** (brain + nerve cord, 10⁸ synapses): *Drosophila melanogaster*. Nothing else is remotely that close.

Mouse cortex cubes, zebrafish larval volumes, and human tractography are **not** in this class.

## Drosophila — best mapped *brain*, not yet a whole animal

| Dataset | Sex | What | Neurons | Open dump |
|---------|-----|------|--------:|-----------|
| FlyWire FAFB v783 | female | whole brain | 139,255 | GitHub annotations + Zenodo connections (v630 public used here) |
| BANC | female | brain **+** VNC, intact neck | **175,401** neurons (glia dropped) | Bates et al. *Nature* 2026; live `scripts/banc_connectome.py` |
| Male CNS v1.0 | male | brain **+** VNC | **166,700** | neuPrint `male-cns:v1.0`; Berg et al. *Cell* 2026-09-03 |
| MANC / FANC | male / female | VNC only | ~16k / sparse | neuPrint / GitHub |
| Larval CNS | — | complete first-instar brain | **2,952** | Winding et al. 2023; live `scripts/larva_connectome.py` |

What we have on disk: female brain (v630), **BANC v888** (175,401 neurons, 13.5 M edges, intact neck), Male CNS v1.0 (165,122 traced, brain+VNC, 25.6 M edges), **and the first-instar larval brain** (2,952 cells, 110,677 edges). Walking seed hits **vnc_motor** on both sexes. Larval mechanosensory seed hits **DN-VNC**. Cell → gene → product Cα: `data/organism_product_join.json`. Homologs of those proteins in bee / mosquito / beetle: `data/homolog_correspondence.json` (protein only — not a connectome).

Genome / genetics: FlyBase complete; ~14k protein-coding genes; UniProt proteome `UP000000803`. **304 Drosophilidae genomes** are annotated (Zenodo 2025) — genomes, not brains.

## Whole-body maps (the only path to “full organism”)

*C. elegans* is the only animal where **every cell is named, the lineage is complete, the genome is complete, and neurons synapse onto named muscles**. That is closer to a molecular recreation than the fly brain dump, at 300 neurons instead of 140k.

Cook 2019 hermaphrodite chemical graph (on `D:\FlyWire_Connectome\C_elegans`, not git):

| Class | n |
|-------|--:|
| Sensory neurons | 83 |
| Interneurons | 81 |
| Motor neurons | 108 |
| Body-wall muscles | 95 |
| Pharynx | 50 |
| Other end organs | 21 |
| Sex-specific | 16 |
| **Total nodes** | **454** |
| Chemical edges | 4,879 |

Boot: `python scripts/worm_connectome.py --boot --seed sensory`

First live boot (Cook 2019 SI5, 453 cells, 4,879 chemical edges, **956 NMJs** onto 95 body-wall muscles; Netzschleuder CSV had listed the muscle nodes and dropped every NMJ — ironed against SI5):

| Hop | Sensory | Interneuron | Motor | Body-wall muscle | Peak cell |
|----:|--------:|------------:|------:|-----------------:|-----------|
| 0 | 83 | 0 | 0 | 0 | ALML (touch) |
| 1 | — | **15.80** | 5.12 | **1.02** | **AVA** (backward command) |
| 2 | — | 11.77 | 9.88 | 2.32 | AVA |
| 3 | — | 12.22 | **16.47** | 6.91 | **RMD** (head motor) |
| 4 | — | 9.12 | 14.80 | **8.13** | RMD |

Sensory in → AVA command interneuron → motor + muscle. Same residual law as the fly. Source: `data/worm_connectome_boot.json`. GABA count is **26** (Pereira DD1–6 / VD1–13 match SI5 `DD01`/`VD01`; earlier boot missed the zero-pad and counted 11).

### C. elegans male — live

`python scripts/worm_connectome.py --boot --sex both --seed sensory`

Cook 2019 SI5 sheet `male chemical`. 575 cells, 5,306 chemical edges, **913 NMJs** onto the same 95 body-wall muscles. SI5 lumped `dBWM`/`vBWM` under MOTOR NEURONS (no group header) — relabeled by the measured Cook names, not a guessed type. Extra **138 sex-specific cells** (rays, CEM, hook, spicules, sex muscles).

| Hop | Herm motor / muscle / sex | Male motor / muscle / sex | Male peak |
|----:|--------------------------:|--------------------------:|-----------|
| 1 | 5.12 / 1.02 / 0.07 | 6.95 / 1.53 / **7.74** | **AIAR** (PVX sex-specific #2; AVA still high) |
| 2 | 9.88 / 2.30 / 0.07 | 10.16 / 2.23 / 6.85 | RIAR |
| 3 | 16.46 / 6.86 / 0.05 | 10.06 / 2.94 / 5.65 | RMDDR |
| 4 | 14.65 / 7.44 / 0.03 | 8.06 / 3.70 / 5.68 | RMDDR |

Same sensory seed (n=83). Male sensory hops dump mass onto the sex circuit; herm does not. Body-wall muscle still lights. Sex-specific seed (n=138) stays in sex cells (hop-1 sex 17.35, motor 1.43) — sex program ≠ locomotion, same split as fly courtship vs walking. Source: `data/worm_male_connectome_boot.json`, `data/worm_sex_compare.json`.

### Ciona (chordate sibling) — live

`python scripts/ciona_connectome.py --boot --seed sensory`

Ryan et al. 2016. 205 nodes, 2,903 edges, contact-depth weights. **NT sign is not annotated** — unsigned residual, no invented GABA.

| Hop | Motor | Muscle | Peak |
|----:|------:|-------:|------|
| 1 | 0 | 0 | Em2 |
| 2 | **1.40** | 0.32 | **MGIN1L** (motor-ganglion interneuron) |
| 3–4 | 2.43–2.49 | 0.33 | MGIN1R |
| 11 | 3.40 | 0.55 | MGIN1R |

Photoreceptor/palp seed → MGIN command → motor / tail muscle. Same residual as fly and worm. Source: `data/ciona_connectome_boot.json`.

### Platynereis (segmented annelid) — live

`python scripts/platynereis_connectome.py --boot --seed sensory`

Verasztó et al. *eLife* 2025. 1,720 cells, 8,451 chemical edges. Unsigned (transmitter mostly unannotated). Effectors include **ciliary bands** (swimming), not only muscle.

| Hop | Motor | Muscle | Cilia | Peak |
|----:|------:|-------:|------:|------|
| 0 | 0 | 0 | 0 | PRC (photoreceptor) |
| 1 | 3.05 | 0.03 | 0 | **IN1** (visual interneuron) |
| 2 | **7.77** | 1.88 | 8.45 | IN1 |
| 4 | 4.05 | 3.41 | **19.09** | **prototroch** (ciliary band) |
| 11 | 2.71 | 1.90 | 18.85 | prototroch |

Sensory → IN1 → motor / muscle / cilia. Same residual. Source: `data/platynereis_connectome_boot.json`. All three whole-body maps now boot.

## Not in this class (do not treat as fly-equivalents)

| Resource | Why it is not a fly-class map |
|----------|-------------------------------|
| Mouse MICrONS / cubic-mm cortex | A cube of cortex, not a mouse |
| Zebrahub / larval zebrafish EM | Partial; Biohub tracking is back-burner (`docs/BIOHUB_FREEZE.md`) |
| Human / macaque “connectomes” | MRI tracts or sparse tracing, not every synapse |
| Honey bee, ant, mosquito, locust | Excellent genomes and behavior; **no** complete EM connectome |
| *Mnemiopsis leidyi* | Aboral-organ / nerve-net piece (~1,000 cells), not the whole ctenophore |
| 304 Drosophilidae genomes | Sequence only |

Thousands of species have a genome. A handful have a cell atlas. **Three** have a whole-body synapse map. **One** has a complete complex brain.

## Full recreation under FSOT (honest stack)

```text
measured genome          → trinary codon / AA map (already live)
measured protein homolog → FSOT product Cα (0.13 Å freeze)
measured cell + synapses → residual hops on the graph (fly brain live; worm whole-animal live)
measured behavior        → observer seed (fly walking + rest/odor + courtship + aggression + sleep live)
```

Missing for a *Drosophila* whole-animal molecular recreation (need these, do not invent them):

1. **VNC** — Male CNS v1.0 **and BANC v888** both boot (female brain+cord, intact neck).
2. **Muscles, gut, cuticle** — not in any fly EM CNS dump.
3. **Per-cell transcriptome joined to `root_id`** — Fly Cell Atlas exists; it is not yet wired to FlyWire IDs here.
4. **Named proteins on named cells** — FlyBase / UniProt join, then product Cα only where a measured homolog exists.
5. **Fly observers** — walking, rest/odor, courtship, aggression, and sleep are live (`scripts/fly_behavior.py`, `scripts/fly_odor.py`, `scripts/fly_courtship.py`, `scripts/fly_aggression.py`, `scripts/fly_sleep.py`).
6. **Other species** — worm (both sexes), Ciona, and Platynereis whole-body maps boot. Bee / mosquito / beetle: protein homologs only.

Genetics join for the walking program (measured, not guessed):

| Gene | FlyBase | UniProt | Sits on |
|------|---------|---------|---------|
| *iav* (TRPV) | FBgn0032043 | Q9W3W0 | Johnston’s organ |
| *nan* (TRPV) | FBgn0036414 | Q9VUD5 | Johnston’s organ |
| *nompC* (TRPN) | FBgn0026324 | Q7KIQ2 | mechanosensory transduction |
| *Gad1* | FBgn0004516 | P20228 | GABA synthesis (inhibitory residual) |

Fold those with the **existing protein product** (measured homologs). Do not MDS a 13 Å “fly protein brain.”

## BANC female brain+cord — live

`python scripts/banc_connectome.py`

Bates et al. *Nature* 2026. Files on `D:\FlyWire_Connectome\banc` (GCS `compiled_data/banc_888/`, Dataverse 10.7910/DVN/7WTH1N). 175,401 neurons after dropping glia/trachea, 13.5 M edges, 21,300 GABA (predicted transmitter). `vnc_motor` = measured `super_class=motor` ∩ `region=ventral_nerve_cord`.

Hop-2 `vnc_motor` (same split as Male CNS):

| Seed | n | hop-2 vnc_motor | hop-1 peak |
|------|--:|----------------:|------------|
| vnc_sensory | 7,845 | **10.16** | AN05B009 (ascending) |
| chordotonal | 2,136 | **5.01** | INXXX007 (VNC intrinsic) |
| JO (`type_prefix jo`) | 1,198 | **6.64** (hop-3 **17.35**) | **DNg29** (same DN as male walking) |
| olfactory ORN | 3,007 | **0.0008** | **il3LN6** (AL local hub) |

Female brain+cord independently reproduces the rule: which sensory class is on decides whether cord motor lights. Olfactory stays in the antennal lobe. Source: `data/banc_connectome_boot.json`.

## What predicts *up* (and what does not)

The experiment — map another insect from genetics without that species’ EM — is **product Cα of measured homologs**, not an invented graph.

| Predicts up | Does not |
|-------------|----------|
| UniProt OrthoDB or UniRef50 member of a residual-mass protein (Gad1, nompC, iav, nan, ChAT, VGlut, Mhc, unc-25, mec-4, myo-3) | A bee / mosquito / beetle connectome |
| Same residual law on a **measured** graph (BANC, Male CNS, larva, worm both sexes, Ciona, Platynereis) | Synapses inferred from a genome or cell atlas |
| Sensory-class split (mechanosensory/JO → motor; olfactory local) as a **hypothesis to test** when a map exists | Seeding an unmeasured mosquito brain with fly types |
| GABA inhibitory **only** where the dump annotates GABA | Invented transmitter signs |

`python scripts/homolog_correspondence.py`

**23 measured homologs, 7 misses.** Bee and beetle have OrthoDB members for most residual-mass genes. Anopheles often lacks the OrthoDB xref; UniRef50 still hits Gad1 / iav / ChAT / VGlut and **misses** nompC / nan / Mhc / mec-4. *iav* and *nan* collapse to the same bee TRPV (paralogs in one OrthoDB group). Bee nompC folds on 5VKQ at **81%** identity (same template as fly nompC). *mec-4* has an OrthoDB member in bee/beetle but **no measured structure map** — product stays `no_measured_map`, no invented 3-D. That is correspondence, not a wiring diagram. Source: `data/homolog_correspondence.json`. Sequences on `D:\FlyWire_Connectome\homologs` (not git).

## Anti-goals

- Do not call a genome a connectome.
- Do not mix v783 annotations onto v630 root IDs.
- Do not claim a mouse or human whole-brain synapse map exists.
- Do not skip the VNC and call the fly brain a whole animal.
- Predict *up* only where measured homologs exist — same rule as Cα.
- Do not invent a honey-bee or mosquito synapse graph from fly hops.
