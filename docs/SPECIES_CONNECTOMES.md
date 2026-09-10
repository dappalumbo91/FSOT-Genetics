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
| BANC | female | brain **+** VNC, intact neck | ~142k | Codex; Bates et al. *Nature* 2026 |
| Male CNS v1.0 | male | brain **+** VNC | **166,700** | neuPrint `male-cns:v1.0`; Berg et al. *Cell* 2026-09-03 |
| MANC / FANC | male / female | VNC only | ~16k / sparse | neuPrint / GitHub |
| Larval CNS | — | complete first-instar CNS | ~3,000 | Winding et al. 2023 |

What we have on disk: **female brain** (v630, 128k neurons, 3.79 M edges). Walking prediction is descending/DN mass. Leg motor neurons live in the VNC — that is why Male CNS / BANC is the next measured graph, not a trained RNN.

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

Sensory in → AVA command interneuron → motor + muscle. Same residual law as the fly. Source: `data/worm_connectome_boot.json`.

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
measured behavior        → observer seed (fly walking live)
```

Missing for a *Drosophila* whole-animal molecular recreation (need these, do not invent them):

1. **VNC** — Male CNS / BANC (published; not yet on `D:`).
2. **Muscles, gut, cuticle** — not in any fly EM CNS dump.
3. **Per-cell transcriptome joined to `root_id`** — Fly Cell Atlas exists; it is not yet wired to FlyWire IDs here.
4. **Named proteins on named cells** — FlyBase / UniProt join, then product Cα only where a measured homolog exists.
5. **Rest / odor / courtship observers** — walking is one program; we do not have a rest trial in the Harvard clip.
6. **Other species** — worm is the next whole-animal boot; Platynereis and Ciona after that; vertebrates only as homologs.

Genetics join for the walking program (measured, not guessed):

| Gene | FlyBase | UniProt | Sits on |
|------|---------|---------|---------|
| *iav* (TRPV) | FBgn0032043 | Q9W3W0 | Johnston’s organ |
| *nan* (TRPV) | FBgn0036414 | Q9VUD5 | Johnston’s organ |
| *nompC* (TRPN) | FBgn0026324 | Q7KIQ2 | mechanosensory transduction |
| *Gad1* | FBgn0004516 | P20228 | GABA synthesis (inhibitory residual) |

Fold those with the **existing protein product** (measured homologs). Do not MDS a 13 Å “fly protein brain.”

## Anti-goals

- Do not call a genome a connectome.
- Do not mix v783 annotations onto v630 root IDs.
- Do not claim a mouse or human whole-brain synapse map exists.
- Do not skip the VNC and call the fly brain a whole animal.
- Predict *up* only where measured homologs exist — same rule as Cα.
