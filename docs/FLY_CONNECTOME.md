# Fly connectome — measured organism graph

**Data is not in this git repo.** It stays on the game drive: `D:\FlyWire_Connectome`.

Same pin `D1D38A`. Same law as Biohub and the protein product: **measured coordinates and measured edges are authority**. Residual scales the interface. We do not invent a 13 Å MDS brain or train a net to hallucinate synapses.

## Why this organism

Adult *Drosophila melanogaster* is the first animal with a complete brain-scale wiring diagram that still fits on a disk:

| Resource | What | Size class |
|----------|------|------------|
| FlyWire FAFB v783 | 139,255 neurons, ~50 M chemical synapses, 8,453 cell types | annotations TSV (small); proofread connections ~0.85 GB; full synapses 9.5 GB |
| Schlegel et al. 2024 | superclass, hemilineage, neurotransmitter, soma xyz, VFB/FBbt | GitHub `flyconnectome/flywire_annotations` |
| Male CNS v1.0 | **165,122 traced** (brain + VNC), 25.6 M edges | Berg et al. *Cell* 2026-09-03; on `D:\FlyWire_Connectome\male_cns` |
| BANC | female brain + cord, intact neck | Bates et al. *Nature* 2026; Codex |

Paper: Dorkenwald et al., *Nature* **634**, 124–138 (2024). Annotations: Schlegel et al., *Nature* **634**, 139–152 (2024). Portal: [codex.flywire.ai](https://codex.flywire.ai/) (sign-in). Open dumps: GitHub annotations + [Zenodo 10676866](https://zenodo.org/records/10676866).

Voxel of the EM volume: **4 × 4 × 40 nm** (anchor/soma columns in the annotation TSV).

## What we read first

`python scripts/fly_connectome.py --inventory`

First live read (2026-09-07), table on `D:\FlyWire_Connectome` (31.7 MB TSV, not git):

| Item | Number |
|------|-------:|
| Neurons in annotation dump | **139,248** |
| With cell type | **137,720** |
| With soma (x,y,z) | **118,104** |
| Superclass optic / central / sensory | 77,541 / 32,383 / 16,907 |
| Flow intrinsic / afferent / efferent | 118,497 / 19,262 / 1,489 |
| Top NT ACh / Glu / GABA | 86,193 / 24,875 / 19,171 |
| Sexually dimorphic + female-specific | 652 + 270 |
| *fru* / *dsx* / coexpress | 3,174 / 54 / 80 |
| Soma span | **817 × 369 × 278 µm** |

Source: `data/fly_connectome_inventory.json`. Voxel 4 × 4 × 40 nm. Paper count 139,255; dump is 7 rows short (non-neuronal / unreleased).

## Boot (measured cascade, not a mind)

`python scripts/fly_connectome.py --boot --seed sensory`

Proofread connections (~852 MB, Zenodo 10676866) on `D:\FlyWire_Connectome`. Edge weight = **measured synapse count**. GABA outgoing is inhibitory (predicted transmitter on the presynaptic cell). Hop count = leftover φ⁵. Each hop is rescaled to the observer max (half-max analog). Biochemistry residual scales the interface.

First boot (v630 public connections, 3.79 M edges, 127,979 neurons; v783 Zenodo 504’d):

`python scripts/fly_connectome.py --boot --seed sensory`

| Hop | Sensory mass | Central | Descending | Motor |
|----:|-------------:|--------:|-----------:|------:|
| 0 | 9,708 | 0 | 0 | 0 |
| 1 | 0.79 | 49.2 | **6.11** | **0.56** |
| 2 | 1.56 | 179 | **2.98** | **1.35** (motor peak) |
| 3–4 | leftover | olfactory LN hubs (lLN1 / lLN2) | falling | falling |

Seed is every `super_class=sensory` cell (includes photoreceptors R1–6). Hop 1–2 already put mass on **descending** and **motor** neurons — the measured path from sensors to effectors. Later hops collapse onto dense olfactory local interneurons (hub leftover). GABA edges are inhibitory. Source: `data/fly_connectome_boot.json`.

This is signal on the measured graph. It is **not** a trained RNN, not inner speech, and not a claim that the fly is solving a human puzzle. Flies navigate, court, fight, and walk. Those are the activities to match.

## Behavior video (measured observer, live)

`python scripts/fly_behavior.py`

Harvard Dataverse [doi:10.7910/DVN/BBNPYX](https://doi.org/10.7910/DVN/BBNPYX) tethered walking on a spherical treadmill. Files stay on `D:\FlyWire_Connectome\behavior` (not git). 3-D leg keypoints are the observer — not a trained pose net. cam-0 mp4 is corroboration only.

| Trial | Frames | s | Tarsus MAD+φ on | Paint | Bouts | Whole clip walking |
|-------|-------:|--:|----------------:|------:|------:|:-------------------|
| Fly01_T001 | 334 | 13.36 | 57 / 333 | 0.171 | 29 | yes (energy floor > median / φ⁵) |
| Fly01_T002 | 348 | 13.92 | 61 / 347 | 0.176 | 26 | yes |
| Fly02_T002 | 338 | 13.52 | 54 / 337 | 0.160 | 34 | yes |

Gate = median + φ·MAD on 6-leg tarsus speed. Tighten to φ² if paint > 1/φ (did not fire). Never loosen. These trials have no rest state — they are walking clips; the gate splits high vs low step energy.

cam-0 frame-diff Jaccard vs tarsus gate is ~0.12. One camera is not 6-leg 3-D. Keypoints stay authority.

Dryad Pratt freely-walking CSV and the Y-maze zip on `D:\` are auth stubs (56–92 bytes), not data.

### Residual boot on the walking program

Walking-on seeds **measured** brain afferents (mechanosensory class / Johnston’s organ), not photoreceptors. Olfactory is the contrast program (same graph, different seed). GPU sparse CSR hops on RTX 5070 (~1 s / 11 hops after the graph is in RAM). Residual Biochemistry 1.092. GABA inhibitory. 0 free parameters.

| Program | n seed | Hop 1 DN / desc | Hop 2 DN / desc | Hop 2 motor | Hop 1 peak |
|---------|-------:|----------------:|----------------:|------------:|------------|
| mechanosensory | 2,646 | **17.97** | **23.82 / 23.84** | **9.71** | **DNg15** (descending) |
| JO (`jo-` types) | 1,105 | 12.74 | 16.50 / 16.55 | 0.28 | CB0478 (central) |
| sensory (incl. R1–6) | 9,708 | 6.06 / 6.11 | 2.97 | 1.35 | ALLN (central hub) |
| olfactory | 2,276 | 0.003 | 0.245 | ~0 | ALLN → MBIN |

Mechanosensory / JO put mass on **descending neurons at hop 1–2**; olfactory does not (hop-2 DN mass ~100× smaller). That is the measured brain→cord walking command lighting up when the observer says the legs are moving.

Hop 2 of the mechanosensory seed peaks on a **motor** neuron. v630 is brain-only: those 100 motor cells are not VNC leg MNs. Do not claim a step-cycle CPG. DNs are the honest product.

`n_active` uses 1/φ of the global max — a hub can zero that count while class mass is still the right monitor.

Source: `data/fly_behavior_flow.json`.

### Male CNS — brain + nerve cord (live)

`python scripts/male_cns.py`

Traced neurons **165,122**, **25,563,197** edges, **22,055** GABA. RTX 5070. Walking seed is VNC sensory (leg/body afferents) and mechanosensory class. **vnc_motor** (708 cells) are the leg/body motor neurons.

| Program | n seed | Hop 2 vnc_motor | Hop 2 descending | Hop 1 peak |
|---------|-------:|----------------:|-----------------:|------------|
| vnc_sensory | 6,365 | **10.74** | 11.39 | IN05B011a (VNC intrinsic) |
| mechanosensory | 5,832 | **7.55** | 13.89 | IN01B001 (VNC intrinsic) |
| JO | 672 | **7.77** | **24.21** | **DNg29** (descending) |
| olfactory | 2,639 | 0.001 | 0.57 | il3LN6 (antennal lobe) |

VNC sensory and JO put mass on **vnc_motor**. Olfactory does not. JO hop 1 peaks on descending neuron DNg29, then the cord. Source: `data/male_cns_boot.json`.

Loop that ran: **video + 3-D tarsus → walking is on → seed mechanosensory / JO → residual hops → descending / DN mass**. Iron inaccuracies only against these kinematics. Predict *up* (other insects, then vertebrates) only where measured homologs exist — same product rule as protein Cα.

### Rest / odor observer (live)

`python scripts/fly_odor.py`

Harvard tethered-walk clips never rest. Álvarez-Salvado et al. *eLife* 7:e37815 Video 1 does: four walking flies, **ACV 10% pulse**, green overlay at the top of the frame. MPEG-4 on `D:\FlyWire_Connectome\behavior\odor` (not git). Dryad kinematics zip is 6.94 GB and was not downloaded. GitHub is LabVIEW / MATLAB, not CSV.

Observer (not a trained pose net): green-excess ≥ observer-max / φ = odor on; MAD+φ on chamber frame-diff = motion on; four lanes ≥ observer-max / φ of column mean.

| Epoch | Frames | Energy | Motion *y* (up is smaller) | Seed |
|-------|-------:|-------:|---------------------------:|------|
| rest (pre-odor, still) | 384 | 0.0046 | 356 | none |
| odor on (green overlay) | 302 / 10.07 s | **0.018** (4.0× rest) | **316** | olfactory |
| offset (odor off, moving) | 384 | **0.044** | 319 | none |
| JO / vnc_sensory (walking contrast) | — | — | — | hop-2 vnc_motor **7.77 / 10.74** |

Olfactory hop 1 peaks on **il3LN6** (antennal lobe). Hop-2 vnc_motor ≈ 0. Rest is not a cell class — odor off and motion below MAD+φ means no afferent seed. Offset search is not an invented OFF class. The pulse raises energy and shifts motion upwind; the walking command on this graph is still JO / vnc_sensory, not the olfactory cascade. Source: `data/fly_odor_flow.json`.

### Courtship observer (live)

`python scripts/fly_courtship.py`

Pan et al. *PLOS ONE* 2011 Movie S1: solitary male **UAS-dTrpA1 / fru-GAL4** at 29 °C (wing extension, abdomen bending, copulation attempts). MOV on `D:\FlyWire_Connectome\behavior\courtship` (not git). Leftover dark CC inside the arena (drop components ≥ max/φ) is the fly — not a trained pose net. MAD+φ on frame-diff energy and fly Rg.

| | |
|--|--|
| Frames | 2,877 / 96 s, tracked 2,848 |
| Motion MAD+φ on | 1,225 (paint 0.426) |
| Spread (Rg) on | 707 (paint 0.246) |
| Union courtship-on | 1,647 |

Seed is **measured** Male CNS genetics, not an invented P1 list: `fruDsx` labeled cells, type prefix `pC1_` (not optic LLPC1), `TN1` song motor neurons.

| Program | n seed | Hop 2 vnc_motor | Hop 2 descending | Hop 1 peak |
|---------|-------:|----------------:|-----------------:|------------|
| fru/dsx | 5,012 | 7.26 | **35.21** | TuTuA_2 (central; fru is broad) |
| pC1_ | 148 | 1.28 | **30.74** | SIP106m (brain, then cord) |
| TN1 | 35 | **10.07** | 1.22 | **hg3 MN** (VNC motor) |
| JO (walking) | 672 | **7.77** | 24.21 | DNg29 |
| olfactory | 2,639 | 0.001 | 0.57 | il3LN6 |

TN1 hits **hg3 MN** at hop 1 — song/flight motor on this dump. pC1 stays in the brain at hop 1 then loads descending. fru/dsx is 5,012 cells, not a command singleton. Olfactory still does not light vnc_motor. Source: `data/fly_courtship_flow.json`.

### Aggression observer (live)

`python scripts/fly_aggression.py`

Gao et al. *eLife* 13:RP104212. MPEG-4 on `D:\FlyWire_Connectome\behavior\aggression` (not git).

- **fig4 video 1** — *pC1SS2>CsChrimson*. Red “Light ON” overlay is the published stimulus (same cut as the odor green dot: observer-max / φ). 596 / 1,888 frames, **19.87 s** pulse. Four wells.
- **fig1 video 1** — Canton-S G14 males tussling. 1,487 live frames, MAD+φ motion 286 (paint 0.192).

Seed is **measured** Male CNS genetics: `pC1_` (paper: pC1SS2 promotes tussling), `dsx`, `male-specific` dimorphism. Same `pC1_` cells also seed courtship — not an invented fight class.

| Program | n seed | Hop 2 vnc_motor | Hop 2 descending | Hop 1 peak |
|---------|-------:|----------------:|-----------------:|------------|
| pC1_ | 148 | 1.28 | **30.74** | SIP106m (brain → cord) |
| dsx | 154 | 4.21 | **56.60** | oviIN |
| male-specific | 1,420 | 4.66 | **36.90** | SIP133m |
| TN1 (courtship song) | 35 | **10.07** | 1.22 | hg3 MN |
| JO (walking) | 672 | **7.77** | 24.21 | DNg29 |
| olfactory | 2,639 | 0.001 | 0.57 | il3LN6 |

pC1 / dsx / male-specific load **descending**, not song MNs. TN1 is the courtship-song contrast. Chamber translation energy *falls* during light-on (tussling is in-place grappling, not locomotion). Leftover two-fly distance did not beat the well rim — the overlay is the authority marker. Source: `data/fly_aggression_flow.json`.

### Sleep observer (live)

`python scripts/fly_sleep.py`

TriKinetics DAM2 IR beam counts (rethomics/damr monitor M064). File on `D:\FlyWire_Connectome\behavior\sleep` (not git). 30 live channels, 3,443 one-minute bins, 12:12 LD from the light column.

The field’s **5-minute** sleep cut is not used. Inactive = zero beam counts. Sleep = consecutive zeros longer than MAD+φ of length-weighted bout durations (leftover consolidated immobility).

| | |
|--|--|
| Sleep cut | **170 min** leftover |
| Fraction sleep | 0.172 |
| Night / day | **0.205 / 0.141** (1.45×) |

Seed is **measured** Male CNS types: `ER*` ellipsoid-body rings (includes ER5 / ER3m), `FB*` fan-shaped body, `LNv` clock neurons.

| Program | n seed | Hop 2 vnc_motor | Hop 2 descending | Hop 1 peak |
|---------|-------:|----------------:|-----------------:|------------|
| ER | 282 | 0 | ~0 | (recurrent; GABA can cancel outgoing) |
| FB | 602 | 0.07 | 1.51 | **hDeltaF** |
| LNv | 20 | 0.28 | 4.30 | SMP368 |
| JO (walking) | 672 | **7.77** | 24.21 | DNg29 |
| olfactory | 2,639 | 0.001 | 0.57 | il3LN6 |

Sleep seeds stay in the central complex and do **not** light the walking command. LNv (clock) leaks more toward descending than ER/FB. Source: `data/fly_sleep_flow.json`.

Other animals that are actually mapped: `docs/SPECIES_CONNECTOMES.md`. Only three whole-body synapse maps exist (*C. elegans*, *Ciona* larva, *Platynereis* larva). The worm whole-animal boot is live (`scripts/worm_connectome.py`): sensory → AVA → motor + body-wall muscle. Genetics join for fly walking proteins (`data/fly_walking_product.json`):

| Gene | n | Mode | Homolog |
|------|--:|------|---------|
| *Gad1* | 510 | **product Cα** | 2OKJ human GAD67, 64% id / 94% cov |
| *nan* | 833 | **product Cα** | 9NVN stink-bug Nan-Iav, 74% id / 86% cov |
| *iav* | 1123 | **product Cα** | 9NVP Inactive chain, 81% id / 60% cov (leftover tails) |
| *nompC* | 1619 | **product Cα** | 5VKQ Drosophila NOMPC, 97% id / 92% cov |

First pass missed *iav* and *nompC*: query-coverage 0.65 vetoed a fully-used Inactive chain (cov_t 0.99, cov_q 0.59), and a globular Rg window vetoed elongated 5VKQ. Leftover floor 1/φ² + close-homolog (id ≥ 1/φ) keeps the measured map. Front door now uses product identity cap 1.0 (not the 0.95 H2H handicap).

PDBs on `D:\FlyWire_Connectome\male_cns\product`.

Cell → gene → product on both live graphs (`data/organism_product_join.json`):

| Gene | Organism | Sits on | Template |
|------|----------|---------|----------|
| ChAT | fly | cholinergic DN / motor (ACh edges) | 2FY4 53% / 78% |
| unc-25 | worm | GABA motor DD/VD/RME | 7LZ6 58% / 94% |
| mec-4 | worm | ALML/ALMR (hop-0 peak) | 6VTL 48% / 56% leftover |
| myo-3 | worm | body-wall muscle (hop-4) | 6XE9 61% / 58% leftover |
| VGlut | fly | glutamatergic vnc_motor NMJ | 7T3O 61% / 69% |
| Mhc | fly | muscle (not in CNS dump) | 5W1A 97% / 40% leftover |

Male CNS already carries measured genetics on cells: **fru/dsx** on 5,012 traced neurons, **receptorType** on 752 (`data/male_cns_genetics_on_cells.json`). Fly Cell Atlas has no published `bodyId` join — do not invent one.

Neuron table only (no 9.5 GB synapse dump until the atlas is live):

- `root_id`, soma and backbone `(x,y,z)` nm-scale voxels
- `flow` / `super_class` / `cell_class` / `cell_type`
- `top_nt` (predicted transmitter)
- `fru_dsx`, `dimorphism` (sex-circuit observers)
- Virtual Fly Brain / FBbt IDs → later FlyBase → UniProt → **FSOT product Cα**

That last arrow is the genetics join: the same product freeze (0.13 Å class) on fly proteins that sit on named cells in the connectome.

## How this sits next to Biohub / protein

```text
Drosophila sequence
    → FSOT product Cα   (measured homologs, 0.13 Å class)
    → molecular 3-D

FlyWire neuron (root_id, soma xyz, type, nt)
    → measured cell in the brain graph
    → organism 3-D  (Biochemistry residual on measured synapses)

Biohub / Zebrahub          back burner — docs/BIOHUB_FREEZE.md
```

## Anti-goals

- Do not embed the 139k graph with MDS and call it a fold.
- Do not train a contact net on FlyWire.
- Do not copy 100 TB of EM voxels. We read tables and proofread connections.
- Do not call leftover hub collapse (lLN1 / MBIN) a thought.
- Do not mix v783 annotations onto v630 root IDs.
- Do not download the 9.5 GB raw synapse dump for this loop.
- Male CNS syn-points (12.7 GB) are not required for the weight-graph boot.
