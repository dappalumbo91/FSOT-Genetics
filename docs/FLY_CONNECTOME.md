# Fly connectome — measured organism graph

**Data is not in this git repo.** It stays on the game drive: `D:\FlyWire_Connectome`.

Same pin `D1D38A`. Same law as Biohub and the protein product: **measured coordinates and measured edges are authority**. Residual scales the interface. We do not invent a 13 Å MDS brain or train a net to hallucinate synapses.

## Why this organism

Adult *Drosophila melanogaster* is the first animal with a complete brain-scale wiring diagram that still fits on a disk:

| Resource | What | Size class |
|----------|------|------------|
| FlyWire FAFB v783 | 139,255 neurons, ~50 M chemical synapses, 8,453 cell types | annotations TSV (small); proofread connections ~0.85 GB; full synapses 9.5 GB |
| Schlegel et al. 2024 | superclass, hemilineage, neurotransmitter, soma xyz, VFB/FBbt | GitHub `flyconnectome/flywire_annotations` |
| Male CNS v1.0 | 166,700 neurons, brain + nerve cord | later (neuPrint / Codex) |
| BANC v888 | female brain + cord | later |

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

Next table (when we need edges): Zenodo `proofread_connections_783.feather` (~852 MB) on the same `D:\` tree. Not the 9.5 GB raw synapse list until the atlas+graph is residual-scored.

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
- Do not copy 100 TB of EM voxels. We read tables and, later, proofread connections.
