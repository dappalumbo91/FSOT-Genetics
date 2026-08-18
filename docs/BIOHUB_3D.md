# Biohub / Zebrahub 3-D — measured organism geometry

**Data is not in this git repo.** It stays on the game drive (and the FSOT archive).  
This repo **reads** it.

| Asset | Where | What |
|-------|--------|------|
| Kaggle Biohub train/test | `D:\Kaggle_Biohub_Data` (~85 GB) | 3D+time light-sheet volumes (OME-Zarr) + sparse GEFF tracks |
| Existing competition code | `C:\Users\damia\biohub-fsot-unet` | U-Net + FSOT linker (junction → `D:\`) |
| Zebrahub public tracks | `I:\FSOT-Physical-Archive\05_Zebrahub-Development` | ~46 M cell detections, 5 DaXi embryos |
| Lean scalar panel | FSOT-2.1-Lean `zebrafish_cell_tracking_panel` | 0.022% residual — **not** a 3-D RMSD |

Competition: [Biohub – Cell Tracking During Development](https://www.kaggle.com/competitions/biohub-cell-tracking-during-development).

## What 3-D we can read now

`python scripts/biohub_3d.py`  
`python scripts/biohub_3d.py --inventory`

First live read (2026-08-17):

| Item | Number |
|------|-------:|
| Train GEFFs readable | **199 / 199** |
| Annotated nodes (min / med / max) | 50 / 659 / 1950 |
| Proxy `44b6_0113de3b` | 52 nodes, 50 edges, 0 divisions, t = 0–75 |
| Volume shape | **100 × 64 × 256 × 256** uint16 (T, Z, Y, X) |
| Voxel | 1.625 / 0.40625 / 0.40625 µm |
| Proxy spatial span | 101 × 63 × 73 µm |
| Median parent→child step | **2.88 µm** (measured) |

The volume header is read **without** loading 85 GB of voxels. Cell centers are the 3-D object we map.

## Voxels (what we were not reading)

`python scripts/biohub_3d.py --voxels`

Light-sheet pixels are now sliced in place (one T is ~8 MB). Nuclei are **brightness peaks** in the volume — measured observer, no trained U-Net in this repo.

Proxy `44b6_0113de3b` (100 frames):

| Item | Number |
|------|-------:|
| Frames loaded | **100 / 100** |
| GT cell vs background intensity | **1687 / 155** (~11×) |
| Dense peaks | **8,000** (estimated true cells 25,755) |
| Official 7 µm recall vs sparse GT | 27/52 (**0.52**) |
| 12 µm recall (one nucleus) | 51/52 (**0.98**) |
| Median match when inside 7 µm | **0.57 µm** |

The annotated GEFF point is an *approximate* cell center. We read the brightness maximum of the nucleus. Those two points sit on the same cell; 7 µm is tight for that difference, 12 µm is a nucleus. Source: `data/biohub_3d_voxels.json`.

**Kaggle submit** (U-Net + ILP, public ~0.848 floor) stays in `C:\Users\damia\biohub-fsot-unet`. This repo is the genetics / 3-D *reader*, not a second competition stack.

## How this sits next to the protein product

```text
Danio rerio sequence
    → FSOT product Cα   (measured homologs, 0.13 Å class)
    → molecular 3-D

Biohub / Zebrahub GEFF
    → measured cell (t, z, y, x) µm
    → organism 3-D  (Biochemistry residual on parent→child steps)
```

Same pin `D1D38A`. Same law: **measured coordinates are authority**. Residual scales the interface; it does not invent cell positions or a 13 Å MDS fold.

Sparse GT tracks are a subset of the true cells (`estimated_number_of_nodes` ~ 25k–30k per video). Evaluation in the Kaggle project is edge Jaccard on that subset. Here we only **read** and residual-score the measured graph.

## Next (when we need more)

1. Stream one time-point from the OME-Zarr (pixel 3-D) for a viewer — still no copy of the dump.  
2. Join Zebrahub gene-expression tracks on `I:\` to Danio UniProt → product Cα on the same embryo.  
3. Leave U-Net training in `biohub-fsot-unet`; this repo stays the genetics / 3-D *reader*.
