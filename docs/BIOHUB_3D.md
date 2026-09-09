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

Refinement (FSOT, 0 free params):

- Center = **half-max first moment** of the observed blob (not the brightest voxel).
- Gate = median + φ·MAD (φ²·MAD if that paints > 1/φ of voxels).
- NMS = φ³ µm (φ⁴ merged an annotated cell with an unannotated neighbor).
- Residual second collapse on leftover brightness (7–12 µm ghosts).
- Lineage = Hungarian on the first collapse (intensity identity). Unmatched primaries get a second pass at φ⁵ µm. Isolated residual peaks (farther than NMS from every primary) may meet an unmatched primary; leftover residual–residual tracks stay off the primary map. Halo residual↔primary is not mixed (that stole tracks). Leftover-length dests (> φ⁴) yield to a closer first-pass dest when 2-opt cannot swap (return > leftover).
- Outcome = the parent's **predicted child** lands within 7 µm of the measured next cell. Pairing both GT ends independently was matching a closer ghost that was not the continuation (pair-match still in `link_meta`). Kaggle scores **edge Jaccard**, not follow.

| Video | GT | 7 µm find | Follow 7 / 12 µm | Edge Jaccard (adj.) | Product nodes |
|-------|---:|----------:|-----------------:|--------------------:|--------------:|
| `44b6_0113de3b` (sparse) | 52 | **1.00** | **1.00 / 1.00** | **1.00 (1.00)** | 26,233 |
| `6bba_09961292` (dense) | 1950 | **0.96** | **0.84 / 0.94** | **0.82 (0.83)** | 28,164 |

Find-recall uses all peaks (observer). Product graph: first-collapse + isolated residual; in-shell residual folds into the primary centroid. Dest fill + 2-opt at leftover φ⁵ µm; leftover-length dests then yield to a closer first-pass dest when the 2-opt return exceeds leftover. Dense: TP 1538 / FP 6 / FN 333. Source: `data/biohub_3d_voxels.json`.

The 7 µm official radius is peak/centroid vs annotator center. At one nucleus (12 µm) we recover almost every annotated cell and most of its next frame.

Kaggle is the **reference metric** (7 µm Hungarian, adjusted edge Jaccard). A U-Net notebook is optional and only after it actually runs; v64 v1 died after predict on `import tracksdata` in the notebook kernel. This repo stays the genetics / 3-D reader.

## Disconnect (do not paper over)

The competition detector is almost nothing but free parameters (a trained U-Net). This reader is 0. Those are not the same object.

FSOT can score the **measured** field: voxels, blob centroids, parent→child steps, residual \(r = 1+|S|·P_{\mathrm{NEW}}\). It does not train an appearance model, and it must not invent a nucleus the brightness field did not show. Where a net “knows” a center from labels, this observer only has the half-max first moment of leftover light.

That gap is real. It is not a missing weight, and wrapping the U-Net does not close it inside the law. If the leftover 7–8 µm shell can move inside the official 7 µm ball, it is because the **blob mathematics** still under-reads the nucleus — not because we failed to fit intelligence. Residual scales the interface; it does not become the net.

### Eye / relay (image as measurement)

The math already derives components from **measurements**. It has not been the camera. A trained net can be the eye — cones: photons → an activation field — without being the brain.

```text
light-sheet pixels          measured scene (already in OME-Zarr)
        ↓
ML appearance field         optional apparatus (isolates what MAD+φ does not)
        ↓
FSOT relay                  half-max centroid, φ³ NMS, leftover residual, identity
        ↓
this cell / this protein    product (lineage now; Zebrahub gene → Cα later)
```

The image **is** the measurement. FSOT translates that field into structure. The net may paint leftover brightness the native gate missed; it must not emit the annotator coordinate. If the net outputs xyz, the relay is a pass-through and we are wrapping free parameters again. Claim path stays 0 free params (`docs/DESIGN.md`). U-Net stays lab apparatus, out of the claim, same as a microscope vendor’s ISP.

First local relay (FT detector **field** × native photons, then the same MAD+φ / φ³ NMS / leftover-yield observer):

| Observer | Video | Find 7 | Jaccard (adj.) | TP / FP / FN | Nodes |
|----------|-------|-------:|----------------|--------------|------:|
| Native photons (claim) | dense | 0.96 | **0.82 (0.83)** | 1538 / 6 / 333 | 28,164 |
| Eye re-detect whole volume | dense | 0.97 | 0.84 (0.81) | 1583 / 3 / 288 | 43,642 |
| Native + isolated eye fill (photon rank) | dense | 0.97 | 0.83 (0.83) | 1565 / 5 / 306 | 31,117 |
| Native + isolated eye fill (sigmoid ≥ 1/φ) | dense | 0.97 | 0.84 (0.84) | 1567 / 5 / 304 | 31,050 |
| Native + correspondence-trained eye fill | dense | **0.97** | **0.848 (0.848)** | **1589 / 2 / 282** | **31,117** |
| Native photons (claim) | proxy | 1.00 | **1.00 (1.00)** | 50 / 0 / 0 | 26,233 |
| Eye re-detect / in-shell peak fold | proxy | 1.00 | 0.94 (0.85–0.94) | 47 / 0 / 3 | extra blobs |
| Native + isolated eye fill | proxy | **1.00** | **1.00 (1.00)** | **50 / 0 / 0** | 26,233 |

What failed: re-detecting through the eye (44k nodes, proxy steal); folding discrete in-shell eye peaks into native xyz (proxy 1.00→0.94); re-centroid of the same native blob on photons×eye (dense Jaccard 0.82→0.78).

What held: **native product stays the graph**. Isolated eye leftover farther than NMS, only up to `estimated_number_of_nodes`.

The FT net already knew the leftover: isolated ghosts sit at sigmoid **0.002**; the 24 native-miss GT fills sit at **0.999**. Ranking fill by photons put the silenced field back — that was the regression. Gate and rank by the trained field (≥ 1/φ) : dense Jaccard 0.835. Then fine-tune the eye on other train videos (hold out dense + proxy) with a steal-shell annulus around each labeled cell (NMS..φ⁵, target 0). Correspondence fill: dense **0.848 / 0.848** (TP 1589 / FN 282) at \(T_{\mathrm{true}}\). Proxy is already over the estimate, so it still receives no extra blobs and stays **1.00**. That 0.848 is this train video’s adj edge Jaccard, not public LB Final. Claim path is still photons. Apparatus: `scripts/_eye_relay.py`, `scripts/_eye_correspondence_train.py`.

### Variety panel (existing 199 train GEFFs on `D:\` — no extra download)

Caches stay on the game drive (`_fsot_detect_cache`, `_fsot_eye_cache`). Test is already 4 volumes (full submit set). Native find **0.96 is not universal**.

| Video | GT / \(T_{\mathrm{true}}\) | Native find / Jac (adj.) | Corr-fill find / Jac (adj.) |
|-------|---------------------------:|--------------------------|-----------------------------|
| `44b6_0113de3b` (holdout proxy) | 52 / 25,755 | **1.00 / 1.00** | **1.00 / 1.00** (already over \(T_{\mathrm{true}}\), no add) |
| `6bba_09961292` (holdout dense) | 1950 / 31,117 | 0.96 / 0.82 (0.83) | **0.97 / 0.848** |
| `44b6_18ced818` (density) | 100 / 78,644 | 0.96 / 0.59 (0.60) | 0.96 / **0.63** |
| `44b6_0b24845f` (low contrast) | 51 / 32,795 | **0.10 / 0.06** | **0.78 / 0.47** |
| `6bba_05db0fb1` | 1229 / 69,800 | 0.55 / 0.41 (0.42) | **0.86 / 0.64** |
| `6bba_48816121` | 935 / 23,965 | 0.95 / 0.75 (0.75) | 0.95 / 0.75 (native already over \(T_{\mathrm{true}}\)) |

`0b24845f`: background median 1120, MAD-gate 2432, GT intensity 1337 — **cell below the photon gate**. Paint fraction 0.6% (too *little*; the existing φ² tightening only fires when we paint too much). The photon observer is structurally blind here. The eye is the cone layer for that video. Source: `data/biohub_eye_variety.json`.

**105 / 197** train videos have MAD-gate paint fraction < 1/φ³ (low contrast). Oversampling those in correspondence train (hold out dense+proxy) raised `0b24845f` fill Jaccard **0.47 → 0.67** (find 0.78 → 0.90) but dropped `6bba_05db0fb1` **0.64 → 0.61** and held-out dense **0.848 → 0.839**. Low-contrast weights stay a specialist for the extreme gate-miss (`0b248`, paint 0.6%, GT below thr). They do **not** replace steal-shell, and paint-fraction < 1/φ³ is **not** a blanket switch. Proxy stays 1.00 either way (over \(T_{\mathrm{true}}\)).

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

Dense leftover after leftover-yield: **333 FN** on the photon claim (Jaccard 0.82). Steal-shell correspondence fill is **282 FN** (Jaccard **0.848** on this train video) without breaking proxy. Nucleus-ball extra train **regressed** (0.848→0.838) — keep the steal-shell weights. Of the 282: **86** still 7–12 µm node misses, **196** first-pass linker collisions (40 true-closer at 3.2 vs 6.0 µm). NMS-floor extra yield +1 TP only. Public Final 0.848 / live top **0.962** is the Kaggle target; 0.848 here is one train video, not public LB. Do not Final a kernel until public > 0.848 and the notebook actually runs.

Re-centroid of the same peaks in a φ⁴ window (nucleus / first-pass scale) **failed**: find 0.96→0.92, product-find 0.95→0.88, Jaccard 0.82→0.67. The 7–8 µm shell is not an under-read of the same blob — a larger first moment merges the neighbor NMS already split. That is the free-parameter disconnect in numbers, not a missing radius.

1. Stream one time-point from the OME-Zarr (pixel 3-D) for a viewer — still no copy of the dump.  
2. Join Zebrahub gene-expression tracks on `I:\` to Danio UniProt → product Cα on the same embryo.  
3. Leave U-Net training in `biohub-fsot-unet`; this repo stays the genetics / 3-D *reader*. If the net is used, it is the eye (activation field), not the reported center.  
4. Competition climb: steal-shell fill is **0.848** on held-out dense, **0.47–0.75** on other densities. Low-contrast videos need the eye (native find 0.10). Do not Final until public > 0.848. No extra Kaggle download — 199 train + 4 test already on `D:\`.
