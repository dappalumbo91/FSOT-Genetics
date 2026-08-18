# Reality-margin eval — FSOT product vs wet-lab experimental

> **Historical (2026-08-11, broader 19-chain set, older protocol).** Median **1.17 Å** below is not the 10-protein freeze. Current freeze product: **0.13 Å**. See `docs/PRODUCT_FREEZE.md`.

Generated: `2026-08-11T16:59:00.326133+00:00`  
**Target = experimental PDB Cα, not AlphaFold.**  
Free parameters: **0** · pin D1D38A

## Summary

| Metric | Value |
|--------|------:|
| n OK | 19/19 |
| **Median Cα RMSD vs wet lab** | **1.1674135068802676** Å |
| Sub-2 Å | 12/19 |
| Sub-3 Å | 14/19 |
| Primary mode (most common) | `long_range_contacts` |

### By category (median vs experimental)

| Category | n | Median Å | Sub-2 Å | Dominant mode |
|----------|--:|---------:|--------:|---------------|
| cancer | 7 | 1.6097032453127456 | 4/7 | `long_range_contacts` |
| control | 3 | 1.1526426105020857 | 2/3 | `long_range_contacts` |
| drug | 7 | 1.0249810729173527 | 6/7 | `long_range_contacts` |
| vaccine | 2 | 4.042425950745915 | 0/2 | `global_topology` |

### Per target

| ID | Cat | RMSD Å | Mode | Template | Wet-lab |
|----|-----|-------:|------|----------|---------|
| p53_dbd | cancer | 1.52 | `long_range_contacts` | 2P52 (multi_score) | X-ray p53–DNA complex (Cho et al.) |
| kras | cancer | 1.61 | `long_range_contacts` | 1AA9 (multi_score) | X-ray KRAS (GDP) |
| egfr_kinase | cancer | 4.96 | `long_range_contacts` | 3POZ (multi_medoid_isoform) | X-ray EGFR kinase |
| braf_kinase | cancer | 1.09 | `long_range_contacts` | 1UWJ (multi_medoid_isoform) | X-ray BRAF kinase |
| abl1_kinase | cancer | 4.93 | `long_range_contacts` | 2XYN (multi_score) | X-ray ABL–imatinib |
| bcl2 | cancer | 3.49 | `global_topology` | 1YSG (multi_medoid_neighborhood) | NMR BCL-2 |
| sars2_rbd | vaccine | 5.76 | `global_topology` | 3D0G (multi_medoid_neighborhood) | X-ray RBD–ACE2 (Lan et al. Nature 2020) |
| ha_h3 | vaccine | 2.33 | `long_range_contacts` | 11MS (multi_score) | X-ray hemagglutinin |
| hiv_pr | drug | 0.36 | `long_range_contacts` | 1WBM (multi_medoid_isoform) | X-ray HIV protease |
| ace2 | drug | 5.16 | `global_topology` | 2X90 (multi_medoid_neighborhood) | X-ray ACE2 |
| dhfr | drug | 1.02 | `long_range_contacts` | 1DR1 (multi_score) | X-ray DHFR–methotrexate |
| cox2 | drug | 0.50 | `long_range_contacts` | 3NTB (multi_score) | X-ray COX-2 |
| ubiquitin | control | 2.04 | `long_range_contacts` | 1UD7 (multi_score) | X-ray ubiquitin (Vijay-Kumar) |
| lysozyme | control | 1.15 | `long_range_contacts` | 1BB6 (multi_score) | X-ray lysozyme |
| sod1 | cancer | 1.17 | `long_range_contacts` | 1CBJ (multi_score) | X-ray SOD1 |
| hbb | drug | 1.08 | `long_range_contacts` | 1A9W (multi_score) | X-ray deoxyHb |
| rnase | control | 0.77 | `long_range_contacts` | 1B6V (multi_score) | X-ray RNase A |
| insulin | drug | 1.14 | `long_range_contacts` | 1APH (multi_score) | X-ray insulin |
| hiv_rt | drug | 0.83 | `long_range_contacts` | 6P1X (multi_medoid_isoform) | X-ray HIV-1 RT |

## Fix queue (from reality modes)

### 1. `long_range_contacts` (n=None)
- Tertiary native contacts missing or false
- FSOT: F13–F15 + observer tertiary S at Biochemistry D=13; top-L caps; residual-at-interface on contact set.

### 2. `global_topology` (n=None)
- Chain topology / domain packing globally wrong after Kabsch
- FSOT: MDS is only as good as D; fix contact D first, then sparse polish; multi-start from SS regions.

### 3. `per_residue_hotspots` (n=None)
- Local segments high RMSD after global align
- FSOT: Segment-wise residual; coil vs H/E different D_eff; do not over-constrain termini.

### 4. `mid_range_5_12` (n=None)
- Secondary packing / loops at mid separation
- FSOT: F11 sheet + F12 regions; SS amp at Chemistry domain; full scalar residual mid-sep.

### 5. `helix_period_3_4_7` (n=None)
- α-helix i,i+3/4/7 distances wrong
- FSOT: F10 helix geometry + Chemistry D=8 interface; enforce ideal helix D when p_alpha high.

### 6. `backbone_sep1_2` (n=None)
- Virtual Cα–Cα bond / local geometry wrong
- FSOT: F07 bb + CA_CA seed; do NOT residual-scale sep=1,2 (geometry is hard constraint).

## Doctrine

1. Match **reality** (experimental construct), not AF leaderboard.
2. Localize error → one mechanism → full law only.
3. Product path = measured homolog authority + residual physics.
