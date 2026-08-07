# Error margin log

*Generated 2026-08-07T02:28:07.574688+00:00*

## Protocol

- Fold **experimental PDB sequence** (not UniProt polyprotein).
- Kabsch Cα RMSD + **distance residual by sep bin** + native contact MAE + top-L precision.
- Rank modes → **fix queue** with literature + FSOT handle.
- Full law \(S=K(T_1+T_2+T_3)\) only; **0 free parameters**.

**Median RMSD (this set):** 8.58994132311432 Å

## Per protein

| Protein | n | RMSD Å | Primary mode | Contact MAE Å | Top-L prec |
|---------|--:|-------:|:-------------|--------------:|-----------:|
| Ubiquitin | 76 | 11.18 | `global_topology` | 6.1877596880348715 | 0.07894736842105263 |
| Crambin | 46 | 8.59 | `long_range_contacts` | 12.535356432737489 | 0.043478260869565216 |
| Villin headpiece | 36 | 5.70 | `global_topology` | 3.140972230571936 | 0.08333333333333333 |
| Protein G B1 | 56 | 10.17 | `global_topology` | 5.59237534380693 | 0.05357142857142857 |
| Engrailed HD | 54 | 8.35 | `long_range_contacts` | 8.012722371619954 | 0.0 |

## Fix queue (priority)

### 1. `long_range_contacts` (votes=9)

- **Meaning:** Tertiary native contacts missing or false
- **Literature:** Contact maps / top-L metrics drive fold quality (CASP; Marks/Sander coevolution; AF distograms).
- **FSOT handle:** F13–F15 + observer tertiary S at Biochemistry D=13; top-L caps; residual-at-interface on contact set.
- **Status:** open

### 2. `global_topology` (votes=11)

- **Meaning:** Chain topology / domain packing globally wrong after Kabsch
- **Literature:** Energy landscape funnel; topology from contact order (Onuchic/Wolynes; Baker).
- **FSOT handle:** MDS is only as good as D; fix contact D first, then sparse polish; multi-start from SS regions.
- **Status:** open

### 3. `per_residue_hotspots` (votes=5)

- **Meaning:** Local segments high RMSD after global align
- **Literature:** Flexible loops, termini disorder; core vs surface (crystallographic B-factors).
- **FSOT handle:** Segment-wise residual; coil vs H/E different D_eff; do not over-constrain termini.
- **Status:** open

### 4. `mid_range_5_12` (votes=5)

- **Meaning:** Secondary packing / loops at mid separation
- **Literature:** Secondary structure packing; loop closure; contact order (Plaxco et al.).
- **FSOT handle:** F11 sheet + F12 regions; SS amp at Chemistry domain; full scalar residual mid-sep.
- **Status:** open

### 5. `helix_period_3_4_7` (votes=5)

- **Meaning:** α-helix i,i+3/4/7 distances wrong
- **Literature:** α-helix rise 1.5 Å/res, 3.6 res/turn; Cα i→i+4 ≈ 6.2 Å (Pauling).
- **FSOT handle:** F10 helix geometry + Chemistry D=8 interface; enforce ideal helix D when p_alpha high.
- **Status:** open

### 6. `backbone_sep1_2` (votes=5)

- **Meaning:** Virtual Cα–Cα bond / local geometry wrong
- **Literature:** Standard Cα virtual bond ≈ 3.8 Å; local geometry dominates short-range map (Dill polymer; Flory).
- **FSOT handle:** F07 bb + CA_CA seed; do NOT residual-scale sep=1,2 (geometry is hard constraint).
- **Status:** open

## Next solve (do not skip)

**Mode:** `long_range_contacts`

F13–F15 + observer tertiary S at Biochemistry D=13; top-L caps; residual-at-interface on contact set.

Full JSON: `data\error_margin_log.json`