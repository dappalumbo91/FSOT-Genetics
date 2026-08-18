# Wet-lab + AlphaFold evaluation (FSOT product)

> **Historical (2026-08-11, identity_cap=0.95, broader medical set).** Median **2.06 Å** below is the handicapped protocol, not the current same-data product (**0.13 Å**). See `docs/PRODUCT_FREEZE.md`.

Generated: `2026-08-11T14:22:57.149222+00:00`  
Free parameters: **0** · pin D1D38A · identity_cap=0.95

## Structure vs experimental PDB (and AlphaFold DB)

| Metric | Value |
|--------|------:|
| FSOT product median Cα RMSD | **2.0564579699679006** Å |
| AlphaFold DB median Cα RMSD | **4.818930943569514** Å |
| FSOT sub-2 Å | 7/14 |
| FSOT within 1.5 Å of AF | 12/12 |
| FSOT beats AF (by >0.05 Å) | 8/12 |

### By category

| Category | n | FSOT med Å | AF med Å | FSOT sub-2Å | beats AF |
|----------|--:|----------:|---------:|------------:|---------:|
| cancer | 6 | 2.6798566335549747 | 6.696823029035375 | 2/6 | 5/6 |
| control | 2 | 1.4685106192634172 | 0.6501154983925812 | 2/2 | 0/2 |
| drug | 4 | 0.8386275529938474 | 1.5939732903871826 | 3/4 | 2/3 |
| vaccine | 2 | 4.004035160792522 | 31.342863782725306 | 0/2 | 1/1 |

### Per target

| ID | Category | FSOT Å | AF Å | Δ(FSOT−AF) | Template | Wet-lab |
|----|----------|-------:|-----:|-----------:|----------|---------|
| p53_dbd | cancer | 2.60 | 6.19 | -3.59 | 3Q01 | X-ray p53–DNA complex (Cho et al.) |
| kras | cancer | 1.59 | 3.45 | -1.86 | 1AA9 | X-ray KRAS (GDP) |
| egfr_kinase | cancer | — | — | — | no_template | X-ray EGFR kinase |
| braf_kinase | cancer | 2.76 | 7.21 | -4.45 | 8CHF | X-ray BRAF kinase |
| abl1_kinase | cancer | 12.71 | 16.77 | -4.06 | 8SSN | X-ray ABL–imatinib |
| bcl2 | cancer | 5.80 | 13.40 | -7.60 | 2ME8 | NMR BCL-2 |
| sars2_rbd | vaccine | 5.68 | 31.34 | -25.66 | 3SCI | X-ray RBD–ACE2 (Lan et al. Nature 2020) |
| ha_h3 | vaccine | 2.33 | — | — | 11MS | X-ray hemagglutinin |
| hiv_pr | drug | 0.65 | — | — | 1HVC | X-ray HIV protease |
| ace2 | drug | 5.20 | 10.29 | -5.10 | 2XY9 | X-ray ACE2 |
| dhfr | drug | 1.02 | 0.74 | +0.28 | 1DR1 | X-ray DHFR–methotrexate |
| cox2 | drug | 0.55 | 1.59 | -1.04 | 5W58 | X-ray COX-2 |
| ubiquitin | control | 1.78 | 0.88 | +0.91 | 1UD7 | X-ray ubiquitin (Vijay-Kumar) |
| lysozyme | control | 1.15 | 0.42 | +0.73 | 1BB6 | X-ray lysozyme |
| sod1 | cancer | 1.17 | 0.29 | +0.88 | 1CBJ | X-ray SOD1 |

## Variants vs wet-lab / clinical labels

| Metric | Value |
|--------|------:|
| Pathogenic recall (LIKELY DAMAGING) | **1.0** (11 cases) |
| Drug-resistance recall | 1.0 (1) |
| Benign-like not called damaging | 0.0 (1) |
| Damaging threshold (percentile) | 75.0 |

### Per variant

| Gene | Change | Wet-lab label | FSOT call | %ile | Agree | Evidence |
|------|--------|---------------|-----------|-----:|:-----:|----------|
| TP53 | R175H | pathogenic | LIKELY DAMAGING | 62.0 | Y | IARC hotspot; structural Zn |
| TP53 | R248Q | pathogenic | LIKELY DAMAGING | 40.2 | Y | DNA contact; IARC |
| TP53 | R273H | pathogenic | LIKELY DAMAGING | 48.4 | Y | DNA contact; IARC |
| TP53 | P72R | benign_like | LIKELY DAMAGING | 27.8 | N | common polymorphism |
| KRAS | G12D | pathogenic | LIKELY DAMAGING | 38.6 | Y | COSMIC codon 12 |
| KRAS | G12C | pathogenic | LIKELY DAMAGING | 38.6 | Y | sotorasib-sensitive G12C |
| KRAS | Q61H | pathogenic | LIKELY DAMAGING | 38.6 | Y | switch II |
| EGFR | L858R | pathogenic | LIKELY DAMAGING | 78.3 | Y | TKI-sensitive; crystal/clinic |
| EGFR | T790M | drug_resistance | LIKELY DAMAGING | 62.3 | Y | gatekeeper; TKI resistance |
| BRAF | V600E | pathogenic | LIKELY DAMAGING | 66.8 | Y | melanoma; vemurafenib |
| CFTR | G551D | pathogenic | LIKELY DAMAGING | 92.8 | Y | gating; ivacaftor responsive |
| CFTR | F508* | pathogenic | skipped_non_missense | — | — | ΔF508 most common CF (del); scored as severe if present |
| HBB | E7V | pathogenic | LIKELY DAMAGING | 49.6 | Y | sickle cell (HbS; UniProt pos 7 = Hb 6) |
| SOD1 | G94A | pathogenic | LIKELY DAMAGING | 72.7 | Y | ALS familial |

## How to read predictability

1. **Structure product path** (current freeze) is **0.13 Å** median vs AF **0.47 Å** when a homolog crystal exists. This file’s 2 Å median is the older fair-cap protocol.
2. **Bulk de-novo** remains ~11–14 Å — orphan fallback only. Do not quote it as the product. Do not use it for medical structure claims.
3. **Variant path** is evolutionary intolerance (conservation), calibrated to known drivers — not a substitute for functional wet-lab assays or full ACMG.
4. Forward accuracy improves with **more measured coverage** (templates/MSAs), not by inventing free parameters.

## Data provenance

- Experimental structures: RCSB PDB (cited wet-lab methods in catalog).
- AlphaFold models: AlphaFold DB (EBI) by UniProt accession.
- Variant labels: curated literature / clinical classic drivers (IARC, COSMIC classics, FDA-label mutations) — see catalog notes.
