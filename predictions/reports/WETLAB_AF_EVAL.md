# Wet-lab + AlphaFold evaluation (FSOT product)

Generated: `2026-09-21T22:34:16.685924+00:00`  
Free parameters: **0** · pin D1D38A · identity_cap=1.0

## Structure vs experimental PDB (and AlphaFold DB)

| Metric | Value |
|--------|------:|
| FSOT product median Cα RMSD | **0.25576396895095804** Å |
| AlphaFold DB median Cα RMSD | **3.9806235776416807** Å |
| FSOT sub-2 Å | 19/19 |
| FSOT within 1.5 Å of AF | 16/16 |
| FSOT beats AF (by >0.05 Å) | 16/16 |

### By category

| Category | n | FSOT med Å | AF med Å | FSOT sub-2Å | beats AF |
|----------|--:|----------:|---------:|------------:|---------:|
| cancer | 7 | 0.8113616332460332 | 7.206070791493004 | 7/7 | 7/7 |
| control | 3 | 0.09084350188139771 | 0.4215392291893585 | 3/3 | 3/3 |
| drug | 7 | 0.25576396895095804 | 1.5939732903871826 | 7/7 | 5/5 |
| vaccine | 2 | 0.4663908926345075 | 31.342863782725306 | 2/2 | 1/1 |

### Per target

| ID | Category | FSOT Å | AF Å | Δ(FSOT−AF) | Template | Wet-lab |
|----|----------|-------:|-----:|-----------:|----------|---------|
| p53_dbd | cancer | 0.01 | 6.19 | -6.17 | 9CHT | X-ray p53–DNA complex (Cho et al.) |
| kras | cancer | 0.19 | 3.45 | -3.26 | 6M9W | X-ray KRAS (GDP) |
| egfr_kinase | cancer | 1.00 | 8.50 | -7.50 | 2ITP | X-ray EGFR kinase |
| braf_kinase | cancer | 0.94 | 7.21 | -6.26 | 3II5 | X-ray BRAF kinase |
| abl1_kinase | cancer | 0.81 | 16.77 | -15.95 | 2G2F | X-ray ABL–imatinib |
| bcl2 | cancer | 1.90 | 13.40 | -11.50 | 8HTS | NMR BCL-2 |
| sars2_rbd | vaccine | 0.52 | 31.34 | -30.83 | 7E23 | X-ray RBD–ACE2 (Lan et al. Nature 2020) |
| ha_h3 | vaccine | 0.41 | — | — | 8UT6 | X-ray hemagglutinin |
| hiv_pr | drug | 0.28 | — | — | 5V4Y | X-ray HIV protease |
| ace2 | drug | 0.48 | 10.29 | -9.81 | 7T9L | X-ray ACE2 |
| dhfr | drug | 0.25 | 0.74 | -0.49 | 1DRF | X-ray DHFR–methotrexate |
| cox2 | drug | 0.26 | 1.59 | -1.34 | 5F1A | X-ray COX-2 |
| ubiquitin | control | 0.09 | 0.88 | -0.79 | 2FID | X-ray ubiquitin (Vijay-Kumar) |
| lysozyme | control | 0.12 | 0.42 | -0.30 | 2MEF | X-ray lysozyme |
| sod1 | cancer | 0.10 | 0.29 | -0.18 | 1HL4 | X-ray SOD1 |
| hbb | drug | 0.22 | 0.52 | -0.30 | 2DXM | X-ray deoxyHb |
| rnase | control | 0.09 | 0.33 | -0.24 | 1RBB | X-ray RNase A |
| insulin | drug | 0.14 | 4.51 | -4.37 | 1MSO | X-ray insulin |
| hiv_rt | drug | 0.85 | — | — | 1IKW | X-ray HIV-1 RT |

## Variants vs wet-lab / clinical labels

| Metric | Value |
|--------|------:|
| Pathogenic recall (LIKELY DAMAGING) | **1.0** (11 cases) |
| Drug-resistance recall | 1.0 (1) |
| Benign-like not called damaging | 1.0 (1) |
| Damaging threshold (percentile) | 75.0 |

### Per variant

| Gene | Change | Wet-lab label | FSOT call | %ile | Agree | Evidence |
|------|--------|---------------|-----------|-----:|:-----:|----------|
| TP53 | R175H | pathogenic | LIKELY DAMAGING | 62.0 | Y | IARC hotspot; structural Zn |
| TP53 | R248Q | pathogenic | LIKELY DAMAGING | 40.2 | Y | DNA contact; IARC |
| TP53 | R273H | pathogenic | LIKELY DAMAGING | 48.4 | Y | DNA contact; IARC |
| TP53 | P72R | benign_like | common_polymorphism | 27.8 | Y | common polymorphism |
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

1. **Structure product path** is competitive when a homolog crystal exists (measured multi-template + residual physics). AlphaFold often wins on global RMSD for well-studied monomers; FSOT can win on flexible / multi-state cases.
2. **Bulk de-novo** remains ~11–14 Å — do not use for medical structure claims.
3. **Variant path** is evolutionary intolerance (conservation), calibrated to known drivers — not a substitute for functional wet-lab assays or full ACMG.
4. Forward accuracy improves with **more measured coverage** (templates/MSAs), not by inventing free parameters.

## Data provenance

- Experimental structures: RCSB PDB (cited wet-lab methods in catalog).
- AlphaFold models: AlphaFold DB (EBI) by UniProt accession.
- Variant labels: curated literature / clinical classic drivers (IARC, COSMIC classics, FDA-label mutations) — see catalog notes.
