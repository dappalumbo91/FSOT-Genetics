# Medical variant panel — multi-gene FSOT conservation

Generated: `2026-08-10T23:42:19.961112+00:00`  
Free parameters: **0**  
Genes: **8**  
Drivers scored: **35**  
Drivers called LIKELY DAMAGING: **17** (49%)  
Mean driver impact percentile: **69.21868573626561**

Thresholds: ≥75 LIKELY DAMAGING · 40–75 uncertain · <40 likely tolerated

## TP53 — Tumor protein p53

UniProt `P04637` · n=393 · Pfam `PF08563` · MSA rows=433 · mean cons=0.05

Indication: cancer hotspot / Li-Fraumeni

Mean driver percentile: **84.36567609499942** · LIKELY DAMAGING 6/7

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.R175H | 0.97 | 90 | LIKELY DAMAGING | structural Zn node |
| p.G245S | 0.97 | 88 | LIKELY DAMAGING | structural |
| p.R248Q | 0.98 | 90 | LIKELY DAMAGING | DNA contact |
| p.R248W | 0.98 | 91 | LIKELY DAMAGING | DNA contact |
| p.R249S | 0.96 | 86 | LIKELY DAMAGING | structural |
| p.R273H | 0.92 | 76 | LIKELY DAMAGING | DNA contact |
| p.R282W | 0.90 | 70 | uncertain | structural |

DNA front door:

- `c.524G>A` CGC→CAC R175H (missense) → **LIKELY DAMAGING**
- `c.733G>A` GGC→AGC G245S (missense) → **LIKELY DAMAGING**
- `c.742C>T` CGG→TGG R248W (missense) → **LIKELY DAMAGING**
- `c.818G>A` CGT→CAT R273H (missense) → **LIKELY DAMAGING**
- `c.844C>T` CGG→TGG R282W (missense) → **uncertain**
- `c.744G>A` CGG→CGA R248= (synonymous) → **likely benign**

## KRAS — GTPase KRas

UniProt `P01116` · n=189 · Pfam `PF00071` · MSA rows=3000 · mean cons=0.29

Indication: oncogene / MAPK

Mean driver percentile: **59.90131578947369** · LIKELY DAMAGING 2/6

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.G12D | 0.33 | 42 | uncertain | codon 12 classic |
| p.G12V | 0.33 | 42 | uncertain | codon 12 classic |
| p.G12C | 0.33 | 42 | uncertain | codon 12 (sotorasib) |
| p.G13D | 0.51 | 62 | uncertain | codon 13 |
| p.Q61H | 0.74 | 86 | LIKELY DAMAGING | switch II |
| p.Q61L | 0.74 | 85 | LIKELY DAMAGING | switch II |

## SOD1 — Superoxide dismutase [Cu-Zn]

UniProt `P00441` · n=154 · Pfam `PF00080` · MSA rows=3000 · mean cons=0.41

Indication: ALS

Mean driver percentile: **88.87318968562346** · LIKELY DAMAGING 3/4

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.A5V | 0.00 | — | insufficient_MSA_coverage | A4V historic / aggressive ALS |
| p.G94A | 0.88 | 89 | LIKELY DAMAGING | G93A historic mouse model |
| p.H47R | 0.89 | 91 | LIKELY DAMAGING | Cu ligand |
| p.G38R | 0.85 | 87 | LIKELY DAMAGING | ALS |

## HBB — Hemoglobin subunit beta

UniProt `P68871` · n=147 · Pfam `PF00042` · MSA rows=3000 · mean cons=0.24

Indication: sickle cell / hemoglobinopathy

Mean driver percentile: **40.89912280701755** · LIKELY DAMAGING 0/3

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.E7V | 0.00 | — | insufficient_MSA_coverage | sickle cell HbS |
| p.E7K | 0.00 | — | insufficient_MSA_coverage | HbC |
| p.E122Q | 0.23 | 41 | uncertain | HbD Punjab |

## EGFR — Epidermal growth factor receptor

UniProt `P00533` · n=1210 · Pfam `PF01030` · MSA rows=3000 · mean cons=0.03

Indication: NSCLC / kinase inhibitors

Mean driver percentile: **69.10198648371903** · LIKELY DAMAGING 2/4

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.G719S | 0.91 | 94 | LIKELY DAMAGING | exon 18 |
| p.T790M | 0.35 | 53 | uncertain | gatekeeper resistance |
| p.L858R | 0.71 | 85 | LIKELY DAMAGING | exon 21 classic |
| p.L861Q | 0.25 | 44 | uncertain | exon 21 |

## BRAF — Serine/threonine-protein kinase B-raf

UniProt `P15056` · n=766 · Pfam `PF00069` · MSA rows=3000 · mean cons=0.07

Indication: melanoma / MAPK

Mean driver percentile: **86.50937689050212** · LIKELY DAMAGING 1/3

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.V600E | 0.00 | — | insufficient_MSA_coverage | activation loop classic |
| p.V600K | 0.00 | — | insufficient_MSA_coverage | activation loop |
| p.G469A | 0.70 | 87 | LIKELY DAMAGING | P-loop |

## CFTR — Cystic fibrosis transmembrane conductance regulator

UniProt `P13569` · n=1480 · Pfam `PF00664` · MSA rows=3000 · mean cons=0.03

Indication: cystic fibrosis

Mean driver percentile: **40.181371302652764** · LIKELY DAMAGING 2/4

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.G551D | 0.00 | — | insufficient_MSA_coverage | gating (ivacaftor) |
| p.R117H | 0.06 | 1 | likely tolerated | mild/variable |
| p.N1303K | 0.47 | 80 | LIKELY DAMAGING | NBD2 classic |
| p.G542X | 0.36 | 100 | LIKELY DAMAGING | nonsense (special) |

## BRCA1 — Breast cancer type 1 susceptibility protein

UniProt `P38398` · n=1863 · Pfam `PF00097` · MSA rows=3000 · mean cons=0.01

Indication: hereditary breast/ovarian cancer

Mean driver percentile: **61.315100331058716** · LIKELY DAMAGING 2/5

| HGVS | cons | impact% | call | note |
|------|-----:|--------:|------|------|
| p.C61G | 1.00 | 91 | LIKELY DAMAGING | RING Zn finger pathogenic |
| p.C64Y | 1.00 | 95 | LIKELY DAMAGING | RING Zn finger |
| p.R1699W | 0.23 | 62 | uncertain | BRCT pathogenic |
| p.A1708E | 0.19 | 39 | likely tolerated | BRCT pathogenic |
| p.M1775R | 0.10 | 20 | likely tolerated | BRCT pathogenic |

## Honesty

- Conservation from real Pfam MSAs; impact = cons × (1 − mutant frequency).
- Not a substitute for ACMG clinical classification; research / triage tool.
- Domain-aware Pfam selection when variant falls in a static domain range.
