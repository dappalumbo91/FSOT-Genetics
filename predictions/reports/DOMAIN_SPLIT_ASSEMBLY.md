# Domain-split FSOT assembly

Generated: `2026-08-10T23:45:42.618344+00:00`  
Free parameters: **0**

| Gene | N | domains | global RMSD | bulk RMSD | best domain RMSD |
|------|--:|--------:|------------:|----------:|-----------------:|
| KRAS | 189 | 1 | 9.09 | 15.71 | 1.65 |
| SOD1 | 154 | 1 | 10.80 | 16.92 | 0.57 |
| HBB | 147 | 1 | 23.36 | 14.83 | 0.31 |
| TP53 | 393 | 4 | 30.83 | 18.74 | 0.57 |
| EGFR | 1210 | 6 | 123.40 | 24.22 | 3.35 |

## Notes

- Per-domain templates beat full-chain bulk on multi-domain targets when homologs exist.
- Global RMSD can stay large if domain–domain orientation is unknown (no joint template).
- Medical use: trust **per-domain** coordinates + confidence; treat inter-domain pose as low-confidence.

### KRAS

- Ras family (PF00071) 5-164: RMSD=1.6480835632505704 n=159
  source=template_physics tmpl=7VVB id=1.0

### SOD1

- Copper/zinc superoxide dismutase (SODC) (PF00080) 15-150: RMSD=0.5726784251337873 n=136
  source=template_msa_fuse tmpl=1HL4 id=1.0

### HBB

- Globin (PF00042) 27-142: RMSD=0.3073663281545889 n=116
  source=template_physics tmpl=1A0U id=1.0

### TP53

- P53 transactivation motif (PF08563) None-None: RMSD=None n=0
- Transactivation domain 2 (PF18521) None-None: RMSD=None n=0
- P53 DNA-binding domain (PF00870) 100-288: RMSD=0.571344075636834 n=182
- P53 tetramerisation motif (PF07710) None-None: RMSD=None n=5
  source=template_physics tmpl=2K8F id=1.0
  source=template tmpl=2B3G id=1.0
  source=template_physics tmpl=1GZH id=1.0
  source=template_physics tmpl=1OLG id=1.0

### EGFR

- Receptor L domain (PF01030) None-None: RMSD=None n=0
- Furin-like cysteine rich region (PF00757) None-None: RMSD=None n=0
- Receptor L domain (PF01030) None-None: RMSD=None n=0
- Growth factor receptor domain IV (PF14843) None-None: RMSD=None n=0
- Epidermal growth factor receptor transmembrane-j (PF21314) None-None: RMSD=None n=0
- Protein tyrosine and serine/threonine kinase (PF07714) 714-966: RMSD=3.3520627222553014 n=249
  source=bulk_single tmpl=None id=None
  source=template_msa_fuse tmpl=1IVO id=1.0
  source=template_msa_fuse tmpl=1IVO id=1.0
  source=template_msa_fuse tmpl=3B2V id=1.0
  source=template_msa_fuse tmpl=2N5S id=1.0
  source=template_msa_fuse tmpl=1M14 id=1.0
