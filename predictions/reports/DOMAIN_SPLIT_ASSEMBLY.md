# Domain-split FSOT assembly

Generated: `2026-08-11T03:57:09.615031+00:00`  
Free parameters: **0**

| Gene | N | domains | global RMSD | bulk RMSD | best domain RMSD |
|------|--:|--------:|------------:|----------:|-----------------:|
| KRAS | 189 | 1 | 6.74 | 15.71 | 1.65 |
| SOD1 | 154 | 1 | 0.29 | 16.92 | 0.30 |
| HBB | 147 | 1 | 0.30 | 14.83 | 0.25 |
| TP53 | 393 | 4 | 10.31 | 18.74 | 6.07 |

## Notes

- Per-domain templates beat full-chain bulk on multi-domain targets when homologs exist.
- Global RMSD can stay large if domain–domain orientation is unknown (no joint template).
- Medical use: trust **per-domain** coordinates + confidence; treat inter-domain pose as low-confidence.

### KRAS

- Ras family (PF00071) 5-164: RMSD=1.650529054703886 n=159
  source=joint_multi_domain_template tmpl=7VVB id=1.0

### SOD1

- Copper/zinc superoxide dismutase (SODC) (PF00080) 15-150: RMSD=0.3005138677830555 n=136
  source=joint_multi_domain_template tmpl=4B3E id=1.0

### HBB

- Globin (PF00042) 27-142: RMSD=0.2541571741138644 n=116
  source=joint_multi_domain_template tmpl=1DXT id=1.0

### TP53

- P53 transactivation motif (PF08563) None-None: RMSD=None n=0
- Transactivation domain 2 (PF18521) None-None: RMSD=None n=0
- P53 DNA-binding domain (PF00870) 100-288: RMSD=6.068131349087482 n=182
- P53 tetramerisation motif (PF07710) None-None: RMSD=None n=5
  source=joint_multi_domain_template tmpl=6XRE id=1.0
  source=joint_multi_domain_template tmpl=6XRE id=1.0
  source=joint_multi_domain_template tmpl=6XRE id=1.0
  source=joint_multi_domain_template tmpl=6XRE id=1.0
