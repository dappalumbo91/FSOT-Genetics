# Solidify data dump — residual + evidence contacts

*Generated 2026-08-07T02:35:38.524592+00:00*

**Median Cα RMSD:** 8.590 Å

| Protein | RMSD | R_g err | Top-L structure | Top-L evidence | Evidence recall | Natives LR |
|---------|-----:|--------:|----------------:|---------------:|----------------:|-----------:|
| Ubiquitin | 11.18 | 3.33 | 0.07894736842105263 | 0.11842105263157894 | 0.06293706293706294 | 143 |
| Crambin | 8.59 | 0.70 | 0.043478260869565216 | 0.021739130434782608 | 0.013888888888888888 | 72 |
| Villin headpiece | 5.70 | 2.11 | 0.08333333333333333 | 0.08333333333333333 | 0.15789473684210525 | 19 |
| Protein G B1 | 10.17 | 1.66 | 0.05357142857142857 | 0.14285714285714285 | 0.07272727272727272 | 110 |
| Engrailed HD | 8.29 | 0.90 | 0.0 | 0.037037037037037035 | 0.07142857142857142 | 28 |

## What the evidence ranker hits vs misses

### Ubiquitin

- mean evidence on hits: `{'hydrophobic': 4.717574892279288, 'salt': 1.4660907964826488, 'disulfide': 0.0, 'register': 0.1111111111111111, 'm_ij': 3.405657232850558, 'polarizability': 0.3379580849381736}`
- top predicted (native?):
  - ER 17-53 score=11.597 native=False
  - II 2-60 score=11.347 native=False
  - II 12-60 score=11.296 native=False
  - II 2-43 score=11.246 native=False
  - II 22-60 score=11.219 native=False
  - II 2-35 score=11.165 native=False
- missed natives (best ranks still outside top-L):
  - LV 14-25 rank=87 score=9.154 ev={'m_ij': 1.2763869713473026, 'hydrophobic': 6.477261708828154, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.9999999999999999, 'polarizability': 0.4001536612713431, 'total': 9.153802341446799}
  - FI 44-60 rank=110 score=8.410 ev={'m_ij': 2.2663627898014145, 'hydrophobic': 5.743323886170517, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 8.409840337243274}
  - FI 3-12 rank=115 score=8.256 ev={'m_ij': 2.1123760852713467, 'hydrophobic': 5.743323886170517, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 8.255853632713206}
  - IL 22-55 rank=137 score=8.199 ev={'m_ij': 1.0946166597445492, 'hydrophobic': 6.704672541775057, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 8.19944286279095}
  - IL 2-66 rank=141 score=8.198 ev={'m_ij': 1.0936657767378095, 'hydrophobic': 6.704672541775057, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 8.19849197978421}

### Crambin

- mean evidence on hits: `{'hydrophobic': 4.054689706125153, 'salt': 0.0, 'disulfide': 3.6302882838135426, 'register': 0.9999999999999999, 'm_ij': 10.80877177215399, 'polarizability': 0.6593675647413986}`
- top predicted (native?):
  - CC 15-25 score=20.153 native=True
  - CC 31-39 score=19.331 native=False
  - CC 3-15 score=16.201 native=False
  - II 6-24 score=11.221 native=False
  - CC 2-15 score=10.771 native=False
  - VI 14-24 score=10.416 native=False
- missed natives (best ranks still outside top-L):
  - CI 2-34 rank=47 score=7.287 ev={'m_ij': 1.3398782469971546, 'hydrophobic': 5.433343050173111, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.5136617030252888, 'total': 7.286883000195554}
  - CI 2-33 rank=48 score=7.285 ev={'m_ij': 1.337734937916714, 'hydrophobic': 5.433343050173111, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.5136617030252888, 'total': 7.284739691115114}
  - CI 2-32 rank=50 score=7.282 ev={'m_ij': 1.3354403849586087, 'hydrophobic': 5.433343050173111, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.5136617030252888, 'total': 7.282445138157009}
  - CI 3-32 rank=52 score=7.280 ev={'m_ij': 1.3329823429958565, 'hydrophobic': 5.433343050173111, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.5136617030252888, 'total': 7.279987096194256}
  - FC 12-31 rank=64 score=7.043 ev={'m_ij': 1.2436040731914733, 'hydrophobic': 4.28601618288327, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.9999999999999999, 'polarizability': 0.5136617030252888, 'total': 7.043281959100032}

### Villin headpiece

- mean evidence on hits: `{'hydrophobic': 3.8210185291136995, 'salt': 1.4660907964826488, 'disulfide': 0.0, 'register': 0.0, 'm_ij': 1.3412960872188933, 'polarizability': 0.33795808493817364}`
- top predicted (native?):
  - AL 16-28 score=8.287 native=False
  - FL 10-20 score=8.153 native=False
  - FL 10-22 score=8.139 native=False
  - LV 1-9 score=7.812 native=False
  - VL 9-34 score=7.806 native=False
  - VL 9-20 score=7.806 native=False
- missed natives (best ranks still outside top-L):
  - LF 1-35 rank=37 score=6.390 ev={'m_ij': 0.701036542826107, 'hydrophobic': 5.288886556512465, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.390076760609915}
  - FF 10-17 rank=38 score=6.283 ev={'m_ij': 1.35190036660039, 'hydrophobic': 4.530540201926432, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.282594229798165}
  - FQ 17-25 rank=106 score=1.953 ev={'m_ij': 0.6831463382562837, 'hydrophobic': 0.0, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.9999999999999999, 'polarizability': 0.27012302744893857, 'total': 1.9532693657052222}
  - AQ 18-25 rank=110 score=1.932 ev={'m_ij': 0.6616505391260155, 'hydrophobic': 0.0, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.9999999999999999, 'polarizability': 0.27012302744893857, 'total': 1.9317735665749538}
  - FR 6-14 rank=132 score=1.154 ev={'m_ij': 0.8378780922368003, 'hydrophobic': 0.0, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.3163728419681587, 'total': 1.154250934204959}

### Protein G B1

- mean evidence on hits: `{'hydrophobic': 5.685176491534343, 'salt': 0.5497840486809933, 'disulfide': 0.0, 'register': 0.12499999999999999, 'm_ij': 1.539121460147647, 'polarizability': 0.3768303201464045}`
- top predicted (native?):
  - VV 20-53 score=10.312 native=False
  - IV 5-53 score=9.792 native=True
  - IV 5-38 score=9.721 native=False
  - IV 5-28 score=9.631 native=False
  - IV 5-20 score=9.500 native=False
  - LV 11-53 score=9.086 native=False
- missed natives (best ranks still outside top-L):
  - FF 29-51 rank=58 score=6.583 ev={'m_ij': 1.65191710508232, 'hydrophobic': 4.530540201926432, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.582610968280095}
  - LF 4-29 rank=63 score=6.477 ev={'m_ij': 0.7877514039443224, 'hydrophobic': 5.288886556512465, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.476791621728131}
  - LF 4-51 rank=65 score=6.447 ev={'m_ij': 0.7575721868936773, 'hydrophobic': 5.288886556512465, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.446612404677485}
  - MV 0-20 rank=70 score=6.328 ev={'m_ij': 1.2435030856825144, 'hydrophobic': 4.570958272617417, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.5136617030252888, 'total': 6.32812306132522}
  - KE 3-14 rank=77 score=6.095 ev={'m_ij': 1.4832383241746798, 'hydrophobic': 0.0, 'salt': 4.398272389447946, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.2135669322718347, 'total': 6.095077645894461}

### Engrailed HD

- mean evidence on hits: `{'hydrophobic': 6.704672541775057, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'm_ij': 0.9758772784458095, 'polarizability': 0.4001536612713431}`
- top predicted (native?):
  - ER 8-50 score=10.269 native=False
  - RE 0-25 score=10.043 native=False
  - ER 19-50 score=10.038 native=False
  - II 42-53 score=9.999 native=False
  - RE 2-25 score=9.986 native=False
  - ER 25-50 score=9.944 native=False
- missed natives (best ranks still outside top-L):
  - RE 28-39 rank=60 score=7.885 ev={'m_ij': 3.2727056811020443, 'hydrophobic': 0.0, 'salt': 4.398272389447946, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.2135669322718347, 'total': 7.884545002821826}
  - FF 17-46 rank=89 score=6.583 ev={'m_ij': 1.6522037732695902, 'hydrophobic': 4.530540201926432, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.582897636467365}
  - LF 23-46 rank=97 score=6.474 ev={'m_ij': 0.784791424300783, 'hydrophobic': 5.288886556512465, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.473831642084591}
  - FL 5-37 rank=100 score=6.458 ev={'m_ij': 0.7690647168026518, 'hydrophobic': 5.288886556512465, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.45810493458646}
  - LF 13-46 rank=101 score=6.457 ev={'m_ij': 0.7675825844135733, 'hydrophobic': 5.288886556512465, 'salt': 0.0, 'disulfide': 0.0, 'register': 0.0, 'polarizability': 0.4001536612713431, 'total': 6.456622802197382}

