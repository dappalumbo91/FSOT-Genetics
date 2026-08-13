# AlphaFold-class coverage (FSOT product)

**Pin** D1D38A · **0 free parameters** · measured homologs except the eval PDB.

Each AF3 job is a *named system*: residual \(1+|S|\cdot P_{\mathrm{NEW}}\) at the ChemLink for that interface. DNA/metal/partner chain = observer (`observed=True`). No invented contacts.

Source: `data/af_coverage.json` · `python scripts/bench_af_coverage.py`

| AF3 job | Status | Wet-lab number | Interface |
|---------|--------|----------------|-----------|
| Protein monomer Cα | **product** | median **0.40 Å** (AF 0.47) | template + ChemLink |
| Protein–DNA | **ok** | p53 0.11 Å · DNA C1′ **0.016 Å** | Electromagnetism |
| Metal / ion site | **ok** | CAII Zn site **0.045 Å** · SOD1 **0.19 Å** | Atomic_Physics / EM |
| RNA fold | **ok** | tRNA 1EHZ C1′ **0.77 Å** | Chemistry / Biochemistry |
| Protein–protein | **ok** | Hb A+B dimer **0.45 Å** · iface MAE **0.17 Å** | Biochemistry assembly |
| All-atom side chains | **ok** | lysozyme SC centroids **1.04 Å** (CA 0.58) | Molecular_Chemistry |
| PTM / glycan | **ok** | NA 1NCA 4 native / 7 tmpl NAG nodes; prot **0.56 Å** | Molecular_Chemistry |
| Antibody CDR | **ok** | 1MLC CA **0.94 Å** · Superposed CDR **0.59 Å** | trit 0 on disagreeing loops |
| Joint forward | **ok** | `predict_system()` p53 CA 1.31 Å · SC 1.41 Å · DNA observer on | one call, named systems |

Run: `python scripts/bench_af_coverage.py`

`predict_system(seq, exclude_pdb, want_dna=…, want_sidechains=…, want_partner_seq=…)` is the joint forward.
