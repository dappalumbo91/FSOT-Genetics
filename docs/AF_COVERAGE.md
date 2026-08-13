# AlphaFold-class coverage (FSOT product)

**Pin** D1D38A · **0 free parameters** · measured homologs except the eval PDB.

Each AF3 job is a *named system*: residual \(1+|S|\cdot P_{\mathrm{NEW}}\) at the ChemLink for that interface. DNA/metal/partner/ligand = observer (`observed=True`). No invented contacts.

Source: `data/af_coverage.json` · `python scripts/bench_af_coverage.py`

| AF3 job | Status | Wet-lab number | Interface |
|---------|--------|----------------|-----------|
| Protein monomer Cα | **product** | median **0.40 Å** (AF 0.47) | template + ChemLink |
| Protein–DNA | **ok** | p53 0.11 Å · DNA C1′ **0.016 Å** | Electromagnetism |
| Metal / ion site | **ok** | CAII Zn site **0.045 Å** · SOD1 **0.19 Å** | Atomic_Physics / EM |
| RNA fold | **ok** | tRNA 1EHZ C1′ **0.77 Å** | Chemistry / Biochemistry |
| Protein–protein | **ok** | Hb A+B dimer **0.45 Å** · iface MAE **0.17 Å** | Biochemistry assembly |
| All-atom side chains | **ok** | lysozyme SC centroids **1.04 Å** · heavy **1.88 Å** · bb N/C/O **0.77 Å** (CA 0.58) | Molecular_Chemistry |
| PTM / glycan | **ok** | NA 1NCA 4 native / 10 tmpl NAG nodes; prot **0.55 Å** | Molecular_Chemistry |
| PTM / phospho | **ok** | PKA 1ATP 2 SEP/TPO nodes; prot **0.77 Å** | Molecular_Chemistry |
| Antibody CDR | **ok** | 1MLC CA **0.94 Å** · Superposed CDR **0.59 Å** | trit 0 on disagreeing loops |
| Antibody H+L | **ok** | 1MLC pair **1.03 Å** · iface MAE **0.40 Å** | Biochemistry assembly |
| Protein–RNA | **ok** | U1A 1URN prot **1.60 Å** · RNA seed C1′ **0.28 Å** (9 nt) | Electromagnetism |
| Ligand | **ok** | trypsin 3PTB BEN site **0.24 Å** (prot 0.77) | Molecular_Chemistry |
| Joint forward | **ok** | `predict_system()` p53 CA **0.39 Å** · SC 1.36 Å · DNA C1′ **0.016 Å** | one call, apparatus min |

Run: `python scripts/bench_af_coverage.py`

`predict_system(seq, exclude_pdb, want_dna=…, want_sidechains=…, want_partner_seq=…)` is the joint forward.

## Depth notes

- Side-chain **heavy atoms** sit in the product N–CA–C residue frame (not just centroids).
- Superposed CDRs are reported, not blended, when the primary is a near-self crystal (1MLB). Consensus rebuild is only for remote Fabs (`identity < 1/φ`).
- Short RNA uses a **contiguous seed** (longest exact substring). Needleman–Wunsch on 17-nt hairpins invents the wrong register; the shared AUUGCAC loop is 1.00 Å raw, 0.28 Å after the 1M5K seed.
- Ligand springs are first-shell (`e+φ ≈ 4.3 Å`), not the 8.5 Å contact envelope.
- Joint applies DNA-observer springs per collapse and scores the apparatus (min over `state_reps`).
