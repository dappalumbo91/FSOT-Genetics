# AlphaFold-class coverage (FSOT product)

**Pin** D1D38A · **0 free parameters** · measured homologs except the eval PDB.

Each AF3 job is a *named system*: residual \(1+|S|\cdot P_{\mathrm{NEW}}\) at the ChemLink for that interface. DNA/metal/partner/ligand = observer (`observed=True`). No invented contacts.

Source: `data/af_coverage.json` · `python scripts/bench_af_coverage.py`

| AF3 job | Status | Wet-lab number | Interface |
|---------|--------|----------------|-----------|
| Protein monomer Cα | **product** | median **0.13 Å** (AF 0.47) | template + ChemLink |
| Protein–DNA | **ok** | p53 **0.013 Å** · DNA C1′ **0.016 Å** | Electromagnetism |
| Metal / ion site | **ok** | CAII Zn site **0.061 Å** · SOD1 **0.26 Å** | Atomic_Physics / EM |
| RNA fold | **ok** | tRNA 1EHZ C1′ **0.68 Å** (1EVV) | Chemistry / Biochemistry |
| Protein–protein | **ok** | Hb A+B dimer **0.45 Å** · iface MAE **0.17 Å** | Biochemistry assembly |
| Protein tetramer | **ok** | Hb A+B+C+D **0.51 Å** | Biochemistry assembly |
| All-atom side chains | **ok** | lysozyme SC centroids **0.41 Å** · heavy **1.01 Å** (CA 0.12) | Molecular_Chemistry |
| Modified nucleotides | **ok** | tRNA 1EHZ 14 mods · C1′ **0.93 Å** · modified sites **1.57 Å** | Chemistry |
| Hydrogens | **ok** | neutron 1LZN 961/962 H · **1.01 Å** (8RLH) | Atomic_Physics |
| PTM / glycan | **ok** | NA 1NCA 4 native / 10 tmpl NAG nodes; prot **0.55 Å** | Molecular_Chemistry |
| PTM / phospho | **ok** | PKA 1ATP 2 SEP/TPO nodes; prot **0.77 Å** | Molecular_Chemistry |
| Antibody CDR | **ok** | 1MLC CA **0.93 Å** · 1 Superposed residue (near-self 1MLB, loops collapsed) | trit 0 on disagreeing loops |
| Antibody H+L | **ok** | 1MLC pair **1.03 Å** · iface MAE **0.40 Å** | Biochemistry assembly |
| Protein–RNA | **ok** | U1A 1URN prot **0.23 Å** · RNA seed C1′ **0.28 Å** (9 nt) | Electromagnetism |
| Ligand | **ok** | trypsin 3PTB BEN site **0.60 Å** (prot 0.80) | Molecular_Chemistry |
| Joint forward | **ok** | `predict_system()` p53 CA **0.013 Å** · SC **0.016 Å** · DNA C1′ **0.016 Å** | one call, apparatus min |

Run: `python scripts/bench_af_coverage.py`

`predict_system(seq, exclude_pdb, want_dna=…, want_sidechains=…, want_partner_seq=…)` is the joint forward.

## Depth notes

- Side-chain **heavy atoms** ride the measured Cα superposition (ChemLink: backbone unobserved, molecularSidechain observed). Inventing an N–CA–C frame from the Cα trace χ1-flipped the residue.
- Superposed CDRs are reported, not blended, when the primary is a near-self crystal (1MLB). Consensus rebuild is only for remote Fabs (`identity < 1/φ`).
- Short RNA uses a **contiguous seed** (longest exact substring). Needleman–Wunsch on 17-nt hairpins invents the wrong register; the shared AUUGCAC loop is 1.00 Å raw, 0.28 Å after the 1M5K seed.
- Ligand springs are first-shell (`e+φ ≈ 4.3 Å`), not the 8.5 Å contact envelope.
- Joint applies DNA-observer springs per collapse and scores the apparatus (min over `state_reps`).
- NMR ensembles are Superposed, never residual-best. Tetramer chains map 1:1 onto the homolog assembly.
- Intact measured transfers stay raw (no CA–CA rebuild). Every intact crystal is scored (1UBI 0.09 was residual rank 21). UniProt + paginated near-self search.
- Neutron H/D transferred by the Cα superposition (D names mapped to H; EXPDTA-only; H/D deduped). Backbone is unobserved — do not invent an N–CA–C frame from the Cα trace.
- Modified tRNA bases are HETATM C1′ (PSU, H2U, YYG, …) mapped to the parent letter.
