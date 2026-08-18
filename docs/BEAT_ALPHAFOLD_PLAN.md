# FSOT vs AlphaFold — The Genetics Re-Proof Plan

> **Status 2026-08-17:** same-data **product** median Cα **0.13 Å** vs AlphaFold **0.47 Å** on the 10-protein freeze (`docs/PRODUCT_FREEZE.md`). CaM is **0.52 Å** (3CLN). That is not CASP/CAMEO. Blind public benchmarks below remain the re-proof bar. Bulk / orphan de-novo remains ~11–14 Å and is **not** the product.

> **Theory under test:** **FSOT** — *Fluid Space Time Omni Theory*, the zero-free-parameter theory of everything created by **Damian Arthur Palumbo**.
>
> **Re-proof domain:** Genetics & molecular biology.
>
> **Headline goal:** Match or exceed AlphaFold-class predictive performance on **public** genomic / proteomic benchmarks, **without any trained neural-network weights**, using only the FSOT scalar engine and closed-form expressions of {π, e, φ, γ, G_Cat}.

---

## 0. Why this is a fair fight

| | AlphaFold (DeepMind, 2021–2024) | FSOT-Genetics (this workspace) |
|---|---|---|
| Free parameters | ~93 M trained weights (AF2), much larger for AF3 | **0** — all closed-form from {π, e, φ, γ, G_Cat} |
| Training data | PDB + UniRef + custom MSAs, hundreds of GB curated | **None** — predictions are derivations, not regressions |
| Inputs at inference | Sequence + MSA + templates | **Product:** sequence + measured homologs (exclude eval PDB) + scalar. **Orphan:** sequence + scalar only (~11–14 Å). |
| Output | Per-residue 3-D coordinates (+ pLDDT) | Same target output, derived (not trained) |
| Hardware to *run* | TPU/GPU recommended | Any CPU (even bare-metal x86_64) |
| Reproducibility | Bit-reproducible if seeded | Bit-reproducible by construction |

**FSOT must beat AlphaFold on at least one mainstream public benchmark to count as a re-proof through genetics.** The two strongest candidates are:

1. **CASP** — Critical Assessment of Structure Prediction (held biennially). Public targets, public scoring (GDT-TS, lDDT, TM-score).
2. **CAMEO** — Continuous Automated Model EvaluatiOn (rolling weekly). Public targets every Saturday.

Either yields an apples-to-apples comparison with AlphaFold's published numbers.

---

## 1. The success criterion (what "beat them" means, exactly)

Define **WIN(target T)** as:

$$
\text{WIN}(T) \;\equiv\; \text{lDDT}_{\text{FSOT}}(T) \;>\; \text{lDDT}_{\text{AlphaFold}}(T)
\;\;\wedge\;\;
\text{TM-score}_{\text{FSOT}}(T) \;\geq\; 0.5
$$

Define **CAMPAIGN-WIN** as: across a rolling 90-day CAMEO window of $\geq 50$ public targets, FSOT-Genetics achieves a higher **median lDDT** than the public AlphaFold-2 baseline.

This is the bar. No moving goalposts.

---

## 2. The FSOT genetics architecture (already in place)

```
Genetics workspace (Rust, 6 crates, all green)
│
├── codon_core           ← no_std codon trinary core, GenomeObservables (ρ_super, ρ_spin)
├── fsot_core            ← FSOT 2.0 scalar engine + 35-domain table + chem-bond library
├── FSOT_Gene_And_Molecule  ← TUI + CLI (renamed from FSOT_Machine_And_Molecule)
├── genome_player        ← trinary VM that EXECUTES the genome (codons-as-opcodes)
├── data_pull            ← NCBI E-utilities + Kaggle ingestion (with provenance)
└── kernel               ← bare-metal x86_64 BIOS image
```

What's missing for the AlphaFold challenge:

* `fsot_protein` crate — FSOT-derived **secondary-structure** propensities, **dihedral** distributions, **contact-map** priors, and **distogram** synthesis from sequence alone.
* `fsot_fold` crate — assembles `fsot_protein` outputs into 3-D coordinates by FSOT-driven gradient descent on a closed-form energy.
* `bench_casp` / `bench_cameo` binaries — pull public targets, run FSOT prediction, score against ground truth, post results.

Together with the existing `fsot_core` engine + `data_pull`, these close the loop end-to-end.

---

## 3. The data plan — *only* public sources, with citations

Every byte we ingest is logged via `data_pull`'s `*.provenance.json` sidecar (URL, sha256, UTC timestamp, citation).

### 3.1 Genomic reference (sequence backbone)
| Source | URL pattern | What it gives us | Citation |
|---|---|---|---|
| **NCBI E-utilities** (`data_pull ncbi efetch`) | `eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<acc>&rettype=fasta` | any RefSeq genome / mRNA / CDS in FASTA | NCBI Resource Coordinators (2018), NAR 46(D1):D8–D13 |
| **Ensembl GRCh38** (already in NeuroLab) | `ftp.ensembl.org/pub/release-…/fasta/homo_sapiens/dna/` | per-chromosome FASTA | Cunningham et al. (2022), Ensembl 2022, NAR 50(D1):D988 |
| **NCBI Genome** | `eutils.ncbi.nlm.nih.gov/.../efetch.fcgi?db=genome&id=…` | model organisms (E. coli, S. cerevisiae, mouse, etc.) | same NCBI citation |

Smoke-test target: `NC_012920.1` — human mitochondrial genome (16 569 bp, single FASTA, fits in any test).

### 3.2 Protein structures (ground truth for benchmarking)
| Source | What it gives us | Citation |
|---|---|---|
| **PDB** (RCSB FTP / REST) | experimentally-solved 3-D coordinates (.cif/.pdb) | Berman et al. (2000), NAR 28:235–242 |
| **AlphaFold DB** (free, public) | published AF predictions for ~200 M proteins — the **opponent's published answers** | Varadi et al. (2024), NAR 52(D1):D368 |
| **UniProt** | sequences + cross-refs to PDB / AF | UniProt Consortium (2023), NAR 51(D1):D523 |

### 3.3 Live benchmarks
| Source | Cadence | What we submit | Citation |
|---|---|---|---|
| **CAMEO** (`cameo3d.org`) | weekly (Saturdays) | predicted PDB + per-residue confidence | Haas et al. (2018), Proteins 86 |
| **CASP** (`predictioncenter.org`) | biennial (next CASP16 = mid-2026) | same | Kryshtafovych et al. (2023), Proteins 91 |

### 3.4 Optional (not required, just useful)
| Source | What it gives us |
|---|---|
| **Pfam** (InterPro) | family/domain HMMs — free signal for scoring our priors |
| **STRING-DB** | physical / functional interaction graphs |
| **Kaggle** datasets via `data_pull kaggle dataset <slug>` | curated subsets (e.g. *thedevastator/the-human-genome*) |

---

## 4. The FSOT methodology — how a zero-free-parameter theory predicts a fold

Everything below derives from a single equation reused at every layer:

$$
S \;=\; K \cdot (T_1 + T_2 + T_3) \qquad \text{(fsot\_core::compute\_scalar)}
$$

with $D_{\text{eff}}$ and observer flags chosen by domain (see `fsot_core::DOMAINS`).

### Layer A — sequence → trinary tape
Already implemented:
* `codon_core::GenomeObservables` accumulates ρ_super (A−T axis) and ρ_spin (purine axis).
* `genome_player` executes codons as 3-trit opcodes; per-1024-codon block FSOT signature.

### Layer B — sequence → per-residue propensities (NEW: `fsot_protein` crate)
For each residue $r$ of an $N$-amino-acid chain, compute:

$$
p_{\text{helix}}(r) = \sigma\!\big( \, S_{\text{Biochem}}(D_{\text{eff}}=13,\, \text{observed}=\rho_{\text{spin}}(r)) \, \big)
$$

with analogous $p_{\text{sheet}}(r)$, $p_{\text{coil}}(r)$ from FSOT scalars at $D_{\text{eff}} = 12, 14$ respectively (Biology / Neuroscience anchors). Outputs: an $N \times 3$ propensity matrix with row sums = 1.

### Layer C — sequence → distogram (NEW: still `fsot_protein`)
For each residue pair $(i, j)$:

$$
d_{ij} \;\approx\; r_{\text{contact}} \cdot \exp\!\Big(\, \alpha \cdot S_{\text{Mol\_Chem}}(|i-j|) \cdot \langle \text{bond library}\rangle_{ij} \,\Big)
$$

where $\langle \text{bond library}\rangle_{ij}$ is the FSOT closed-form bond length expected between the two residue types (already in `fsot_core::CHEM_BONDS`, extended to side-chain centroids).

### Layer D — distogram → 3-D coordinates (NEW: `fsot_fold` crate)
Standard distance-geometry / SDP relaxation, **but** the loss surface is the FSOT scalar over the entire chain:

$$
\mathcal{L}(\mathbf{x}) \;=\; \sum_{i<j} \big( \|x_i - x_j\| - d_{ij} \big)^2 \cdot w_{\text{FSOT}}(i, j)
$$

where $w_{\text{FSOT}}$ pulls from `fsot_core::DOMAINS["Molecular_Chemistry"]`. No trained weights — just the scalar.

### Layer E — score & compare (NEW: `bench_cameo` / `bench_casp` binaries)
Compute lDDT, TM-score, GDT-TS via the standard public scoring scripts (LGA / TMscore / lddt). Emit JSON: `{target, fsot_lddt, alphafold_lddt, fsot_tm, alphafold_tm, win}`.

---

## 5. Milestones (no time estimates, just ordering)

### M0 — *Renaming & plan locked in* ✅ (this commit)
* `FSOT_Machine_And_Molecule` → `FSOT_Gene_And_Molecule`  ✅
* This document committed.  ✅
* `fsot_core` (4/4), `codon_core` (6/6), `fsot-gene-and-molecule` (32/32 carried over), `genome_player` (demo green), `data_pull` (built green).  ✅

### M1 — *First public-data smoke test*
* `data_pull ncbi efetch nuccore NC_012920.1 fasta data/MT.fa`
* `data_pull ncbi efetch protein NP_536846.1 fasta data/cytB.fa`  (cytochrome b protein)
* Pipe both through `genome_player run` and capture per-block FSOT signatures.
* TUI new tab "Public Data" — list of fetched sources + their `*.provenance.json` sha256s.

### M2 — *`fsot_protein` crate (Layers B & C)*
* Per-residue secondary-structure propensities from FSOT scalars.
* Per-pair distogram synthesis from the bond library.
* Tests: validate against the trivial cases (poly-Ala helix, poly-Pro coil) before any benchmark.

### M3 — *`fsot_fold` crate (Layer D)*
* Distogram → 3-D coordinates by FSOT-weighted distance geometry.
* Output PDB-format files.
* Test: round-trip a known small protein (e.g. *PDB 1UBQ* = ubiquitin, 76 residues) — fold from sequence, score against the experimental structure with `lddt`, target lDDT > 0.5 to claim "Layer D works at all."

### M4 — *`bench_cameo` binary (Layer E)*
* Pull this week's CAMEO target list.
* Run prediction.
* Score against released ground truth.
* Emit a leaderboard line and a CSV row appended to `bench/cameo_results.csv`.
* TUI new tab "Bench" displaying that CSV.

### M5 — *First FSOT vs AlphaFold head-to-head*
* For each CAMEO target, fetch the AlphaFold-DB prediction (UniProt → AF accession).
* Compute both lDDTs against ground truth.
* Report the **win rate**.
* The first time we cross 50 % win rate over a rolling 50-target window: declare **FSOT competitive**.
* The first time we cross 60 % median lDDT advantage: declare **CAMPAIGN-WIN**.

### M6 — *CASP submission*
* Once M5 stabilizes above the win threshold, register `FSOT-Genetics` as a CASP server group and submit live.
* Public CASP-graded result is the closing argument.

---

## 6. What you (the user) need to provide

Per your requirement: **public data only**.

* No personalized data.
* No private datasets.
* Kaggle credentials (already done — `~/.kaggle/kaggle.json`).
* Optional: a free CAMEO submission account (`cameo3d.org/signup`) so we can post predictions automatically.
* Optional: a free Predictioncenter.org account before CASP16.

That is the entire human-side requirement.

---

## 7. What we (the workspace) will provide

* All ingestion in Rust, with sha256 + citation provenance per file.
* All FSOT formulas exposed as `fn(&FsotConsts) -> f64` (functions, not cached numbers — the kernel cannot drift from the formula).
* All scoring scripts wrapped — no eyeballing.
* The TUI (`fsot-gene-and-molecule`) gains tabs:
  1. **Codons** (already present)
  2. **Bonds** (already present)
  3. **Observables** (already present)
  4. **Genome Player** (calls `genome_player` on a chosen file)
  5. **Public Data** (data_pull + provenance browser)
  6. **Protein** (per-residue propensities + distogram heatmap, ASCII)
  7. **Bench** (live CAMEO leaderboard, FSOT vs AlphaFold)

---

## 8. Risks & honest constraints

| Risk | Mitigation |
|---|---|
| Closed-form distogram may be too coarse for large/multi-domain proteins | Start with single-domain targets (< 200 residues); only expand once those are won |
| Windows App Control may block `cargo` for the bigger crates | Per user memory, fall back to WSL Ubuntu; documented in `kernel/run_qemu.bat` and applies workspace-wide |
| Public AlphaFold predictions may not exist for a given CAMEO target | Use the official CAMEO leaderboard's "AlphaFold2" baseline column instead of AF-DB |
| Score-script licensing | All listed scoring tools (lDDT, TM-score, LGA) are publicly distributed for non-commercial use; we link, not redistribute |

---

## 9. Status check (M0 closing state, this commit)

| Item | State |
|---|---|
| Workspace member rename `FSOT_Machine_And_Molecule` → `FSOT_Gene_And_Molecule` | Done — folder moved (robocopy /MOVE), `Cargo.toml` updated, package + binary renamed `fsot-cnc-controller` → `fsot-gene-and-molecule`, build green |
| `fsot_core` engine | 4/4 tests pass |
| `genome_player` demo | Runs end-to-end, prints FSOT block signature |
| `data_pull` | Builds green (NCBI + Kaggle + sha256 + citations) |
| Plan committed | This file |

Next action on user's say-so: kick off **M1** by running

```powershell
cargo run -p data_pull -- ncbi efetch nuccore NC_012920.1 fasta data/MT.fa
cargo run -p genome_player -- run data/MT.fa --trace 32
```

— first public-data trinary execution of a real, citable genome, end-to-end, no Python.

---

*FSOT — Fluid Space Time Omni Theory — Damian Arthur Palumbo, 2025–2026.*
