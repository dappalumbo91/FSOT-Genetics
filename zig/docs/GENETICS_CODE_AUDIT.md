# Genetics-as-code audit — biologically accurate brain

**Date:** 2026-08-01  
**Doctrine:** Build a biologically accurate brain whose structure and variance come from **FSOT trinary genetics** (64-codon map → ORFs → expression → phenotype → \(W\)), not free FI tables or open polish.

**Archive-solved math (use first):** [`ARCHIVE_SOLVED_MATH_FOR_MIND.md`](ARCHIVE_SOLVED_MATH_FOR_MIND.md) — pin D1D38A, biological brain derivations, codon cert, Allen SDK cohort notes.

**Codon authority in-repo:** [`data/64_codon_trinary_map.txt`](../data/64_codon_trinary_map.txt) · twin `src/codon.zig`

---

## 1. Correct pipeline (what “following the role” means)

```text
FSOT seeds (π, e, φ, γ, G) + pin D1D38A
        │
        ▼
64-codon PRIMARY / SECONDARY trinary map
        │
        ▼
class ORFs (SCN / KCN / CACNA / LEAK per Pyr·PV·SST·VIP)
        │
        ▼
codon_fixed.geneExpression  (φ^spin · e^{|q|/(π n)} · …)
        │
        ▼
genotype_fixed.phenotypeFromGenes  → d_eff, fire, ref, AHP, fi_stim
        │
        ▼
brain_fixed.applyPhenotype + genetic_fixed.fsotPairWeight → W
        │
        ▼
network step / wet cascade / STDP
        │
        ▼
Allen / scalpel / CSV  =  READOUT of genetic population
```

**Allen does not design cells.** Genetics designs cells; Allen scores them.

---

## 2. What already follows the role (KEEP)

| Component | Path | Status |
|-----------|------|--------|
| `codon.zig` / `codon_fixed.zig` | Full 64-codon primary law + expression | **Aligned** |
| `genotype_fixed.buildCellTypeGenotype` | Class ORFs → genes → phenotype | **Aligned** |
| `brain_fixed.initSeeded` | Every unit genotype + `applyPhenotype` + `wireGenetic` | **Aligned** |
| `genetic_fixed.fsotPairWeight` / `wireFromGenotypesF` | Trinary spin/charge → \(W\) | **Aligned** |
| `organism_fixed` | Uses `BrainF` (genetic) | **Aligned** |
| `stdp_fixed` / pathways | Plasticity on genetic \(W\) | **Aligned** (process model) |
| Intel / claim / compose / think | Organism on genetic brain | **Aligned** (behavior layer) |

Core mind body is **genetic code first**. That is the architecture we keep.

---

## 3. Violations / drift (FIX)

| Site | Problem | Severity | Fix direction |
|------|---------|----------|----------------|
| `bio_probe_fixed.defaultBioParams` (old) | Free tables + ad-hoc diversity | **Critical** | → `fillFromGenetics` (**done**) |
| `analyticalLockBioMatch` (old) | Homogenized all cells to Allen mean | **Critical** | Soft blend only; keep genetic identity (**done**) |
| `scalpel_rate_fixed.seedClassParams` (old) | Hand-written FI tables | **Critical** | → `paramsFromCellType` genetics (**done**) |
| `allen_dist_fixed.mapSpecimen` | Allen row → knobs without ORF path | **High** | Bridge only; migrate samples to genetic diversity |
| `allen_class_dist_fixed` | Same bridge + class free nudge | **High** | Seed from class ORFs; refine ORFs for PV rate |
| `applyClassNudge` free multipliers | Extra non-ORF class law on top of class ORFs | **Medium** | Fold into ORF/expression differences only |
| `expressionBias` tables | Hardcoded scales not codon-derived | **Medium** | Eventually seed-fold or ORF-only |
| Dual f64 path (`brain.zig` / `genotype.zig`) | Lab twin; must not become mind authority | **Low** | Keep as parity lab only |

---

## 4. How Allen variance must work under this role

| Wrong | Right |
|-------|--------|
| Free `UnitParamsF` polish until CSV KS closes | Population of **mutated ORFs / class genotypes** → FI → score vs CSV |
| Every cell forced to mean ISI | Genetic diversity produces variance; mean/sd/KS are **readouts** |
| Scalpel hand-tunes PV ref=6 | PV **class ORFs** + expression already set short refractory |

Refinement loop:

```text
open Allen residual
  → change class ORF / expression / phenotype law (genetics)
  → re-express population
  → re-measure
  → not: invent fi_stim outside gene product
```

---

## 5. Gate inventory (who uses what)

| Gate | Must source params from | Notes |
|------|-------------------------|--------|
| `fixed` brain structure | `brain_fixed` genotypes | Already |
| `fixed` Allen pop FI | `fillFromGenetics` + soft lock | Genetics base |
| `scalpel` class rates | `paramsFromCellType` | Genetics base |
| `allen-dist` / class-dist | Prefer genetics; mapSpecimen bridge debt | High priority next |
| `stress` organism | genetic brain | Already |

---

## 6. Immediate code actions (this pass)

1. In-repo `data/64_codon_trinary_map.txt` (authority file next to `codon.zig`).  
2. `bio_probe_fixed.fillFromGenetics` / `paramsFromCellType`.  
3. Scalpel seeds from genotypes.  
4. Soft genetic-preserving Allen blend (no full homogenize).  
5. This audit document.

## 7. Next passes (ordered)

1. **Cre-class FI** entirely from class ORFs; drop free PV stim hacks.  
2. **Allen specimen samples** → optional: pair specimen with nearest genetic phenotype, or replace sample diversity with `mutateOrf` populations sized to match CSV CV.  
3. **Eliminate double class law** (`applyClassNudge` vs class ORFs) by folding into ORF set.  
4. **Wet cascade** channel counts from genotype `n_channels` only.  
5. Document every new module: “knobs from genetics Y/N”.

---

## 8. Review verdict

| Layer | Genetics-as-code? |
|-------|-------------------|
| Mind body (`brain_fixed` / organism / \(W\)) | **Yes** |
| Wet / STDP / intel schedules | **Mostly** (process on genetic lattice) |
| Allen / scalpel / dist lab | **Was drifting free-param; re-anchored to genetics this pass** |
| Full CSV variance from ORF diversity alone | **Work in progress** — correct target |

**Bottom line:** The living brain path was already genetic. Lab FI gates had drifted into free tables. We re-anchor those gates to codon→phenotype and treat further Allen accuracy as **genetic refinement**, not polish theater.
