# Mechanism gap map — where FSOT fails and what to solve

**Source:** `data/wetlab_af_eval.json` (2026-08-11) + product freeze + error-margin doctrine.  
**Rule:** fix **one mechanism at a time** under \(S=K(T_1+T_2+T_3)\), pin `D1D38A`, **0 free parameters**.  
**Attitude:** failures are **diagnostic**, not embarrassment — they tell us which FSOT interface is wrong.

---

## Snapshot (what worked / what didn’t)

| Regime | FSOT median | Note |
|--------|------------:|------|
| Drug targets (good templates) | **0.84 Å** | Mechanism OK: measured Cα + residual physics |
| Controls (Ubq, lysozyme) | 1.47 Å | OK but AF better (~0.65) — **precision ceiling**, not collapse |
| Cancer (mixed) | 2.68 Å | Wins vs AF on several; still not sub-2 on hard kinases |
| Vaccine antigens | ~4 Å | Domain/template confusion |
| **EGFR kinase** | *no template* | **Coverage failure** |
| **ABL1 kinase** | **12.7 Å** | **Wrong measured map** (high id, wrong pose) |
| **BCL-2** | 5.8 Å | Remote homolog / conformation class |
| **SARS-CoV-2 RBD** | 5.7 Å | Wrong family template (3SCI) + domain scope |
| Variants pathogenic | **11/11** | Conservation mechanism works |
| Variants benign-like (P72R) | **false +** | Specificity mechanism missing |

---

## Failure modes → FSOT mechanism → solve handle

Ranked by medical impact × solvability under pure FSOT.

### M1 — Wrong homolog authority (high identity, wrong structure)

| | |
|--|--|
| **Evidence** | ABL1: id≈0.95, cov≈0.88 → **12.7 Å**; BCL-2: id≈0.67 → **5.8 Å**; SARS RBD → template **3SCI** (wrong family) at 5.7 Å |
| **What failed** | Selection still trusts **seq score** (id×cov) over **structural majority of measured homologs** |
| **Not the residual law** | Residual physics cannot fix a wrong measured scaffold |
| **FSOT handle** | Density-weighted multi-template among fair candidates (native-free); Pfam **domain** accession filter; ChemLink D_eff for kinase vs antibody vs viral RBD; reject templates that fail `model_is_sane` + Rg vs `target_rg_fsot` harder |
| **Success metric** | ABL1 product **&lt; 3 Å**; RBD product **&lt; 2.5 Å** with a true sarbecovirus RBD template |
| **Status** | **OPEN — top priority** |

### M2 — Template search coverage holes

| | |
|--|--|
| **Evidence** | EGFR kinase **no_template** (307 aa chain 2ITX) |
| **What failed** | Homolog/Pfam pool empty or filters too strict for that polymer entity |
| **FSOT handle** | Same multi-entity UniProt/Pfam path as p53 fix; InterPro domain range → search **domain sequence** not full EGFR; RCSB sequence search with domain-only query; kinase Pfam PF07714 structure list |
| **Success metric** | EGFR status=ok and product **&lt; 2.5 Å** |
| **Status** | **OPEN** |

### M3 — Domain vs full-chain evaluation mismatch

| | |
|--|--|
| **Evidence** | AF medians inflated when UniProt is multi-domain (spike, ACE2, full kinases) while wet-lab PDB is a **domain** |
| **What failed** | Not always “FSOT better than AF” — sometimes **unfair AF full-length vs domain crystal** |
| **FSOT handle** | Domain-split (already `domain_split_assemble.py`): evaluate product **on wet-lab chain sequence**; AF align only on that span; report domain-scoped RMSD as medical truth |
| **Success metric** | Every wet-lab case has `eval_scope: domain|full` explicit; no silent full-chain AF |
| **Status** | **OPEN (evaluation integrity + product path)** |

### M4 — Conformation / ligand / assembly state

| | |
|--|--|
| **Evidence** | p53 2.60 vs wet-lab DNA complex; BRAF 2.76; kinases open/closed; BCL-2 NMR vs crystal |
| **What failed** | Single template (or multi-fill near high-id) freezes **one** biological state |
| **FSOT handle** | Multi-template **ensemble** residual: weight templates by Biochemistry residual + ChemLink class; optional soft termini (seed-closed); do **not** average open+closed indiscriminately — cluster by structural density first (M1) |
| **Success metric** | p53/BRAF product **&lt; 2.0 Å** without tanking RNase/CaM class winners |
| **Status** | **OPEN** (depends on M1) |

### M5 — Precision gap on already-correct topology

| | |
|--|--|
| **Evidence** | Ubiquitin 1.78 vs AF 0.88; lysozyme 1.15 vs 0.42; DHFR 1.02 vs 0.74 |
| **What failed** | Topology OK; local geometry / packing / flexible termini lag AF’s learned polish |
| **FSOT handle** | Residual-weighted physics **channels** already on (bond/clash/anchor); next: per-pair ChemLink residual on **SS / salt / hbond only** (not shotgun springs); MSA packing polish only in near-contact envelope (already); C-term soft rebuild only when multi-template variance high |
| **Success metric** | Controls median **≤ 1.0 Å** without median product &gt; 1.16 on freeze H2H |
| **Status** | **OPEN — second wave** (do not thrash before M1–M2) |

### M6 — Bulk de-novo information wall

| | |
|--|--|
| **Evidence** | Bulk ~11–14 Å; error-margin log: `long_range_contacts` + `global_topology` dominate |
| **What failed** | Pairwise F15/scalar does not determine full distance field (contacts underdetermine structure) |
| **FSOT handle** | Many-body / fluid geometry distance law (research); until then **product = measured authority only**; bulk is orphan fallback with honest ceiling |
| **Success metric** | Never market bulk as medical structure; orphan path documents ceiling |
| **Status** | **ACKNOWLEDGED ceiling** (not this sprint’s product fix) |

### M7 — Variant specificity (false positives)

| | |
|--|--|
| **Evidence** | TP53 **P72R** (common polymorphism) → LIKELY DAMAGING; pathogenic recall 100% |
| **What failed** | Conservation alone cannot separate polymorphic conserved sites from drivers |
| **FSOT handle** | Population allele frequency as **data** (gnomAD/ClinVar AF when available) — not a trained weight; trinary codon class for synonymous; dual-gate: absolute cons × f_mut **and** “common poly demotion” when AF_pop &gt; seed threshold (e.g. 1/φ³); benign control panel expansion |
| **Success metric** | P72R not LIKELY DAMAGING; pathogenic recall ≥ 0.9 on expanded set |
| **Status** | **OPEN — medical specificity** |

### M8 — Evaluation / AF comparison fairness

| | |
|--|--|
| **Evidence** | SARS AF 31 Å vs domain wet-lab — numbers scare but mix scopes |
| **FSOT handle** | Always report (a) domain-scoped RMSD (b) full-chain RMSD separately; AF model cropped to wet-lab residues after NW align |
| **Success metric** | Report tables never mix scopes without a flag |
| **Status** | **OPEN — metrics hygiene** |

---

## What is *not* broken (do not “fix”)

| Mechanism | Evidence | Keep |
|-----------|----------|------|
| Residual law \(1+\|S\|P_{\mathrm{NEW}}\) | Zig≡Python; drug median 0.84 Å | Keep pin domains |
| Multi-template fill (φ³ / φ⁶) | Freeze product 1.16 Å H2H | Keep |
| Measured template authority | Beats bulk everywhere medical | Never replace with pure bulk |
| Pathogenic conservation signal | 11/11 drivers | Keep absolute + percentile dual gate |
| Zero free parameters | All ship gates | No learned RMSD polish |

---

## Solve queue (FSOT-only, one mechanism at a time)

```text
[1] M1  Wrong homolog authority     → density / domain-Pfam template select
[2] M2  Coverage holes (EGFR)       → domain-scoped RCSB + PF07714
[3] M3+M8 Eval scope honesty        → domain-scoped AF + product metrics
[4] M7  Variant specificity         → pop-AF data demotion + benign panel
[5] M4  Conformation classes        → multi-template density clusters
[6] M5  Precision polish            → ChemLink residual SS only, gated
[7] M6  Bulk wall                   → research many-body; not product
```

**Hard gate for any ship:** product freeze H2H median **≤ 1.16 Å** and wet-lab ABL1 **improves** if M1 ships; no silent regressions on RNase/CaM-class winners.

---

## Mechanism ↔ FSOT domain table (interfaces)

| Mechanism | Primary domain S | ChemLink / D_eff | Residual use |
|-----------|------------------|------------------|--------------|
| Backbone bond | Physical_Chemistry | D=8 backbone | residual on bond channel |
| Local steric | Chemistry | D=8 clash | residual on clash |
| Fold observation / anchor | Biochemistry | D=13 | residual strengthens template fidelity |
| Kinase / catalytic map | Molecular_Chemistry + Biochemistry | D=9 / 13 | template class, not free D search |
| Variant intolerance | Biology / Biochemistry | conservation data | not residual invent |
| Wrong template | — | data selection | **fix selection, not S** |

Doctrine: **wrong interface → worse residual.** If RMSD rises after a residual change, **revert** — math applied at wrong place.

---

## Immediate next experiment (single mechanism)

**M1 only:** on ABL1, BCL-2, SARS RBD, BRAF —

1. Collect full fair homolog pool (no early-exit).  
2. Native-free **density** pick among top-k by score **and** among high-coverage Pfam-domain PDBs.  
3. Compare product RMSD vs current greedy multi-fill.  
4. Ship only if wet-lab medians improve **and** freeze H2H ≤ 1.16.

Do **not** combine with M5 polish in the same commit.

---

## Predictability going forward

| If we solve… | Expected predictability gain |
|--------------|------------------------------|
| M1+M2 | Kinases / antigens become **deployable** structure products |
| M3+M8 | Honest AF comparisons; stop false “wins/losses” |
| M7 | Medical variant **specificity** usable beyond hotspot recall |
| M5 | Close remaining ~0.5–1 Å on controls |
| M6 | Only path to orphan folds — long research |

**Bottom line:** we are not “bad at FSOT.” We are **under-using measured structural authority selection** and **over-trusting sequence identity**, while residual physics and conservation already work where the map is right.
