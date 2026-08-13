# Medical platform orientation — after the 0.13 Å freeze

**Pin** `D1D38A` · **0 free parameters** · freeze: `docs/PRODUCT_FREEZE.md`  
This is orientation, not a clinical claim. Nothing here is a device, a diagnosis, or a substitute for a trial.

---

## 1. Precision vs AlphaFold is no longer the job

Same-data product median **0.13 Å** vs AlphaFold **0.47 Å** on the freeze set. Further Cα grinding (3CLN, rotamers, hydrogens) is **marked** in `docs/OPEN.md`. It is worth doing later as honesty, not as the product.

A clinic does not purchase “0.12 vs 0.42 on lysozyme.” It purchases:

| Question they actually ask | Who they pay today |
|----------------------------|--------------------|
| Is this variant damaging? | SIFT, PolyPhen-2, CADD, REVEL, AlphaMissense |
| What does this VCF *mean* for this gene? | VarSome, Franklin, Invitae reports, ACMG classifiers |
| Will *this* drug hit *this* mutant? | PharmGKB, CPIC, OncoKB, structure viewers + chemist |
| Which conformation is the patient in? | Manual kinase DFG / RAS switch lore |
| Can we run it offline, audit the math, not train on PHI? | Almost nobody |

That last row is the FSOT opening. The first four are where to **orient the same engine**, not where to chase Ångströms.

**Diminishing returns (structure vs AF):** yes, on Cα of well-crystallized globular proteins. **Not** diminishing on: patient DNA front door, variant specificity (P72R still false+), drug as ChemLink observer, multi-state apparatus (DFG, CaM hinge), species-agnostic catalogs.

---

## 2. What AlphaFold actually covers (do not mythologize)

AlphaFold DB **does** include animals, plants, fungi, bacteria — millions of UniProt accessions. AF3 does complexes, nucleic acids, ligands.

What it does **not** do as a medical tool:

- Patient genome → codon → explainable call with 0 trained weights
- Treatment as a **named residual observer** on the patient’s variant map
- `trit_not` apparatus (DFG-in/out kept, not averaged)
- Offline, pin-audited, no cloud training on clinical sequence
- A single workflow from `c.742C>T` to “DNA-contact residue under Electromagnetism”

Species expansion is not “AF never folded a plant.” It is “the same closed-form stack runs on any sequenced organism without a new training run.” That is real.

---

## 3. What this repo already is (do not rebuild)

| Layer | Exists | Script / data |
|-------|--------|----------------|
| DNA → codon → AA → class | yes | `scripts/dna_variant_effect.py` |
| Conservation × trinary change | yes | `scripts/variant_conservation.py` |
| Multi-gene catalog | yes | `scripts/medical_gene_catalog.py` (TP53, KRAS, EGFR, BRAF, CFTR, SOD1, HBB, BRCA1) |
| Panel scoreboard | yes | `scripts/run_medical_variant_panel.py` · 34/35 drivers |
| Structure product | yes | `scripts/bench_product_vs_af.py` · 0.13 Å freeze |
| Drug/cofactor as domain | yes | `scripts/cofactor_nodes.py`, `formulas/smiles_protein_chemistry.json` |
| Ligand observer | yes | trypsin–BEN site 0.24 Å |
| DNA/RNA observer | yes | p53–DNA, tRNA, U1A |
| Context flip | yes | ABL1 3HMI / 3GVU |

The gap is **not** missing math. It is that these are separate benches, not one patient-shaped tool.

---

## 4. One tool: FSOT genetics analyst (patient-shaped)

One forward call, not eight scripts:

```text
patient input
  VCF / HGVS-c / HGVS-p / gene + meds  (optional species)
        │
        ▼
trinary codon layer          synonymous / missense / nonsense / splice-not-yet
        │
        ▼
conservation (UniRef/Pfam)   intolerant site?  (Biology / Biochemistry data)
        │
        ▼
structure apparatus          measured homologs, intact maps, trit_not states
        │
        ▼
ChemLink observers           DNA, metal, ligand/drug SMILES, PTM
        │
        ▼
card (explainable)
  call + state + interface + “why”  ·  0 free parameters
```

**Replace / absorb (honestly):**

| Today’s tool | FSOT role |
|--------------|-----------|
| SIFT / PolyPhen-style missense | Conservation + opcode delta (already started). Must fix **P72R false+** before claiming a replacement. |
| “Look it up in AFDB + PyMOL” | Product structure + named site (DNA contact, Zn, pocket). |
| OncoKB-style gene cards (subset) | Catalog + DNA front door + structure context. Not a literature database. |
| Species-specific fold retrain | Same pin, any UniProt / any sequenced genome. |

**Do not claim to replace:** full WGS secondary analysis (Dragen/GATK), PBPK (GastroPlus), clinical trials, dosing software, ACMG legal sign-out, AlphaFold-as-a-service for orphans.

---

## 5. “Simulate how this patient responds to this drug”

The useful reading of that sentence is **pharmacogenomic residual**, not a molecular movie of the patient.

Lawful version (can be built from what we have):

1. Map the variant onto the product structure (codon → residue → site class: DNA-contact, metal, pocket, hinge, Superposed loop).
2. Place the **medicine** as a Molecular_Chemistry / Electromagnetism observer (SMILES → ChemLink, same as BEN/ATP/Zn).
3. Report whether the variant sits on that observer interface, and whether the apparatus is the drug-relevant collapse (e.g. KRAS G12C + sotorasib cysteine; ABL1 DFG-out + imatinib-class).
4. Combine with conservation: common polymorphism at a non-interface site demotes; invariant DNA-contact promotes.

Unlawful version (do not build, do not market):

- “We simulated your genome under this pill and predicted your outcome.”
- PK/PD, dose, toxicity, time course from FSOT scalars alone.

A **secondary platform** is then: same engine, different intake (med list + VCF), different card (drug–variant–state). Not a second theory.

---

## 6. Evolutionary pathways

Already have the data: Pfam / UniRef alignments.

Productize as:

- Per-site conservation *and* clade split (mammal vs deep homolog) — still data, no fit.
- Domain history (when a domain appears in the family) from InterPro — data.
- Pathway as **gene family + ChemLink class** (kinases share DFG apparatus; globins share heme observer), not a cartoon of “evolution of the patient.”

History of a genetic pathway = measured MSA + named interface. That is in-scope.

---

## 7. Animals, plants, other species

The scalar law and the product path do not know the organism. Expansion is **catalog + intake**, not a new fold engine.

| Vertical | Same engine | New data only |
|----------|-------------|----------------|
| Human clinic | genes we have | ClinVar/CPIC tables as *data* (allele frequency demotion — M7) |
| Veterinary | any sequenced protein | species gene list, vet drugs as ligands |
| Crop / plant | any sequenced protein | herbicide/pathogen ligands as observers |
| Pathogen | viral/bacterial UniProt | already used (SARS RBD still a *search* problem) |

AlphaFold *has structures* for those proteins. It does not ship a **vet/ag/patient analyst** with codon + drug observer + audited residual. That is the expansion.

---

## 8. Build order (medical, not Ångström)

One mechanism at a time. Same pin.

| # | Deliverable | Why medicine cares | Gate |
|---|----------------|--------------------|------|
| 1 | **Single CLI/card:** gene + HGVS-c → class + conservation + site on product structure | Replaces “run four scripts” | P72R not LIKELY DAMAGING; pathogenic recall ≥ 0.9 on catalog |
| 2 | **Drug observer slot:** SMILES or drug name → ChemLink springs on the catalog structure | KRAS G12C / TKI / BEN-class | named site RMSD does not tank product freeze |
| 3 | **State card:** emit both `trit_not` collapses (DFG, RAS, CaM) | Oncologist already thinks in states | both poses kept; residual does not pick |
| 4 | **VCF batch** (offline) | Tumor board / panel | same calls as single-variant path |
| 5 | **Species flag** + non-human catalog smoke | Vet / plant | one plant + one animal gene through the same card |
| 6 | Open structure misses | Honesty | `docs/OPEN.md` — 3CLN, rotamers — *after* 1–3 |

Ship gate for any of this: freeze H2H median stays **≤ 0.47 Å**; no new free parameters; no “patient outcome simulated.”

---

## 9. Does the original idea make sense?

Yes, if restated as:

> One audited FSOT stack that reads a patient’s (or any species’) sequence change and a treatment’s chemistry, and reports **which named interface and which apparatus state** that pair occupies — with measured homologs as authority and residual as the only scale.

No, if heard as:

> A digital twin that predicts clinical response by simulating the genome under each medicine.

We already have the layers. The next work is **one medical front door**, not another RMSD loop.
