# Domain interface (D_eff) for protein fold — multi-scale FSOT

## The mistake we may have been making

v7 protein derivations route **all** amplitudes through:

| Role | Domain | D_eff (pin table) | δψ (pin) | Notes |
|------|--------|-------------------|----------|-------|
| chem_amp | Molecular_Chemistry | 9 | 0.5 | OK for residue chemistry |
| region_amp | Biochemistry | 13 | **0.35** | Derivations text said δψ=0.1 — **doc/table mismatch** |
| long-range gate | `ceil(η·13)` | 7 | — | Hard-coded Biochemistry D |

Folding is **not** one interface:

| Physical stage | Typical FSOT domain | D_eff | Why |
|----------------|---------------------|------|-----|
| Covalent / virtual bond | Physical_Chemistry | 8 | backbone geometry |
| H-bond secondary structure | Chemistry | 8 | local polar chemistry |
| Side-chain / residue pairs | Molecular_Chemistry | 9 | molecular interactions |
| Hydrophobic packing / tertiary | Biochemistry **or** Condensed_Matter | 13 / 14 | macromolecule / dense pack |
| Cellular / evolutionary context | Biology | 12 | sequence context (optional) |

Literature (energy-landscape theory) says the *free-energy* funnel can be described with a **few reaction coordinates**, while configuration space is high-dimensional. That is **not** the same object as FSOT `D_eff`. We do **not** free-fit D to RMSD. We **route terms** to named domains from the pin table.

## Lawful multi-interface

### Claim default (`legacy_v7`) — matches protein derivations

```text
F07 backbone          → seed only
F03–F06 + F10–F11     → |S(Molecular_Chemistry)| · P_NEW   # D=9
F13 region pairs      → |S(Biochemistry)| · P_NEW · C_EFF   # D=13
long_range_gate       → ceil(ETA_EFF · 13) = 7
```

### Experimental multi-scale (`multi_scale_v9`) — theory ladder

```text
F07 backbone          → seed only
F03–F06 chemistry     → |S(Molecular_Chemistry)| · P_NEW   # D=9
F10–F11 SS bonuses    → scaled by |S(Chemistry)|           # D=8 H-bond
F13 region pairs      → |S(Biochemistry)| · P_NEW · C_EFF   # D=13
packing top-L         → Condensed_Matter influence          # D=14
```

**Probe (offline 1UBQ/1CRN/…):** `legacy_v7` median ~8.3 Å beat `multi_scale_v9` ~9.4 Å on that set.  
We do **not** auto-switch the claim default to the lowest RMSD (that would be free fitting). Multi-scale stays a diagnostic routing.

All domains are **pre-registered** in `vendor/fsot_compute.py`. No continuous D search.

## Residual-at-interface doctrine

If the wrong domain is used, residual vs experiment grows.  
Probe (`scripts/run_deff_interface_probe.py`) **reports** RMSD/contact under lawful routings for diagnosis — it does **not** auto-pick the winner as a trained dial. Default routing stays theory-first (table above).

## Chemical connection → D_eff (v15 — primary pair routing)

Pair residual now uses **connecting chemical system**, not separation alone:

| Chem link | Domain | D_eff | Observed |
|-----------|--------|------:|:--------:|
| Backbone Cα geometry | Physical_Chemistry | 8 | no |
| Disulfide C–C | Atomic_Physics | 7 | yes |
| Salt bridge (opp. charge) | Electromagnetism | 9 | yes |
| Hydrophobic core (KD>0 both, long) | Condensed_Matter | 14 | yes |
| H-bond secondary (α/β) | Chemistry | 8 | yes |
| Mid-range sidechain | Molecular_Chemistry | 9 | domain default |
| Tertiary long-range | Biochemistry | 13 | yes |

`δψ` / `δθ` come from that domain’s pin row, then trinary-modulated.  
`S = K(T1+T2+T3)` and residual `(1+|S|·P_NEW)` are evaluated **at that interface**.

## Pin table |S| (D1D38A, this machine)

Run: `python scripts/domain_interface.py`

## Gap to close next

1. Bare-metal trit path for pair products (Python trit VM → Zig).
2. Improve F15→D **after** chem-link D_eff residual is honest (error log).
3. Optional Biology D=12 for cellular-context packs only.
