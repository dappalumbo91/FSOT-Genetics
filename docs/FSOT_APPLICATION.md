# How FSOT is applied (structure product) — do not invert

**Pin:** `D1D38A` · **Law:** \(S = K(T_1+T_2+T_3)\) · **0 free parameters**

## Residual law (archive)

\[
r = 1 + |S_{\mathrm{domain}}| \cdot P_{\mathrm{NEW}}
\]

| Channel | Named domain | Use |
|---------|--------------|-----|
| bond | Physical_Chemistry | residual-weight CA–CA springs |
| clash | Chemistry | residual-weight steric |
| fold / anchor | Biochemistry | residual-weight fidelity + Rg observation |

**Direction:** residual **scales the force/energy of the correct interface**.  
It does **not** invent new ranking free-forms (length×score, medoid shotguns, etc.).

## Correct order of operations

```text
1. MEASURED data eligibility
   - homolog PDB Cα transfer
   - fair: identity ∈ [MIN, CAP], coverage ≥ MIN, exclude self PDB
   - optional isoform expand CAP→0.99 if pool starved (still exclude self)

2. DATA authority when alignment is strong (id×cov ≥ 1/φ)
   - primary template = best id×cov among *crystals* (NMR = Superposed)
   - every intact crystal is a state_rep (residual does not drop 1UBI)
   - residual does NOT re-pick conformational state (`trit_not` stays)

3. RESIDUAL-at-interface when alignment is remote/moderate
   - E = r_bond·Σ(L−CA_CA)² + r_clash·clashes + r_fold·(Rg−target)²
   - rank measured maps by E (lower = better under law)
   - multi-fill weights ∝ r_fold / E
   - never residual-best an NMR ensemble (2LGF E=0.08 / 14 Å)

4. PRODUCT physics (only if transfer is bond-broken)
   - intact: mean (L−CA_CA)² ≤ 1/φ² → keep the measured map
   - broken: fuse_relax residual-weighted bond/clash/anchor
```

## Wrong applications (caused regressions)

| Mistake | Why wrong |
|---------|-----------|
| `score = id × cov × length_sim` | free geometric invent, not residual law |
| residual rank over **full** pool | residual at wrong interface → CaM apo/holo flip |
| residual override of high-id primary | invents state against measured sequence authority |
| residual-best among observations of one collapse | dropped 1UBI 0.09 Å (rank 21) for 2C7M 0.86 Å |
| bond-idealize an intact crystal | 1EXR 0.80 → 1.16 Å (wrong Physical_Chemistry) |
| blend / discard context-flips | DFG-in and DFG-out are `trit_not` of one apparatus (0 = Superposed). Residual must not pick between them. |
| medoid-all / soft disagree switches | geometric shotgun, not \(S=K(T_1+T_2+T_3)\) |

## Ship gate

Same-data product median **≤ 0.47 Å** (AlphaFold median on the freeze set).  
Current freeze: **0.13 Å** (`docs/PRODUCT_FREEZE.md`). Guards (RNase, CaM) must not jump above 3 Å.  
If residual change raises median → **revert** (math at wrong interface).
