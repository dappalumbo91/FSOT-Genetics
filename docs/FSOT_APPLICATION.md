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
   - primary template = best id×cov (measured)
   - multi-fill score-powered measured Cα among top-k
   - residual does NOT re-pick conformational state over clear sequence homolog

3. RESIDUAL-at-interface when alignment is remote/moderate
   - E = r_bond·Σ(L−CA_CA)² + r_clash·clashes + r_fold·(Rg−target)²
   - rank measured maps by E (lower = better under law)
   - multi-fill weights ∝ r_fold / E

4. PRODUCT physics (always)
   - fuse_relax: residual-weighted bond/clash/anchor on the measured map
```

## Wrong applications (caused regressions)

| Mistake | Why wrong |
|---------|-----------|
| `score = id × cov × length_sim` | free geometric invent, not residual law |
| residual rank over **full** pool | residual at wrong interface → CaM apo/holo flip |
| residual override of high-id primary | invents state against measured sequence authority |
| medoid-all / soft disagree switches | geometric shotgun, not \(S=K(T_1+T_2+T_3)\) |

## Ship gate

Product H2H median **≤ 1.16 Å** (freeze) and guards (RNase, CaM) must not tank.  
If residual change raises median → **revert** (math at wrong interface).
