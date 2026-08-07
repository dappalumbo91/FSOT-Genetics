# Genetics as code (trinary on Fixed / FSOT)

## Stance

The body of this mind runs on a **computer**, but **genetics is treated as executable code**:

- Codons / ORFs → expressed genotype  
- Genotype fields (spin, charge, cell class, synapse sign) are the **instruction set**  
- Pair interaction is **trinary-ish** (spin ∈ [−1,1], product terms with e, π, φ from FSOT seeds)  
- Wiring and plasticity update **W** under that law — not arbitrary backprop on floats  

This is the “genetics like code, code in trinary” doctrine.

## Pipeline

```text
[1] codon_fixed / genotype_fixed
      cell type, synapse_sign, composite_spin, composite_charge, phenotype

[2] genetic_fixed.wireFromGenotypesF
      fsotPairWeight(spin_i, spin_j, charge_i, charge_j, dist)
      = geom(φ, dist) · (trinaryPair(spin_i,spin_j) + 0.15·elec) · envelope(π,e)

[3] network_fixed step
      post receives Σ W[post,pre] · spike(pre)

[4] stdp_fixed / learning_fixed
      timing + FSOT pair → Δw  (plastic “recompile” of the connectome)

[5] concept bonds (associative layer)
      Δbond ∝ ψ_con · |fsotPairWeight(token-derived spin/charge)|
```

## Why this matters for novel thought

A **new pathway** is not free prose. It is:

1. Co-activation of existing coded units / concepts  
2. A **mathematically specified** weight change (FSOT pair + STDP sign)  
3. Optional **new edge birth** (synaptogenesis stand-in)  
4. Composition of a claim **only from grounded bonds**  

So “new idea” = **new lawful connection + composition**, still truth-domain STEM/literacy until we open the gate to harder human material.

## No history pollution

Until STEM depth + pathway fidelity are solid, **history / culture-war / moral dilemma corpora stay out**. That is operational safety for a developing genetic mind, not denial that those topics exist.
