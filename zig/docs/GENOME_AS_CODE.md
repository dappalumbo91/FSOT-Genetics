# Genome as code — cellular expand (Python lab → Zig body)

**Thesis:** Treat **DNA codons, chemical bonds, and molecular programs** as the real “source code” of the neural body.  
**Python** is the **lab / immune-like prototyping layer** (test a patch safely).  
**Zig** is the **body / genome expression layer** (once the patch is accurate).  
**FSOT math** (`S = K(T1+T2+T3)`, archive D1D38A) is the **physics of the cell**, not a free-fit optimizer.

This is **not** “train an LLM on GitHub.”  
It is **cellular software engineering**: ORF → phenotype → dynamics → verify vs wet-lab → promote to bare metal.

---

## 1. Biological analogy

| Biology | FSOT-2.1-Neural |
|---------|-----------------|
| DNA / codon table | [`data/64_codon_trinary_map.txt`](../data/64_codon_trinary_map.txt) + class ORFs in `genotype.zig` · audit [`GENETICS_CODE_AUDIT.md`](GENETICS_CODE_AUDIT.md) |
| Transcription / translation | Codon → AA/process → ion-channel phenotype |
| Membrane dynamics | `FSOTNeuronBatch` / Zig `neuron.zig` step |
| Synaptic proteins | Genetic \(W\) from trinary spins + φ geometry |
| Immune / proofreading | Wet-lab battery + stress + Python parity before Zig promote |
| Cell division / expand | Grow population / regions under same genotype laws |
| Communication | Machine words + future NLP/Shakespeare bridges |

---

## 2. Two coding apparatuses, one genome

```text
  [Genome / chemistry]     codon map + gene ORFs + FSOT seeds
            │
            ▼
  [Python lab]             prototype patch · scalpel · SME · unit tests
            │  PASS wet-lab gates
            ▼
  [Zig body]               promote step / pack / inject · QEMU / host
            │
            ▼
  [Measure]                Allen rates · SME · pin · parity ΔS
```

| Layer | Language | Role |
|-------|----------|------|
| Lab | **Python** | Fast iteration, data locks, console, battery |
| Body | **Zig** | Trinary substrate, freestanding path, MachineFrame inject |
| Law | **FSOT archive** | Scalar + folds — zero free parameters |

---

## 3. Patch protocol (how “code learning” works here)

1. **Propose** a cellular change in Python (e.g. ORF tweak, region drive, encode schedule).  
2. **Express** through codon/genotype APIs — not arbitrary weight matrices.  
3. **Verify** with `run_wetlab_accuracy_battery.py` + relevant probe.  
4. **Promote** only the verified kernel/ABI pieces to Zig (`frame_inject`, neuron step, …).  
5. **Parity** Python ↔ Zig (`scripts/parity_zig_neuron.py`).  

Scaffold module: `fsot_nuron/cellular_expand.py`  
CLI: `python run_cellular_expand.py --check`

---

## 4. What we do *not* do

- Fit free parameters to chase leaderboards  
- Replace FSOT with gradient soup  
- Treat Morse as the genome  
- Claim DNA “writes Zig” without wet-lab gates  

---

## 5. Commands

```powershell
python run_wetlab_accuracy_battery.py
python run_cellular_expand.py --check
python run_stress_suite.py
python scripts/parity_zig_neuron.py
```
