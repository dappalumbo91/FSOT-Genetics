# Trinary bare metal — porting FSOT down to the substrate

**Critique accepted:** the permanent body is **not** “logical trits on forever-binary silicon.”  
The point of porting downward is a **trinary computational substrate** at bare metal (firmware / process image / future custom or multi-level hardware), with FSOT codons, spins, and neural states as first-class \(\{-1,0,+1\}\) values.

Python (and any temporary binary host) only **simulates** that substrate for research. Embodiment languages (Zig / Rust / Ada) implement the **trinary machine model** and eventually sit as close to metal as the platform allows.

---

## 1. Why trinary (not binary encoding as the end state)

| FSOT structure | Native trinary meaning |
|----------------|------------------------|
| Codon primary map | base → trit; codon → trit³ |
| Membrane trinary | damped / superposed / emergent |
| Genetic spin | \(\tau \in [-1,1]\) collapsed to trit classes |
| Synaptic polarity | −1 / 0 / +1 conductance class |
| SR-ITE continuous flow | phase_diff Spin_Up / Spin_Down / superposed (archive Zig) |

Binary is a **host contingency** while we develop. **Target:** trinary all the way down — memory cells, ops, and interconnect semantics.

Industry has pursued multi-level / non-binary computation (including large-vendor research programs). We do not need their stack; we need a **clear FSOT trinary machine model** we can implement and harden.

---

## 2. Machine model (specification)

### 2.1 Trit

\[
\mathbb{T} = \{-1,\, 0,\, +1\}
\]

| Symbol | Name | Neural role |
|--------|------|-------------|
| \(+1\) | up / emergent / purine class | fire / SCN-favoring / E drive |
| \(0\) | superposed / quiet | refractory grey / secondary codon GC |
| \(-1\) | down / damped / pyrimidine class | silence / I polarity / damp |

### 2.2 Registers & memory

- **Trit cell:** one value in \(\mathbb{T}\) (bare-metal layout: own packing format, not “just an `i8` forever”).  
- **Trit word:** \(w = (t_{k-1},\ldots,t_0)\) for codon-sized \(k=3\) or bus-width \(k=27,81,\ldots\) (powers aligned to codon geometry preferred).  
- **State slab:** per-unit fields stored as trinary where discrete; continuous \(S\) may stay fixed-point until pure discrete collapse paths exist.

### 2.3 Core ops (must be native in body kernel)

| Op | Definition | Use |
|----|------------|-----|
| `neg(t)` | \(-t\) | polarity flip |
| `abs0(t)` | \(0\) if \(t=0\) else \(1\) as class | activity mask |
| `consensus(a,b)` | \(a\) if \(a=b\) else \(0\) | agreement / gate |
| `sum_sat(a,b)` | clamp \(a+b\) to \(\mathbb{T}\) | local field |
| `pair(a,b)` | \(a\cdot b\) in \(\mathbb{T}\) | trinary pair term of \(W\) |
| `primary(base)` | A,G→+1; C,T→−1 | genetic law |
| `from_S(s, ℓ, h)` | piecewise thresholds | membrane trit |

Reference implementation: `fsot_nuron/trinary_substrate.py` (host oracle).  
Bare-metal ports must match these tables bit-for-trit.

### 2.4 Packing (transition formats)

Until physical ternary RAM exists on the host machine, the **firmware ABI** still uses a **defined trinary packing** (not ad-hoc floats):

| Scheme | Layout | Role |
|--------|--------|------|
| **T1** | 2 bits / trit: `00=0, 01=+1, 11=−1` (spare `10` illegal) | wire format |
| **T3** | 6 bits / codon (3×T1) | genetic bus |
| **T27** | 54 bits / 27-trit word (padded to 64-bit container *only as carrier*) | bulk state |

**Important:** containers may ride on binary buses *as transport* during development; the **semantics and ALU are trinary**. The roadmap end-state is hardware whose cells and ops are multi-level/trinary, not a permanent confession that “everything is binary.”

---

## 3. Port ladder (downward)

```text
[0] Python torch host          — research accuracy, Allen, thesis
[1] Lean formal panel          — codon / fold / signs (math)
[2] Trinary oracle module      — pure trit ops + packing (this doc)
[3] Zig|Rust|Ada kernel        — no_std / freestanding step in trit words
[4] Bare-metal / QEMU image    — archive-style boot; state on disk substrate
[5] Multi-level / custom HW    — physical trinary or multi-level cells when available
```

SR-ITE archive already runs Zig continuous flow + virtual substrate binary *blob* with trinary *logic*. We extend that lineage: **trinary is the law of the body**, not a paint job on a binary neural net.

---

## 4. ABI sketch (`fsot_trit_abi`)

```text
struct TritWord { u64 pack; u8 n_trits; }   // carrier only
struct UnitState {
  Trit  membrane;      // from_S
  Trit  spin_class;    // genetic composite quantized
  i16   refractory;    // may stay integer count
  i16   S_fp;          // fixed-point S until discrete-only path
}
fn step(units: []UnitState, W_trit: TritMatrix, ext: []Trit) void
fn inject_sensory(region, TritWord features) void
fn inject_metrics(MetricTrits) void   // subconscious plant
```

Python host implements the same contracts for parity tests.

---

## 5. Sensory + subconscious stay trinary at the edge

- Vision/audio features **quantize to trit streams** before cortex inject (not float forever).  
- System metrics → trit interoception word → thalamus.  
- Human-like + computer-native senses share one **trinary sensory bus**.

---

## 6. Non-negotiables

1. Codon primary law unchanged (Lean + map file).  
2. Trit ops table is the portability contract.  
3. No free parameters on scalar path.  
4. Binary may carry; it must not **define** the architecture.  
5. Thesis / FORMULAS updated when trit ISA changes.

---

## 7. Commands

```powershell
python -c "from fsot_nuron.trinary_substrate import self_test; print(self_test())"
python run_brain_design.py --profile ai_efficient --sensory
```
