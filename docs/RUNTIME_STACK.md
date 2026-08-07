# Runtime stack — trinary law vs bare metal

## Short answer

| Layer | What runs today | Trinary? |
|-------|-----------------|----------|
| **Law** | FSOT seeds + trit ontology `{-1,0,+1}` | Yes (as *syntax*) |
| **Python fold** | `float64` / NumPy evaluating those laws | **Emulated** trinary, not trit silicon |
| **Zig twin** | `zig/src/*` — same pair law; fixed-point path in neuron-zig | Closer to bare metal; **not** the H2H runner yet |
| **Bare-metal trit machine** | Planned / partial in neuron-zig `trit.zig` packing | **Not** driving structure H2H today |

So: we are **not** currently dropping the full fold to a trit CPU and executing opcodes in hardware.  
Python is the **lab / fast formula runner**. Trinary is the **code syntax and interaction law**. Zig is the **body twin** for parity. Bare metal is the **direction**, not the current scoreboard path.

```text
  [1] Codon / AA trinary SYNTAX     ← real {-1,0,+1} words (unique 7-trit opcodes)
  [2] Pair law (φ, π, e, γ)         ← closed form (Zig genetic_pair + Python twin)
  [3] Distogram → MDS → Cα          ← Python NumPy  (this is the H2H engine)
  [4] Optional Zig / fixed / trit   ← parity & future bare-metal path
```

## Why Python first

- Iterate F15 / D_eff interfaces in seconds on the Omen.
- Same residual-gate discipline as Lean (verify_cross).
- Promote kernels to Zig/trit only after law is stable (genome-as-code doctrine).

## What “running trinary” will mean when bare metal is on

1. AA opcodes as packed trit words (`trit.zig` T1 packing).
2. Pair products as trit algebra where possible; seed floats only at geometry boundary.
3. Zig freestanding or QEMU path as the claim runtime; Python stays the lab.

Until then: **honest scoreboard = Python formula path under D1D38A**.

## D_eff interfaces (separate issue)

`D_eff` is the **domain dimensional interface** in the 35-domain table — not the PCA dimension of a folding landscape.  
Protein fold spans **multiple** FSOT domains (chemistry → molecular → biochemistry).  
See `docs/DOMAIN_INTERFACE_FOLD.md`.
