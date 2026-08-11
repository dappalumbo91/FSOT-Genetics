# Bare-metal genetics — template program

**You hate Python. Good.** Python froze the scoreboard. Metal runs the system.

This doc is the **template** for compiling and running the genetics product on
real substrate (host first, then QEMU freestanding), using the same patterns as:

| Reference | Path / repo | What we steal |
|-----------|-------------|----------------|
| **fsot-neuron-zig** | Multiboot kernel, COM1 serial, `run_qemu.ps1` | Boot + gate printout |
| **Fsot trinary / fsot_os** | QEMU q35, FSOTB cell, trit ALU | Trinary cell + kernel loop |
| **FSOT-2.1-Lean** | Pin `D1D38A`, scalar law | Math authority |
| **FSOT-Genetics zig/** | codon / trit / scalar | Already in-tree twins |

## Port ladder (genetics)

```text
[0] Python scripts/*          — FREEZE metrics (PRODUCT_FREEZE.md)
[1] Lean / formulas           — law + codon map authority
[2] zig/ trit + codon + scalar — in-tree twins (host tests)
[3] zig product residual      — residual law on force channels
[4] freestanding Multiboot    — QEMU -kernel, serial gates
[5] (optional) fsot_os cell   — FSOTB program that calls genetics ops
[6] multi-level HW            — physical trinary when available
```

Binary may **carry** trits (T1 packing). It must not **define** the architecture.
See `zig/docs/TRINARY_BARE_METAL.md`.

## Commands (this repo)

```powershell
cd zig

# Host (fast, no QEMU) — product residual + codon + scalar
zig build host

# Freestanding Multiboot kernel
zig build kernel

# QEMU serial gate (needs qemu-system-x86_64)
.\run_qemu.ps1
```

**Parity gate (Zig host vs Python pin oracle — required before claim):**

```powershell
python scripts/parity_zig_python.py
# → PARITY_GATE PASS  (data/parity_zig_python.json)
```

Compares residual channels, domain S, seeds, ATG codon, DNA→AA fragment,
and one residual-physics bond step. Tolerances: seeds/residual ≤ 5e−6 abs.

Expected serial markers (kernel):

```text
FSOT_GENETICS_KERNEL
FSOT_TRIT PASS
FSOT_CODON PASS
FSOT_SCALAR PASS
FSOT_RESIDUAL PASS
FSOT_PRODUCT_CELL PASS
FSOT_STAGE_GENETICS_OK
```

## Product cell (what “run real stuff” means)

A **genetics product cell** is a freestanding unit that:

1. Loads a DNA/AA sequence (ROM table or host-injected buffer).
2. Maps DNA → trinary codons → AA (`codon.zig`).
3. Evaluates domain scalars \(S = K(T_1+T_2+T_3)\) for Physical_Chemistry /
   Chemistry / Biochemistry (`scalar.zig` + pin domains).
4. Forms residual factors \(r = 1 + |S|\cdot P_{\mathrm{NEW}}\).
5. (Host path later) applies residual-weighted bond/clash/anchor physics to
   measured template Cα buffers — same law as Python `msa_template_fuse.py`.
6. Emits serial / FSOTB-visible gates. **No free parameters.**

Template Cα coordinates stay **measured data** (homolog PDBs baked as tables
or streamed). The cell does not invent bulk folds as the product path.

## Wiring to fsot_os (optional later)

From `Fsot trinary/fsot_os`:

- Compile a small **FSOTB** program (see `fsotb_asm`) that:
  - `SPAWN` genetics residual check
  - `JOIN` and report residual word on serial
- Mount genetics ROM as ramdisk blob (same idea as brain safetensors).
- Kernel already has `trinary.rs`, `fsot_scalar.rs` — keep **one pin**.

Do **not** reimplement the scalar with different constants.

## Rust path (optional twin)

`crates/fsot_core` is already `no_std`-capable. Product residual can live as
`crates/fsot_protein` residual module and link into a freestanding kernel the
same way `fsot_os` does. Zig is the first genetics freestanding target because
neuron-zig already proved the QEMU Multiboot path.

## Non-negotiables

1. Pin `D1D38A` / seed table matches Lean + `vendor/fsot_compute.py`.
2. Residual law only: \(1+|S|\cdot P_{\mathrm{NEW}}\) on named domains.
3. Measured structure authority for product coordinates.
4. Gate printout on serial; host test must pass before kernel claim.
5. Python is oracle for RMSD benches — not the shipping runtime.
