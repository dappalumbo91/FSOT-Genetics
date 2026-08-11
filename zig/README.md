# Zig genetics product runtime (bare metal)

**Shipping runtime lives here.** Python under `../scripts/` is the research
oracle (RMSD benches). This tree is host + freestanding QEMU — same pin
`D1D38A`, residual law, trinary codon front door.

Pattern refs: **fsot-neuron-zig** Multiboot kernel · **Fsot trinary/fsot_os** cell ·  
**FSOT-2.1-Lean** math authority. Roadmap: `../docs/BARE_METAL_GENETICS_ROADMAP.md`.

## Build / run

```powershell
# Host self-test (required gate)
zig build host

# Freestanding Multiboot kernel
zig build kernel

# QEMU serial gate (if qemu-system-x86_64 installed)
.\run_qemu.ps1
```

## Layout

| File | Role |
|------|------|
| `src/seeds.zig` | π, e, φ, γ, derived seeds |
| `src/trit.zig` | Trit ontology {-1,0,+1}, packing |
| `src/codon.zig` | 64-codon PRIMARY/SECONDARY + DNA→AA |
| `src/scalar.zig` | \(S=K(T_1+T_2+T_3)\) |
| `src/product.zig` | **Product residual** \(1+\|S\|\cdot P_{\mathrm{NEW}}\) on pin domains |
| `src/serial.zig` | COM1 UART freestanding console |
| `src/main_host.zig` | Host product gate |
| `src/main_kernel.zig` | Multiboot kernel product cell |
| `src/genetic_pair.zig` | Pair geometry (protein-ready) |
| `linker.ld` / `run_qemu.ps1` | QEMU Multiboot path (neuron-zig twin) |

**Law:** 0 free parameters. Measured homolog Cα remains authority for coordinates;
this cell runs the **law + codon + residual** on metal. Full mind body:
https://github.com/dappalumbo91/fsot-neuron-zig
