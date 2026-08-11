# Related embodiments (same FSOT pin D1D38A)

| Project | URL / path | Role |
|---------|------------|------|
| **FSOT-Genetics (this repo)** | https://github.com/dappalumbo91/FSOT-Genetics | Genetics product: residual structure + **Zig bare metal cell** |
| **FSOT-2.1-Lean** | https://github.com/dappalumbo91/FSOT-2.1-Lean | Scalar law hub, Lean, multi-domain gates |
| **fsot-neuron-zig** | https://github.com/dappalumbo91/fsot-neuron-zig | Multiboot QEMU mind kernel pattern (boot twin) |
| **Fsot trinary / fsot_os** | desktop `Fsot trinary/fsot_os` | FSOTB cell, QEMU OS, trit ALU (optional cell host) |
| **FSOT-2.1-Neural** | https://github.com/dappalumbo91/FSOT-2.1-Neural | Allen wet-lab / neural monorepo |

**Runtime ladder:** Python research freeze → `zig build host` → `zig build kernel` → QEMU → (optional) FSOTB cell on fsot_os.

Do not treat this repo as a neural fold network.  
Do not evaluate without the D1D38A pin.  
Do not claim product RMSD moves without re-running the freeze gate (`docs/PRODUCT_FREEZE.md`).
