# Zig runtime vs Python oracle — parity

**Gate script:** `scripts/parity_zig_python.py`  
**Artifact:** `data/parity_zig_python.json`  
**Last result:** **PASS** (all keys within tolerance)

## What runs where

| Layer | Role |
|-------|------|
| **Python** `scripts/` + `vendor/fsot_compute.py` | Research oracle, RMSD benches, RCSB I/O |
| **Zig host** `zig build host` | Shipping residual/codon/scalar on Windows |
| **Zig kernel** QEMU Multiboot | Same law freestanding (serial gates) |

Python remains the **metric authority** for product RMSD (freeze 1.16 Å).  
Zig is the **runtime authority** for the pin residual cell — it must match the oracle.

## Compared quantities

| Key | Python | Zig | Status |
|-----|-------:|----:|:------:|
| `P_NEW` | 0.300302276670 | 0.300302276670 | PASS |
| `r_bond` (Physical_Chemistry) | 1.100304945540 | 1.100304945540 | PASS |
| `r_clash` (Chemistry) | 1.122488513106 | 1.122488513106 | PASS |
| `r_anchor` (Biochemistry) | 1.091958932350 | 1.091958932350 | PASS |
| domain \|S\| values | pin table | pin table | PASS (~1e−13) |
| `multi_top_k` / `multi_power` | 4 / φ⁶ | 4 / φ⁶ | PASS |
| ATG codon trit | [+1,−1,+1] AA=M | same | PASS |
| DNA→AA fragment | `MALWMRLLPLLALLALWGPDPA` | same | PASS |
| 1-step residual physics bond | 3.98239512079 Å | same | PASS |
| neuro sample S | 1.56177176444 | same | PASS |

## QEMU (freestanding)

```
FSOT_GENETICS_KERNEL pin=D1D38A
FSOT_TRIT PASS
FSOT_CODON PASS
FSOT_SCALAR PASS
FSOT_RESIDUAL PASS  r_bond=1.100 r_clash=1.122 r_anchor=1.091
FSOT_PRODUCT_CELL PASS aa_len=22
FSOT_STAGE_GENETICS_OK
```

Kernel prints 3 decimals on serial; full f64 parity is the **host** gate.

## What is *not* in Zig yet (honest)

- RCSB multi-template fetch / ensemble Cα build (still Python)
- Full product fold RMSD vs AF bench (Python freeze)
- MSA packing fuse

Those stay on the research path until measured tables are baked into the cell.
The **law** (scalar + residual + codon) is not screwed up.

## Re-run

```powershell
python scripts/parity_zig_python.py
cd zig; zig build host; .\run_qemu.ps1
```
