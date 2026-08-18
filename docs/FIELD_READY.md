# Field / production readiness — FSOT-Genetics

**Audience:** you, a collaborator, or a lab machine that needs to *see* and *run* the system without spelunking 80 scripts.

## Design choice (UI)

You do not need a React app or a design agency.

| Approach | Why |
|----------|-----|
| **Static field console** (`field/console.html`) | One page, dark scientific UI, tables + 3D, no npm |
| **CLI front doors** | Real work: predict, parity, medical panel, QEMU |
| **Zig host / kernel** | Production *runtime* law; not the GUI |

GUI shows **frozen gates**. CLI/Zig **does** the science.

---

## Production checklist

### A. Correctness gates (must be green)

| Gate | Command | Pass means |
|------|---------|------------|
| Cross-verify | `python scripts/verify_cross.py` | Lean-style pin consistency |
| Zig ≡ Python | `python scripts/parity_zig_python.py` | Residual/codon law not drifted |
| Zig host | `cd zig && zig build host` | `FSOT_STAGE_GENETICS_OK` |
| QEMU (if installed) | `cd zig && .\run_qemu.ps1` | Serial `FSOT_STAGE_GENETICS_OK` |
| Product freeze | `data/product_vs_alphafold.json` | Median ≤ **0.47 Å** AF (freeze **0.13 Å**; see PRODUCT_FREEZE) |

### B. Field surfaces (what a human sees)

| Surface | How to open | Shows |
|---------|-------------|--------|
| **Field console** | `python scripts/build_field_console.py --open` | Scoreboard, H2H table, variants, 3D, commands |
| Or serve | `python scripts/serve_field_console.py` | Same over `http://127.0.0.1:8765` |
| Medical CLI | `python scripts/fsot_predict.py --id 1UBQ` | Structure + regime |
| Variant panel | `python scripts/run_medical_variant_panel.py` | Driver calls |

### C. Packaging / ops

1. **One-button field pack** — `python scripts/build_field_pack.py`  
   → `dist/FSOT-Genetics-field-<sha>.zip` (console + host exe + pin + freeze + data)  
2. **Version stamp** — console shows git SHA + dirty flag + pin prefix  
3. **Offline** — tables + host gate offline; 3Dmol CDN needs network once  
4. **Network honesty** — full predict/template still needs RCSB/UniRef; pack does not  
5. **Claim boundaries** — PRODUCT_FREEZE + console honest-limits section  

### D. Visual polish principles (so UI does not hurt)

1. **Three numbers first** — product median, parity, medical drivers  
2. **Tables over charts** until data is boringly stable  
3. **Commands copy-pasteable** on the same page as results  
4. **Dark, high-contrast, monospace for science values**  
5. **Never present bulk as the product** — orphan ceiling ~11–14 Å is real; the product is **0.13 Å**. Hiding either number is a marketing failure.  

---

## Recommended daily field ritual

```powershell
python scripts/run_field_stress_suite.py   # full lab stress (required before ship)
python scripts/build_field_console.py --open
python scripts/build_field_pack.py         # USB zip
```

## What “production ready” means here

| Ready | Not ready (yet) |
|-------|-----------------|
| Auditable 0-param law on metal | Hospital EMR integration |
| Reproducible H2H + medical panel | Full ACMG clinical report |
| Visual console for demos / your desk | Polished SaaS multi-tenant UI |
| Clear oversell boundaries | Claiming AF-beating orphans, or quoting 15 Å as if it were the product |

Ship the **console + gates**. Iterate UI only when a real user (you) hits friction.
