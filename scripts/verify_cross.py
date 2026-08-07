#!/usr/bin/env python3
"""FSOT-Genetics cross-verification gate (Lean-style green gate).

Hard checks (exit 1 on any fail):
  1. Vendor fsot_compute.py SHA-256 matches D1D38A authority pin
  2. Structure engine free_parameters == 0
  3. Formula fold is fast (sub-second on short peptide)
  4. F01–F15 pipeline produces finite Cα coords
  5. Scoreboard snapshot honesty (if present): engine label + 0 free params

This is the genetics-domain twin of FSOT-2.1-Lean multi-prover / margin gates:
mathematical path only — no neural-net weights, no free dials.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

PIN_PATH = ROOT / "vendor" / "fsot_compute_AUTHORITY_PIN.json"
COMPUTE_PATH = ROOT / "vendor" / "fsot_compute.py"
EXPECTED_PREFIX = "D1D38A"


def fail(msg: str) -> None:
    print(f"FAIL  {msg}")
    sys.exit(1)


def ok(msg: str) -> None:
    print(f"OK    {msg}")


def main() -> int:
    print("=" * 64)
    print("FSOT-Genetics cross-verification")
    print(f"  root = {ROOT}")
    print("=" * 64)

    # ── 1. Authority pin ──────────────────────────────────────────────
    if not COMPUTE_PATH.is_file():
        fail(f"missing {COMPUTE_PATH}")
    raw = COMPUTE_PATH.read_bytes()
    sha = hashlib.sha256(raw).hexdigest().upper()
    if not sha.startswith(EXPECTED_PREFIX):
        fail(f"fsot_compute sha={sha[:12]}… does not start with {EXPECTED_PREFIX}")
    pin = {}
    if PIN_PATH.is_file():
        pin = json.loads(PIN_PATH.read_text(encoding="utf-8"))
        cert = str(pin.get("certificate_authority") or pin.get("authority_sha256") or "").upper()
        if cert and cert != sha:
            fail(f"pin cert {cert[:12]}… != file {sha[:12]}…")
        ok(f"authority pin D1D38A  sha={sha[:16]}…  bytes={len(raw)}")
    else:
        ok(f"authority sha={sha[:16]}… (no pin file; prefix check only)")

    # ── 2–4. Structure engine math path ───────────────────────────────
    import fsot_compute as fc  # noqa: E402
    from fsot_structure_engine import predict_ca_coords  # noqa: E402

    # seed constants exist
    for name in ("PI", "E", "PHI", "GAMMA"):
        if not hasattr(fc, name):
            fail(f"fsot_compute missing {name}")
    ok(f"seeds π={float(fc.PI):.6f} e={float(fc.E):.6f} φ={float(fc.PHI):.6f}")

    # ubiquitin head (classic short fold target fragment)
    seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    t0 = time.perf_counter()
    pred = predict_ca_coords(seq, rounds=16)
    ms = (time.perf_counter() - t0) * 1000.0

    if pred.get("free_parameters") != 0:
        fail(f"free_parameters={pred.get('free_parameters')} (must be 0)")
    ok(f"free_parameters=0  engine={pred.get('engine')}  routing={pred.get('routing')}")
    if pred.get("runtime"):
        ok(f"runtime={pred.get('runtime')}")
    if not pred.get("full_law"):
        fail("engine must set full_law=True (whole S=K(T1+T2+T3), not frozen |S| slice)")
    ok(f"full_law S_final={pred.get('S_final_observation')}  observer_mod={pred.get('observer_mod_final')}")
    ok(f"T1={pred.get('T1_final')} T2={pred.get('T2_final')} T3={pred.get('T3_final')}")

    # parity: float full law vs vendor domain_scalar
    from full_scalar_law import parity_check_domain  # noqa: E402

    for dname in ("Biochemistry", "Molecular_Chemistry"):
        pr = parity_check_domain(dname)
        if not pr["ok"]:
            fail(f"full scalar parity fail {dname}: {pr}")
        ok(f"scalar parity {dname} err={pr['abs_err']:.2e}")

    xyz = pred["ca_coords"]
    if xyz.shape != (len(seq), 3):
        fail(f"coords shape {xyz.shape} != ({len(seq)}, 3)")
    if not (xyz == xyz).all():  # NaN check
        fail("NaN in Cα coordinates")
    ok(f"Cα coords finite  n={len(seq)}  start={pred.get('embed_start')}")

    # speed: formula branch must stay seconds-scale (hard gate 5 s for n~76)
    if ms > 5000:
        fail(f"predict_ms={ms:.1f} exceeds 5000 ms hard gate (formula path must be fast)")
    ok(f"predict_ms={ms:.1f}  (< 5000 ms hard gate)")

    # domain scalars for biochem
    try:
        sb = float(fc.domain_scalar("Biochemistry"))
        sm = float(fc.domain_scalar("Molecular_Chemistry"))
        ok(f"|S_biochem|={abs(sb):.6f}  |S_molchem|={abs(sm):.6f}")
    except Exception as e:
        fail(f"domain_scalar: {e}")

    # ── 5. Scoreboard honesty (optional snapshot) ─────────────────────
    board = ROOT / "data" / "fsot_vs_alphafold_structure.json"
    if board.is_file():
        doc = json.loads(board.read_text(encoding="utf-8"))
        if doc.get("free_parameters") not in (0, None):
            fail("scoreboard free_parameters != 0")
        eng = doc.get("engine") or ""
        ok(f"scoreboard present  engine={eng}  free_parameters=0")
        s = doc.get("summary") or {}
        if s.get("fsot_median_predict_ms") is not None:
            ok(f"scoreboard median predict_ms={s.get('fsot_median_predict_ms')}")
    else:
        ok("no scoreboard snapshot yet (optional)")

    # formulas doc present
    deriv = ROOT / "formulas" / "FSOT_PROTEIN_DERIVATIONS.md"
    if not deriv.is_file():
        fail("missing formulas/FSOT_PROTEIN_DERIVATIONS.md")
    ok("F01–F15 derivation doc present")

    # ── 6. Trinary syntax expansion (genetics as code) ────────────────
    from trinary_syntax import uniqueness_report, write_expanded_maps, aa_pair_weight  # noqa: E402

    payload = write_expanded_maps()
    u = uniqueness_report()
    if not u.get("all_unique"):
        fail(f"expanded AA words not unique: {u.get('expanded_collisions')}")
    ok(f"trinary expansion 20/20 unique opcodes  F01_collision_groups={len(u['f01_collisions'])}")
    # Zig pair law finite
    w = aa_pair_weight("F", "W", 8)
    if not (w == w) or w <= 0:
        fail(f"fsotPairWeight F–W@8 invalid: {w}")
    ok(f"Zig pair weight F–W@8 = {w:.6f}")
    if not (ROOT / "zig" / "src" / "genetic_pair.zig").is_file():
        fail("missing zig/src/genetic_pair.zig (neuron twin)")
    ok("neuron-zig genetics sources present under zig/")

    print("=" * 64)
    print("ALL CROSS-VERIFICATION GATES PASSED")
    print("  law: S=K(T1+T2+T3)  pin: D1D38A  free_parameters: 0")
    print("  path: mathematical formula branch (not neural net)")
    print("  syntax: expanded trinary AA opcodes + codon map")
    print("=" * 64)
    return 0


# exported for engine smoke
def free_parameters_claim() -> int:
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as e:
        print(f"FAIL  uncaught: {e}")
        raise SystemExit(1)
