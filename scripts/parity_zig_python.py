#!/usr/bin/env python3
"""Zig product host vs Python pin oracle — hard parity gate.

Fails if residual law, seeds, codon ATG, DNA→AA fragment, or one-step
residual physics diverge beyond float tolerances.

Usage (from repo root or scripts/):
  python scripts/parity_zig_python.py
"""
from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from full_scalar_law import P_NEW, residual_scale  # noqa: E402
from residual_physics_refine import residual_physics_relax  # noqa: E402

ZIG = ROOT / "zig"
TOL_STRICT = 1e-9  # seeds / residual must match f64 path tightly
TOL_S = 5e-6  # domain S: Zig f64 vs Python mpmath→float
TOL_PHYS = 1e-9

DNA = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCC"

# Standard genetic code (match codon.zig)
_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def dna_to_aa(dna: str) -> str:
    out = []
    for i in range(0, len(dna) - 2, 3):
        aa = _CODON.get(dna[i : i + 3].upper(), "X")
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


def primary_atg() -> tuple[int, int, int]:
    # A,G → +1; C,T → −1
    def b(x: str) -> int:
        return 1 if x in "AGag" else (-1 if x in "CTct" else 0)

    return b("A"), b("T"), b("G")


def python_oracle() -> dict[str, float | int | str]:
    phi = float(fc.PHI)
    domains = {
        "Physical_Chemistry": abs(float(fc.domain_scalar("Physical_Chemistry"))),
        "Chemistry": abs(float(fc.domain_scalar("Chemistry"))),
        "Biochemistry": abs(float(fc.domain_scalar("Biochemistry"))),
    }
    # signed S for report
    s_pc = float(fc.domain_scalar("Physical_Chemistry"))
    s_ch = float(fc.domain_scalar("Chemistry"))
    s_bc = float(fc.domain_scalar("Biochemistry"))
    r_bond = residual_scale(s_pc)
    r_clash = residual_scale(s_ch)
    r_anchor = residual_scale(s_bc)

    # one-step residual physics (same as Zig residualPhysicsParitySample)
    X0 = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]])
    # residual_physics_relax does many iters — for parity use one manual step
    CA_CA = 3.8
    lr = 0.08
    w_anchor = 0.05 * r_anchor
    X = X0.copy()
    G = w_anchor * (X - X0)
    d = X[1:] - X[:-1]
    L = np.linalg.norm(d, axis=1) + 1e-9
    f = ((L - CA_CA) / L)[:, None] * d * r_bond
    G[:-1] -= f
    G[1:] += f
    X = X - lr * G
    X = X - X.mean(axis=0)
    bond_len = float(np.linalg.norm(X[1] - X[0]))
    end_x = float(X[2, 0])

    atg = primary_atg()
    aa = dna_to_aa(DNA)

    # neuro scalar via vendor (same defaults as Zig computeNeuro)
    # computeNeuro(delta_psi=0.7, hits=1, rho=1)
    si = fc.ScalarInput(
        N=fc.mpf(4),
        P=fc.mpf(3),
        D_eff=fc.mpf(13),
        recent_hits=fc.mpf(1),
        delta_psi=fc.mpf("0.7"),
        delta_theta=fc.mpf(1),
        rho=fc.mpf(1),
        scale=fc.mpf(1),
        amplitude=fc.mpf(1),
        trend_bias=fc.mpf(0),
        observed=True,
    )
    neuro_S = float(fc.compute_scalar(si))
    # clamp like Zig
    neuro_S = max(-3.0, min(3.0, neuro_S))

    return {
        "P_NEW": float(P_NEW),
        "K": float(fc.K),
        "PHI": phi,
        "S_Physical_Chemistry": s_pc,
        "S_Chemistry": s_ch,
        "S_Biochemistry": s_bc,
        "r_bond": r_bond,
        "r_clash": r_clash,
        "r_anchor": r_anchor,
        "multi_top_k": max(2, int(round(phi**3))),
        "multi_power": phi**6,
        "codon_atg0": atg[0],
        "codon_atg1": atg[1],
        "codon_atg2": atg[2],
        "aa_atg": "M",
        "aa_translate": aa,
        "aa_len": len(aa),
        "phys_bond_len": bond_len,
        "phys_end_x": end_x,
        "neuro_S": neuro_S,
        "trit_pair": -1,
        "trit_consensus": 1,
    }


def run_zig_host() -> dict[str, str]:
    # rebuild host so PARITY lines are current
    subprocess.run(
        ["zig", "build", "host"],
        cwd=ZIG,
        check=True,
        capture_output=True,
        text=True,
    )
    # run installed binary (build host already ran it; re-run for capture)
    exe = ZIG / "zig-out" / "bin" / "fsot_genetics_host.exe"
    if not exe.exists():
        exe = ZIG / "zig-out" / "bin" / "fsot_genetics_host"
    proc = subprocess.run([str(exe)], capture_output=True, text=True, check=True)
    text = proc.stdout + proc.stderr
    kv: dict[str, str] = {}
    for line in text.splitlines():
        if line.startswith("PARITY "):
            body = line[len("PARITY ") :]
            # multi key=value on one line
            for part in body.split():
                if "=" in part:
                    k, v = part.split("=", 1)
                    kv[k] = v
    return kv


def cmp_float(name: str, py: float, zig_s: str, tol: float, rows: list) -> bool:
    try:
        zg = float(zig_s)
    except ValueError:
        rows.append((name, py, zig_s, "PARSE", False))
        return False
    ok = abs(py - zg) <= tol or (
        abs(py) > 1e-12 and abs(py - zg) / abs(py) <= tol
    )
    # for values near 1, absolute tol preferred
    ok = abs(py - zg) <= tol
    rows.append((name, py, zg, abs(py - zg), ok))
    return ok


def main() -> int:
    print("=== PYTHON ORACLE ===")
    py = python_oracle()
    for k in sorted(py.keys()):
        print(f"  {k}={py[k]}")

    print("\n=== ZIG HOST ===")
    try:
        zg = run_zig_host()
    except subprocess.CalledProcessError as e:
        print("zig build/run failed:", e)
        if e.stdout:
            print(e.stdout)
        if e.stderr:
            print(e.stderr)
        return 2
    for k in sorted(zg.keys()):
        print(f"  {k}={zg[k]}")

    print("\n=== COMPARE ===")
    rows: list = []
    all_ok = True

    float_checks = [
        ("P_NEW", TOL_STRICT),
        ("K", TOL_S),  # f64 seed table vs mpmath
        ("PHI", TOL_STRICT),
        ("S_Physical_Chemistry", TOL_S),
        ("S_Chemistry", TOL_S),
        ("S_Biochemistry", TOL_S),
        ("r_bond", TOL_S),
        ("r_clash", TOL_S),
        ("r_anchor", TOL_S),
        ("multi_power", TOL_S),
        ("phys_bond_len", TOL_PHYS),
        ("phys_end_x", TOL_PHYS),
        ("neuro_S", 1e-4),  # clamp + f64 vs mpmath
    ]
    for name, tol in float_checks:
        if name not in zg:
            print(f"  MISSING zig key: {name}")
            all_ok = False
            continue
        if not cmp_float(name, float(py[name]), zg[name], tol, rows):
            all_ok = False

    # ints
    for name in ("multi_top_k", "codon_atg0", "codon_atg1", "codon_atg2", "aa_len", "trit_pair", "trit_consensus"):
        if name not in zg:
            print(f"  MISSING zig key: {name}")
            all_ok = False
            continue
        pv, zv = int(py[name]), int(float(zg[name]))
        ok = pv == zv
        rows.append((name, pv, zv, 0, ok))
        if not ok:
            all_ok = False

    # strings
    for name in ("aa_atg", "aa_translate"):
        if name not in zg:
            print(f"  MISSING zig key: {name}")
            all_ok = False
            continue
        pv, zv = str(py[name]), zg[name]
        ok = pv == zv
        rows.append((name, pv, zv, 0, ok))
        if not ok:
            all_ok = False

    print(f"{'key':<28}{'python':>18}{'zig':>18}{'|d|':>14}  ok")
    print("-" * 82)
    for name, a, b, d, ok in rows:
        mark = "PASS" if ok else "FAIL"
        if isinstance(a, float) or (isinstance(a, (int, float)) and not isinstance(a, bool)):
            try:
                print(f"{name:<28}{float(a):18.12g}{float(b):18.12g}{float(d):14.3e}  {mark}")
            except (TypeError, ValueError):
                print(f"{name:<28}{str(a):>18}{str(b):>18}{str(d):>14}  {mark}")
        else:
            print(f"{name:<28}{str(a):>18}{str(b):>18}{'':>14}  {mark}")

    print("-" * 82)
    if all_ok:
        print("PARITY_GATE PASS — Zig host matches Python pin oracle")
        out = ROOT / "data" / "parity_zig_python.json"
        import json
        from datetime import datetime, timezone

        out.write_text(
            json.dumps(
                {
                    "generated_at": datetime.now(timezone.utc).isoformat(),
                    "status": "PASS",
                    "python": {k: (float(v) if isinstance(v, float) else v) for k, v in py.items()},
                    "zig": zg,
                    "rows": [
                        {
                            "key": n,
                            "python": a if not isinstance(a, float) else a,
                            "zig": b if not isinstance(b, float) else b,
                            "ok": ok,
                        }
                        for n, a, b, d, ok in rows
                    ],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"Wrote {out}")
        return 0
    print("PARITY_GATE FAIL — fix Zig/Python before shipping bare metal")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
