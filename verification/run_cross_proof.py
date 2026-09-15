#!/usr/bin/env python3
"""FSOT-Genetics multi-prover gauntlet.

Independent of FSOT-2.1-Lean hub overall_ok. Every layer re-runs on THIS tree.

  python verification/run_cross_proof.py

Layers (use if present; missing optional tools SKIP, present-and-fail FAIL):
  python verify_cross · catalog obligations · SMT python · Z3 · Lean lake
  Coq · Isabelle · F* · Rust kernel · TLA+ TLC
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verification"))
REPORT = ROOT / "data" / "cross_proof_report.json"
HUB = Path(r"C:\Users\damia\Desktop\FSOT-2.1-Lean")
ISA_HOME = Path(r"C:\Users\damia\Desktop\Isabelle2025-2")
FSTAR_HOME = Path(os.environ.get("FSTAR_HOME") or r"I:\FSOT-Physical-Archive\07_Portable-Toolchain\fstar")


def run(cmd: list[str], cwd: Path | None = None, timeout: int = 180, env: dict | None = None) -> tuple[int, str]:
    try:
        p = subprocess.run(
            cmd,
            cwd=str(cwd or ROOT),
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
        out = (p.stdout or "") + (p.stderr or "")
        return p.returncode, out[-6000:]
    except FileNotFoundError as e:
        return 127, str(e)
    except subprocess.TimeoutExpired:
        return 124, "timeout"


def which(names: tuple[str, ...], extra: list[Path] | None = None) -> str | None:
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    for path in extra or []:
        if path.is_file():
            return str(path)
    return None


def resolve_isabelle() -> dict | None:
    p = shutil.which("isabelle")
    if p:
        return {"mode": "posix", "tool": p}
    bash = ISA_HOME / "contrib" / "cygwin" / "bin" / "bash.exe"
    sh = ISA_HOME / "bin" / "isabelle"
    if bash.is_file() and sh.is_file():
        return {"mode": "cygwin", "tool": str(sh), "bash": str(bash), "home": str(ISA_HOME)}
    return None


def cygpath(win: Path) -> str:
    resolved = win.resolve()
    drive = resolved.drive.rstrip(":").lower()
    tail = resolved.as_posix().split(":", 1)[-1]
    return f"/cygdrive/{drive}{tail}"


def ensure_tla_jar() -> Path | None:
    local = ROOT / "verification" / "tla" / "tla2tools.jar"
    if local.is_file():
        return local
    for cand in (
        HUB / "tools" / "tla" / "tla2tools.jar",
        Path(r"C:\Users\damia\Desktop\Circuit\tools\tla\tla2tools.jar"),
    ):
        if cand.is_file():
            shutil.copy2(cand, local)
            return local
    return None


def main() -> int:
    layers: list[dict] = []

    def record(name: str, required: bool, rc: int, detail: str) -> None:
        if rc == 127 and not required:
            status = "SKIP"
            ok = True
        elif rc == 0:
            status = "PASS"
            ok = True
        else:
            status = "FAIL"
            ok = False
        layers.append(
            {
                "name": name,
                "required": required,
                "ok": ok,
                "status": status,
                "returncode": rc,
                "detail": detail[-1500:],
            }
        )
        print(f"[{status}] {name}")
        if detail and status == "FAIL":
            print(detail[-1200:])

    # 0. Export catalog from live JSON
    rc, out = run([sys.executable, str(ROOT / "verification" / "export_obligations.py")], timeout=60)
    record("export_obligations", True, rc, out)

    # 1. Engine gate
    rc, out = run([sys.executable, str(ROOT / "scripts" / "verify_cross.py")], timeout=3600)
    record("python_verify_cross", True, rc, out)

    # 2. SMT python replay
    rc, out = run([sys.executable, str(ROOT / "verification" / "smt_replay.py")])
    record("smt_python", True, rc, out)

    # 3. Z3
    z3 = which(
        ("z3", "z3.exe"),
        [HUB / "tools" / "z3" / "z3-4.15.4-x64-win" / "bin" / "z3.exe"],
    )
    if z3:
        rc, out = run([z3, str(ROOT / "verification" / "smt" / "genetics_bounds.smt2")])
        ok = rc == 0 and "sat" in out and "unsat" not in out.split()[:3]
        record("smt_z3", True, 0 if ok else 1, out)
    else:
        record("smt_z3", False, 127, "z3 not found")

    # 4. Lean
    lake = which(("lake", "lake.exe"))
    if lake:
        rc, out = run([lake, "build"], timeout=3600)
        record("lean_lake_build", True, rc, out)
    else:
        record("lean_lake_build", False, 127, "lake not on PATH")

    # 5. Coq / Rocq
    coqc = which(("coqc", "coqc.exe", "rocqc", "rocqc.exe"))
    if coqc:
        rc, out = run([coqc, "-q", "GeneticsSpine.v"], cwd=ROOT / "verification" / "coq", timeout=180)
        record("coq_genetics_spine", True, rc, out)
    else:
        record("coq_genetics_spine", False, 127, "coqc not on PATH")

    # 6. Isabelle
    isa = resolve_isabelle()
    if isa:
        thy = ROOT / "verification" / "isabelle"
        if isa["mode"] == "cygwin":
            cmd = f"cd '{cygpath(Path(isa['home']))}' && bin/isabelle build -D '{cygpath(thy)}' FSOTGenetics"
            rc, out = run([isa["bash"], "--login", "-c", cmd], timeout=900)
        else:
            rc, out = run([isa["tool"], "build", "-D", str(thy), "FSOTGenetics"], timeout=900)
        record("isabelle_genetics_spine", True, rc, out)
    else:
        record("isabelle_genetics_spine", False, 127, "isabelle not found")

    # 7. F*
    fstar = which(
        ("fstar", "fstar.exe"),
        [FSTAR_HOME / "bin" / "fstar.exe"],
    )
    if fstar:
        env = os.environ.copy()
        env["FSTAR_HOME"] = str(FSTAR_HOME)
        fstar_z3 = which(
            ("z3-4.13.3.exe", "z3-4.13.3"),
            [
                FSTAR_HOME / "bin" / "z3-4.13.3.exe",
                FSTAR_HOME / "bin" / "z3.exe",
            ],
        )
        cmd = [fstar]
        if fstar_z3:
            cmd += ["--smt", fstar_z3]
        elif z3:
            cmd += ["--smt", z3, "--z3version", "4.15.4"]
        cmd.append("FSOTGenetics.fst")
        rc, out = run(cmd, cwd=ROOT / "verification" / "fstar", timeout=300, env=env)
        record("fstar_genetics", True, rc, out)
    else:
        record("fstar_genetics", False, 127, "fstar not found")

    # 8. Rust kernel
    cargo = which(("cargo", "cargo.exe"))
    if cargo:
        rc, out = run(
            [cargo, "run", "--quiet", "--release"],
            cwd=ROOT / "verification" / "rust_genetics_kernel",
            timeout=180,
        )
        ok = rc == 0 and "FSOT_GENETICS_RUST_OK" in out
        record("rust_genetics_kernel", True, 0 if ok else rc or 1, out)
    else:
        record("rust_genetics_kernel", False, 127, "cargo not on PATH")

    # 9. TLA+ TLC
    jar = ensure_tla_jar()
    java = which(("java", "java.exe"))
    if java and jar:
        rc, out = run(
            [java, "-cp", str(jar), "tlc2.TLC", "-config", "GeneticsSolve.cfg", "GeneticsSolve.tla"],
            cwd=ROOT / "verification" / "tla",
            timeout=180,
        )
        ok = rc == 0 and "No error has been found" in out
        record("tla_tlc", True, 0 if ok else 1, out)
    else:
        record("tla_tlc", False, 127, "java or tla2tools.jar missing")

    required_fail = [L for L in layers if L["required"] and L["status"] == "FAIL"]
    doc = {
        "application": "FSOT-Genetics",
        "inherits_hub_overall_ok": False,
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3)",
        "layers": layers,
        "required_fail_count": len(required_fail),
        "overall_ok": len(required_fail) == 0,
        "frameworks_passed": [L["name"] for L in layers if L["status"] == "PASS"],
        "frameworks_skipped": [L["name"] for L in layers if L["status"] == "SKIP"],
    }
    REPORT.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    print("=" * 64)
    print(f"overall_ok={doc['overall_ok']}  required_fail={len(required_fail)}  layers={len(layers)}")
    print("passed:", ", ".join(doc["frameworks_passed"]) or "-")
    print("skipped:", ", ".join(doc["frameworks_skipped"]) or "-")
    print("wrote", REPORT)
    return 0 if doc["overall_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
