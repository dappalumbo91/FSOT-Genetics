#!/usr/bin/env python3
"""One-button field pack: console + Zig host + pin + freeze + LICENSE → zip.

USB / demo ready. No full repo required on the receiving machine for *viewing*
the console or running the residual host gate.

Usage:
  python scripts/build_field_pack.py
  python scripts/build_field_pack.py --skip-zig   # HTML + docs only
  python scripts/build_field_pack.py --skip-parity

Outputs:
  dist/FSOT-Genetics-field-<sha>.zip
  dist/FSOT-Genetics-field-<sha>/   (unpacked mirror)
  field/MANIFEST.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_field_console import git_stamp, main as build_console  # noqa: E402

DIST = ROOT / "dist"
ZIG = ROOT / "zig"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def build_zig_host() -> Path | None:
    zig = shutil.which("zig")
    if not zig:
        print("WARN: zig not on PATH — pack without host binary")
        return None
    r = subprocess.run(
        [zig, "build", "host"],
        cwd=ZIG,
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        print("WARN: zig build host failed:")
        print(r.stderr or r.stdout)
        return None
    for name in ("fsot_genetics_host.exe", "fsot_genetics_host"):
        p = ZIG / "zig-out" / "bin" / name
        if p.exists():
            return p
    return None


def run_parity() -> dict:
    try:
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "parity_zig_python.py")],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=180,
        )
        status = "PASS" if r.returncode == 0 and "PARITY_GATE PASS" in (r.stdout + r.stderr) else "FAIL"
        return {
            "status": status,
            "returncode": r.returncode,
            "tail": (r.stdout + r.stderr)[-800:],
        }
    except Exception as e:
        return {"status": "ERROR", "error": str(e)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-zig", action="store_true")
    ap.add_argument("--skip-parity", action="store_true")
    args = ap.parse_args(argv)

    stamp = git_stamp()
    sha = str(stamp.get("sha") or "nogit")
    dirty = bool(stamp.get("dirty"))
    pack_name = f"FSOT-Genetics-field-{sha}" + ("-dirty" if dirty else "")
    pack_dir = DIST / pack_name
    zip_path = DIST / f"{pack_name}.zip"

    if pack_dir.exists():
        shutil.rmtree(pack_dir)
    pack_dir.mkdir(parents=True)

    # 1) Rebuild console with SHA stamp
    print("=== build field console ===")
    rc = build_console([])
    if rc != 0:
        print("console build failed", file=sys.stderr)
        return rc

    # 2) Optional parity
    parity_info: dict = {"status": "SKIPPED"}
    if not args.skip_parity:
        print("=== parity gate ===")
        parity_info = run_parity()
        print(f"  parity={parity_info.get('status')}")

    # 3) Zig host
    host_path: Path | None = None
    if not args.skip_zig:
        print("=== zig build host ===")
        host_path = build_zig_host()
        if host_path:
            print(f"  host={host_path}")

    # 4) Assemble pack tree
    files_meta: list[dict] = []

    def add(src: Path, rel: str, required: bool = True) -> None:
        if not src.exists():
            if required:
                raise FileNotFoundError(f"missing required pack file: {src}")
            print(f"  skip missing {src}")
            return
        dst = pack_dir / rel
        _copy(src, dst)
        files_meta.append(
            {
                "path": rel.replace("\\", "/"),
                "sha256": _sha256(src),
                "bytes": src.stat().st_size,
            }
        )

    add(ROOT / "field" / "console.html", "console.html")
    add(ROOT / "field" / "console_data.json", "console_data.json", required=False)
    add(ROOT / "docs" / "PRODUCT_FREEZE.md", "docs/PRODUCT_FREEZE.md")
    add(ROOT / "docs" / "FIELD_READY.md", "docs/FIELD_READY.md")
    add(ROOT / "docs" / "PARITY_ZIG_PYTHON.md", "docs/PARITY_ZIG_PYTHON.md", required=False)
    add(ROOT / "docs" / "BARE_METAL_GENETICS_ROADMAP.md", "docs/BARE_METAL_GENETICS_ROADMAP.md", required=False)
    add(ROOT / "LICENSE", "LICENSE")
    add(ROOT / "vendor" / "fsot_compute_AUTHORITY_PIN.json", "pin/fsot_compute_AUTHORITY_PIN.json")
    add(ROOT / "data" / "product_vs_alphafold.json", "data/product_vs_alphafold.json")
    add(ROOT / "data" / "parity_zig_python.json", "data/parity_zig_python.json", required=False)
    add(ROOT / "data" / "medical_variant_panel.json", "data/medical_variant_panel.json", required=False)
    add(ROOT / "predictions" / "ubq_fsot.pdb", "data/ubq_fsot.pdb", required=False)
    add(ROOT / "data" / "pdb_samples" / "1UBQ.pdb", "data/1UBQ.pdb", required=False)

    if host_path:
        rel_host = f"bin/{host_path.name}"
        add(host_path, rel_host)
        # small runner
        if host_path.suffix.lower() == ".exe":
            (pack_dir / "bin" / "run_host_gate.cmd").write_text(
                "@echo off\r\n"
                "cd /d %~dp0\r\n"
                "fsot_genetics_host.exe\r\n"
                "echo.\r\n"
                "echo Exit %ERRORLEVEL% — expect FSOT_STAGE_GENETICS_OK above\r\n"
                "pause\r\n",
                encoding="utf-8",
            )
        else:
            (pack_dir / "bin" / "run_host_gate.sh").write_text(
                "#!/usr/bin/env bash\n"
                "cd \"$(dirname \"$0\")\"\n"
                "./fsot_genetics_host\n",
                encoding="utf-8",
            )

    readme = f"""# FSOT-Genetics Field Pack

**Git:** `{sha}` ({stamp.get("sha_full")}){"  **DIRTY worktree**" if dirty else ""}  
**Branch:** `{stamp.get("branch")}`  
**Built:** {datetime.now(timezone.utc).isoformat()}  
**Pin:** D1D38A  
**Free parameters:** 0  

## What's inside

| Path | Purpose |
|------|---------|
| `console.html` | Visual field console (open in browser) |
| `bin/fsot_genetics_host*` | Residual / codon / scalar host gate |
| `bin/run_host_gate.cmd` | Double-click host gate (Windows) |
| `data/` | Frozen product, parity, medical snapshots |
| `docs/PRODUCT_FREEZE.md` | Scoreboard freeze |
| `docs/FIELD_READY.md` | Ops checklist |
| `pin/fsot_compute_AUTHORITY_PIN.json` | Authority pin certificate |
| `LICENSE` | License |

## Quick start (offline-ish)

1. Open `console.html` in Chrome/Edge/Firefox.  
   - Tables work fully offline.  
   - 3D viewer needs network once for 3Dmol CDN (or skip 3D).
2. Run host gate:
   - Windows: double-click `bin/run_host_gate.cmd`  
   - Expect: `FSOT_STAGE_GENETICS_OK` and residual r_bond≈1.100
3. Read claim boundaries in `docs/PRODUCT_FREEZE.md` — do not oversell bulk orphans.

## What this pack is NOT

- Full RCSB multi-template prediction pipeline (needs repo + network).
- AlphaFold replacement for de-novo folds.
- Clinical ACMG report generator.

## Rebuild from full repo

```powershell
python scripts/build_field_pack.py
```

Repo: https://github.com/dappalumbo91/FSOT-Genetics  
"""
    (pack_dir / "README.md").write_text(readme, encoding="utf-8")
    files_meta.append(
        {
            "path": "README.md",
            "sha256": _sha256(pack_dir / "README.md"),
            "bytes": (pack_dir / "README.md").stat().st_size,
        }
    )

    product = {}
    pjson = ROOT / "data" / "product_vs_alphafold.json"
    if pjson.exists():
        product = json.loads(pjson.read_text(encoding="utf-8")).get("summary") or {}

    manifest = {
        "pack_name": pack_name,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git": stamp,
        "pin": "D1D38A",
        "free_parameters": 0,
        "parity": parity_info,
        "product_summary": product,
        "files": files_meta,
        "how_to_open_console": "Open console.html in a browser",
        "how_to_run_host": "bin/run_host_gate.cmd or bin/fsot_genetics_host*",
    }
    man_path = pack_dir / "MANIFEST.json"
    man_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    # also drop a copy under field/ for repo visibility
    (ROOT / "field" / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    # 5) Zip
    print("=== zip ===")
    DIST.mkdir(parents=True, exist_ok=True)
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for f in pack_dir.rglob("*"):
            if f.is_file():
                zf.write(f, arcname=str(f.relative_to(pack_dir.parent)))
    zip_sha = _sha256(zip_path)

    print(f"Pack dir: {pack_dir}")
    print(f"Zip:      {zip_path}")
    print(f"Zip SHA256: {zip_sha}")
    print(f"Files: {len(files_meta)}")
    print(f"Parity: {parity_info.get('status')}")
    if dirty:
        print("NOTE: worktree dirty — pack name includes -dirty")
    if parity_info.get("status") == "FAIL":
        print("WARN: parity FAIL — pack still written; do not ship as green", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
