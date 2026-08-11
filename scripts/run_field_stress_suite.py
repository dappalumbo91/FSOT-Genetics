#!/usr/bin/env python3
"""FSOT-Genetics field stress suite — lab-style end-to-end gate.

Runs every production surface we ship:
  1. Authority pin + cross-verify
  2. Residual law / domain S (Python pin)
  3. Zig host residual gate (if zig present)
  4. Zig ≡ Python parity harness
  5. Product freeze scoreboard (median ≤ 1.16 Å, 10/10 within 1.5 of AF)
  6. Medical variant panel integrity (drivers, free_params=0)
  7. Lab structure assay: ubiquitin Cα product vs experimental 1UBQ
  8. Residual physics one-step geometry sanity
  9. DNA→AA + codon primary law
 10. Field console build + HTML contract (git stamp, hero data, no broken refs)
 11. Field pack zip integrity (MANIFEST hashes)
 12. QEMU kernel optional (if qemu present)

Exit 0 only if all required gates pass.
Write: data/field_stress_suite.json

Usage:
  python scripts/run_field_stress_suite.py
  python scripts/run_field_stress_suite.py --skip-qemu --skip-pack
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

OUT = ROOT / "data" / "field_stress_suite.json"
# Product freeze gate (PRODUCT_FREEZE.md)
PRODUCT_MEDIAN_MAX = 1.165  # allow float noise around 1.16
PRODUCT_WITHIN_1P5_MIN = 10
UBQ_PRODUCT_MAX_A = 2.5  # lab assay soft ceiling (product freeze ~1.78)


class Suite:
    def __init__(self) -> None:
        self.rows: list[dict] = []
        self.t0 = time.perf_counter()

    def check(self, name: str, ok: bool, detail: str = "", *, required: bool = True) -> bool:
        status = "PASS" if ok else ("FAIL" if required else "WARN")
        self.rows.append(
            {
                "name": name,
                "status": status,
                "required": required,
                "detail": detail[:2000],
            }
        )
        mark = {"PASS": "OK  ", "FAIL": "FAIL", "WARN": "WARN"}[status]
        print(f"  {mark}  {name}: {detail}")
        return ok

    def failed_required(self) -> list[dict]:
        return [r for r in self.rows if r["status"] == "FAIL" and r["required"]]


def _run(cmd: list[str], *, cwd: Path | None = None, timeout: int = 300) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=cwd or ROOT,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )


def test_cross_verify(s: Suite) -> None:
    r = _run([sys.executable, str(ROOT / "scripts" / "verify_cross.py")], timeout=120)
    s.check(
        "cross_verify_pin",
        r.returncode == 0,
        f"exit={r.returncode} tail={(r.stdout + r.stderr)[-300:]}",
    )


def test_residual_python(s: Suite) -> None:
    import fsot_compute as fc
    from full_scalar_law import residual_scale
    from residual_physics_refine import residuals_report

    rep = residuals_report()
    r_bond = rep["Physical_Chemistry_residual"]
    r_clash = rep["Chemistry_residual"]
    r_anchor = rep["Biochemistry_residual"]
    s.check("residual_r_bond_band", 1.05 < r_bond < 1.20, f"r_bond={r_bond:.12f}")
    s.check("residual_r_clash_band", 1.05 < r_clash < 1.25, f"r_clash={r_clash:.12f}")
    s.check("residual_r_anchor_band", 1.05 < r_anchor < 1.20, f"r_anchor={r_anchor:.12f}")
    # recompute from domain_scalar must match report
    s_pc = float(fc.domain_scalar("Physical_Chemistry"))
    s.check(
        "residual_recompute_match",
        abs(residual_scale(s_pc) - r_bond) < 1e-12,
        f"S={s_pc:.12f}",
    )
    s.check("p_new_pin", abs(float(fc.P_NEW) - 0.30030227667037146) < 1e-9, f"P_NEW={float(fc.P_NEW)}")


def test_zig_host(s: Suite) -> None:
    zig = shutil.which("zig")
    if not zig:
        s.check("zig_host", False, "zig not on PATH", required=False)
        return
    r = _run([zig, "build", "host"], cwd=ROOT / "zig", timeout=180)
    out = r.stdout + r.stderr
    s.check(
        "zig_host_gate",
        r.returncode == 0 and "FSOT_STAGE_GENETICS_OK" in out,
        f"exit={r.returncode} has_ok={'FSOT_STAGE_GENETICS_OK' in out}",
    )
    # residual lines present
    s.check(
        "zig_host_residuals_printed",
        "r_bond=" in out and "1.100" in out,
        "residual printout",
    )


def test_parity(s: Suite) -> None:
    r = _run([sys.executable, str(ROOT / "scripts" / "parity_zig_python.py")], timeout=180)
    out = r.stdout + r.stderr
    s.check(
        "parity_zig_python",
        r.returncode == 0 and "PARITY_GATE PASS" in out,
        f"exit={r.returncode}",
    )
    pj = ROOT / "data" / "parity_zig_python.json"
    if pj.exists():
        data = json.loads(pj.read_text(encoding="utf-8"))
        s.check("parity_json_status", data.get("status") == "PASS", str(data.get("status")))


def test_product_freeze(s: Suite) -> None:
    path = ROOT / "data" / "product_vs_alphafold.json"
    if not path.exists():
        s.check("product_freeze_file", False, "missing product_vs_alphafold.json")
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    sm = data.get("summary") or {}
    med = sm.get("fsot_product_median_A")
    within = sm.get("product_within_1p5A_of_alphafold")
    n = sm.get("n")
    free = sm.get("free_parameters")
    s.check("product_n_10", n == 10, f"n={n}")
    s.check("product_free_params_0", free == 0, f"free={free}")
    s.check(
        "product_median_gate",
        med is not None and med <= PRODUCT_MEDIAN_MAX,
        f"median={med} max={PRODUCT_MEDIAN_MAX}",
    )
    s.check(
        "product_within_1p5_of_af",
        within is not None and within >= PRODUCT_WITHIN_1P5_MIN,
        f"within={within}",
    )
    # every row has product rmsd
    rows = data.get("results") or []
    ok_rows = all(r.get("fsot_product_rmsd_A") is not None for r in rows)
    s.check("product_all_rows_scored", ok_rows and len(rows) == 10, f"rows={len(rows)}")


def test_medical_panel(s: Suite) -> None:
    path = ROOT / "data" / "medical_variant_panel.json"
    if not path.exists():
        s.check("medical_panel_file", False, "missing medical_variant_panel.json")
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    sm = data.get("summary") or {}
    s.check("medical_free_params", data.get("free_parameters") == 0, "free_parameters")
    s.check(
        "medical_driver_recall",
        (sm.get("n_likely_damaging") or 0) >= 30 and (sm.get("n_drivers") or 0) >= 30,
        f"{sm.get('n_likely_damaging')}/{sm.get('n_drivers')}",
    )
    genes = data.get("genes") or []
    s.check("medical_genes_present", len(genes) >= 6, f"n_genes={len(genes)}")
    # TP53 R175H must exist and be damaging-class
    found = False
    for g in genes:
        if g.get("symbol") != "TP53":
            continue
        for d in g.get("drivers") or []:
            if d.get("pos") == 175 and d.get("wt") == "R" and d.get("mut") == "H":
                found = True
                call = str(d.get("call") or "")
                s.check(
                    "lab_tp53_r175h_call",
                    "DAMAGING" in call.upper(),
                    f"call={call}",
                )
    s.check("lab_tp53_r175h_present", found, "R175H in panel")


def test_lab_ubiquitin_assay(s: Suite) -> None:
    """Classic lab control: model ubiquitin vs experimental 1UBQ Cα."""
    from run_fsot_vs_alphafold_structure import fetch_pdb, kabsch_rmsd
    from run_rcsb_template_holdout import best_template
    from msa_template_fuse import fuse_predict

    cache = Path.home() / ".cache" / "fsot-genetics" / "stress_lab"
    cache.mkdir(parents=True, exist_ok=True)
    hit = fetch_pdb("1UBQ", "A", cache)
    if not hit:
        s.check("lab_ubq_fetch", False, "could not fetch 1UBQ")
        return
    seq, nat = hit
    s.check("lab_ubq_length", len(seq) == 76, f"n={len(seq)}")
    t0 = time.perf_counter()
    tmpl = best_template(seq, "1UBQ", identity_cap=0.95)
    s.check("lab_ubq_template_found", tmpl is not None, f"tmpl={tmpl and tmpl.get('pdb_id')}")
    if tmpl is None:
        return
    prod = fuse_predict(seq, tmpl["model"], None)
    X = prod["ca_coords"]
    rmsd = float(kabsch_rmsd(X, nat))
    dt = time.perf_counter() - t0
    s.check(
        "lab_ubq_product_rmsd",
        rmsd < UBQ_PRODUCT_MAX_A,
        f"rmsd={rmsd:.3f} A  template={tmpl.get('pdb_id')}  regime={prod.get('regime')}  {dt:.1f}s",
    )
    s.check("lab_ubq_free_params", prod.get("free_parameters") == 0, str(prod.get("free_parameters")))
    # residual keys present on product engine
    res = prod.get("residual") or {}
    s.check(
        "lab_ubq_residual_channels",
        "Physical_Chemistry" in res and "Chemistry" in res,
        str(res),
    )


def test_residual_physics_step(s: Suite) -> None:
    from residual_physics_refine import residual_physics_relax, residuals_report

    X0 = np.array([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]])
    X = residual_physics_relax(X0, iters=1)
    bl = float(np.linalg.norm(X[1] - X[0]))
    # after one residual-weighted step, bond should move toward CA_CA=3.8 from 4.0
    s.check("physics_bond_moves_toward_ideal", 3.7 < bl < 4.0, f"bond={bl:.6f}")
    s.check("physics_finite", np.isfinite(X).all(), "coords finite")
    rep = residuals_report()
    s.check("physics_report_keys", "P_NEW" in rep and "Physical_Chemistry_residual" in rep, "keys")


def test_codon_dna(s: Suite) -> None:
    def base_p(b: str) -> int:
        return 1 if b in "AGag" else (-1 if b in "CTctUu" else 0)

    atg = (base_p("A"), base_p("T"), base_p("G"))
    s.check("codon_atg_primary", atg == (1, -1, 1), str(atg))

    # DNA fragment used in zig product cell
    dna = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCC"
    table = {
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
    aa = []
    for i in range(0, len(dna) - 2, 3):
        c = table.get(dna[i : i + 3], "X")
        if c == "*":
            break
        aa.append(c)
    aa_s = "".join(aa)
    s.check(
        "dna_translate_insulin_fragment",
        aa_s == "MALWMRLLPLLALLALWGPDPA" and len(aa_s) == 22,
        aa_s,
    )


def test_field_console(s: Suite) -> None:
    r = _run([sys.executable, str(ROOT / "scripts" / "build_field_console.py")], timeout=60)
    s.check("console_build", r.returncode == 0, f"exit={r.returncode}")
    html_path = ROOT / "field" / "console.html"
    data_path = ROOT / "field" / "console_data.json"
    s.check("console_html_exists", html_path.exists(), str(html_path))
    s.check("console_data_exists", data_path.exists(), str(data_path))
    if not html_path.exists():
        return
    html = html_path.read_text(encoding="utf-8")
    s.check("ui_has_title", "FSOT-Genetics" in html and "Field Console" in html, "title")
    s.check("ui_has_3dmol_or_fallback", "3Dmol" in html or "viewer" in html, "viewer")
    s.check("ui_has_parity_badge", "badge-parity" in html, "parity badge")
    s.check("ui_has_git_badge", "badge-git" in html and "git-line" in html, "git stamp")
    s.check("ui_has_product_table", "prod-rows" in html, "product table")
    s.check("ui_has_variant_table", "var-rows" in html, "variant table")
    s.check("ui_has_honest_limits", "Honest limits" in html or "do not oversell" in html, "limits")
    s.check("ui_charset_utf8", "charset=utf-8" in html.lower() or "charset=\"utf-8\"" in html.lower(), "utf8")
    # DATA embed must parse
    m = re.search(r"const DATA = (\{.*?\});\s*\nconst PDB", html, re.S)
    if not m:
        # fallback looser
        m = re.search(r"const DATA = (\{.+\});", html, re.S)
    if m:
        try:
            data = json.loads(m.group(1))
            s.check("ui_data_json_parse", True, f"keys={list(data.keys())[:8]}")
            s.check(
                "ui_data_git_sha",
                bool((data.get("git") or {}).get("sha")),
                str((data.get("git") or {}).get("sha")),
            )
            s.check(
                "ui_data_product_median",
                (data.get("product") or {}).get("summary", {}).get("fsot_product_median_A") is not None,
                "median present",
            )
            s.check(
                "ui_data_parity_pass",
                (data.get("parity") or {}).get("status") == "PASS",
                str((data.get("parity") or {}).get("status")),
            )
        except json.JSONDecodeError as e:
            s.check("ui_data_json_parse", False, str(e))
    else:
        s.check("ui_data_json_parse", False, "DATA blob not found")

    # console_data sidecar
    if data_path.exists():
        cd = json.loads(data_path.read_text(encoding="utf-8"))
        s.check("console_data_pin", cd.get("pin") == "D1D38A", str(cd.get("pin")))


def test_field_pack(s: Suite, *, skip: bool) -> None:
    if skip:
        s.check("field_pack", True, "skipped", required=False)
        return
    r = _run(
        [sys.executable, str(ROOT / "scripts" / "build_field_pack.py"), "--skip-parity"],
        timeout=300,
    )
    s.check("field_pack_build", r.returncode == 0, f"exit={r.returncode}")
    dist = ROOT / "dist"
    zips = sorted(dist.glob("FSOT-Genetics-field-*.zip"), key=lambda p: p.stat().st_mtime, reverse=True)
    s.check("field_pack_zip_exists", bool(zips), f"count={len(zips)}")
    if not zips:
        return
    zpath = zips[0]
    with zipfile.ZipFile(zpath, "r") as zf:
        names = set(zf.namelist())
        # zip entries are pack_name/console.html
        has_console = any(n.endswith("console.html") for n in names)
        has_manifest = any(n.endswith("MANIFEST.json") for n in names)
        has_license = any(n.endswith("LICENSE") for n in names)
        has_host = any("fsot_genetics_host" in n for n in names)
        has_pin = any("AUTHORITY_PIN" in n for n in names)
        s.check("pack_contains_console", has_console, "console.html")
        s.check("pack_contains_manifest", has_manifest, "MANIFEST.json")
        s.check("pack_contains_license", has_license, "LICENSE")
        s.check("pack_contains_host", has_host, "host binary", required=False)
        s.check("pack_contains_pin", has_pin, "authority pin")
        # verify manifest hashes against zip members
        man_name = next(n for n in names if n.endswith("MANIFEST.json"))
        man = json.loads(zf.read(man_name))
        ok_hash = True
        bad = []
        for ent in man.get("files") or []:
            rel = ent["path"]
            # find in zip
            match = next((n for n in names if n.endswith("/" + rel) or n.endswith(rel)), None)
            if not match:
                # README etc at root of pack
                match = next((n for n in names if n.replace("\\", "/").endswith(rel)), None)
            if not match:
                continue
            raw = zf.read(match)
            h = hashlib.sha256(raw).hexdigest()
            if h != ent.get("sha256"):
                ok_hash = False
                bad.append(rel)
        s.check("pack_manifest_hashes", ok_hash, f"bad={bad[:5]}")


def test_qemu(s: Suite, *, skip: bool) -> None:
    if skip:
        s.check("qemu_kernel", True, "skipped", required=False)
        return
    qemu = shutil.which("qemu-system-x86_64")
    if not qemu and Path(r"C:\Program Files\qemu\qemu-system-x86_64.exe").exists():
        qemu = r"C:\Program Files\qemu\qemu-system-x86_64.exe"
    if not qemu:
        s.check("qemu_kernel", False, "qemu not found", required=False)
        return
    r = _run(
        ["powershell", "-NoProfile", "-File", str(ROOT / "zig" / "run_qemu.ps1")],
        cwd=ROOT / "zig",
        timeout=120,
    )
    out = r.stdout + r.stderr
    s.check(
        "qemu_kernel_gate",
        r.returncode == 0 and "FSOT_STAGE_GENETICS_OK" in out,
        f"exit={r.returncode}",
        required=False,  # soft: environment-dependent
    )


def test_no_junk_probes(s: Suite) -> None:
    junk = [
        "scripts/_probe_multifill.py",
        "scripts/_probe_coverage_select.py",
        "scripts/bulk_refine_v2.py",
        "scripts/_audit_templates.py",
    ]
    present = [j for j in junk if (ROOT / j).exists()]
    s.check("repo_no_junk_probes", not present, f"still present: {present}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-qemu", action="store_true")
    ap.add_argument("--skip-pack", action="store_true")
    ap.add_argument("--skip-lab-net", action="store_true", help="Skip 1UBQ network assay")
    args = ap.parse_args(argv)

    print("=" * 64)
    print("FSOT-Genetics FIELD STRESS SUITE")
    print(f"  root = {ROOT}")
    print("=" * 64)
    s = Suite()

    print("\n[1] Authority / cross-verify")
    test_cross_verify(s)

    print("\n[2] Residual law (Python)")
    test_residual_python(s)

    print("\n[3] Zig host")
    test_zig_host(s)

    print("\n[4] Zig vs Python parity")
    test_parity(s)

    print("\n[5] Product freeze scoreboard")
    test_product_freeze(s)

    print("\n[6] Medical variant panel")
    test_medical_panel(s)

    print("\n[7] Lab assay: ubiquitin product vs 1UBQ")
    if args.skip_lab_net:
        s.check("lab_ubq_product_rmsd", True, "skipped", required=False)
    else:
        try:
            test_lab_ubiquitin_assay(s)
        except Exception as e:
            s.check("lab_ubq_product_rmsd", False, f"exception: {e}")

    print("\n[8] Residual physics step")
    test_residual_physics_step(s)

    print("\n[9] Codon / DNA translate")
    test_codon_dna(s)

    print("\n[10] Field console UI contract")
    test_field_console(s)

    print("\n[11] Field pack")
    test_field_pack(s, skip=args.skip_pack)

    print("\n[12] QEMU kernel (optional)")
    test_qemu(s, skip=args.skip_qemu)

    print("\n[13] Repo hygiene")
    test_no_junk_probes(s)

    elapsed = time.perf_counter() - s.t0
    fails = s.failed_required()
    n_pass = sum(1 for r in s.rows if r["status"] == "PASS")
    n_fail = sum(1 for r in s.rows if r["status"] == "FAIL")
    n_warn = sum(1 for r in s.rows if r["status"] == "WARN")

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "elapsed_s": elapsed,
        "summary": {
            "n": len(s.rows),
            "pass": n_pass,
            "fail": n_fail,
            "warn": n_warn,
            "required_fails": len(fails),
            "overall": "PASS" if not fails else "FAIL",
        },
        "gates": s.rows,
        "product_median_max_A": PRODUCT_MEDIAN_MAX,
        "ubq_product_max_A": UBQ_PRODUCT_MAX_A,
    }
    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("\n" + "=" * 64)
    print(
        f"RESULT {report['summary']['overall']}  "
        f"pass={n_pass} fail={n_fail} warn={n_warn}  {elapsed:.1f}s"
    )
    if fails:
        print("Required failures:")
        for f in fails:
            print(f"  - {f['name']}: {f['detail']}")
    print(f"Wrote {OUT}")
    print("=" * 64)
    return 0 if not fails else 1


if __name__ == "__main__":
    raise SystemExit(main())
