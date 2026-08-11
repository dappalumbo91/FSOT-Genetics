#!/usr/bin/env python3
"""Build a self-contained field console (HTML) from frozen data artifacts.

No framework. No backend. Open field/console.html in a browser.
Optionally serve: python scripts/serve_field_console.py

Usage:
  python scripts/build_field_console.py
  python scripts/build_field_console.py --open
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import webbrowser
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "field"
OUT_HTML = OUT_DIR / "console.html"


def _load(path: Path) -> dict | list | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _read_text(path: Path, max_chars: int = 400_000) -> str | None:
    if not path.exists():
        return None
    t = path.read_text(encoding="utf-8", errors="replace")
    return t[:max_chars]


def git_stamp() -> dict[str, str | bool]:
    """Return short/full SHA, dirty flag, branch (best-effort)."""
    out: dict[str, str | bool] = {
        "sha": "unknown",
        "sha_full": "unknown",
        "branch": "unknown",
        "dirty": False,
    }

    def _run(args: list[str]) -> str | None:
        try:
            r = subprocess.run(
                args,
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
            if r.returncode != 0:
                return None
            return (r.stdout or "").strip()
        except Exception:
            return None

    full = _run(["git", "rev-parse", "HEAD"])
    short = _run(["git", "rev-parse", "--short", "HEAD"])
    branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
    dirty_s = _run(["git", "status", "--porcelain"])
    if full:
        out["sha_full"] = full
    if short:
        out["sha"] = short
    if branch:
        out["branch"] = branch
    if dirty_s is not None:
        out["dirty"] = bool(dirty_s.strip())
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--open", action="store_true", help="Open in default browser")
    args = ap.parse_args(argv)

    stamp = git_stamp()

    product = _load(ROOT / "data" / "product_vs_alphafold.json") or {}
    parity = _load(ROOT / "data" / "parity_zig_python.json") or {}
    medical = _load(ROOT / "data" / "medical_variant_panel.json") or {}
    stress = _load(ROOT / "data" / "medical_stress_suite.json") or {}
    freeze_note = _read_text(ROOT / "docs" / "PRODUCT_FREEZE.md") or ""

    # Prefer generated ubiquitin model; fall back to sample PDB
    pdb_text = (
        _read_text(ROOT / "predictions" / "ubq_fsot.pdb")
        or _read_text(ROOT / "data" / "pdb_samples" / "1UBQ.pdb")
        or ""
    )

    ps = product.get("summary") or {}
    ms = medical.get("summary") or {}
    pr = (parity.get("python") or {}) if isinstance(parity, dict) else {}
    results = product.get("results") or []
    genes = medical.get("genes") or []

    # Flatten drivers for table
    drivers = []
    for g in genes:
        for d in g.get("drivers") or []:
            drivers.append(
                {
                    "gene": g.get("symbol"),
                    "pos": d.get("pos"),
                    "change": f"{d.get('wt','')}{d.get('pos','')}{d.get('mut','')}",
                    "call": d.get("call"),
                    "percentile": d.get("impact_percentile"),
                    "conservation": d.get("conservation"),
                }
            )

    pin_meta = _load(ROOT / "vendor" / "fsot_compute_AUTHORITY_PIN.json") or {}
    pin_sha = (pin_meta.get("authority_sha256") or "D1D38A")[:12]

    payload = {
        "built_at": datetime.now(timezone.utc).isoformat(),
        "pin": "D1D38A",
        "pin_sha256_prefix": pin_sha,
        "git": stamp,
        "free_parameters": 0,
        "product": {
            "summary": ps,
            "results": results,
            "generated_at": product.get("generated_at"),
        },
        "parity": {
            "status": parity.get("status"),
            "generated_at": parity.get("generated_at"),
            "r_bond": pr.get("r_bond"),
            "r_clash": pr.get("r_clash"),
            "r_anchor": pr.get("r_anchor"),
            "P_NEW": pr.get("P_NEW"),
        },
        "medical": {
            "summary": ms,
            "drivers": drivers[:80],
            "generated_at": medical.get("generated_at"),
            "n_genes": ms.get("n_genes"),
        },
        "stress_summary": (stress.get("summary") if isinstance(stress, dict) else None),
        "freeze_excerpt": freeze_note[:1200],
        "pdb_name": "ubq_fsot.pdb" if (ROOT / "predictions" / "ubq_fsot.pdb").exists() else "1UBQ.pdb",
    }

    # Escape for embedding in JS
    data_js = json.dumps(payload, ensure_ascii=True)
    pdb_js = json.dumps(pdb_text)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>FSOT-Genetics Field Console</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  :root {{
    --bg: #0b0f14;
    --panel: #121a22;
    --panel2: #18222d;
    --border: #243140;
    --text: #e7eef6;
    --muted: #8fa3b8;
    --accent: #3d9cf0;
    --good: #3ecf8e;
    --warn: #e6b84d;
    --bad: #e85d5d;
    --mono: "Cascadia Code", "SF Mono", Consolas, monospace;
    --sans: "Segoe UI", system-ui, -apple-system, sans-serif;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; background: var(--bg); color: var(--text);
    font-family: var(--sans); line-height: 1.45;
  }}
  header {{
    padding: 1.25rem 1.5rem;
    border-bottom: 1px solid var(--border);
    background: linear-gradient(180deg, #101820 0%, var(--bg) 100%);
    display: flex; flex-wrap: wrap; gap: 1rem; align-items: flex-end;
    justify-content: space-between;
  }}
  header h1 {{ margin: 0; font-size: 1.35rem; font-weight: 650; letter-spacing: 0.02em; }}
  header .sub {{ color: var(--muted); font-size: 0.9rem; margin-top: 0.25rem; }}
  .badge {{
    display: inline-flex; align-items: center; gap: 0.35rem;
    padding: 0.35rem 0.7rem; border-radius: 999px;
    font-size: 0.8rem; font-weight: 600; border: 1px solid var(--border);
    background: var(--panel);
  }}
  .badge.good {{ color: var(--good); border-color: #245c44; }}
  .badge.bad {{ color: var(--bad); border-color: #5c2424; }}
  .badge.neutral {{ color: var(--accent); }}
  main {{
    max-width: 1200px; margin: 0 auto; padding: 1.25rem 1.5rem 3rem;
    display: grid; gap: 1rem;
  }}
  .hero {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 0.75rem;
  }}
  .metric {{
    background: var(--panel); border: 1px solid var(--border); border-radius: 12px;
    padding: 1rem 1.1rem;
  }}
  .metric .label {{ color: var(--muted); font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.06em; }}
  .metric .value {{ font-size: 1.75rem; font-weight: 700; margin-top: 0.2rem; font-variant-numeric: tabular-nums; }}
  .metric .hint {{ color: var(--muted); font-size: 0.8rem; margin-top: 0.2rem; }}
  .grid2 {{
    display: grid; grid-template-columns: 1.1fr 0.9fr; gap: 1rem;
  }}
  @media (max-width: 900px) {{ .grid2 {{ grid-template-columns: 1fr; }} }}
  section.card {{
    background: var(--panel); border: 1px solid var(--border); border-radius: 12px;
    padding: 1rem 1.15rem; overflow: hidden;
  }}
  section.card h2 {{
    margin: 0 0 0.75rem; font-size: 1rem; font-weight: 650;
    display: flex; align-items: center; justify-content: space-between; gap: 0.5rem;
  }}
  table {{
    width: 100%; border-collapse: collapse; font-size: 0.88rem;
  }}
  th, td {{
    text-align: left; padding: 0.45rem 0.5rem; border-bottom: 1px solid var(--border);
    font-variant-numeric: tabular-nums;
  }}
  th {{ color: var(--muted); font-weight: 600; font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.04em; }}
  tr:hover td {{ background: var(--panel2); }}
  .call-ok {{ color: var(--good); font-weight: 600; }}
  .call-warn {{ color: var(--warn); font-weight: 600; }}
  .call-bad {{ color: var(--bad); font-weight: 600; }}
  #viewer {{
    width: 100%; height: 380px; border-radius: 8px; background: #0a1016;
    border: 1px solid var(--border); position: relative;
  }}
  pre.cmd {{
    background: #0a1016; border: 1px solid var(--border); border-radius: 8px;
    padding: 0.85rem 1rem; overflow-x: auto; font-family: var(--mono); font-size: 0.82rem;
    color: #c5d6e8; margin: 0.4rem 0 0.8rem;
  }}
  .tabs {{ display: flex; gap: 0.4rem; flex-wrap: wrap; margin-bottom: 0.75rem; }}
  .tab {{
    background: var(--panel2); border: 1px solid var(--border); color: var(--muted);
    border-radius: 8px; padding: 0.35rem 0.7rem; cursor: pointer; font-size: 0.85rem;
  }}
  .tab.active {{ color: var(--text); border-color: var(--accent); background: #1a2a3a; }}
  .muted {{ color: var(--muted); font-size: 0.85rem; }}
  footer {{
    max-width: 1200px; margin: 0 auto; padding: 0 1.5rem 2rem;
    color: var(--muted); font-size: 0.8rem;
  }}
  a {{ color: var(--accent); }}
  .barwrap {{ display: flex; align-items: center; gap: 0.5rem; }}
  .bar {{
    height: 8px; border-radius: 4px; background: #1e2a36; flex: 1; overflow: hidden;
  }}
  .bar > i {{ display: block; height: 100%; background: var(--accent); }}
  .bar > i.good {{ background: var(--good); }}
  .bar > i.warn {{ background: var(--warn); }}
</style>
</head>
<body>
<header>
  <div>
    <h1>FSOT-Genetics · Field Console</h1>
    <div class="sub">Zero free parameters · pin <code id="pin">D1D38A</code> · law S = K(T1+T2+T3)</div>
    <div class="sub" id="git-line" style="margin-top:0.35rem;font-family:var(--mono);font-size:0.8rem;"></div>
  </div>
  <div style="display:flex; gap:0.5rem; flex-wrap:wrap;">
    <span class="badge" id="badge-parity">parity …</span>
    <span class="badge neutral" id="badge-git">git …</span>
    <span class="badge neutral" id="badge-params">0 free params</span>
    <span class="badge neutral" id="badge-built">built …</span>
  </div>
</header>

<main>
  <div class="hero" id="hero"></div>

  <div class="grid2">
    <section class="card">
      <h2>Structure product vs AlphaFold <span class="muted" id="prod-n"></span></h2>
      <div style="overflow-x:auto;">
        <table>
          <thead>
            <tr><th>Protein</th><th>AF Å</th><th>Product Å</th><th>Δ vs AF</th><th>Regime</th></tr>
          </thead>
          <tbody id="prod-rows"></tbody>
        </table>
      </div>
      <p class="muted" style="margin-top:0.75rem;">Product = multi-template measured Cα + residual-weighted physics. Bulk orphan path is separate (~11–14 Å ceiling).</p>
    </section>

    <section class="card">
      <h2>Structure viewer <span class="muted" id="pdb-label"></span></h2>
      <div id="viewer"></div>
      <p class="muted" style="margin-top:0.6rem;">Cα ribbon · rotate drag · scroll zoom. Offline if 3Dmol CDN blocked: still use tables.</p>
    </section>
  </div>

  <section class="card">
    <h2>Medical variant panel <span class="muted" id="med-sum"></span></h2>
    <div style="overflow-x:auto; max-height: 360px; overflow-y:auto;">
      <table>
        <thead>
          <tr><th>Gene</th><th>Variant</th><th>Call</th><th>Impact %ile</th><th>Conservation</th></tr>
        </thead>
        <tbody id="var-rows"></tbody>
      </table>
    </div>
  </section>

  <section class="card">
    <h2>Pin residual (Zig ≡ Python)</h2>
    <div class="hero" id="residual-hero" style="margin-bottom:0.5rem;"></div>
    <p class="muted">Runtime cell must match research oracle. Re-check: <code>python scripts/parity_zig_python.py</code></p>
  </section>

  <section class="card">
    <h2>Field commands</h2>
    <div class="tabs" id="cmd-tabs">
      <button class="tab active" data-tab="predict">Predict</button>
      <button class="tab" data-tab="parity">Parity / metal</button>
      <button class="tab" data-tab="medical">Medical</button>
      <button class="tab" data-tab="console">This console</button>
    </div>
    <div id="cmd-body"></div>
  </section>

  <section class="card">
    <h2>Honest limits (do not oversell)</h2>
    <ul class="muted" style="margin:0; padding-left:1.2rem;">
      <li>Product path needs homolog coverage — not a trained AlphaFold replacement for orphans.</li>
      <li>Bulk single-sequence median ~11–14 Å is a known information ceiling.</li>
      <li>Variant calls are conservation × substitution specificity (0 trained weights), not a full clinical ACMG engine.</li>
      <li>Visual console is a read-only field dashboard; production runs stay CLI / Zig host / QEMU gates.</li>
    </ul>
  </section>
</main>

<footer>
  Built from local <code>data/*.json</code> · regenerate with <code>python scripts/build_field_console.py</code><br/>
  FSOT-Genetics · github.com/dappalumbo91/FSOT-Genetics · pin D1D38A ·
  <span id="foot-git"></span>
</footer>

<script>
const DATA = {data_js};
const PDB = {pdb_js};

function fmt(x, d=2) {{
  if (x === null || x === undefined || Number.isNaN(x)) return "—";
  return Number(x).toFixed(d);
}}

function callClass(c) {{
  if (!c) return "";
  const u = String(c).toUpperCase();
  if (u.includes("DAMAGING") && !u.includes("UNLIKELY")) return "call-ok";
  if (u.includes("UNCERTAIN") || u.includes("BENIGN")) return "call-warn";
  return "";
}}

// Hero metrics
(function() {{
  const s = DATA.product.summary || {{}};
  const m = DATA.medical.summary || {{}};
  const p = DATA.parity || {{}};
  const hero = document.getElementById("hero");
  const items = [
    {{ label: "Product median", value: fmt(s.fsot_product_median_A) + " Å", hint: "vs AF " + fmt(s.alphafold_median_A) + " Å" }},
    {{ label: "Within 1.5 Å of AF", value: (s.product_within_1p5A_of_alphafold ?? "—") + "/" + (s.n ?? "—"), hint: "sub-2 Å: " + (s.product_sub2A ?? "—") }},
    {{ label: "Variant drivers", value: (m.n_likely_damaging ?? "—") + "/" + (m.n_drivers ?? "—"), hint: "LIKELY DAMAGING · " + fmt((m.driver_recall_at_75pct||0)*100, 0) + "% recall" }},
    {{ label: "Parity gate", value: p.status || "—", hint: "Zig host vs Python pin" }},
  ];
  hero.innerHTML = items.map(it => `
    <div class="metric">
      <div class="label">${{it.label}}</div>
      <div class="value">${{it.value}}</div>
      <div class="hint">${{it.hint}}</div>
    </div>`).join("");

  const b = document.getElementById("badge-parity");
  if ((p.status || "").toUpperCase() === "PASS") {{
    b.textContent = "parity PASS";
    b.className = "badge good";
  }} else {{
    b.textContent = "parity " + (p.status || "unknown");
    b.className = "badge bad";
  }}
  document.getElementById("badge-built").textContent = "built " + (DATA.built_at || "").slice(0,19).replace("T"," ") + "Z";
  const g = DATA.git || {{}};
  const dirty = g.dirty ? " dirty" : "";
  const sha = g.sha || "unknown";
  const branch = g.branch || "?";
  const gitBadge = document.getElementById("badge-git");
  gitBadge.textContent = "git " + sha + dirty;
  if (g.dirty) gitBadge.className = "badge bad";
  else gitBadge.className = "badge neutral";
  const gitLine = document.getElementById("git-line");
  if (gitLine) {{
    gitLine.textContent =
      "branch " + branch + " · " + (g.sha_full || sha) + (g.dirty ? " · WORKTREE DIRTY" : " · clean") +
      " · pin " + (DATA.pin_sha256_prefix || "D1D38A") + "…";
  }}
  const foot = document.getElementById("foot-git");
  if (foot) foot.textContent = "git " + sha + (g.dirty ? " (dirty)" : "");
}})();

// Product table
(function() {{
  const rows = DATA.product.results || [];
  document.getElementById("prod-n").textContent = rows.length ? (rows.length + " proteins") : "";
  const tb = document.getElementById("prod-rows");
  tb.innerHTML = rows.map(r => {{
    const af = r.alphafold_rmsd_A, pr = r.fsot_product_rmsd_A;
    const d = (af != null && pr != null) ? (pr - af) : null;
    return `<tr>
      <td>${{r.name || "—"}}</td>
      <td>${{fmt(af)}}</td>
      <td><strong>${{fmt(pr)}}</strong></td>
      <td>${{d==null ? "—" : ((d>=0?"+":"") + fmt(d))}}</td>
      <td class="muted">${{r.regime || "—"}}</td>
    </tr>`;
  }}).join("") || `<tr><td colspan="5" class="muted">No product_vs_alphafold.json — run bench_product_vs_af.py</td></tr>`;
}})();

// Variants
(function() {{
  const m = DATA.medical.summary || {{}};
  document.getElementById("med-sum").textContent =
    (m.n_likely_damaging!=null) ? (`${{m.n_likely_damaging}}/${{m.n_drivers}} LIKELY DAMAGING · ${{m.n_genes}} genes`) : "";
  const tb = document.getElementById("var-rows");
  const rows = DATA.medical.drivers || [];
  tb.innerHTML = rows.map(d => {{
    const pct = d.percentile;
    const w = Math.max(0, Math.min(100, pct || 0));
    return `<tr>
      <td>${{d.gene || "—"}}</td>
      <td><code>${{d.change || "—"}}</code></td>
      <td class="${{callClass(d.call)}}">${{d.call || "—"}}</td>
      <td>
        <div class="barwrap">
          <span style="width:3rem">${{fmt(pct,1)}}</span>
          <div class="bar"><i class="${{w>=75?"good":"warn"}}" style="width:${{w}}%"></i></div>
        </div>
      </td>
      <td>${{fmt((d.conservation||0)*100,1)}}%</td>
    </tr>`;
  }}).join("") || `<tr><td colspan="5" class="muted">No medical_variant_panel.json</td></tr>`;
}})();

// Residuals
(function() {{
  const p = DATA.parity || {{}};
  const el = document.getElementById("residual-hero");
  const items = [
    {{ label: "r_bond", value: fmt(p.r_bond, 6), hint: "Physical_Chemistry" }},
    {{ label: "r_clash", value: fmt(p.r_clash, 6), hint: "Chemistry" }},
    {{ label: "r_anchor", value: fmt(p.r_anchor, 6), hint: "Biochemistry" }},
    {{ label: "P_NEW", value: fmt(p.P_NEW, 6), hint: "pin factor" }},
  ];
  el.innerHTML = items.map(it => `
    <div class="metric">
      <div class="label">${{it.label}}</div>
      <div class="value" style="font-size:1.25rem">${{it.value}}</div>
      <div class="hint">${{it.hint}}</div>
    </div>`).join("");
}})();

// Commands
const CMDS = {{
  predict: `python scripts/fsot_predict.py --id 1UBQ --pdb-out model.pdb
python scripts/bench_product_vs_af.py`,
  parity: `python scripts/parity_zig_python.py
cd zig
zig build host
.\\run_qemu.ps1`,
  medical: `python scripts/run_medical_variant_panel.py
python scripts/dna_variant_effect.py --help
python scripts/run_medical_stress_suite.py`,
  console: `python scripts/build_field_console.py --open
python scripts/serve_field_console.py`,
}};
(function() {{
  const body = document.getElementById("cmd-body");
  function show(tab) {{
    body.innerHTML = `<pre class="cmd">${{CMDS[tab] || ""}}</pre>`;
    document.querySelectorAll(".tab").forEach(t => t.classList.toggle("active", t.dataset.tab === tab));
  }}
  document.querySelectorAll(".tab").forEach(t => t.addEventListener("click", () => show(t.dataset.tab)));
  show("predict");
}})();

// 3D viewer
(function() {{
  document.getElementById("pdb-label").textContent = DATA.pdb_name || "";
  if (!PDB || !window.$3Dmol) {{
    document.getElementById("viewer").innerHTML =
      "<div style='padding:1rem;color:#8fa3b8'>3D viewer unavailable (CDN or empty PDB). Tables still work.</div>";
    return;
  }}
  const element = document.getElementById("viewer");
  const config = {{ backgroundColor: "0a1016" }};
  const viewer = $3Dmol.createViewer(element, config);
  viewer.addModel(PDB, "pdb");
  viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }}, stick: {{ hidden: true }} }});
  // Cα-only models look better as spheres+line
  viewer.setStyle({{ atom: "CA" }}, {{ sphere: {{ radius: 0.6, color: "spectrum" }}, line: {{ color: "white" }} }});
  viewer.zoomTo();
  viewer.render();
}})();
</script>
</body>
</html>
"""

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(html, encoding="utf-8")
    # sidecar JSON for tools
    (OUT_DIR / "console_data.json").write_text(
        json.dumps(payload, indent=2)[:50_000], encoding="utf-8"
    )
    print(f"Wrote {OUT_HTML}")
    print(f"  product median={ps.get('fsot_product_median_A')}")
    print(f"  parity={parity.get('status')}")
    print(f"  medical drivers={ms.get('n_likely_damaging')}/{ms.get('n_drivers')}")
    if args.open:
        webbrowser.open(OUT_HTML.resolve().as_uri())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
