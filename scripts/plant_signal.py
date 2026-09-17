#!/usr/bin/env python3
"""Arabidopsis information graph — measured PPIs, not a connectome.

Plants do not have a fly-class EM wiring diagram. Information moves on
experimental protein interactions (IntAct, taxid 3702). Same residual law
as the animal graphs: seed a stimulus class → hops → does the effector
class light?

  python scripts/plant_signal.py

MITAB cache on D:\\FlyWire_Connectome\\plants\\intact (not git).
"""
from __future__ import annotations

import json
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402

_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))
_UA = "FSOT-Genetics plant signal (mailto:local)"

CACHE = Path(r"D:\FlyWire_Connectome\plants\intact\arabidopsis_intact.mitab25.txt")
OUT_GIT = ROOT / "data" / "plant_signal_boot.json"
PSICQUIC = (
    "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/"
    "current/search/query/species:3702"
)

# Reviewed Arabidopsis accessions that sit on measured signaling jobs.
# Not invented. UniProt Swiss-Prot.
CLASSES: dict[str, dict[str, str]] = {
    # light receptors
    "O48963": {"symbol": "PHOT1", "class": "photoreceptor", "job": "blue-light phototropism"},
    "P93025": {"symbol": "PHOT2", "class": "photoreceptor", "job": "blue-light phototropism"},
    "P14713": {"symbol": "PHYB", "class": "photoreceptor", "job": "red-light phytochrome"},
    "P42497": {"symbol": "PHYA", "class": "photoreceptor", "job": "far-red phytochrome"},
    "Q43125": {"symbol": "CRY1", "class": "photoreceptor", "job": "blue-light cryptochrome"},
    "Q96524": {"symbol": "CRY2", "class": "photoreceptor", "job": "blue-light cryptochrome"},
    "Q9FN03": {"symbol": "UVR8", "class": "photoreceptor", "job": "UV-B receptor"},
    # ABA / drought → stomata
    "Q940H6": {"symbol": "OST1", "class": "aba_kinase", "job": "ABA SnRK2 stomatal closure"},
    "Q8VZS9": {"symbol": "PYR1", "class": "aba_receptor", "job": "ABA receptor PYR/PYL"},
    "Q8S8E3": {"symbol": "PYL4", "class": "aba_receptor", "job": "ABA receptor PYR/PYL"},
    "P49597": {"symbol": "ABI1", "class": "aba_phosphatase", "job": "PP2C ABA negative"},
    "P25042": {"symbol": "ABI2", "class": "aba_phosphatase", "job": "PP2C ABA negative"},
    # stomatal effectors
    "Q9FLV9": {"symbol": "SLAC1", "class": "stomata_channel", "job": "guard-cell anion channel"},
    "Q39153": {"symbol": "KAT1", "class": "stomata_channel", "job": "guard-cell K+ inward"},
    "Q94KI8": {"symbol": "GORK", "class": "stomata_channel", "job": "guard-cell K+ outward"},
    # auxin / growth
    "Q9C6B8": {"symbol": "PIN1", "class": "auxin_transport", "job": "auxin efflux"},
    "Q570C0": {"symbol": "TIR1", "class": "auxin_receptor", "job": "auxin F-box receptor"},
    "P93024": {"symbol": "ARF5", "class": "auxin_tf", "job": "auxin response factor MP"},
    # light transcriptional
    "O24646": {"symbol": "HY5", "class": "light_tf", "job": "bZIP light transcription"},
    "Q495N3": {"symbol": "PIF3", "class": "light_tf", "job": "phytochrome interacting factor"},
    # already-folded metabolic (should NOT be the stomatal job)
    "O03042": {"symbol": "rbcL", "class": "calvin", "job": "carbon fixation"},
    "P83755": {"symbol": "psbA", "class": "psii", "job": "PSII D1"},
}

SEEDS = {
    "photoreceptor": ["photoreceptor"],
    "aba": ["aba_kinase", "aba_receptor"],
    "auxin": ["auxin_transport", "auxin_receptor"],
    "calvin": ["calvin", "psii"],
}


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=120) as fh:
        return fh.read()


def ensure_mitab() -> Path:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    if CACHE.exists() and CACHE.stat().st_size > 1_000_000:
        return CACHE
    print("  download IntAct Arabidopsis MITAB (species:3702)", flush=True)
    count = int(_get(PSICQUIC + "?format=count").decode().strip())
    print(f"  intact count={count}", flush=True)
    page = 2000
    chunks: list[bytes] = []
    for start in range(0, count, page):
        url = f"{PSICQUIC}?format=tab25&firstResult={start}&maxResults={page}"
        raw = _get(url)
        chunks.append(raw if raw.endswith(b"\n") else raw + b"\n")
        print(f"  fetched {start + page if start + page < count else count}/{count}", flush=True)
    CACHE.write_bytes(b"".join(chunks))
    print(f"  wrote {CACHE} ({CACHE.stat().st_size} bytes)", flush=True)
    return CACHE


_ACC = re.compile(r"uniprotkb:([A-Z0-9]{6,10})(?:-\d+)?")
_GENE = re.compile(r"uniprotkb:([A-Za-z0-9_.-]+)\(gene name\)")
_TAX = re.compile(r"taxid:(\d+)")


def _accs(field: str) -> list[str]:
    return _ACC.findall(field or "")


def _gene(field: str) -> str:
    m = _GENE.search(field or "")
    return m.group(1) if m else ""


def _tax(field: str) -> str:
    m = _TAX.search(field or "")
    return m.group(1) if m else ""


def load_graph(path: Path) -> dict[str, Any]:
    edges: Counter[tuple[str, str]] = Counter()
    genes: dict[str, str] = {}
    n_rows = 0
    n_keep = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        n_rows += 1
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        if _tax(cols[9]) != "3702" or _tax(cols[10]) != "3702":
            continue
        itype = (cols[11] or "").lower()
        if "genetic" in itype:
            continue
        a_list = _accs(cols[0]) or _accs(cols[2])
        b_list = _accs(cols[1]) or _accs(cols[3])
        if not a_list or not b_list:
            continue
        a, b = a_list[0], b_list[0]
        ga, gb = _gene(cols[4]), _gene(cols[5])
        if ga:
            genes[a] = ga
        if gb:
            genes[b] = gb
        if a == b:
            continue
        key = (a, b) if a < b else (b, a)
        edges[key] += 1
        n_keep += 1
    nodes = sorted({n for e in edges for n in e})
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    row, col, dat = [], [], []
    for (a, b), w in edges.items():
        i, j = idx[a], idx[b]
        row += [i, j]
        col += [j, i]
        dat += [float(w), float(w)]
    W = sparse.csr_matrix((dat, (row, col)), shape=(n, n), dtype=np.float64)
    return {
        "nodes": nodes,
        "idx": idx,
        "W": W,
        "genes": genes,
        "n_rows": n_rows,
        "n_physical_undirected": len(edges),
        "n_report": n_keep,
        "n_nodes": n,
        "n_directed_entries": int(W.nnz),
    }


def _info(acc: str, genes: dict[str, str]) -> dict[str, str]:
    if acc in CLASSES:
        return CLASSES[acc]
    g = (genes.get(acc) or "").upper()
    for info in CLASSES.values():
        if info["symbol"].upper() == g:
            return info
    return {}


def _mass_by_class(a: np.ndarray, nodes: list[str], genes: dict[str, str]) -> dict[str, float]:
    pos = np.maximum(a, 0.0)
    out: dict[str, float] = defaultdict(float)
    for i, acc in enumerate(nodes):
        cls = _info(acc, genes).get("class") or "unlabeled"
        out[cls] += float(pos[i])
    return dict(out)


def _top(a: np.ndarray, nodes: list[str], genes: dict[str, str], k: int = 8) -> list[dict[str, Any]]:
    pos = np.maximum(a, 0.0)
    order = np.argsort(-pos)[:k]
    rows = []
    for i in order:
        if pos[int(i)] <= 0:
            continue
        acc = nodes[int(i)]
        info = _info(acc, genes)
        rows.append(
            {
                "uniprot": acc,
                "gene": info.get("symbol") or genes.get(acc) or "",
                "class": info.get("class") or "unlabeled",
                "a": float(pos[int(i)]),
            }
        )
    return rows


def hop_seed(graph: dict[str, Any], class_names: list[str], *, seed: str) -> dict[str, Any]:
    nodes: list[str] = graph["nodes"]
    idx: dict[str, int] = graph["idx"]
    W = graph["W"]
    genes: dict[str, str] = graph["genes"]
    want = set(class_names)
    seed_i = []
    seen: set[int] = set()
    for acc, i in idx.items():
        if _info(acc, genes).get("class") in want and i not in seen:
            seed_i.append(i)
            seen.add(i)
    a = np.zeros(len(nodes), dtype=np.float64)
    if seed_i:
        a[seed_i] = 1.0 / len(seed_i)
    hops = int(round(_PHI ** 5))
    keep = {0, 1, 2, 3, hops}
    trace = []

    def snap(step: int, vec: np.ndarray) -> dict[str, Any]:
        pos = np.maximum(vec, 0.0)
        mass = _mass_by_class(vec, nodes, genes)
        return {
            "hop": step,
            "l1": float(pos.sum()),
            "n_active": int((pos > 1.0 / _PHI).sum()),
            "stomata_channel": float(mass.get("stomata_channel") or 0.0),
            "light_tf": float(mass.get("light_tf") or 0.0),
            "auxin_transport": float(mass.get("auxin_transport") or 0.0),
            "photoreceptor": float(mass.get("photoreceptor") or 0.0),
            "aba_kinase": float(mass.get("aba_kinase") or 0.0),
            "calvin": float(mass.get("calvin") or 0.0),
            "mass_by_class": dict(sorted(mass.items(), key=lambda kv: -kv[1])[:12]),
            "top": _top(vec, nodes, genes),
        }

    trace.append(snap(0, a))
    for h in range(1, hops + 1):
        a = _R_BIO * W.dot(a)
        mx = float(np.max(np.abs(a))) + 1e-12
        a = a / mx
        if h in keep:
            s = snap(h, a)
            print(
                f"  {seed} hop {h} active={s['n_active']} "
                f"stomata={s['stomata_channel']:.4f} light_tf={s['light_tf']:.4f} "
                f"top={s['top'][0]['gene'] if s['top'] else '-'}",
                flush=True,
            )
            trace.append(s)
    present = [nodes[i] for i in seed_i]
    present_sym = [(_info(a, genes).get("symbol") or genes.get(a) or a) for a in present]
    wanted_sym = [info["symbol"] for info in CLASSES.values() if info["class"] in class_names]
    return {
        "seed": seed,
        "n_seed": len(seed_i),
        "seed_accessions": present,
        "seed_on_graph": present_sym,
        "seed_missing": [s for s in wanted_sym if s not in present_sym],
        "trace": trace,
    }


def main() -> int:
    path = ensure_mitab()
    print("  build graph", flush=True)
    g = load_graph(path)
    print(
        f"  nodes={g['n_nodes']} undirected={g['n_physical_undirected']} "
        f"reports={g['n_report']} residual={_R_BIO:.6f}",
        flush=True,
    )
    runs = {}
    for name, classes in SEEDS.items():
        print(f"== seed {name}", flush=True)
        runs[name] = hop_seed(g, classes, seed=name)

    def hop_field(seed: str, hop: int, field: str) -> float:
        tr = next((t for t in runs[seed]["trace"] if t["hop"] == hop), {})
        return float(tr.get(field) or 0.0)

    compare = {
        f"hop_{h}": {
            seed: {
                "stomata_channel": hop_field(seed, h, "stomata_channel"),
                "light_tf": hop_field(seed, h, "light_tf"),
                "auxin_transport": hop_field(seed, h, "auxin_transport"),
                "n_active": int(hop_field(seed, h, "n_active")),
                "top": ((next((t for t in runs[seed]["trace"] if t["hop"] == h), {}) or {}).get("top") or [None])[0],
            }
            for seed in SEEDS
        }
        for h in (1, 2, 3)
    }
    join = {
        "product": "Arabidopsis IntAct physical PPI → residual hops (not a connectome)",
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3); r=1+|S|·P_NEW Biochemistry on measured PPI edges",
        "not": (
            "Not a plant connectome. Edges are experimental IntAct physical "
            "interactions (taxid 3702), not synapses or invented plasmodesmata."
        ),
        "authority": "IntAct PSICQUIC species:3702 MITAB 2.5; physical only; both ends Arabidopsis",
        "residual_Biochemistry": _R_BIO,
        "n_nodes": g["n_nodes"],
        "n_undirected_edges": g["n_physical_undirected"],
        "n_intact_rows": g["n_rows"],
        "n_physical_reports": g["n_report"],
        "programs": list(runs.values()),
        "compare": compare,
        "cache": str(path),
    }
    OUT_GIT.write_text(json.dumps(join, indent=2), encoding="utf-8")
    print(f"  wrote {OUT_GIT}", flush=True)
    h2 = compare["hop_2"]
    print(
        f"  hop-2 stomata  aba={h2['aba']['stomata_channel']:.4f} "
        f"light={h2['photoreceptor']['stomata_channel']:.4f} "
        f"calvin={h2['calvin']['stomata_channel']:.4f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
