#!/usr/bin/env python3
"""Read the FlyWire adult Drosophila connectome (measured authority).

Annotations live on the game drive, not in git:

  D:\\FlyWire_Connectome

  python scripts/fly_connectome.py --inventory

Same pin as protein product and Biohub: measured soma/anchor coordinates
and measured types/transmitters. Residual does not invent synapses.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402

FLY_ROOT = Path(r"D:\FlyWire_Connectome")
ANN_URL = (
    "https://raw.githubusercontent.com/flyconnectome/flywire_annotations/"
    "main/supplemental_files/Supplemental_file1_neuron_annotations.tsv"
)
ANN_NAME = "Supplemental_file1_neuron_annotations.tsv"
CONN_NAME = "proofread_connections_783.feather"
CONN_URL = (
    "https://zenodo.org/records/10676866/files/proofread_connections_783.feather"
    "?download=1"
)
# FlyWire FAFB EM voxel (Schlegel et al. 2024)
SCALE_XY_NM = 4.0
SCALE_Z_NM = 40.0

_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))


def ensure_annotations(dest_dir: Path = FLY_ROOT) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / ANN_NAME
    if path.exists() and path.stat().st_size > 1_000_000:
        return path
    print(f"  downloading {ANN_URL}", flush=True)
    urllib.request.urlretrieve(ANN_URL, path)
    print(f"  wrote {path} ({path.stat().st_size} bytes)", flush=True)
    return path


def _float(v: str) -> float | None:
    v = (v or "").strip()
    if not v or v.lower() in {"na", "nan", "none"}:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def read_annotations(path: Path) -> dict[str, Any]:
    """Measured neuron atlas: types, transmitters, soma xyz (voxels)."""
    import csv

    with path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    n = len(rows)
    super_c = Counter((r.get("super_class") or "").strip() or "unlabeled" for r in rows)
    flow = Counter((r.get("flow") or "").strip() or "unlabeled" for r in rows)
    nt = Counter((r.get("top_nt") or "").strip() or "unlabeled" for r in rows)
    side = Counter((r.get("side") or "").strip() or "unlabeled" for r in rows)
    dimorph = Counter((r.get("dimorphism") or "").strip() or "unlabeled" for r in rows)
    fru = Counter((r.get("fru_dsx") or "").strip() or "unlabeled" for r in rows)
    n_type = sum(1 for r in rows if (r.get("cell_type") or "").strip())
    n_soma = 0
    soma = []
    for r in rows:
        x, y, z = _float(r.get("soma_x")), _float(r.get("soma_y")), _float(r.get("soma_z"))
        if x is None or y is None or z is None:
            continue
        n_soma += 1
        soma.append((x * SCALE_XY_NM, y * SCALE_XY_NM, z * SCALE_Z_NM))
    span = None
    if soma:
        a = np.asarray(soma, dtype=np.float64)
        span = (a.max(axis=0) - a.min(axis=0)).tolist()
    return {
        "path": str(path),
        "n_neurons": n,
        "n_with_cell_type": n_type,
        "n_with_soma_xyz": n_soma,
        "super_class": dict(super_c.most_common(24)),
        "flow": dict(flow.most_common()),
        "top_nt": dict(nt.most_common()),
        "side": dict(side.most_common()),
        "dimorphism": dict(dimorph.most_common(12)),
        "fru_dsx": dict(fru.most_common(12)),
        "soma_span_nm": span,
        "voxel_nm": [SCALE_XY_NM, SCALE_XY_NM, SCALE_Z_NM],
        "authority": "FlyWire FAFB v783 + Schlegel et al. 2024 annotations",
        "free_parameters": 0,
        "phi": _PHI,
        "residual_Biochemistry": _R_BIO,
    }


def ensure_connections(dest_dir: Path = FLY_ROOT) -> Path:
    """Prefer local v630 dump if Zenodo v783 is not on disk (504s are common)."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    v783 = dest_dir / CONN_NAME
    if v783.exists() and v783.stat().st_size > 10_000_000:
        return v783
    v630 = dest_dir / "v630" / "connections.csv.gz"
    if v630.exists() and v630.stat().st_size > 1_000_000:
        return v630
    print(f"  downloading {CONN_URL}", flush=True)
    req = urllib.request.Request(CONN_URL, headers={"User-Agent": "FSOT-Genetics/fly"})
    with urllib.request.urlopen(req, timeout=600) as resp, v783.open("wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            print(f"  ... {v783.stat().st_size / 1e6:.1f} MB", flush=True)
    print(f"  wrote {v783} ({v783.stat().st_size} bytes)", flush=True)
    return v783


def _load_meta(ann_path: Path) -> dict[str, dict[str, str]]:
    import csv

    meta: dict[str, dict[str, str]] = {}
    with ann_path.open("r", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rid = (r.get("root_id") or "").strip()
            if not rid:
                continue
            meta[rid] = {
                "super_class": (r.get("super_class") or "").strip(),
                "cell_class": (r.get("cell_class") or "").strip(),
                "cell_type": (r.get("cell_type") or "").strip(),
                "top_nt": (r.get("top_nt") or "").strip().lower(),
                "flow": (r.get("flow") or "").strip(),
            }
    return meta


def _edge_columns(df) -> tuple[str, str, str]:
    cols = {c.lower(): c for c in df.columns}
    pre = cols.get("pre_root_id") or cols.get("pre") or cols.get("pre_pt_root_id")
    post = cols.get("post_root_id") or cols.get("post") or cols.get("post_pt_root_id")
    w = (
        cols.get("syn_count")
        or cols.get("n_synapses")
        or cols.get("synapses")
        or cols.get("weight")
        or cols.get("count")
    )
    if not pre or not post:
        raise KeyError(f"need pre/post columns, got {list(df.columns)}")
    if not w:
        raise KeyError(f"need synapse-count column, got {list(df.columns)}")
    return pre, post, w


def boot_activity(
    ann_path: Path,
    conn_path: Path,
    *,
    seed: str = "sensory",
    hops: int | None = None,
) -> dict[str, Any]:
    """Residual-scaled cascade on the measured synapse graph.

    Seed = all neurons whose super_class (or cell_class/type) matches *seed*.
    Edge weight = measured synapse count. GABA outgoing is inhibitory
    (measured transmitter). Hop count defaults to leftover φ⁵. Amplitude
    is rescaled to the observer max each hop (half-max analog). Not a
    trained RNN and not a thought.
    """
    import pandas as pd
    from scipy.sparse import csr_matrix

    if hops is None:
        hops = int(round(_PHI ** 5))
    print(f"  reading {conn_path}", flush=True)
    if str(conn_path).endswith(".feather"):
        df = pd.read_feather(conn_path)
        meta = _load_meta(ann_path)
    else:
        df = pd.read_csv(conn_path)
        cls_path = conn_path.parent / "classification.csv.gz"
        meta = {}
        if cls_path.exists():
            cl = pd.read_csv(cls_path)
            for _, r in cl.iterrows():
                rid = str(int(r["root_id"])) if pd.notna(r["root_id"]) else ""
                if not rid:
                    continue
                meta[rid] = {
                    "super_class": str(r.get("super_class") or ""),
                    "cell_class": str(r.get("class") or ""),
                    "cell_type": str(r.get("cell_type") or ""),
                    "top_nt": "",
                    "flow": str(r.get("flow") or ""),
                }
        else:
            meta = _load_meta(ann_path)
    pre_c, post_c, w_c = _edge_columns(df)
    pre = df[pre_c].astype(str)
    post = df[post_c].astype(str)
    # root ids may be ints
    pre = pre.str.replace(r"\.0$", "", regex=True)
    post = post.str.replace(r"\.0$", "", regex=True)
    ids = sorted(set(meta) | set(pre) | set(post))
    idx = {rid: i for i, rid in enumerate(ids)}
    n = len(ids)
    seed_l = seed.lower()
    seed_i = [
        idx[rid]
        for rid, m in meta.items()
        if rid in idx
        and (
            seed_l in (m.get("super_class") or "").lower()
            or seed_l in (m.get("cell_class") or "").lower()
            or seed_l in (m.get("cell_type") or "").lower()
            or seed_l in (m.get("flow") or "").lower()
        )
    ]
    wt = df[w_c].astype(np.float64)
    ok = pre.isin(idx) & post.isin(idx)
    pre_i = pre[ok].map(idx).to_numpy()
    post_i = post[ok].map(idx).to_numpy()
    w = wt[ok].to_numpy()
    # Inhibitory: edge nt_type if present, else presynaptic top_nt.
    sign = np.ones(len(pre_i), dtype=np.float64)
    if "nt_type" in df.columns:
        nt = df.loc[ok, "nt_type"].astype(str).str.upper().to_numpy()
        sign = np.where(nt == "GABA", -1.0, 1.0)
    else:
        for k, pi in enumerate(pre_i):
            rid = ids[int(pi)]
            if (meta.get(rid) or {}).get("top_nt") == "gaba":
                sign[k] = -1.0
    W = csr_matrix((w * sign, (post_i, pre_i)), shape=(n, n))
    a = np.zeros(n, dtype=np.float64)
    a[seed_i] = 1.0
    trace = []
    targets = ("motor", "descending", "endocrine", "sensory")

    def snapshot(step: int) -> dict[str, Any]:
        pos = np.maximum(a, 0.0)
        by = Counter()
        for rid in ids:
            m = meta.get(rid) or {}
            sc = (m.get("super_class") or "unlabeled") or "unlabeled"
            by[sc] += float(pos[idx[rid]])
        top = np.argsort(-pos)[:8]
        return {
            "hop": step,
            "l1": float(pos.sum()),
            "n_active": int((pos > 1.0 / _PHI).sum()),
            "mass_by_super_class": dict(by.most_common(12)),
            "target_mass": {
                t: float(by.get(t, 0.0)) for t in targets
            },
            "top": [
                {
                    "root_id": ids[int(i)],
                    "super_class": (meta.get(ids[int(i)]) or {}).get("super_class", ""),
                    "cell_type": (meta.get(ids[int(i)]) or {}).get("cell_type", ""),
                    "a": float(pos[int(i)]),
                }
                for i in top
                if pos[int(i)] > 0
            ],
        }

    trace.append(snapshot(0))
    for h in range(1, hops + 1):
        a = _R_BIO * W.dot(a)
        mx = float(np.max(np.abs(a))) + 1e-12
        a = a / mx
        if h in {1, 2, 3, hops} or h == int(round(_PHI ** 3)):
            snap = snapshot(h)
            trace.append(snap)
            print(
                f"  hop {h} active={snap['n_active']} "
                f"motor={snap['target_mass']['motor']:.4f} "
                f"desc={snap['target_mass']['descending']:.4f}",
                flush=True,
            )
    return {
        "seed": seed,
        "n_seed": len(seed_i),
        "n_neurons": n,
        "n_edges_used": int(ok.sum()),
        "hops": hops,
        "residual_Biochemistry": _R_BIO,
        "gaba_inhibitory": True,
        "authority": (
            "FlyWire v630 connections (Murthy/Seung public dump) when v783 "
            "Zenodo is unavailable; GABA sign from edge nt_type"
        ),
        "free_parameters": 0,
        "trace": trace,
        "note": (
            "Cascade on measured synapse counts. GABA outgoing negative. "
            "Not a trained dynamics model and not a thought."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", action="store_true")
    ap.add_argument(
        "--boot",
        action="store_true",
        help="Residual-scaled sensory→motor cascade on proofread connections.",
    )
    ap.add_argument("--seed", default="sensory", help="super_class / cell_class substring")
    ap.add_argument("--root", default=str(FLY_ROOT))
    args = ap.parse_args(argv)
    root = Path(args.root)
    path = ensure_annotations(root)
    if args.boot:
        conn = ensure_connections(root)
        run = boot_activity(path, conn, seed=args.seed)
        print(json.dumps(run, indent=2))
        out = ROOT / "data" / "fly_connectome_boot.json"
        out.write_text(json.dumps(run, indent=2), encoding="utf-8")
        print(f"  wrote {out}", flush=True)
        return 0
    inv = read_annotations(path)
    print(json.dumps(inv, indent=2))
    out = ROOT / "data" / "fly_connectome_inventory.json"
    out.write_text(json.dumps(inv, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
