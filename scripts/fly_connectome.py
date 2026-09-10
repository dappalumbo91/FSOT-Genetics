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
                "cell_class": (r.get("cell_class") or r.get("class") or "").strip(),
                "cell_type": (r.get("cell_type") or "").strip(),
                "top_nt": (r.get("top_nt") or "").strip().lower(),
                "flow": (r.get("flow") or "").strip(),
                "sub_class": (r.get("sub_class") or "").strip(),
                "side": (r.get("side") or "").strip(),
            }
    return meta


def _meta_from_classification(cls_path: Path) -> dict[str, dict[str, str]]:
    import pandas as pd

    cl = pd.read_csv(cls_path)
    meta: dict[str, dict[str, str]] = {}
    rids = cl["root_id"].to_numpy()

    def col(name: str) -> list[str]:
        if name not in cl.columns:
            return [""] * len(cl)
        return cl[name].fillna("").astype(str).tolist()

    super_c = col("super_class")
    klass = col("class")
    cell_t = col("cell_type")
    flow = col("flow")
    sub = col("sub_class")
    side = col("side")
    for i, rid_v in enumerate(rids):
        if rid_v is None or (isinstance(rid_v, float) and rid_v != rid_v):
            continue
        try:
            rid = str(int(rid_v))
        except (TypeError, ValueError):
            rid = str(rid_v).replace(".0", "")
        if not rid:
            continue
        meta[rid] = {
            "super_class": super_c[i],
            "cell_class": klass[i],
            "cell_type": cell_t[i],
            "top_nt": "",
            "flow": flow[i],
            "sub_class": sub[i],
            "side": side[i],
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


_GRAPH_CACHE: dict[str, Any] | None = None


def cuda_name() -> str | None:
    try:
        import torch

        if torch.cuda.is_available():
            return str(torch.cuda.get_device_name(0))
    except Exception:
        return None
    return None


def load_boot_graph(ann_path: Path, conn_path: Path) -> dict[str, Any]:
    """Measured synapse CSR + neuron meta. Cached per process."""
    global _GRAPH_CACHE
    key = f"{ann_path}|{conn_path}"
    if _GRAPH_CACHE is not None and _GRAPH_CACHE.get("key") == key:
        return _GRAPH_CACHE

    import pandas as pd
    from scipy.sparse import csr_matrix

    print(f"  reading {conn_path}", flush=True)
    if str(conn_path).endswith(".feather"):
        df = pd.read_feather(conn_path)
        meta = _load_meta(ann_path)
        authority = (
            "FlyWire proofread connections v783 + Schlegel annotations; "
            "GABA sign from edge nt_type"
        )
    else:
        df = pd.read_csv(conn_path)
        cls_path = conn_path.parent / "classification.csv.gz"
        if cls_path.exists():
            meta = _meta_from_classification(cls_path)
        else:
            meta = _load_meta(ann_path)
        authority = (
            "FlyWire v630 connections (Murthy/Seung public dump) when v783 "
            "Zenodo is unavailable; GABA sign from edge nt_type"
        )
    pre_c, post_c, w_c = _edge_columns(df)
    pre = df[pre_c].astype(str).str.replace(r"\.0$", "", regex=True)
    post = df[post_c].astype(str).str.replace(r"\.0$", "", regex=True)
    ids = sorted(set(meta) | set(pre) | set(post))
    idx = {rid: i for i, rid in enumerate(ids)}
    n = len(ids)
    wt = df[w_c].astype(np.float64)
    ok = pre.isin(idx) & post.isin(idx)
    pre_i = pre[ok].map(idx).to_numpy()
    post_i = post[ok].map(idx).to_numpy()
    w = wt[ok].to_numpy()
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
    graph = {
        "key": key,
        "ids": ids,
        "idx": idx,
        "meta": meta,
        "W": W,
        "n": n,
        "n_edges": int(ok.sum()),
        "authority": authority,
    }
    _GRAPH_CACHE = graph
    return graph


def seed_indices(graph: dict[str, Any], seed: str, *, how: str = "substring") -> list[int]:
    """Select measured neurons. *how*: substring | class | type_prefix."""
    seed_l = seed.lower()
    idx = graph["idx"]
    meta = graph["meta"]
    out: list[int] = []
    for rid, m in meta.items():
        if rid not in idx:
            continue
        sc = (m.get("super_class") or "").lower()
        cc = (m.get("cell_class") or "").lower()
        ct = (m.get("cell_type") or "").lower()
        fl = (m.get("flow") or "").lower()
        sub = (m.get("sub_class") or "").lower()
        hit = False
        if how == "class":
            hit = sc == seed_l or cc == seed_l
        elif how == "type_prefix":
            hit = ct.startswith(seed_l)
        else:
            hit = (
                seed_l in sc
                or seed_l in cc
                or seed_l in ct
                or seed_l in fl
                or seed_l in sub
            )
        if hit:
            out.append(idx[rid])
    return out


def _snapshot(
    a: np.ndarray,
    *,
    step: int,
    ids: list[str],
    idx: dict[str, int],
    meta: dict[str, dict[str, str]],
) -> dict[str, Any]:
    pos = np.maximum(a, 0.0)
    by_sc: Counter = Counter()
    by_cl: Counter = Counter()
    for rid in ids:
        m = meta.get(rid) or {}
        p = float(pos[idx[rid]])
        by_sc[(m.get("super_class") or "unlabeled") or "unlabeled"] += p
        by_cl[(m.get("cell_class") or "unlabeled") or "unlabeled"] += p
    top = np.argsort(-pos)[:8]

    def _sum(*keys: str) -> float:
        return float(sum(by_sc.get(k, 0.0) for k in keys))

    return {
        "hop": step,
        "l1": float(pos.sum()),
        "n_active": int((pos > 1.0 / _PHI).sum()),
        "mass_by_super_class": dict(by_sc.most_common(12)),
        "mass_by_class": dict(by_cl.most_common(12)),
        "target_mass": {
            "motor": _sum("motor", "vnc_motor", "cb_motor"),
            "descending": _sum("descending", "descending_neuron"),
            "endocrine": _sum("endocrine", "cb_endocrine", "vnc_endocrine"),
            "sensory": _sum(
                "sensory",
                "vnc_sensory",
                "cb_sensory",
                "ol_sensory",
                "sensory_ascending",
            ),
            "vnc_motor": _sum("vnc_motor"),
            "cb_motor": _sum("cb_motor"),
            "vnc_sensory": _sum("vnc_sensory"),
        },
        "class_mass": {
            "DN": float(by_cl.get("DN", 0.0)),
            "mechanosensory": float(
                by_cl.get("mechanosensory", 0.0)
                + by_cl.get("mechanosensory_tactile", 0.0)
                + by_cl.get("mechanosensory_proprioceptive", 0.0)
            ),
            "olfactory": float(by_cl.get("olfactory", 0.0)),
            "motor": float(by_cl.get("motor", 0.0)),
        },
        "top": [
            {
                "root_id": ids[int(i)],
                "super_class": (meta.get(ids[int(i)]) or {}).get("super_class", ""),
                "cell_class": (meta.get(ids[int(i)]) or {}).get("cell_class", ""),
                "cell_type": (meta.get(ids[int(i)]) or {}).get("cell_type", ""),
                "a": float(pos[int(i)]),
            }
            for i in top
            if pos[int(i)] > 0
        ],
    }


def _hop_cpu(W, a: np.ndarray) -> np.ndarray:
    a = _R_BIO * W.dot(a)
    mx = float(np.max(np.abs(a))) + 1e-12
    return a / mx


def _hop_gpu_prepare(W):
    import torch

    device = torch.device("cuda")
    Wc = W.tocsr()
    Wc.sort_indices()
    crow = torch.from_numpy(Wc.indptr.astype(np.int64)).to(device)
    col = torch.from_numpy(Wc.indices.astype(np.int64)).to(device)
    val = torch.from_numpy((Wc.data * _R_BIO).astype(np.float64)).to(device)
    try:
        torch.sparse.check_sparse_tensor_invariants(False)
    except Exception:
        pass
    Wt = torch.sparse_csr_tensor(crow, col, val, size=Wc.shape, device=device)
    return torch, device, Wt


def residual_cascade(
    graph: dict[str, Any],
    seed_i: list[int],
    *,
    seed: str,
    how: str = "substring",
    hops: int | None = None,
    gpu: bool = True,
) -> dict[str, Any]:
    """Residual hops on a loaded measured graph. 0 free parameters."""
    import time

    if hops is None:
        hops = int(round(_PHI ** 5))
    ids = graph["ids"]
    idx = graph["idx"]
    meta = graph["meta"]
    W = graph["W"]
    n = graph["n"]
    a = np.zeros(n, dtype=np.float64)
    if seed_i:
        a[np.asarray(seed_i, dtype=np.int64)] = 1.0
    keep = {1, 2, 3, hops, int(round(_PHI ** 3))}
    trace = [_snapshot(a, step=0, ids=ids, idx=idx, meta=meta)]
    device_name = None
    used_gpu = False
    t0 = time.perf_counter()
    if gpu and cuda_name() and seed_i:
        try:
            torch, device, Wt = _hop_gpu_prepare(W)
            at = torch.from_numpy(a).to(device)
            device_name = str(torch.cuda.get_device_name(0))
            used_gpu = True
            def _spmv(vec):
                try:
                    return torch.mv(Wt, vec)
                except RuntimeError:
                    return torch.sparse.mm(Wt, vec.unsqueeze(1)).squeeze(1)

            for h in range(1, hops + 1):
                at = _spmv(at)
                at = at / (at.abs().max() + 1e-12)
                if h in keep:
                    a = at.detach().cpu().numpy()
                    snap = _snapshot(a, step=h, ids=ids, idx=idx, meta=meta)
                    trace.append(snap)
                    tm = snap["target_mass"]
                    print(
                        f"  hop {h} active={snap['n_active']} "
                        f"motor={tm['motor']:.4f} "
                        f"vnc_motor={tm.get('vnc_motor', 0.0):.4f} "
                        f"desc={tm['descending']:.4f}",
                        flush=True,
                    )
        except Exception as exc:
            print(f"  GPU hop failed ({exc}); CPU residual", flush=True)
            used_gpu = False
            device_name = None
            a = np.zeros(n, dtype=np.float64)
            if seed_i:
                a[np.asarray(seed_i, dtype=np.int64)] = 1.0
            trace = [_snapshot(a, step=0, ids=ids, idx=idx, meta=meta)]
    if not used_gpu:
        for h in range(1, hops + 1):
            a = _hop_cpu(W, a)
            if h in keep:
                snap = _snapshot(a, step=h, ids=ids, idx=idx, meta=meta)
                trace.append(snap)
                tm = snap["target_mass"]
                print(
                    f"  hop {h} active={snap['n_active']} "
                    f"motor={tm['motor']:.4f} "
                    f"vnc_motor={tm.get('vnc_motor', 0.0):.4f} "
                    f"desc={tm['descending']:.4f}",
                    flush=True,
                )
    elapsed = time.perf_counter() - t0
    return {
        "seed": seed,
        "seed_how": how,
        "n_seed": len(seed_i),
        "n_neurons": n,
        "n_edges_used": graph["n_edges"],
        "hops": hops,
        "residual_Biochemistry": _R_BIO,
        "gaba_inhibitory": True,
        "gpu": used_gpu,
        "device": device_name,
        "elapsed_s": elapsed,
        "authority": graph["authority"],
        "free_parameters": 0,
        "trace": trace,
        "note": (
            "Cascade on measured synapse counts. GABA outgoing negative. "
            "Not a trained dynamics model and not a thought. "
            "v630 is brain-only: motor super_class is not VNC leg MNs; "
            "descending / DN is the brain→cord walking command."
        ),
    }


def boot_activity(
    ann_path: Path,
    conn_path: Path,
    *,
    seed: str = "sensory",
    hops: int | None = None,
    how: str = "substring",
    gpu: bool = True,
) -> dict[str, Any]:
    """Residual-scaled cascade on the measured synapse graph.

    Seed = measured neurons matching *seed* (see *how*). Edge weight =
    measured synapse count. GABA outgoing is inhibitory (measured
    transmitter). Hop count defaults to leftover φ⁵. Amplitude is
    rescaled to the observer max each hop (half-max analog). Not a
    trained RNN and not a thought.
    """
    graph = load_boot_graph(ann_path, conn_path)
    seed_i = seed_indices(graph, seed, how=how)
    return residual_cascade(
        graph, seed_i, seed=seed, how=how, hops=hops, gpu=gpu
    )


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
    ap.add_argument(
        "--how",
        default="substring",
        choices=("substring", "class", "type_prefix"),
        help="How to match --seed against measured labels.",
    )
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    ap.add_argument("--root", default=str(FLY_ROOT))
    args = ap.parse_args(argv)
    root = Path(args.root)
    path = ensure_annotations(root)
    if args.boot:
        conn = ensure_connections(root)
        run = boot_activity(path, conn, seed=args.seed, how=args.how, gpu=args.gpu)
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
