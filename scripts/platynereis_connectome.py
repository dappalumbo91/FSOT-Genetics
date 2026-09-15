#!/usr/bin/env python3
"""Platynereis dumerilii 3-day larva — measured whole-body graph.

Verasztó et al. eLife 2025. Files on D:\\FlyWire_Connectome\\Platynereis
(GitHub JekelyLab/Platynereis_3D_connectome_2024).

Cell-level chemical adjacency (synapse counts). Transmitter is mostly
unannotated — unsigned residual, no invented GABA. Sensory seed →
motor / muscle / ciliary effectors.

Same pin as fly, worm, Ciona. 0 free parameters.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402

PLAT = Path(r"D:\FlyWire_Connectome\Platynereis")
_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))


def _cell_class(name: str) -> str:
    u = name.upper()
    if (
        u.startswith("PRC")
        or u.startswith("SN")
        or "CHAEMECH" in u
        or u.startswith("ANTPU")
        or "NUCH" in u
    ):
        return "sensory"
    if u.startswith("MN") or u.startswith("MS"):
        return "motor"
    if u.startswith("MUS"):
        return "muscle"
    if u.startswith("IN"):
        return "interneuron"
    if "PROTOTROCH" in u or "PARATROCH" in u:
        return "ciliated"
    if "GLAND" in u or "COVERCELL" in u:
        return "gland"
    if "PIGMENT" in u or "CRESCENT" in u:
        return "pigment"
    return "other"


def _read_graph() -> dict[str, Any]:
    import pandas as pd
    from scipy.sparse import csr_matrix

    adj = pd.read_csv(PLAT / "full_connectome_adjacency_matrix.csv")
    names = adj["Neurons"].astype(str).tolist()
    types = [_cell_class(nm) for nm in names]
    n = len(names)
    idx = {nm: i for i, nm in enumerate(names)}
    # second 'fragment' column was renamed fragment.1
    colmap = {}
    for c in adj.columns[1:]:
        key = str(c)
        if key.endswith(".1") and key[:-2] in idx and key not in idx:
            key = key[:-2]
        colmap[str(c)] = key
    src: list[int] = []
    tgt: list[int] = []
    wts: list[float] = []
    mat = adj.iloc[:, 1:]
    for j, col in enumerate(mat.columns):
        post = colmap.get(str(col), str(col))
        if post not in idx:
            continue
        pj = idx[post]
        colv = mat.iloc[:, j].to_numpy()
        hits = np.nonzero(colv)[0]
        for i in hits:
            wt = float(colv[i])
            if wt == 0:
                continue
            src.append(int(i))
            tgt.append(pj)
            wts.append(wt)
    W = csr_matrix(
        (np.asarray(wts, dtype=np.float64), (np.asarray(tgt), np.asarray(src))),
        shape=(n, n),
    )
    return {
        "names": names,
        "types": types,
        "W": W,
        "n": n,
        "n_edges": int(len(src)),
        "type_counts": dict(Counter(types)),
    }


def boot_activity(seed: str = "sensory", hops: int | None = None) -> dict[str, Any]:
    if hops is None:
        hops = int(round(_PHI ** 5))
    g = _read_graph()
    seed_l = seed.lower()
    seed_i = [i for i, t in enumerate(g["types"]) if seed_l in t.lower()]
    n = g["n"]
    a = np.zeros(n, dtype=np.float64)
    if seed_i:
        a[np.asarray(seed_i, dtype=np.int64)] = 1.0
    keep = {1, 2, 3, hops, int(round(_PHI ** 3))}
    names, types, W = g["names"], g["types"], g["W"]

    def snap(step: int) -> dict[str, Any]:
        pos = np.maximum(a, 0.0)
        by = Counter()
        for i, t in enumerate(types):
            by[t] += float(pos[i])
        top = np.argsort(-pos)[:8]
        return {
            "hop": step,
            "l1": float(pos.sum()),
            "n_active": int((pos > 1.0 / _PHI).sum()),
            "mass_by_type": dict(by.most_common()),
            "target_mass": {
                "sensory": float(by.get("sensory", 0.0)),
                "motor": float(by.get("motor", 0.0)),
                "muscle": float(by.get("muscle", 0.0)),
                "ciliated": float(by.get("ciliated", 0.0)),
                "interneuron": float(by.get("interneuron", 0.0)),
            },
            "top": [
                {
                    "name": names[int(i)],
                    "type": types[int(i)],
                    "a": float(pos[int(i)]),
                }
                for i in top
                if pos[int(i)] > 0
            ],
        }

    trace = [snap(0)]
    print(
        f"  Platynereis n={n} edges={g['n_edges']} seed={seed} "
        f"n_seed={len(seed_i)} UNSIGNED (NT mostly unannotated)",
        flush=True,
    )
    for h in range(1, hops + 1):
        a = _R_BIO * W.dot(a)
        mx = float(np.max(np.abs(a))) + 1e-12
        a = a / mx
        if h in keep:
            s = snap(h)
            trace.append(s)
            tm = s["target_mass"]
            print(
                f"  hop {h} active={s['n_active']} "
                f"motor={tm['motor']:.4f} muscle={tm['muscle']:.4f} "
                f"cilia={tm['ciliated']:.4f}",
                flush=True,
            )
    return {
        "organism": "Platynereis dumerilii 3-day larva",
        "seed": seed,
        "n_seed": len(seed_i),
        "n_cells": n,
        "n_edges": g["n_edges"],
        "type_counts": g["type_counts"],
        "hops": hops,
        "residual_Biochemistry": _R_BIO,
        "gaba_inhibitory": False,
        "authority": (
            "Verasztó et al. eLife 2025; "
            "JekelyLab/Platynereis_3D_connectome_2024 "
            "full_connectome_adjacency_matrix.csv; synapse counts; "
            "unsigned residual — transmitter mostly unannotated"
        ),
        "free_parameters": 0,
        "trace": trace,
        "note": (
            "Segmented annelid whole body: neurons, muscle, ciliary bands. "
            "Sensory seed → motor / muscle / cilia. Not a trained dynamics model."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--boot", action="store_true")
    ap.add_argument("--seed", default="sensory")
    args = ap.parse_args(argv)
    if not args.boot:
        g = _read_graph()
        inv = {
            "n_cells": g["n"],
            "n_edges": g["n_edges"],
            "type_counts": g["type_counts"],
            "free_parameters": 0,
        }
        print(json.dumps(inv, indent=2))
        out = ROOT / "data" / "platynereis_connectome_inventory.json"
        out.write_text(json.dumps(inv, indent=2), encoding="utf-8")
        print(f"  wrote {out}", flush=True)
        return 0
    run = boot_activity(seed=args.seed)
    print(json.dumps({k: run[k] for k in run if k != "trace"}, indent=2))
    out = ROOT / "data" / "platynereis_connectome_boot.json"
    out.write_text(json.dumps(run, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
