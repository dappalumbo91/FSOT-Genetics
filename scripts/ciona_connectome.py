#!/usr/bin/env python3
"""Ciona intestinalis larval CNS — measured chordate graph.

Ryan, Lu & Meinertzhagen 2016 (eLife). Netzschleuder dump on
D:\\FlyWire_Connectome\\Ciona. 205 nodes (177 CNS + muscle / periphery),
2,903 directed edges, weight = contact depth (µm).

Synaptic sign is NOT annotated in this dump. Residual hops are unsigned
(no invented GABA). Photoreceptor seed → motor / tail muscle.

Same pin, 0 free parameters. Not a thought.
"""
from __future__ import annotations

import json
import re
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

CIONA = Path(r"D:\FlyWire_Connectome\Ciona\graph")
_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))


def _cell_class(name: str) -> str:
    s = name.strip()
    sl = s.lower()
    if sl.startswith("pr") or sl.startswith("pns") or sl.startswith("aten"):
        return "sensory"
    if sl.startswith("ant") or re.fullmatch(r"pn[a-z]", sl):
        return "sensory"
    if sl.startswith("lens"):
        return "sensory"
    if sl.startswith("mn"):
        return "motor"
    if sl.startswith("midtail"):
        return "muscle"
    if sl.startswith("coronet"):
        return "coronet"
    if "in" in sl or sl.startswith("amg") or sl.startswith("neck"):
        return "interneuron"
    if sl.startswith("ddn"):
        return "interneuron"
    if sl.isdigit():
        return "unlabeled"
    return "other"


def _read_graph() -> dict[str, Any]:
    import pandas as pd
    from scipy.sparse import csr_matrix

    nodes = pd.read_csv(CIONA / "nodes.csv")
    edges = pd.read_csv(CIONA / "edges.csv")
    nodes.columns = [c.strip().lstrip("#").strip() for c in nodes.columns]
    edges.columns = [c.strip().lstrip("#").strip() for c in edges.columns]
    names = nodes["name"].astype(str).tolist()
    types = [_cell_class(nm) for nm in names]
    n = len(names)
    src = edges["source"].to_numpy(dtype=np.int64)
    tgt = edges["target"].to_numpy(dtype=np.int64)
    w = edges["depth"].to_numpy(dtype=np.float64)
    W = csr_matrix((w, (tgt, src)), shape=(n, n))
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
                "interneuron": float(by.get("interneuron", 0.0)),
                "coronet": float(by.get("coronet", 0.0)),
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
        f"  Ciona n={n} edges={g['n_edges']} seed={seed} n_seed={len(seed_i)} "
        f"UNSIGNED (NT not annotated)",
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
                f"motor={tm['motor']:.4f} muscle={tm['muscle']:.4f}",
                flush=True,
            )
    return {
        "organism": "Ciona intestinalis tadpole larva",
        "seed": seed,
        "n_seed": len(seed_i),
        "n_cells": n,
        "n_edges": g["n_edges"],
        "type_counts": g["type_counts"],
        "hops": hops,
        "residual_Biochemistry": _R_BIO,
        "gaba_inhibitory": False,
        "authority": (
            "Ryan, Lu & Meinertzhagen 2016 via Netzschleuder; "
            "weight = contact depth µm; synaptic sign not annotated — "
            "unsigned residual, no invented GABA"
        ),
        "free_parameters": 0,
        "trace": trace,
        "note": (
            "Chordate sibling of vertebrates. Photoreceptor seed → motor / "
            "tail muscle. Not a trained dynamics model."
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
        out = ROOT / "data" / "ciona_connectome_inventory.json"
        out.write_text(json.dumps(inv, indent=2), encoding="utf-8")
        print(f"  wrote {out}", flush=True)
        return 0
    run = boot_activity(seed=args.seed)
    print(json.dumps({k: run[k] for k in run if k != "trace"}, indent=2))
    out = ROOT / "data" / "ciona_connectome_boot.json"
    out.write_text(json.dumps(run, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
