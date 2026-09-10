#!/usr/bin/env python3
"""C. elegans whole-animal connectome — measured authority.

Cook et al. 2019 hermaphrodite chemical graph lives on the game drive:

  D:\\FlyWire_Connectome\\C_elegans

This is the only adult animal with named neurons AND named muscles on
the same synaptic graph. Seed sensory → residual hops → motor + body-wall
muscle mass. GABA names are the WormAtlas / Pereira set (measured), not
a trained classifier.

Same pin as the fly boot. 0 free parameters. Not a thought.
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

WORM_ROOT = Path(r"D:\FlyWire_Connectome\C_elegans")
CHEM = WORM_ROOT / "hermaphrodite_chemical_corrected"
SI5 = WORM_ROOT / "SI5.xlsx"

_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))

# Measured GABAergic neurons (WormAtlas / Pereira et al. 2015).
_GABA = frozenset(
    {
        "AVL",
        "DVB",
        "RIS",
        "RMED",
        "RMEV",
        "RMEL",
        "RMER",
        *(f"DD{i}" for i in range(1, 7)),
        *(f"VD{i}" for i in range(1, 14)),
    }
)


def _read_si5(path: Path = SI5) -> dict[str, Any] | None:
    """Cook 2019 SI5 hermaphrodite chemical, including NMJs onto muscles.

    Netzschleuder CSV lists muscle nodes but drops every edge into them.
    SI5 is the measured adjacency (pre = rows, post = columns).
    """
    if not path.exists() or path.stat().st_size < 100_000:
        return None
    import pandas as pd
    from scipy.sparse import csr_matrix

    df = pd.read_excel(path, sheet_name="hermaphrodite chemical", header=None)
    groups: list[str] = []
    cur = "UNLABELED"
    for c in range(df.shape[1]):
        v = df.iloc[0, c]
        if pd.notna(v):
            cur = str(v).strip().upper()
        groups.append(cur)
    post_names = [None if pd.isna(x) else str(x).strip() for x in df.iloc[2].tolist()]
    pre_rows: list[tuple[int, str]] = []
    for r in range(3, df.shape[0]):
        v = df.iloc[r, 2]
        if pd.notna(v):
            pre_rows.append((r, str(v).strip()))
    type_of: dict[str, str] = {}
    for c, nm in enumerate(post_names):
        if nm:
            type_of[nm] = groups[c]
    for r, nm in pre_rows:
        type_of.setdefault(nm, "UNLABELED")
    names = sorted(type_of)
    idx = {nm: i for i, nm in enumerate(names)}
    types = [type_of[nm] for nm in names]
    src: list[int] = []
    tgt: list[int] = []
    wts: list[float] = []
    for r, pre in pre_rows:
        pi = idx[pre]
        for c, post in enumerate(post_names):
            if not post or c < 3:
                continue
            v = df.iloc[r, c]
            if pd.isna(v):
                continue
            try:
                wt = float(v)
            except (TypeError, ValueError):
                continue
            if wt == 0.0:
                continue
            src.append(pi)
            tgt.append(idx[post])
            wts.append(wt)
    sign = np.ones(len(src), dtype=np.float64)
    for k, si in enumerate(src):
        if names[int(si)] in _GABA:
            sign[k] = -1.0
    n = len(names)
    W = csr_matrix(
        (np.asarray(wts) * sign, (np.asarray(tgt), np.asarray(src))),
        shape=(n, n),
    )
    n_nmj = sum(
        1
        for t, s in zip(tgt, src)
        if types[int(t)] == "BODYWALL MUSCLES"
    )
    return {
        "names": names,
        "types": types,
        "W": W,
        "n": n,
        "n_edges": int(len(src)),
        "n_nmj_to_bwm": int(n_nmj),
        "type_counts": dict(Counter(types)),
        "n_gaba": int(sum(1 for nm in names if nm in _GABA)),
        "source": str(path),
    }


def _read_graph(chem_dir: Path = CHEM) -> dict[str, Any]:
    si5 = _read_si5()
    if si5 is not None:
        return si5
    import pandas as pd
    from scipy.sparse import csr_matrix

    nodes = pd.read_csv(chem_dir / "nodes.csv")
    edges = pd.read_csv(chem_dir / "edges.csv")
    nodes.columns = [c.strip().lstrip("#").strip() for c in nodes.columns]
    edges.columns = [c.strip().lstrip("#").strip() for c in edges.columns]
    names = nodes["name"].astype(str).tolist()
    types = (
        nodes["node_type"].astype(str).str.strip().str.upper().tolist()
        if "node_type" in nodes.columns
        else ["UNLABELED"] * len(nodes)
    )
    n = len(names)
    src = edges["source"].to_numpy(dtype=np.int64)
    tgt = edges["target"].to_numpy(dtype=np.int64)
    w = edges["connectivity"].to_numpy(dtype=np.float64)
    sign = np.ones(len(src), dtype=np.float64)
    for k, si in enumerate(src):
        if names[int(si)] in _GABA:
            sign[k] = -1.0
    W = csr_matrix((w * sign, (tgt, src)), shape=(n, n))
    return {
        "names": names,
        "types": types,
        "W": W,
        "n": n,
        "n_edges": int(len(src)),
        "n_nmj_to_bwm": 0,
        "type_counts": dict(Counter(types)),
        "n_gaba": int(sum(1 for nm in names if nm in _GABA)),
        "source": str(chem_dir),
    }


def _seed_i(graph: dict[str, Any], seed: str) -> list[int]:
    seed_l = seed.lower()
    out = []
    for i, t in enumerate(graph["types"]):
        if seed_l in t.lower() or seed_l == graph["names"][i].lower():
            out.append(i)
    return out


def boot_activity(seed: str = "sensory", hops: int | None = None) -> dict[str, Any]:
    if hops is None:
        hops = int(round(_PHI ** 5))
    g = _read_graph()
    seed_i = _seed_i(g, seed)
    n = g["n"]
    a = np.zeros(n, dtype=np.float64)
    if seed_i:
        a[np.asarray(seed_i, dtype=np.int64)] = 1.0
    keep = {1, 2, 3, hops, int(round(_PHI ** 3))}
    names = g["names"]
    types = g["types"]
    W = g["W"]

    def snap(step: int) -> dict[str, Any]:
        pos = np.maximum(a, 0.0)
        by = Counter()
        for i, t in enumerate(types):
            by[t] += float(pos[i])
        top = np.argsort(-pos)[:8]
        muscle = float(by.get("BODYWALL MUSCLES", 0.0))
        motor = float(by.get("MOTOR NEURONS", 0.0))
        sensory = float(by.get("SENSORY NEURONS", 0.0))
        return {
            "hop": step,
            "l1": float(pos.sum()),
            "n_active": int((pos > 1.0 / _PHI).sum()),
            "mass_by_type": dict(by.most_common()),
            "target_mass": {
                "sensory": sensory,
                "motor": motor,
                "bodywall_muscle": muscle,
                "interneuron": float(by.get("INTERNEURONS", 0.0)),
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
        f"  C. elegans hermaphrodite chemical n={n} edges={g['n_edges']} "
        f"NMJ={g.get('n_nmj_to_bwm', 0)} seed={seed} n_seed={len(seed_i)} "
        f"GABA={g['n_gaba']}",
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
                f"motor={tm['motor']:.4f} muscle={tm['bodywall_muscle']:.4f}",
                flush=True,
            )
    return {
        "organism": "Caenorhabditis elegans hermaphrodite",
        "seed": seed,
        "n_seed": len(seed_i),
        "n_cells": n,
        "n_edges": g["n_edges"],
        "n_nmj_to_bwm": g.get("n_nmj_to_bwm", 0),
        "type_counts": g["type_counts"],
        "n_gaba": g["n_gaba"],
        "graph_source": g.get("source"),
        "hops": hops,
        "residual_Biochemistry": _R_BIO,
        "gaba_inhibitory": True,
        "authority": (
            "Cook et al. 2019 SI5 hermaphrodite chemical adjacency "
            "(OpenWorm ConnectomeToolbox copy); GABA names "
            "WormAtlas/Pereira 2015; NMJs onto body-wall muscles included. "
            "Netzschleuder CSV listed muscle nodes but dropped every NMJ."
        ),
        "free_parameters": 0,
        "trace": trace,
        "note": (
            "Whole-animal graph: neurons + muscles. Not a trained dynamics "
            "model and not a thought. This is the organism-complete analog "
            "of the fly brain boot."
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
            "n_gaba": g["n_gaba"],
            "path": str(CHEM),
            "free_parameters": 0,
        }
        print(json.dumps(inv, indent=2))
        out = ROOT / "data" / "worm_connectome_inventory.json"
        out.write_text(json.dumps(inv, indent=2), encoding="utf-8")
        print(f"  wrote {out}", flush=True)
        return 0
    run = boot_activity(seed=args.seed)
    print(json.dumps({k: run[k] for k in run if k != "trace"}, indent=2))
    out = ROOT / "data" / "worm_connectome_boot.json"
    out.write_text(json.dumps(run, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
