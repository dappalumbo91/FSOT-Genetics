#!/usr/bin/env python3
"""Drosophila first-instar larva brain — measured CNS graph.

Winding et al. Science 2023. Files on D:\\FlyWire_Connectome\\larva
(GitHub brain-networks/larval-drosophila-connectome Supplementary-Data-S1).

~2,952 neurons, all-all synapse counts. Cell types from the paper
annotations (sensory, DN-VNC, DN-SEZ, KC, …). Transmitter is not in
this dump — unsigned residual, no invented GABA.

Same pin as adult fly, worm, Ciona, Platynereis. 0 free parameters.
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

LARVA = Path(r"D:\FlyWire_Connectome\larva\Supplementary-Data-S1")
_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))

PROGRAMS = [
    {"name": "sensory", "seed": "sensory", "how": "celltype"},
    {"name": "olfactory", "seed": "olfactory", "how": "extra"},
    {"name": "mechano", "seed": "mechano", "how": "extra"},
    {"name": "gustatory", "seed": "gustatory", "how": "extra"},
]


def _read_graph() -> dict[str, Any]:
    import pandas as pd
    from scipy.sparse import csr_matrix

    ann = pd.read_csv(LARVA / "annotations.csv")
    meta: dict[str, dict[str, str]] = {}
    for _, row in ann.iterrows():
        ct = str(row.get("celltype") or "")
        extra = str(row.get("additional_annotations") or "")
        for col in ("left_id", "right_id"):
            rid = str(row[col]).strip()
            if not rid or rid.lower() == "no pair":
                continue
            meta[rid] = {"celltype": ct, "extra": extra}

    print(f"  reading all-all matrix", flush=True)
    adj = pd.read_csv(LARVA / "all-all_connectivity_matrix.csv")
    names = [str(x) for x in adj.iloc[:, 0].tolist()]
    n = len(names)
    idx = {nm: i for i, nm in enumerate(names)}
    types = []
    extras = []
    for nm in names:
        m = meta.get(nm) or {}
        types.append(m.get("celltype") or "unlabeled")
        extras.append(m.get("extra") or "")
    mat = adj.iloc[:, 1:]
    col_ids = [str(c) for c in mat.columns]
    post_of_col = np.array(
        [idx[c] if c in idx else -1 for c in col_ids], dtype=np.int64
    )
    arr = mat.to_numpy(dtype=np.float64)
    pre_i, col_j = np.nonzero(arr)
    post_i = post_of_col[col_j]
    ok = post_i >= 0
    pre_i = pre_i[ok]
    post_i = post_i[ok]
    wts = arr[pre_i, col_j[ok]]
    W = csr_matrix((wts, (post_i, pre_i)), shape=(n, n))
    n_typed = sum(1 for t in types if t != "unlabeled")
    print(
        f"  n={n} edges={int(ok.sum())} typed={n_typed}",
        flush=True,
    )
    return {
        "names": names,
        "types": types,
        "extras": extras,
        "W": W,
        "n": n,
        "n_edges": int(ok.sum()),
        "type_counts": dict(Counter(types).most_common()),
        "n_typed": n_typed,
    }


def _seed_i(g: dict[str, Any], seed: str, how: str) -> list[int]:
    seed_l = seed.lower()
    out: list[int] = []
    for i, (t, e) in enumerate(zip(g["types"], g["extras"])):
        hit = False
        if how == "celltype":
            hit = t.lower() == seed_l
        elif how == "extra":
            hit = seed_l in e.lower()
        else:
            hit = seed_l in t.lower() or seed_l in e.lower()
        if hit:
            out.append(i)
    return out


def boot_activity(
    g: dict[str, Any],
    seed: str,
    how: str,
    hops: int | None = None,
) -> dict[str, Any]:
    if hops is None:
        hops = int(round(_PHI ** 5))
    seed_i = _seed_i(g, seed, how)
    n = g["n"]
    a = np.zeros(n, dtype=np.float64)
    if seed_i:
        a[np.asarray(seed_i, dtype=np.int64)] = 1.0
    keep = {1, 2, 3, hops, int(round(_PHI ** 3))}
    names, types, extras, W = g["names"], g["types"], g["extras"], g["W"]

    def snap(step: int) -> dict[str, Any]:
        pos = np.maximum(a, 0.0)
        by = Counter()
        for i, t in enumerate(types):
            by[t] += float(pos[i])
        top = np.argsort(-pos)[:8]

        def _sum(*keys: str) -> float:
            return float(sum(by.get(k, 0.0) for k in keys))

        return {
            "hop": step,
            "l1": float(pos.sum()),
            "n_active": int((pos > 1.0 / _PHI).sum()),
            "mass_by_type": dict(by.most_common(12)),
            "target_mass": {
                "sensory": _sum("sensory"),
                "DN_VNC": _sum("DN-VNC"),
                "DN_SEZ": _sum("DN-SEZ"),
                "pre_DN_VNC": _sum("pre-DN-VNC"),
                "ascending": _sum("ascending"),
                "KC": _sum("KC"),
            },
            "top": [
                {
                    "id": names[int(i)],
                    "celltype": types[int(i)],
                    "extra": extras[int(i)],
                    "a": float(pos[int(i)]),
                }
                for i in top
                if pos[int(i)] > 0
            ],
        }

    trace = [snap(0)]
    print(
        f"  larva n={n} edges={g['n_edges']} seed={seed} how={how} "
        f"n_seed={len(seed_i)} UNSIGNED",
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
            top0 = (s["top"] or [{}])[0]
            print(
                f"  hop {h} active={s['n_active']} "
                f"DN-VNC={tm['DN_VNC']:.4f} DN-SEZ={tm['DN_SEZ']:.4f} "
                f"top={top0.get('celltype')}",
                flush=True,
            )
    return {
        "seed": seed,
        "how": how,
        "n_seed": len(seed_i),
        "hops": hops,
        "trace": trace,
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--boot", action="store_true")
    ap.add_argument("--seed", default="sensory")
    args = ap.parse_args(argv)
    if not (LARVA / "all-all_connectivity_matrix.csv").exists():
        raise SystemExit(f"missing {LARVA}")
    g = _read_graph()
    if not args.boot:
        inv = {
            "n_cells": g["n"],
            "n_edges": g["n_edges"],
            "n_typed": g["n_typed"],
            "type_counts": g["type_counts"],
            "free_parameters": 0,
        }
        print(json.dumps(inv, indent=2))
        out = ROOT / "data" / "larva_connectome_inventory.json"
        out.write_text(json.dumps(inv, indent=2), encoding="utf-8")
        print(f"  wrote {out}", flush=True)
        return 0
    programs = []
    for spec in PROGRAMS:
        run = boot_activity(g, spec["seed"], spec["how"])
        flow = []
        for snap in run["trace"]:
            tm = snap.get("target_mass") or {}
            top = (snap.get("top") or [{}])[0]
            flow.append(
                {
                    "hop": snap.get("hop"),
                    "l1": snap.get("l1"),
                    "n_active": snap.get("n_active"),
                    "DN_VNC": tm.get("DN_VNC"),
                    "DN_SEZ": tm.get("DN_SEZ"),
                    "pre_DN_VNC": tm.get("pre_DN_VNC"),
                    "sensory": tm.get("sensory"),
                    "ascending": tm.get("ascending"),
                    "top": {
                        "celltype": top.get("celltype"),
                        "extra": top.get("extra"),
                        "a": top.get("a"),
                    },
                }
            )
        programs.append(
            {
                "name": spec["name"],
                "seed": spec["seed"],
                "how": spec["how"],
                "n_seed": run["n_seed"],
                "flow": flow,
            }
        )
    compare = {}
    for hop in (1, 2, 3):
        row = {}
        for p in programs:
            hit = next((f for f in p["flow"] if f.get("hop") == hop), None)
            if hit:
                row[p["name"]] = {
                    "DN_VNC": hit.get("DN_VNC"),
                    "DN_SEZ": hit.get("DN_SEZ"),
                    "n_active": hit.get("n_active"),
                    "top": hit.get("top"),
                }
        compare[f"hop_{hop}"] = row
    report = {
        "organism": "Drosophila melanogaster first-instar larva",
        "n_cells": g["n"],
        "n_edges": g["n_edges"],
        "n_typed": g["n_typed"],
        "type_counts": g["type_counts"],
        "residual_Biochemistry": _R_BIO,
        "gaba_inhibitory": False,
        "free_parameters": 0,
        "phi": _PHI,
        "authority": (
            "Winding et al. Science 2023; all-all synapse counts + cell-type "
            "annotations (Supplementary-Data-S1). Unsigned — NT not in this dump."
        ),
        "programs": programs,
        "compare": compare,
        "note": (
            "Larval brain (not the adult Male CNS). DN-VNC is the brain→cord "
            "crawling command analog of adult DNg. Sensory vs olfactory vs "
            "mechanosensory on the same measured graph. Not a trained RNN."
        ),
    }
    out = ROOT / "data" / "larva_connectome_boot.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(compare, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
