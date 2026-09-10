#!/usr/bin/env python3
"""Male CNS (brain + VNC) residual boot — measured graph.

Berg et al. Cell 2026. Files on D:\\FlyWire_Connectome\\male_cns (not git):

  body-annotations-male-cns-v1.0-minconf-0.5.feather
  body-neurotransmitters-male-cns-v1.0.feather
  connectome-weights-male-cns-v1.0-minconf-0.5.feather

Traced neurons only. GABA from consensus_nt. Walking seed is
vnc_sensory (leg/body afferents in the cord) and mechanosensory
class. Report vnc_motor — those are the leg motor neurons.

Same residual as FlyWire v630. 0 free parameters.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

from fly_connectome import (  # noqa: E402
    FLY_ROOT,
    residual_cascade,
    seed_indices,
)

MALE = FLY_ROOT / "male_cns"
ANN = MALE / "body-annotations-male-cns-v1.0-minconf-0.5.feather"
NT = MALE / "body-neurotransmitters-male-cns-v1.0.feather"
WTS = MALE / "connectome-weights-male-cns-v1.0-minconf-0.5.feather"

PROGRAMS = [
    {"name": "vnc_sensory", "seed": "vnc_sensory", "how": "class"},
    {"name": "mechanosensory", "seed": "mechanosensory", "how": "substring"},
    {"name": "JO", "seed": "jo", "how": "type_prefix"},
    {"name": "olfactory", "seed": "olfactory", "how": "class"},
]


def load_male_graph() -> dict[str, Any]:
    import pandas as pd
    from scipy.sparse import csr_matrix

    print(f"  reading {ANN.name}", flush=True)
    ann = pd.read_feather(ANN)
    traced = ann[ann["status"] == "Traced"].copy()
    print(f"  traced {len(traced)} / {len(ann)}", flush=True)
    print(f"  reading {NT.name}", flush=True)
    nt = pd.read_feather(NT)
    nt = nt.drop_duplicates("body")
    nt_map = dict(zip(nt["body"].to_numpy(), nt["consensus_nt"].fillna("").astype(str).str.lower()))
    print(f"  reading {WTS.name}", flush=True)
    wdf = pd.read_feather(WTS)
    traced_ids = set(traced["bodyId"].to_numpy().tolist())
    ok = wdf["body_pre"].isin(traced_ids) & wdf["body_post"].isin(traced_ids)
    wdf = wdf.loc[ok]
    print(f"  traced edges {len(wdf)}", flush=True)

    def _s(col: str) -> list[str]:
        if col not in traced.columns:
            return [""] * len(traced)
        return traced[col].fillna("").astype(str).replace({"nan": ""}).tolist()

    bodies = [int(b) for b in traced["bodyId"].to_numpy()]
    ids = [str(b) for b in bodies]
    idx = {rid: i for i, rid in enumerate(ids)}
    super_c = _s("superclass")
    klass = _s("class")
    types = _s("type")
    sub = _s("subclass")
    side = _s("somaSide")
    meta: dict[str, dict[str, str]] = {}
    for i, bid in enumerate(bodies):
        rid = ids[i]
        meta[rid] = {
            "super_class": super_c[i],
            "cell_class": klass[i],
            "cell_type": types[i],
            "top_nt": (nt_map.get(bid) or "").lower(),
            "flow": "",
            "sub_class": sub[i],
            "side": side[i],
        }

    pre = wdf["body_pre"].map(lambda x: str(int(x)))
    post = wdf["body_post"].map(lambda x: str(int(x)))
    wt = wdf["weight"].to_numpy(dtype=np.float64)
    pre_i = pre.map(idx).to_numpy(dtype=np.int64)
    post_i = post.map(idx).to_numpy(dtype=np.int64)
    is_gaba = np.zeros(len(ids), dtype=bool)
    for rid, m in meta.items():
        if m.get("top_nt") == "gaba":
            is_gaba[idx[rid]] = True
    sign = np.where(is_gaba[pre_i], -1.0, 1.0)
    n = len(ids)
    W = csr_matrix((wt * sign, (post_i, pre_i)), shape=(n, n))
    n_gaba = int(is_gaba.sum())
    print(f"  GABA presynaptic cells {n_gaba}", flush=True)
    return {
        "key": "male-cns:v1.0",
        "ids": ids,
        "idx": idx,
        "meta": meta,
        "W": W,
        "n": n,
        "n_edges": int(len(pre_i)),
        "n_gaba": n_gaba,
        "authority": (
            "Male CNS v1.0 Berg et al. Cell 2026; traced neurons only; "
            "GABA from consensus_nt; weights = measured synapse counts"
        ),
    }


def _trace_row(run: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for snap in run.get("trace") or []:
        tm = snap.get("target_mass") or {}
        top = (snap.get("top") or [{}])[0]
        rows.append(
            {
                "hop": snap.get("hop"),
                "l1": snap.get("l1"),
                "n_active": snap.get("n_active"),
                "motor": tm.get("motor"),
                "vnc_motor": tm.get("vnc_motor"),
                "cb_motor": tm.get("cb_motor"),
                "descending": tm.get("descending"),
                "vnc_sensory": tm.get("vnc_sensory"),
                "top": {
                    "cell_type": top.get("cell_type"),
                    "super_class": top.get("super_class"),
                    "a": top.get("a"),
                },
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    args = ap.parse_args(argv)
    for p in (ANN, NT, WTS):
        if not p.exists():
            raise SystemExit(f"missing {p}")
    graph = load_male_graph()
    programs = []
    for spec in PROGRAMS:
        print(
            f"== {spec['name']} seed={spec['seed']} how={spec['how']}",
            flush=True,
        )
        seed_i = seed_indices(graph, spec["seed"], how=spec["how"])
        run = residual_cascade(
            graph,
            seed_i,
            seed=spec["seed"],
            how=spec["how"],
            gpu=args.gpu,
        )
        programs.append(
            {
                "name": spec["name"],
                "seed": spec["seed"],
                "how": spec["how"],
                "n_seed": run["n_seed"],
                "gpu": run.get("gpu"),
                "device": run.get("device"),
                "elapsed_s": run.get("elapsed_s"),
                "flow": _trace_row(run),
            }
        )
    compare = {}
    for hop in (1, 2, 3):
        row = {}
        for p in programs:
            hit = next((f for f in p["flow"] if f.get("hop") == hop), None)
            if hit:
                row[p["name"]] = {
                    "vnc_motor": hit.get("vnc_motor"),
                    "cb_motor": hit.get("cb_motor"),
                    "motor": hit.get("motor"),
                    "descending": hit.get("descending"),
                    "n_active": hit.get("n_active"),
                    "top": hit.get("top"),
                }
        compare[f"hop_{hop}"] = row
    report = {
        "dataset": "male-cns:v1.0",
        "n_neurons": graph["n"],
        "n_edges": graph["n_edges"],
        "n_gaba": graph["n_gaba"],
        "authority": graph["authority"],
        "free_parameters": 0,
        "programs": programs,
        "compare": compare,
        "note": (
            "Traced Male CNS brain+VNC. vnc_motor are leg/body motor neurons. "
            "Walking observer seeds vnc_sensory / mechanosensory. "
            "GABA inhibitory. Not a trained RNN."
        ),
    }
    out = ROOT / "data" / "male_cns_boot.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(compare, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
