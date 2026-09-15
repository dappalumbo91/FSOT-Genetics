#!/usr/bin/env python3
"""BANC female brain + VNC residual boot — measured graph.

Bates et al. Nature 2026. Files on D:\\FlyWire_Connectome\\banc (not git):

  banc_888_meta.feather
  banc_888_edgelist_simple_v3.feather

Harvard Dataverse 10.7910/DVN/7WTH1N; GCS
gs://lee-lab_brain-and-nerve-cord-fly-connectome/compiled_data/banc_888/

Female CNS, intact neck. GABA from neurotransmitter_predicted (same
authority as FlyWire / Male CNS consensus_nt). Super_class motor in
the VNC is reported as vnc_motor so the walking split is comparable
to Male CNS — region × class, both measured, not invented.

Drop glia / trachea / not_a_neuron. 0 free parameters.
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

from fly_connectome import residual_cascade, seed_indices  # noqa: E402

BANC = Path(r"D:\FlyWire_Connectome\banc")
META = BANC / "banc_888_meta.feather"
EDGES = BANC / "banc_888_edgelist_simple_v3.feather"

_DROP = frozenset({"glia", "not_a_neuron", "trachea"})

PROGRAMS = [
    {"name": "vnc_sensory", "seed": "vnc_sensory", "how": "class"},
    {"name": "chordotonal", "seed": "chordotonal", "how": "substring"},
    {"name": "JO", "seed": "jo", "how": "type_prefix"},
    {"name": "olfactory", "seed": "olfactory_receptor_neuron", "how": "class"},
]


def _map_super(sc: str, region: str) -> str:
    sc = (sc or "").strip().lower()
    rg = (region or "").strip().lower()
    if sc == "motor":
        if "ventral" in rg:
            return "vnc_motor"
        if "central" in rg:
            return "cb_motor"
        return "motor"
    if sc == "sensory":
        if "ventral" in rg:
            return "vnc_sensory"
        return "sensory"
    return sc


def load_banc_graph() -> dict[str, Any]:
    import pandas as pd
    from scipy.sparse import csr_matrix

    print(f"  reading {META.name}", flush=True)
    meta_df = pd.read_feather(META)
    sc_raw = meta_df["super_class"].fillna("").astype(str)
    keep = ~sc_raw.str.lower().isin(_DROP)
    meta_df = meta_df.loc[keep].copy()
    print(f"  cells after drop glia/trachea/not_a_neuron {len(meta_df)}", flush=True)

    def _s(col: str) -> list[str]:
        if col not in meta_df.columns:
            return [""] * len(meta_df)
        return (
            meta_df[col]
            .fillna("")
            .astype(str)
            .replace({"nan": "", "None": ""})
            .tolist()
        )

    ids = [str(x) for x in meta_df["banc_888_id"].tolist()]
    idx = {rid: i for i, rid in enumerate(ids)}
    super_raw = _s("super_class")
    region = _s("region")
    mapped = [_map_super(s, r) for s, r in zip(super_raw, region)]
    klass = _s("cell_class")
    types = _s("cell_type")
    sub = _s("cell_sub_class")
    flow = _s("flow")
    side = _s("side")
    nt = [x.lower() for x in _s("neurotransmitter_predicted")]
    dim = _s("sexually_dimorphic")
    meta: dict[str, dict[str, str]] = {}
    for i, rid in enumerate(ids):
        meta[rid] = {
            "super_class": mapped[i],
            "banc_super_class": super_raw[i],
            "region": region[i],
            "cell_class": klass[i],
            "cell_type": types[i],
            "top_nt": nt[i],
            "flow": flow[i],
            "sub_class": sub[i],
            "side": side[i],
            "fru_dsx": "",
            "dimorphism": dim[i],
        }

    print(f"  reading {EDGES.name}", flush=True)
    wdf = pd.read_feather(EDGES)
    pre = wdf["pre"].astype(str)
    post = wdf["post"].astype(str)
    ok = pre.isin(idx) & post.isin(idx)
    wdf = wdf.loc[ok]
    print(f"  neuron-neuron edges {len(wdf)}", flush=True)

    pre_i = pre.loc[ok].map(idx).to_numpy(dtype=np.int64)
    post_i = post.loc[ok].map(idx).to_numpy(dtype=np.int64)
    wt = wdf["count"].to_numpy(dtype=np.float64)
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
        "key": "banc:v888",
        "ids": ids,
        "idx": idx,
        "meta": meta,
        "W": W,
        "n": n,
        "n_edges": int(len(pre_i)),
        "n_gaba": n_gaba,
        "authority": (
            "BANC v888 Bates et al. Nature 2026; female brain+VNC; "
            "GABA from neurotransmitter_predicted; weights = measured "
            "synapse counts (edgelist_simple_v3); glia/trachea dropped; "
            "vnc_motor = super_class motor ∩ region VNC"
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
                    "cell_class": top.get("cell_class"),
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
    for p in (META, EDGES):
        if not p.exists() or p.stat().st_size < 1_000_000:
            raise SystemExit(f"missing {p}")
    graph = load_banc_graph()
    programs = []
    for spec in PROGRAMS:
        print(
            f"== {spec['name']} seed={spec['seed']} how={spec['how']}",
            flush=True,
        )
        seed_i = seed_indices(graph, spec["seed"], how=spec["how"])
        print(f"  n_seed={len(seed_i)}", flush=True)
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
        "dataset": "banc:v888",
        "n_neurons": graph["n"],
        "n_edges": graph["n_edges"],
        "n_gaba": graph["n_gaba"],
        "authority": graph["authority"],
        "free_parameters": 0,
        "programs": programs,
        "compare": compare,
        "note": (
            "Female BANC brain+VNC, intact neck. vnc_motor = motor in the "
            "cord. Walking observer analog is vnc_sensory / chordotonal / JO. "
            "Olfactory is the local-hub contrast. GABA inhibitory. "
            "Not a trained RNN and not an invented bee brain."
        ),
    }
    out = ROOT / "data" / "banc_connectome_boot.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(compare, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
