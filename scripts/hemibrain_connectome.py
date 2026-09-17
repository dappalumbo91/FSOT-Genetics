#!/usr/bin/env python3
"""Janelia hemibrain residual boot — measured central-brain graph.

Scheffer et al. eLife 2020. Live graph: neuPrint hemibrain:v1.2.1 typed neurons
(including cropped Leaves so JO / ORN afferents are present). Compact GCS dump
v1.2 drops cropped cells — JO count there is 0, so it is the wrong graph for
the JO vs olfactory question.

No VNC in this volume. Report descending, not vnc_motor.
DNg29 is absent (whole-CNS type; not in the hemibrain cut).
Neuron properties have no predictedNt — unsigned residual, no invented GABA.

Data cache on D:\\FlyWire_Connectome\\hemibrain (not git). 0 free parameters.
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

from fly_connectome import residual_cascade, seed_indices  # noqa: E402

CACHE = Path(r"D:\FlyWire_Connectome\hemibrain")
NEURON_FEA = CACHE / "neuprint_v1_2_1_typed_neurons.feather"
EDGE_FEA = CACHE / "neuprint_v1_2_1_typed_edges.feather"
NEUPRINT = "https://neuprint.janelia.org/api/custom/custom"
DATASET = "hemibrain:v1.2.1"
_DN = re.compile(r"^DN[a-z]")

PROGRAMS = [
    {"name": "JO", "seed": "JO", "how": "type_prefix"},
    {"name": "olfactory", "seed": "ORN", "how": "type_prefix"},
]


def _cypher(q: str, timeout: int = 180) -> dict[str, Any]:
    req = urllib.request.Request(
        NEUPRINT,
        data=json.dumps({"cypher": q, "dataset": DATASET}).encode(),
        headers={
            "User-Agent": "FSOT-fly-pack",
            "Accept": "application/json",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    last: Exception | None = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.loads(fh.read().decode())
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            last = e
            time.sleep(2.0 * (attempt + 1))
    raise SystemExit(f"neuPrint failed: {last}")


def _family(typ: str) -> tuple[str, str]:
    t = (typ or "").strip()
    if t.startswith("JO"):
        return "sensory", "mechanosensory"
    if t.startswith("ORN"):
        return "sensory", "olfactory"
    # Giant Fiber is the named hemibrain descending escape neuron (not a DN* type string).
    if t == "Giant Fiber" or t.startswith("GF"):
        return "descending", "DN"
    if _DN.match(t):
        return "descending", "DN"
    if t.startswith("KC"):
        return "", "kenyon"
    return "", ""


def fetch_typed_graph() -> None:
    import pandas as pd

    CACHE.mkdir(parents=True, exist_ok=True)
    print("  neuPrint typed neurons", flush=True)
    raw = _cypher(
        "MATCH (n:Neuron) WHERE n.type IS NOT NULL "
        "RETURN n.bodyId AS bodyId, n.type AS type, n.instance AS instance, "
        "n.status AS status, n.statusLabel AS statusLabel, n.cropped AS cropped"
    )
    cols = raw.get("columns") or []
    neurons = pd.DataFrame(raw.get("data") or [], columns=cols)
    neurons.to_feather(NEURON_FEA)
    print(f"  wrote {NEURON_FEA.name} n={len(neurons)}", flush=True)
    ids = [int(x) for x in neurons["bodyId"].tolist()]
    rows: list[tuple[int, int, int]] = []
    batch = 200
    n_batch = (len(ids) + batch - 1) // batch
    for b, i in enumerate(range(0, len(ids), batch), start=1):
        chunk = ids[i : i + batch]
        id_list = ",".join(str(x) for x in chunk)
        q = (
            "MATCH (a:Neuron)-[c:ConnectsTo]->(b:Neuron) "
            f"WHERE a.bodyId IN [{id_list}] AND b.type IS NOT NULL "
            "RETURN a.bodyId AS pre, b.bodyId AS post, c.weight AS weight"
        )
        got = _cypher(q, timeout=240)
        for pre, post, w in got.get("data") or []:
            rows.append((int(pre), int(post), int(w)))
        if b == 1 or b % 10 == 0 or b == n_batch:
            print(f"  edges batch {b}/{n_batch} cumulative={len(rows)}", flush=True)
    edges = pd.DataFrame(rows, columns=["pre", "post", "weight"])
    edges.to_feather(EDGE_FEA)
    print(f"  wrote {EDGE_FEA.name} n={len(edges)}", flush=True)


def load_hemibrain_graph() -> dict[str, Any]:
    import pandas as pd
    from scipy.sparse import csr_matrix

    if not NEURON_FEA.is_file() or not EDGE_FEA.is_file():
        print("  cache miss — fetching typed graph from neuPrint", flush=True)
        fetch_typed_graph()
    print(f"  reading {NEURON_FEA.name}", flush=True)
    neurons = pd.read_feather(NEURON_FEA)
    print(f"  reading {EDGE_FEA.name}", flush=True)
    edges = pd.read_feather(EDGE_FEA)
    bodies = [int(b) for b in neurons["bodyId"].to_numpy()]
    ids = [str(b) for b in bodies]
    idx = {rid: i for i, rid in enumerate(ids)}
    types = neurons["type"].fillna("").astype(str).tolist()
    inst = neurons["instance"].fillna("").astype(str).tolist()
    lab = (
        neurons["statusLabel"].fillna("").astype(str).tolist()
        if "statusLabel" in neurons.columns
        else [""] * len(ids)
    )
    meta: dict[str, dict[str, str]] = {}
    n_dn = 0
    n_jo = 0
    n_orn = 0
    n_kc = 0
    for i, rid in enumerate(ids):
        sc, cc = _family(types[i])
        if sc == "descending":
            n_dn += 1
        if types[i].startswith("JO"):
            n_jo += 1
        if types[i].startswith("ORN"):
            n_orn += 1
        if cc == "kenyon":
            n_kc += 1
        meta[rid] = {
            "super_class": sc,
            "cell_class": cc,
            "cell_type": types[i],
            "top_nt": "",
            "flow": "",
            "sub_class": lab[i],
            "side": "",
            "instance": inst[i],
        }
    pre = edges["pre"].map(lambda x: str(int(x)))
    post = edges["post"].map(lambda x: str(int(x)))
    wt = edges["weight"].to_numpy(dtype=np.float64)
    ok = pre.isin(idx) & post.isin(idx)
    pre_i = pre[ok].map(idx).to_numpy(dtype=np.int64)
    post_i = post[ok].map(idx).to_numpy(dtype=np.int64)
    w = wt[ok]
    n = len(ids)
    W = csr_matrix((w, (post_i, pre_i)), shape=(n, n))
    print(
        f"  typed {n}  edges {int(ok.sum())}  JO={n_jo} ORN={n_orn} "
        f"DN[a-z]={n_dn} KC={n_kc}  unsigned (no NT)",
        flush=True,
    )
    return {
        "key": DATASET,
        "ids": ids,
        "idx": idx,
        "meta": meta,
        "W": W,
        "n": n,
        "n_edges": int(ok.sum()),
        "n_gaba": 0,
        "n_jo": n_jo,
        "n_orn": n_orn,
        "n_dn": n_dn,
        "n_kc": n_kc,
        "authority": (
            "neuPrint hemibrain:v1.2.1 typed neurons including cropped Leaves; "
            "weight = measured synapse count; unsigned residual (no predictedNt); "
            "DNg29 not in this volume; no VNC"
        ),
    }


def _trace_row(run: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for snap in run.get("trace") or []:
        tm = snap.get("target_mass") or {}
        cm = snap.get("class_mass") or {}
        top = (snap.get("top") or [{}])[0]
        rows.append(
            {
                "hop": snap.get("hop"),
                "l1": snap.get("l1"),
                "n_active": snap.get("n_active"),
                "motor": tm.get("motor"),
                "vnc_motor": tm.get("vnc_motor"),
                "descending": tm.get("descending"),
                "sensory": tm.get("sensory"),
                "DN": cm.get("DN"),
                "olfactory": cm.get("olfactory"),
                "mechanosensory": cm.get("mechanosensory"),
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
    ap.add_argument("--fetch", action="store_true", help="refresh neuPrint cache")
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    args = ap.parse_args(argv)
    if args.fetch and NEURON_FEA.exists():
        NEURON_FEA.unlink()
        if EDGE_FEA.exists():
            EDGE_FEA.unlink()
    graph = load_hemibrain_graph()
    dng29 = sum(
        1 for m in graph["meta"].values() if (m.get("cell_type") or "") == "DNg29"
    )
    programs = []
    for spec in PROGRAMS:
        print(f"== {spec['name']} seed={spec['seed']} how={spec['how']}", flush=True)
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
    compare: dict[str, Any] = {}
    for hop in (1, 2, 3):
        row = {}
        for p in programs:
            hit = next((f for f in p["flow"] if f.get("hop") == hop), None)
            if hit:
                row[p["name"]] = {
                    "descending": hit.get("descending"),
                    "vnc_motor": hit.get("vnc_motor"),
                    "motor": hit.get("motor"),
                    "DN": hit.get("DN"),
                    "n_active": hit.get("n_active"),
                    "top": hit.get("top"),
                }
        compare[f"hop_{hop}"] = row
    report = {
        "dataset": DATASET,
        "n_neurons": graph["n"],
        "n_edges": graph["n_edges"],
        "n_gaba": 0,
        "n_jo": graph["n_jo"],
        "n_orn": graph["n_orn"],
        "n_dn": graph["n_dn"],
        "n_kc": graph["n_kc"],
        "n_dng29": dng29,
        "authority": graph["authority"],
        "free_parameters": 0,
        "programs": programs,
        "compare": compare,
        "note": (
            "Hemibrain is a partial central brain — no VNC, so vnc_motor stays 0. "
            "JO vs olfactory is scored on descending (DN[a-z] types; clock DN1/2/3 excluded). "
            "DNg29 is not in this volume. Unsigned residual: neuPrint Neuron has no predictedNt. "
            "Compact GCS dump dropped cropped JO/ORN; this boot uses typed neuPrint including Leaves."
        ),
    }
    out = ROOT / "data" / "hemibrain_connectome_boot.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(compare, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
