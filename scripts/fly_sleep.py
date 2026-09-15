#!/usr/bin/env python3
"""Measured fly sleep observer → FSOT residual boot.

TriKinetics DAM2 beam-break counts (rethomics/damr example monitor M064).
File on D:\\FlyWire_Connectome\\behavior\\sleep (not git). 12:12 LD from the
light column. 2-minute bins.

Observer (not the field's free 5-minute cut):
  inactive = zero beam counts (measured)
  bout lengths of consecutive zeros
  MAD+φ on bout lengths; long leftover = sleep
  night vs day from the light column

Seed is measured Male CNS types, not invented sleep cells:
  ER* ellipsoid-body rings (R2/R5/ER3m), FB* fan-shaped body, LNv clock.

Not a trained RNN and not a thought. 0 free parameters.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

from fly_behavior import mad_phi_gate  # noqa: E402
from fly_connectome import (  # noqa: E402
    FLY_ROOT,
    _PHI,
    _R_BIO,
    cuda_name,
    residual_cascade,
    seed_indices,
)
from male_cns import load_male_graph  # noqa: E402

SLEEP = FLY_ROOT / "behavior" / "sleep"
DAM_NAME = "rethomics_M064.txt"
DAM_URL = "https://raw.githubusercontent.com/rethomics/damr/master/inst/extdata/M064.txt"
SOURCE = "https://github.com/rethomics/damr (DAM2 monitor M064)"

MALE_BOOT = ROOT / "data" / "male_cns_boot.json"
COURT_BOOT = ROOT / "data" / "fly_courtship_flow.json"
AGGR_BOOT = ROOT / "data" / "fly_aggression_flow.json"

PROGRAMS = [
    {"name": "ER", "seed": "er", "how": "type_prefix"},
    {"name": "FB", "seed": "fb", "how": "type_prefix"},
    {"name": "LNv", "seed": "lnv", "how": "substring"},
]


def ensure_dam() -> Path:
    SLEEP.mkdir(parents=True, exist_ok=True)
    path = SLEEP / DAM_NAME
    if path.exists() and path.stat().st_size > 10_000:
        return path
    print(f"  downloading {DAM_URL}", flush=True)
    urllib.request.urlretrieve(DAM_URL, path)
    print(f"  wrote {path} ({path.stat().st_size} bytes)", flush=True)
    return path


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    w = mask.astype(np.int8)
    edges = np.diff(np.concatenate(([0], w, [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    return [(int(s), int(e)) for s, e in zip(starts, ends)]


def load_dam(path: Path) -> dict[str, Any]:
    rows = [ln.split("\t") for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    ok = [r for r in rows if r[3] == "1" and len(r) >= 42]
    times = [datetime.strptime(r[1] + " " + r[2], "%d %b %y %H:%M:%S") for r in ok]
    light = np.array([int(r[9]) for r in ok], dtype=np.int8)
    ch = np.array([[int(x) for x in r[10:42]] for r in ok], dtype=np.float64)
    if len(times) > 1:
        dts = np.diff([t.timestamp() for t in times]) / 60.0
        dt = float(np.median(dts))
    else:
        dt = 1.0
    alive = (ch > 0).mean(axis=0) > 1.0 / (_PHI ** 5)
    return {
        "n_bins": int(ch.shape[0]),
        "n_channels": int(ch.shape[1]),
        "n_alive": int(alive.sum()),
        "bin_min": float(dt),
        "t0": times[0].isoformat(sep=" "),
        "t1": times[-1].isoformat(sep=" "),
        "light": light,
        "activity": ch,
        "alive": alive,
        "hours": np.array([t.hour + t.minute / 60.0 for t in times], dtype=np.float64),
    }


def observe(dam: dict[str, Any]) -> dict[str, Any]:
    ch = dam["activity"]
    light = dam["light"].astype(bool)
    alive = dam["alive"]
    dt = dam["bin_min"]
    flies = np.where(alive)[0]
    all_bouts: list[int] = []
    per_fly = []
    for i in flies:
        act = ch[:, i]
        z = act == 0
        bouts = _runs(z)
        lengths = [e - s for s, e in bouts]
        all_bouts.extend(lengths)
        per_fly.append(
            {
                "channel": int(i),
                "mean_counts": float(act.mean()),
                "frac_zero": float(z.mean()),
                "n_bouts": int(len(lengths)),
                "max_bout_bins": int(max(lengths) if lengths else 0),
            }
        )
    bl = np.asarray(all_bouts, dtype=np.float64) if all_bouts else np.zeros(1)
    # Weight by bout length: half of inactive *time* sits in leftover-long bouts.
    # Unweighted MAD+φ on 1-bin pauses has MAD≈0 and paints everything.
    weighted = np.repeat(bl, np.maximum(bl.astype(np.int64), 1)) if bl.size else np.zeros(1)
    gate = mad_phi_gate(weighted)
    sleep_cut = gate["threshold"]
    # sleep bins: inside a zero-run longer than MAD+φ leftover
    n = ch.shape[0]
    sleep = np.zeros((n, len(flies)), dtype=bool)
    for k, i in enumerate(flies):
        for s, e in _runs(ch[:, i] == 0):
            if (e - s) >= sleep_cut:
                sleep[s:e, k] = True
    night = ~light
    day = light
    frac_sleep = float(sleep.mean())
    frac_night = float(sleep[night].mean()) if night.any() else None
    frac_day = float(sleep[day].mean()) if day.any() else None
    # mean activity night vs day (alive)
    act_a = ch[:, alive]
    return {
        "source": SOURCE,
        "n_bins": dam["n_bins"],
        "bin_min": dt,
        "span_h": dam["n_bins"] * dt / 60.0,
        "t0": dam["t0"],
        "t1": dam["t1"],
        "n_alive": dam["n_alive"],
        "n_channels": dam["n_channels"],
        "frac_light": float(light.mean()),
        "bout_gate": {k: v for k, v in gate.items() if k != "bout_lengths"},
        "sleep_cut_bins": float(sleep_cut),
        "sleep_cut_min": float(sleep_cut * dt),
        "n_inactivity_bouts": int(bl.size),
        "frac_sleep": frac_sleep,
        "frac_sleep_night": frac_night,
        "frac_sleep_day": frac_day,
        "night_over_day": (
            None if not frac_day else (frac_night / frac_day if frac_night is not None else None)
        ),
        "mean_counts_night": float(act_a[night].mean()) if night.any() else None,
        "mean_counts_day": float(act_a[day].mean()) if day.any() else None,
        "per_fly": per_fly,
        "seed_map": {
            "sleep": {
                "seed": "ER / FB",
                "why": "long leftover inactivity; ellipsoid-body rings and fan-shaped body",
            },
            "wake": {
                "seed": "LNv / JO",
                "why": "beam counts above zero; clock LNv and walking JO as contrast",
            },
        },
        "note": (
            "Inactive = zero IR beam counts. Sleep = consecutive zeros longer "
            "than MAD+φ of bout lengths — not the field's free 5-minute cut. "
            "Light column is 12:12 LD."
        ),
        "series": {
            "light": light.astype(int).tolist(),
            "mean_counts": act_a.mean(axis=1).tolist(),
            "frac_sleep": sleep.mean(axis=1).tolist(),
        },
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
                "top": {
                    "cell_type": top.get("cell_type"),
                    "super_class": top.get("super_class"),
                    "a": top.get("a"),
                },
            }
        )
    return rows


def _slim_compare(compare: dict[str, Any]) -> dict[str, Any]:
    keep = {}
    for hop, row in compare.items():
        keep[hop] = {}
        for name, hit in row.items():
            keep[hop][name] = {
                k: hit.get(k)
                for k in (
                    "vnc_motor",
                    "cb_motor",
                    "motor",
                    "descending",
                    "DN",
                    "n_active",
                    "top",
                    "l1",
                )
                if k in hit
            }
    return keep


def run_sleep_hops(*, gpu: bool) -> dict[str, Any]:
    print("== Male CNS sleep seeds", flush=True)
    graph = load_male_graph()
    programs = []
    for spec in PROGRAMS:
        print(f"== {spec['name']} seed={spec['seed']} how={spec['how']}", flush=True)
        seed_i = seed_indices(graph, spec["seed"], how=spec["how"])
        print(f"  n_seed={len(seed_i)}", flush=True)
        run = residual_cascade(
            graph, seed_i, seed=spec["seed"], how=spec["how"], gpu=gpu
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
    for hop in (0, 1, 2, 3):
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
    contrast: dict[str, Any] = {}
    if MALE_BOOT.exists():
        male = json.loads(MALE_BOOT.read_text(encoding="utf-8"))
        contrast["walking_odor"] = _slim_compare(male.get("compare") or {})
    if COURT_BOOT.exists():
        court = json.loads(COURT_BOOT.read_text(encoding="utf-8"))
        contrast["courtship"] = _slim_compare(
            (court.get("connectome") or {}).get("compare") or {}
        )
    if AGGR_BOOT.exists():
        aggr = json.loads(AGGR_BOOT.read_text(encoding="utf-8"))
        contrast["aggression"] = _slim_compare(
            (aggr.get("connectome") or {}).get("compare") or {}
        )
    return {
        "dataset": "male-cns:v1.0",
        "n_neurons": graph["n"],
        "n_edges": graph["n_edges"],
        "n_gaba": graph["n_gaba"],
        "authority": graph["authority"],
        "programs": programs,
        "compare": compare,
        "contrast": contrast,
        "gpu_name": cuda_name(),
    }


def _hop_pick(compare: dict[str, Any], hop: str, name: str, key: str) -> float | None:
    row = (compare.get(hop) or {}).get(name) or {}
    v = row.get(key)
    return float(v) if v is not None else None


def loop_table(obs: dict[str, Any], hops: dict[str, Any]) -> dict[str, Any]:
    cmp_ = hops.get("compare") or {}
    walking = (hops.get("contrast") or {}).get("walking_odor") or {}
    return {
        "sleep_on": {
            "observer": "DAM zeros; MAD+φ leftover bout length",
            "sleep_cut_min": obs.get("sleep_cut_min"),
            "frac_sleep": obs.get("frac_sleep"),
            "frac_sleep_night": obs.get("frac_sleep_night"),
            "frac_sleep_day": obs.get("frac_sleep_day"),
            "night_over_day": obs.get("night_over_day"),
            "seed": "ER / FB",
            "male_cns_hop2_descending_ER": _hop_pick(cmp_, "hop_2", "ER", "descending"),
            "male_cns_hop2_vnc_motor_ER": _hop_pick(cmp_, "hop_2", "ER", "vnc_motor"),
            "male_cns_hop1_peak_ER": ((cmp_.get("hop_1") or {}).get("ER") or {}).get("top"),
            "male_cns_hop2_descending_FB": _hop_pick(cmp_, "hop_2", "FB", "descending"),
            "male_cns_hop1_peak_FB": ((cmp_.get("hop_1") or {}).get("FB") or {}).get("top"),
            "male_cns_hop2_descending_LNv": _hop_pick(cmp_, "hop_2", "LNv", "descending"),
            "male_cns_hop1_peak_LNv": ((cmp_.get("hop_1") or {}).get("LNv") or {}).get("top"),
        },
        "wake_contrast_JO": {
            "seed": "JO",
            "male_cns_hop2_vnc_motor_JO": _hop_pick(walking, "hop_2", "JO", "vnc_motor"),
        },
        "odor_contrast_olfactory": {
            "seed": "olfactory",
            "male_cns_hop2_vnc_motor": _hop_pick(walking, "hop_2", "olfactory", "vnc_motor"),
        },
        "honesty": (
            "The field's 5-minute sleep cut is not used. Sleep is consecutive "
            "zero beam counts longer than MAD+φ of bout lengths. ER* and FB* "
            "are measured type prefixes (ellipsoid body / fan-shaped body). "
            "LNv is the clock contrast. Do not invent a sleep neuron class. "
            "GABA only where annotated."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    ap.add_argument("--skip-hops", action="store_true")
    args = ap.parse_args(argv)
    print("FSOT fly sleep observer → residual boot", flush=True)
    print(f"  phi={_PHI} R_BIO={_R_BIO} gpu={cuda_name()}", flush=True)
    path = ensure_dam()
    print(f"  DAM {path}", flush=True)
    dam = load_dam(path)
    print(
        f"  bins={dam['n_bins']} bin_min={dam['bin_min']} alive={dam['n_alive']}/32 "
        f"{dam['t0']} → {dam['t1']}",
        flush=True,
    )
    obs = observe(dam)
    print(
        f"  sleep_cut={obs['sleep_cut_min']:.2f} min "
        f"frac_sleep={obs['frac_sleep']:.3f} "
        f"night={obs['frac_sleep_night']:.3f} day={obs['frac_sleep_day']:.3f} "
        f"night/day={obs['night_over_day']}",
        flush=True,
    )
    hops: dict[str, Any] = {}
    if not args.skip_hops:
        hops = run_sleep_hops(gpu=args.gpu)
    loop = loop_table(obs, hops)
    report = {
        "product": "measured DAM sleep observer + measured Male CNS types",
        "free_parameters": 0,
        "phi": _PHI,
        "residual_Biochemistry": _R_BIO,
        "authority": (
            "TriKinetics DAM2 beam counts (rethomics/damr M064) + Male CNS "
            "v1.0 ER* / FB* / LNv types"
        ),
        "honesty": (
            "Sleep is leftover-long inactivity on measured IR beam counts, "
            "not a 5-minute free cut and not a video pose net. ER/FB/LNv are "
            "named types in the dump. GABA only where annotated. Not a thought."
        ),
        "observer": {k: v for k, v in obs.items() if k != "series"},
        "series": obs.get("series"),
        "connectome": hops,
        "loop": loop,
        "gpu_name": cuda_name(),
    }
    out = ROOT / "data" / "fly_sleep_flow.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(loop, indent=2), flush=True)
    if obs["frac_sleep"] <= 0 or not obs.get("sleep_cut_min"):
        print("  sleep observer did not boot", flush=True)
        return 1
    if not hops.get("programs"):
        print("  hops did not boot", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
