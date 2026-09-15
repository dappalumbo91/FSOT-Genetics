#!/usr/bin/env python3
"""Measured fly courtship observer → FSOT residual boot.

Pan et al. PLOS ONE 10.1371/journal.pone.0021144 Movie S1: solitary
male UAS-dTrpA1 / fru-GAL4 at 29 °C. Wing extension, abdomen bending,
copulation attempts. File on D:\\FlyWire_Connectome\\behavior\\courtship
(not git).

Observer (not a trained pose net):
  circular arena ≥ observer-max / φ
  leftover dark CC (drop components ≥ max/φ) = the fly
  MAD+φ on frame-diff energy and fly Rg = courtship motor on
  whole clip is courtship (like Harvard walking); gate splits high vs low

Seed is measured Male CNS genetics, not invented P1:
  fru/dsx labeled cells, type prefix pC1_ (not LLPC1), TN1 song MNs.

Hops vs JO / olfactory from the live boots. 0 free parameters.
"""
from __future__ import annotations

import json
import sys
import urllib.request
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
from fly_odor import decode_rgb  # noqa: E402
from male_cns import load_male_graph  # noqa: E402

COURT = FLY_ROOT / "behavior" / "courtship"
VIDEO_NAME = "pone.0021144.s006.mov"
VIDEO_URL = (
    "https://journals.plos.org/plosone/article/file"
    "?id=10.1371/journal.pone.0021144.s006&type=supplementary"
)
PLOS = "https://doi.org/10.1371/journal.pone.0021144"
GENOTYPE = "UAS-dTrpA1 / fru-GAL4(B), 29 C"

MALE_BOOT = ROOT / "data" / "male_cns_boot.json"
WALK_BOOT = ROOT / "data" / "fly_behavior_flow.json"
ODOR_BOOT = ROOT / "data" / "fly_odor_flow.json"

PROGRAMS = [
    {"name": "fru_dsx", "seed": "any", "how": "fru_dsx"},
    {"name": "pC1", "seed": "pc1_", "how": "type_prefix"},
    {"name": "TN1", "seed": "tn1", "how": "type_prefix"},
]


def ensure_video() -> Path:
    COURT.mkdir(parents=True, exist_ok=True)
    path = COURT / VIDEO_NAME
    if path.exists() and path.stat().st_size > 100_000:
        return path
    print(f"  downloading {VIDEO_URL}", flush=True)
    urllib.request.urlretrieve(VIDEO_URL, path)
    print(f"  wrote {path} ({path.stat().st_size} bytes)", flush=True)
    return path


def _fly_on_frame(gray: np.ndarray) -> dict[str, float | None]:
    """Leftover dark CC inside the bright disk. Drop CCs ≥ observer-max / φ."""
    from scipy import ndimage

    peak = float(gray.max()) + 1e-12
    arena = gray >= (peak / _PHI)
    if int(arena.sum()) < 50:
        return {"area": None, "rg": None, "y": None, "x": None}
    med = float(np.median(gray[arena]))
    mad = float(np.median(np.abs(gray[arena] - med))) + 1e-12
    thr = med - _PHI * mad
    dark = np.zeros(gray.shape, dtype=bool)
    dark[arena] = gray[arena] <= thr
    paint = float(dark[arena].mean())
    if paint > 1.0 / _PHI:
        thr = med - (_PHI * _PHI) * mad
        dark[arena] = gray[arena] <= thr
    labeled, nlab = ndimage.label(dark)
    if nlab == 0:
        return {"area": None, "rg": None, "y": None, "x": None}
    sizes = ndimage.sum(dark, labeled, index=np.arange(1, nlab + 1))
    mx = float(np.max(sizes))
    leftover = [i + 1 for i, s in enumerate(sizes) if 2.0 <= float(s) < mx / _PHI]
    if not leftover:
        leftover = [int(np.argmax(sizes)) + 1]
    best = max(leftover, key=lambda i: float(sizes[i - 1]))
    ys, xs = np.where(labeled == best)
    if ys.size < 2:
        return {"area": None, "rg": None, "y": None, "x": None}
    cy = float(ys.mean())
    cx = float(xs.mean())
    rg = float(np.sqrt(((ys - cy) ** 2 + (xs - cx) ** 2).mean()))
    return {"area": float(ys.size), "rg": rg, "y": cy, "x": cx}


def observe(vol: np.ndarray, meta: dict[str, Any]) -> dict[str, Any]:
    gray = vol.mean(axis=3)
    n = int(gray.shape[0])
    energy = np.mean(np.abs(np.diff(gray, axis=0)), axis=(1, 2))
    rgs = np.full(n, np.nan)
    areas = np.full(n, np.nan)
    for t in range(n):
        fly = _fly_on_frame(gray[t])
        if fly["rg"] is not None:
            rgs[t] = fly["rg"]
            areas[t] = fly["area"]
    e_full = np.zeros(n, dtype=np.float64)
    e_full[1:] = energy
    gate_e = mad_phi_gate(energy)
    rg_ok = rgs[np.isfinite(rgs)]
    gate_rg = mad_phi_gate(rg_ok) if rg_ok.size > 8 else mad_phi_gate(np.zeros(1))
    moving = np.zeros(n, dtype=bool)
    moving[1:] = energy >= gate_e["threshold"]
    spread = np.zeros(n, dtype=bool)
    spread[np.isfinite(rgs)] = rgs[np.isfinite(rgs)] >= gate_rg["threshold"]
    on = moving | spread
    rest_floor = float(np.median(energy)) / (_PHI ** 5)
    whole_clip_courtship = float(energy.min()) > rest_floor
    fps = float(meta.get("fps") or 0.0) or 30.0
    return {
        "source": PLOS,
        "genotype": GENOTYPE,
        "video": meta,
        "n_tracked": int(np.isfinite(rgs).sum()),
        "motion_gate": {k: v for k, v in gate_e.items() if k != "bout_lengths"},
        "rg_gate": {k: v for k, v in gate_rg.items() if k != "bout_lengths"},
        "n_motion_on": int(moving.sum()),
        "n_spread_on": int(spread.sum()),
        "n_courtship_on": int(on.sum()),
        "paint_motion": float(moving.mean()),
        "paint_spread": float(spread.mean()),
        "paint_courtship": float(on.mean()),
        "energy_min": float(energy.min()),
        "energy_median": float(np.median(energy)),
        "energy_mean": float(energy.mean()),
        "rg_median": float(np.nanmedian(rgs)),
        "rg_mean": float(np.nanmean(rgs)),
        "rest_floor": rest_floor,
        "whole_clip_is_courtship": whole_clip_courtship,
        "duration_s": float(n) / fps if fps else meta.get("duration_s"),
        "seed_map": {
            "courtship_on": {
                "seed": "fru_dsx labeled / pC1_ / TN1",
                "why": "Movie S1 is fru-GAL4 courtship motor; seed measured genetics",
            },
            "low": {
                "seed": None,
                "why": "MAD+φ low energy/Rg inside a courtship clip; not a rest trial",
            },
        },
        "series": {
            "energy": e_full.tolist(),
            "rg": [None if not np.isfinite(v) else float(v) for v in rgs],
            "area": [None if not np.isfinite(v) else float(v) for v in areas],
            "motion_on": moving.astype(bool).tolist(),
            "spread_on": spread.astype(bool).tolist(),
            "courtship_on": on.astype(bool).tolist(),
        },
        "note": (
            "Leftover dark CC inside the arena is the fly (drop CCs ≥ max/φ, "
            "the chamber ring). MAD+φ on energy and Rg. Not a trained pose net. "
            "This clip has no rest state — it is a courtship trial."
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
                "top": {
                    "cell_type": top.get("cell_type"),
                    "super_class": top.get("super_class"),
                    "cell_class": top.get("cell_class"),
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


def run_courtship_hops(*, gpu: bool) -> dict[str, Any]:
    print("== Male CNS courtship seeds", flush=True)
    graph = load_male_graph()
    programs = []
    for spec in PROGRAMS:
        print(f"== {spec['name']} seed={spec['seed']} how={spec['how']}", flush=True)
        seed_i = seed_indices(graph, spec["seed"], how=spec["how"])
        print(f"  n_seed={len(seed_i)}", flush=True)
        run = residual_cascade(
            graph,
            seed_i,
            seed=spec["seed"],
            how=spec["how"],
            gpu=gpu,
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
                    "vnc_motor": hit.get("vnc_motor"),
                    "cb_motor": hit.get("cb_motor"),
                    "motor": hit.get("motor"),
                    "descending": hit.get("descending"),
                    "n_active": hit.get("n_active"),
                    "top": hit.get("top"),
                }
        compare[f"hop_{hop}"] = row
    contrast = {}
    if MALE_BOOT.exists():
        male = json.loads(MALE_BOOT.read_text(encoding="utf-8"))
        contrast["male_cns_walking_odor"] = _slim_compare(male.get("compare") or {})
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
    walking = (hops.get("contrast") or {}).get("male_cns_walking_odor") or {}
    return {
        "courtship_on": {
            "observer": "fru-GAL4 Movie S1; MAD+φ energy/Rg",
            "n_on": obs.get("n_courtship_on"),
            "paint": obs.get("paint_courtship"),
            "whole_clip_is_courtship": obs.get("whole_clip_is_courtship"),
            "seed": "fru_dsx / pC1_ / TN1",
            "male_cns_hop2_vnc_motor_fru_dsx": _hop_pick(cmp_, "hop_2", "fru_dsx", "vnc_motor"),
            "male_cns_hop2_descending_fru_dsx": _hop_pick(cmp_, "hop_2", "fru_dsx", "descending"),
            "male_cns_hop1_peak_fru_dsx": ((cmp_.get("hop_1") or {}).get("fru_dsx") or {}).get("top"),
            "male_cns_hop2_vnc_motor_pC1": _hop_pick(cmp_, "hop_2", "pC1", "vnc_motor"),
            "male_cns_hop1_peak_pC1": ((cmp_.get("hop_1") or {}).get("pC1") or {}).get("top"),
            "male_cns_hop2_vnc_motor_TN1": _hop_pick(cmp_, "hop_2", "TN1", "vnc_motor"),
            "male_cns_hop1_peak_TN1": ((cmp_.get("hop_1") or {}).get("TN1") or {}).get("top"),
        },
        "walking_contrast_JO": {
            "seed": "JO / vnc_sensory",
            "male_cns_hop2_vnc_motor_JO": _hop_pick(walking, "hop_2", "JO", "vnc_motor"),
            "male_cns_hop2_vnc_motor_vnc_sensory": _hop_pick(
                walking, "hop_2", "vnc_sensory", "vnc_motor"
            ),
        },
        "odor_contrast_olfactory": {
            "seed": "olfactory",
            "male_cns_hop2_vnc_motor": _hop_pick(walking, "hop_2", "olfactory", "vnc_motor"),
            "male_cns_hop1_peak": ((walking.get("hop_1") or {}).get("olfactory") or {}).get("top"),
        },
        "honesty": (
            "fru/dsx is measured genetics on 5012 traced cells, not a command "
            "singleton. pC1_ is the type prefix in this dump (not LLPC1). "
            "TN1 is a measured song motor type. Do not invent P1 from "
            "morphology. This clip is courtship throughout; MAD+φ splits "
            "high vs low motor, it does not invent a rest state."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    ap.add_argument("--skip-hops", action="store_true")
    args = ap.parse_args(argv)
    print("FSOT fly courtship observer → residual boot", flush=True)
    print(f"  phi={_PHI} R_BIO={_R_BIO} gpu={cuda_name()}", flush=True)
    video = ensure_video()
    print(f"  video {video}", flush=True)
    vol, meta = decode_rgb(video)
    print(
        f"  frames={meta['n_frames']} {meta['width']}x{meta['height']} "
        f"{meta['duration_s']:.2f}s fps={meta['fps']}",
        flush=True,
    )
    obs = observe(vol, meta)
    print(
        f"  tracked={obs['n_tracked']}/{meta['n_frames']} "
        f"courtship_on={obs['n_courtship_on']} paint={obs['paint_courtship']:.3f} "
        f"whole_clip={obs['whole_clip_is_courtship']}",
        flush=True,
    )
    hops: dict[str, Any] = {}
    if not args.skip_hops:
        hops = run_courtship_hops(gpu=args.gpu)
    loop = loop_table(obs, hops)
    report = {
        "product": "measured courtship observer + measured Male CNS genetics",
        "free_parameters": 0,
        "phi": _PHI,
        "residual_Biochemistry": _R_BIO,
        "authority": (
            "Pan et al. PLOS ONE 2011 Movie S1 (fru-GAL4 dTrpA1 courtship) + "
            "Male CNS v1.0 fruDsx / type annotations"
        ),
        "honesty": (
            "Harvard walking and eLife odor clips are other programs. This "
            "clip is a fru-GAL4 male performing courtship motor (wing, "
            "abdomen). Seed is measured fru/dsx, pC1_ types, and TN1 — not "
            "an invented P1 list. GABA only where annotated. Not a trained "
            "pose net, not a trained RNN, not a thought."
        ),
        "observer": {k: v for k, v in obs.items() if k != "series"},
        "series": obs.get("series"),
        "connectome": hops,
        "loop": loop,
        "gpu_name": cuda_name(),
    }
    out = ROOT / "data" / "fly_courtship_flow.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(loop, indent=2), flush=True)
    if obs["n_tracked"] < meta["n_frames"] / _PHI:
        print("  observer did not track the fly", flush=True)
        return 1
    if not hops.get("programs"):
        print("  hops did not boot", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
