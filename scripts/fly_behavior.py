#!/usr/bin/env python3
"""Measured fly walking observer → FSOT residual boot.

Harvard Dataverse doi:10.7910/DVN/BBNPYX tethered walking 3D keypoints
are the observer (not a trained pose net). MAD+φ gates high-energy
stepping. Walking-on seeds measured mechanosensory / JO class on the
v630 brain graph. Residual hops (GPU sparse matvec if CUDA) report
descending / DN mass — the brain→cord walking command.

v630 is brain-only. The 100 motor neurons are not VNC leg MNs. This
clip is a walking trial (tarsus energy never near rest). Gate splits
high vs low step energy; it does not invent a rest state.

Not a trained RNN and not a thought. 0 free parameters.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

from fly_connectome import (  # noqa: E402
    FLY_ROOT,
    _PHI,
    _R_BIO,
    boot_activity,
    cuda_name,
    ensure_annotations,
    ensure_connections,
)

BEHAVIOR = FLY_ROOT / "behavior"
LEGS = ("R1", "R2", "R3", "L1", "L2", "L3")
DATAVERSE = "doi:10.7910/DVN/BBNPYX"

TRIALS = [
    {
        "id": "Fly01_T001",
        "body": "Fly01_T001_BodyCoords3D.csv",
        "video": "Fly01_T001cam-0.mp4",
    },
    {
        "id": "Fly01_T002",
        "body": "Fly01_T002_BodyCoords3D.csv",
        "video": "Fly01_T002cam-0.mp4",
    },
    {
        "id": "Fly02_T002",
        "body": "Fly02_T002_BodyCoords3D.csv",
        "video": "Fly02_T002cam-0.mp4",
    },
]

# Walking brain afferents vs contrast programs. Exact class / type prefix
# so photoreceptors do not swamp the seed (first sensory boot lesson).
PROGRAMS = [
    {"name": "mechanosensory", "seed": "mechanosensory", "how": "class"},
    {"name": "JO", "seed": "jo-", "how": "type_prefix"},
    {"name": "olfactory", "seed": "olfactory", "how": "class"},
    {"name": "sensory", "seed": "sensory", "how": "class"},
]


def mad_phi_gate(x: np.ndarray) -> dict[str, Any]:
    """Median + φ·MAD. Tighten to φ² if paint > 1/φ. Never loosen."""
    x = np.asarray(x, dtype=np.float64)
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med))) + 1e-12
    thr = med + _PHI * mad
    on = x >= thr
    paint = float(on.mean())
    tightened = False
    if paint > 1.0 / _PHI:
        thr = med + (_PHI * _PHI) * mad
        on = x >= thr
        paint = float(on.mean())
        tightened = True
    w = on.astype(np.int8)
    edges = np.diff(np.concatenate(([0], w, [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    lengths = (ends - starts).tolist()
    return {
        "median": med,
        "mad": mad,
        "threshold": thr,
        "phi": _PHI,
        "tightened_to_phi2": tightened,
        "n_on": int(on.sum()),
        "n": int(on.size),
        "paint_frac": paint,
        "n_bouts": int(len(starts)),
        "bout_lengths": lengths,
        "max_bout": int(max(lengths) if lengths else 0),
        "mean_bout": float(np.mean(lengths) if lengths else 0.0),
    }


def _xyz(df, prefix: str) -> np.ndarray:
    return np.stack(
        [
            df[f"{prefix}_x"].to_numpy(dtype=np.float64),
            df[f"{prefix}_y"].to_numpy(dtype=np.float64),
            df[f"{prefix}_z"].to_numpy(dtype=np.float64),
        ],
        axis=1,
    )


def kinematics(csv_path: Path) -> dict[str, Any]:
    import pandas as pd

    df = pd.read_csv(csv_path)
    tarsus = np.stack([_xyz(df, f"{leg}-Tar") for leg in LEGS], axis=1)
    step = np.linalg.norm(np.diff(tarsus, axis=0), axis=2)
    energy = step.sum(axis=1)
    per_leg = {leg: step[:, i].tolist() for i, leg in enumerate(LEGS)}
    left = step[:, 3:6].sum(axis=1)
    right = step[:, 0:3].sum(axis=1)
    ant_l = np.linalg.norm(np.diff(_xyz(df, "H-lAnt"), axis=0), axis=1)
    ant_r = np.linalg.norm(np.diff(_xyz(df, "H-rAnt"), axis=0), axis=1)
    antenna = ant_l + ant_r
    gate_t = mad_phi_gate(energy)
    gate_a = mad_phi_gate(antenna)
    return {
        "path": str(csv_path),
        "n_frames": int(len(df)),
        "n_steps": int(energy.size),
        "tarsus_energy": {
            "min": float(energy.min()),
            "median": float(np.median(energy)),
            "mean": float(energy.mean()),
            "max": float(energy.max()),
            "p05": float(np.percentile(energy, 5)),
            "p95": float(np.percentile(energy, 95)),
        },
        "per_leg_mean_speed": {
            leg: float(step[:, i].mean()) for i, leg in enumerate(LEGS)
        },
        "left_right_mean": {
            "left": float(left.mean()),
            "right": float(right.mean()),
        },
        "antenna_energy": {
            "min": float(antenna.min()),
            "median": float(np.median(antenna)),
            "mean": float(antenna.mean()),
            "max": float(antenna.max()),
        },
        "tarsus_gate": gate_t,
        "antenna_gate": gate_a,
        "energy_floor_over_median": float(energy.min()) / max(float(np.median(energy)), 1e-12),
        "rest_floor": float(np.median(energy)) / (_PHI ** 5),
        "whole_clip_is_walking": float(energy.min())
        > float(np.median(energy)) / (_PHI ** 5),
        "series": {
            "tarsus_energy": energy.tolist(),
            "antenna_energy": antenna.tolist(),
            "left": left.tolist(),
            "right": right.tolist(),
        },
        "per_leg_series": per_leg,
    }


def _ffprobe(path: Path) -> dict[str, Any]:
    raw = subprocess.check_output(
        [
            "ffprobe",
            "-v",
            "error",
            "-show_entries",
            "format=duration:stream=width,height,nb_frames,avg_frame_rate,codec_name",
            "-of",
            "json",
            str(path),
        ],
        text=True,
    )
    return json.loads(raw)


def video_observer(mp4: Path) -> dict[str, Any]:
    """Frame-diff motion energy via ffmpeg. Corroborates 3D tarsus gate."""
    info = _ffprobe(mp4)
    stream = (info.get("streams") or [{}])[0]
    fmt = info.get("format") or {}
    w = int(stream.get("width") or 0)
    h = int(stream.get("height") or 0)
    n_meta = int(float(stream.get("nb_frames") or 0))
    duration = float(fmt.get("duration") or 0.0)
    fps_s = str(stream.get("avg_frame_rate") or "0/1")
    try:
        num, den = fps_s.split("/")
        fps = float(num) / float(den) if float(den) else 0.0
    except ValueError:
        fps = 0.0
    ffmpeg = shutil.which("ffmpeg")
    out: dict[str, Any] = {
        "path": str(mp4),
        "bytes": int(mp4.stat().st_size),
        "width": w,
        "height": h,
        "n_frames_meta": n_meta,
        "duration_s": duration,
        "fps": fps,
        "codec": stream.get("codec_name"),
    }
    if not ffmpeg:
        out["frame_diff"] = None
        out["note"] = "ffmpeg not on PATH; 3D keypoints remain the observer"
        return out
    tw, th = 448, 270
    proc = subprocess.run(
        [
            ffmpeg,
            "-v",
            "error",
            "-i",
            str(mp4),
            "-vf",
            f"format=gray,scale={tw}:{th}",
            "-f",
            "rawvideo",
            "-pix_fmt",
            "gray",
            "pipe:1",
        ],
        check=True,
        capture_output=True,
    )
    frames = np.frombuffer(proc.stdout, dtype=np.uint8)
    n = len(frames) // (tw * th)
    if n < 2:
        out["frame_diff"] = None
        return out
    vol = frames[: n * tw * th].reshape(n, th, tw).astype(np.float32)
    energy = np.mean(np.abs(np.diff(vol, axis=0)), axis=(1, 2))
    gate = mad_phi_gate(energy)
    out["n_decoded"] = int(n)
    out["frame_diff"] = {
        "min": float(energy.min()),
        "median": float(np.median(energy)),
        "mean": float(energy.mean()),
        "max": float(energy.max()),
        "gate": gate,
        "series": energy.tolist(),
    }
    return out


def _trace_targets(run: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for snap in run.get("trace") or []:
        tm = snap.get("target_mass") or {}
        cm = snap.get("class_mass") or {}
        rows.append(
            {
                "hop": snap.get("hop"),
                "l1": snap.get("l1"),
                "n_active": snap.get("n_active"),
                "motor": tm.get("motor"),
                "descending": tm.get("descending"),
                "sensory": tm.get("sensory"),
                "DN": cm.get("DN"),
                "mechanosensory": cm.get("mechanosensory"),
                "olfactory": cm.get("olfactory"),
                "top": (snap.get("top") or [])[:4],
            }
        )
    return rows


def run_programs(*, gpu: bool = True) -> dict[str, Any]:
    ann = ensure_annotations(FLY_ROOT)
    conn = ensure_connections(FLY_ROOT)
    programs = []
    for spec in PROGRAMS:
        print(f"== program {spec['name']} seed={spec['seed']} how={spec['how']}", flush=True)
        run = boot_activity(
            ann,
            conn,
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
                "flow": _trace_targets(run),
                "authority": run.get("authority"),
            }
        )
    # Hop-1/2 compare: which program puts more mass on descending / DN.
    compare: dict[str, Any] = {}
    for hop in (1, 2, 3):
        row = {}
        for p in programs:
            hit = next((f for f in p["flow"] if f.get("hop") == hop), None)
            if hit:
                row[p["name"]] = {
                    "descending": hit.get("descending"),
                    "DN": hit.get("DN"),
                    "motor": hit.get("motor"),
                    "n_active": hit.get("n_active"),
                    "l1": hit.get("l1"),
                }
        compare[f"hop_{hop}"] = row
    return {
        "programs": programs,
        "compare": compare,
        "gpu_name": cuda_name(),
        "residual_Biochemistry": _R_BIO,
        "phi": _PHI,
    }


def analyze_trials() -> list[dict[str, Any]]:
    rows = []
    for t in TRIALS:
        body = BEHAVIOR / t["body"]
        video = BEHAVIOR / t["video"]
        if not body.exists() or body.stat().st_size < 1000:
            print(f"  skip {t['id']}: missing {body}", flush=True)
            continue
        print(f"== observer {t['id']}", flush=True)
        kin = kinematics(body)
        # drop per-frame series from the printed path; keep in JSON
        vid = None
        if video.exists() and video.stat().st_size > 1000:
            print(f"  video {video.name}", flush=True)
            vid = video_observer(video)
        # Agreement: tarsus gate vs video frame-diff gate, if same length.
        agree = None
        if vid and vid.get("frame_diff"):
            te = np.asarray(kin["series"]["tarsus_energy"], dtype=np.float64)
            ve = np.asarray(vid["frame_diff"]["series"], dtype=np.float64)
            n = min(te.size, ve.size)
            if n > 8:
                t_on = te[:n] >= kin["tarsus_gate"]["threshold"]
                v_on = ve[:n] >= vid["frame_diff"]["gate"]["threshold"]
                agree = {
                    "n": int(n),
                    "both_on": int((t_on & v_on).sum()),
                    "tarsus_only": int((t_on & ~v_on).sum()),
                    "video_only": int((~t_on & v_on).sum()),
                    "jaccard": float(
                        (t_on & v_on).sum() / max(1, (t_on | v_on).sum())
                    ),
                }
        rows.append(
            {
                "id": t["id"],
                "source": DATAVERSE,
                "kinematics": {
                    k: v
                    for k, v in kin.items()
                    if k not in {"series", "per_leg_series"}
                },
                "series": kin["series"],
                "video": {
                    k: v
                    for k, v in (vid or {}).items()
                    if k != "frame_diff"
                }
                | (
                    {
                        "frame_diff": {
                            k: v
                            for k, v in (vid or {}).get("frame_diff", {}).items()
                            if k != "series"
                        },
                        "frame_diff_series": (vid or {})
                        .get("frame_diff", {})
                        .get("series"),
                    }
                    if vid and vid.get("frame_diff")
                    else {}
                ),
                "tarsus_vs_video_gate": agree,
            }
        )
        tg = kin["tarsus_gate"]
        print(
            f"  frames={kin['n_frames']} tarsus_on={tg['n_on']}/{tg['n']} "
            f"paint={tg['paint_frac']:.3f} bouts={tg['n_bouts']} "
            f"walking_clip={kin['whole_clip_is_walking']}",
            flush=True,
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    ap.add_argument("--skip-video", action="store_true")
    args = ap.parse_args(argv)
    BEHAVIOR.mkdir(parents=True, exist_ok=True)
    print("FSOT fly behavior observer → residual boot", flush=True)
    print(f"  phi={_PHI} R_BIO={_R_BIO} gpu={cuda_name()}", flush=True)
    trials = analyze_trials() if not args.skip_video else []
    flow = run_programs(gpu=args.gpu)
    report = {
        "product": "measured walking observer + measured FlyWire graph",
        "free_parameters": 0,
        "phi": _PHI,
        "residual_Biochemistry": _R_BIO,
        "authority": (
            "Harvard Dataverse 10.7910/DVN/BBNPYX 3D keypoints + cam-0 "
            "video; FlyWire v630 connections when v783 Zenodo unavailable"
        ),
        "honesty": (
            "v630 is the brain, not the VNC. Leg campaniform / chordotonal "
            "and leg motor neurons are not in this dump. Walking prediction "
            "is descending / DN mass after a mechanosensory or JO seed. "
            "These clips are tethered walking (energy never near rest); "
            "MAD+φ splits high vs low step energy. GABA inhibitory. "
            "Not a trained pose net, not a trained RNN, not a thought."
        ),
        "observer": trials,
        "connectome": flow,
        "genetics_join": (
            "Named cell types on this graph (JO, DN, olfactory LN) join "
            "later to FlyBase / UniProt → FSOT product Cα on the same pin. "
            "Do not invent contacts."
        ),
    }
    out = ROOT / "data" / "fly_behavior_flow.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(flow.get("compare"), indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
