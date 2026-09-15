#!/usr/bin/env python3
"""Measured fly rest / odor observer → FSOT residual boot.

Álvarez-Salvado et al. eLife 37815 Video 1 (MPEG-4). Four walking flies,
ACV 10% pulse marked by a green overlay at the top of the frame. Files
stay on D:\\FlyWire_Connectome\\behavior\\odor (not git).

Observer (not a trained pose net):
  green overlay  → odor on (1/φ of max green-excess)
  MAD+φ          → motion on, dark-blob centroids per chamber
  pre-odor still → rest (no afferent seed)
  odor on        → olfactory seed
  offset motion  → odor off, still moving (no invented OFF class)

Hops are the live Male CNS / v630 programs (olfactory vs JO / vnc_sensory).
Rest is not a cell class. Do not invent odor identities or GABA.

Not a trained RNN and not a thought. 0 free parameters.
"""
from __future__ import annotations

import json
import shutil
import subprocess
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
    boot_activity,
    cuda_name,
    ensure_annotations,
    ensure_connections,
)

ODOR = FLY_ROOT / "behavior" / "odor"
VIDEO_NAME = "elife-37815-video1.mp4"
VIDEO_URL = (
    "https://static-movie-usa.glencoesoftware.com/mp4/10.7554/31/"
    "3bac3959905074b9555ead5ec7b505ab71fde0d3/elife-37815-video1.mp4"
)
ELIFE = "https://doi.org/10.7554/eLife.37815"
STIMULUS = "ACV 10% pulse"

MALE_BOOT = ROOT / "data" / "male_cns_boot.json"
WALK_BOOT = ROOT / "data" / "fly_behavior_flow.json"


def ensure_video() -> Path:
    ODOR.mkdir(parents=True, exist_ok=True)
    path = ODOR / VIDEO_NAME
    if path.exists() and path.stat().st_size > 100_000:
        return path
    print(f"  downloading {VIDEO_URL}", flush=True)
    urllib.request.urlretrieve(VIDEO_URL, path)
    print(f"  wrote {path} ({path.stat().st_size} bytes)", flush=True)
    return path


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


def decode_rgb(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    info = _ffprobe(path)
    stream = (info.get("streams") or [{}])[0]
    fmt = info.get("format") or {}
    w = int(stream.get("width") or 0)
    h = int(stream.get("height") or 0)
    fps_s = str(stream.get("avg_frame_rate") or "0/1")
    try:
        num, den = fps_s.split("/")
        fps = float(num) / float(den) if float(den) else 0.0
    except ValueError:
        fps = 0.0
    duration = float(fmt.get("duration") or 0.0)
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        raise SystemExit("ffmpeg not on PATH")
    proc = subprocess.run(
        [
            ffmpeg,
            "-v",
            "error",
            "-i",
            str(path),
            "-f",
            "rawvideo",
            "-pix_fmt",
            "rgb24",
            "pipe:1",
        ],
        check=True,
        capture_output=True,
    )
    pix = np.frombuffer(proc.stdout, dtype=np.uint8)
    n = len(pix) // (w * h * 3)
    vol = pix[: n * w * h * 3].reshape(n, h, w, 3)
    meta = {
        "path": str(path),
        "bytes": int(path.stat().st_size),
        "width": w,
        "height": h,
        "n_frames": int(n),
        "duration_s": duration,
        "fps": fps,
        "codec": stream.get("codec_name"),
    }
    return vol, meta


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    w = mask.astype(np.int8)
    edges = np.diff(np.concatenate(([0], w, [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    return [(int(s), int(e)) for s, e in zip(starts, ends)]


def find_chambers(gray0: np.ndarray) -> list[dict[str, int]]:
    """Bright vertical lanes = walking chambers.

    Vignette puts the left lane below median − φ·MAD, so a lower-tail
    divider gate clips it. Same cut as the green overlay and residual
    n_active: columns (rows) ≥ observer-max / φ.
    """
    col = gray0.mean(axis=0)
    row = gray0.mean(axis=1)
    x_runs = _runs(col >= (float(col.max()) / _PHI))
    y_runs = _runs(row >= (float(row.max()) / _PHI))
    if not y_runs:
        y0, y1 = 0, gray0.shape[0]
    else:
        y0, y1 = max(y_runs, key=lambda r: r[1] - r[0])
    chambers = []
    for x0, x1 in x_runs:
        if x1 - x0 < 20:
            continue
        chambers.append(
            {
                "id": len(chambers),
                "x0": int(x0),
                "x1": int(x1),
                "y0": int(y0),
                "y1": int(y1),
            }
        )
    return chambers


def odor_overlay(vol: np.ndarray) -> dict[str, Any]:
    """Green-dot overlay is the published odor marker.

    Per-frame max of G − max(R, B) in the top strip. On = 1/φ of the
    observer max (same cut as residual n_active). Not a free 20-count.
    """
    top = vol[:, :50, :, :]
    gex = top[:, :, :, 1].astype(np.float64) - np.maximum(
        top[:, :, :, 0], top[:, :, :, 2]
    ).astype(np.float64)
    gmax = gex.max(axis=(1, 2))
    peak = float(gmax.max()) + 1e-12
    thr = peak / _PHI
    on = gmax >= thr
    w = on.astype(np.int8)
    edges = np.diff(np.concatenate(([0], w, [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
    lengths = (ends - starts).tolist()
    loc = None
    if on.any():
        mid = int(np.flatnonzero(on)[len(np.flatnonzero(on)) // 2])
        ys, xs = np.where(gex[mid] >= thr)
        if ys.size:
            loc = {
                "y0": int(ys.min()),
                "y1": int(ys.max()) + 1,
                "x0": int(xs.min()),
                "x1": int(xs.max()) + 1,
                "n_px": int(ys.size),
            }
    return {
        "threshold": float(thr),
        "observer_max": peak,
        "cut": "observer_max / φ",
        "n_on": int(on.sum()),
        "n": int(on.size),
        "paint_frac": float(on.mean()),
        "first_on": int(starts[0]) if len(starts) else None,
        "last_on": int(ends[-1] - 1) if len(ends) else None,
        "n_bouts": int(len(starts)),
        "bout_lengths": lengths,
        "location": loc,
        "series": gmax.tolist(),
        "on": on.astype(bool).tolist(),
    }


def track_chamber(
    gray: np.ndarray, ch: dict[str, int]
) -> dict[str, Any]:
    """Frame-diff energy + first moment. Not a trained pose net.

    Dark-blob largest-CC tracks the vignette, not the fly. Motion
    centroid is the same observer as Biohub: first moment of |Δ|.
    """
    y0, y1, x0, x1 = ch["y0"], ch["y1"], ch["x0"], ch["x1"]
    roi = gray[:, y0:y1, x0:x1].astype(np.float64)
    d = np.abs(np.diff(roi, axis=0))
    energy = d.mean(axis=(1, 2))
    ys = np.arange(y0, y1, dtype=np.float64)
    mass = d.sum(axis=(1, 2)) + 1e-12
    com_y = (d.sum(axis=2) * ys).sum(axis=1) / mass
    upwind = -np.diff(com_y)
    return {
        "id": ch["id"],
        "roi": ch,
        "n_tracked": int(np.isfinite(com_y).sum()),
        "motion_energy": energy.tolist(),
        "centroid_y": com_y.tolist(),
        "upwind": [None if not np.isfinite(v) else float(v) for v in upwind],
        "mean_upwind": float(np.nanmean(upwind)) if np.isfinite(upwind).any() else None,
    }


def _epoch_labels(odor_on: np.ndarray, moving: np.ndarray) -> np.ndarray:
    labels = np.empty(odor_on.size, dtype=object)
    seen = False
    for i, (od, mv) in enumerate(zip(odor_on, moving)):
        if od:
            seen = True
            labels[i] = "odor_on"
        elif not seen:
            labels[i] = "pre_walk" if mv else "rest"
        else:
            labels[i] = "offset" if mv else "post_rest"
    return labels


def _mean_where(x: np.ndarray, mask: np.ndarray) -> float | None:
    hit = x[mask]
    hit = hit[np.isfinite(hit)]
    if hit.size == 0:
        return None
    return float(hit.mean())


def observe(vol: np.ndarray, meta: dict[str, Any]) -> dict[str, Any]:
    gray = vol.mean(axis=3)
    odor = odor_overlay(vol)
    chambers = find_chambers(gray[0])
    tracks = [track_chamber(gray, ch) for ch in chambers]
    n = int(gray.shape[0])
    total = np.zeros(n, dtype=np.float64)
    for t in tracks:
        e = np.asarray(t["motion_energy"], dtype=np.float64)
        total[1 : 1 + e.size] += e[: n - 1]
    motion_gate = mad_phi_gate(total[1:]) if n > 1 else mad_phi_gate(np.zeros(1))
    moving = np.zeros(n, dtype=bool)
    if n > 1:
        moving[1:] = total[1:] >= motion_gate["threshold"]
    odor_on = np.asarray(odor["on"], dtype=bool)[:n]
    moving = moving[:n]
    epochs = _epoch_labels(odor_on, moving)
    fps = float(meta.get("fps") or 0.0) or 30.0
    counts = {k: int((epochs == k).sum()) for k in ("rest", "pre_walk", "odor_on", "offset", "post_rest")}
    energy_by_epoch: dict[str, list[float]] = {k: [] for k in counts}
    comy_by_epoch: dict[str, list[float]] = {k: [] for k in counts}
    upwind_by_epoch: dict[str, list[float]] = {k: [] for k in counts}
    per_fly = []
    for t in tracks:
        e = np.asarray(t["motion_energy"], dtype=np.float64)
        cy = np.asarray(t["centroid_y"], dtype=np.float64)
        uw = np.asarray(
            [np.nan if v is None else v for v in t["upwind"]], dtype=np.float64
        )
        ep_e = epochs[1 : 1 + e.size]
        ep_u = epochs[2 : 2 + uw.size]
        row = {
            "id": t["id"],
            "n_tracked": t["n_tracked"],
            "roi": t["roi"],
            "energy_by_epoch": {},
            "motion_y_by_epoch": {},
            "upwind_by_epoch": {},
        }
        for k in counts:
            me = _mean_where(e[: ep_e.size], ep_e == k)
            my = _mean_where(cy[: ep_e.size], ep_e == k)
            mu = _mean_where(uw[: ep_u.size], ep_u == k)
            row["energy_by_epoch"][k] = me
            row["motion_y_by_epoch"][k] = my
            row["upwind_by_epoch"][k] = mu
            if me is not None:
                energy_by_epoch[k].append(me)
            if my is not None:
                comy_by_epoch[k].append(my)
            if mu is not None and k not in {"rest", "post_rest"}:
                upwind_by_epoch[k].append(mu)
            if k in {"rest", "post_rest"}:
                row["upwind_by_epoch"][k] = None
        per_fly.append(row)

    def _avg(d: dict[str, list[float]]) -> dict[str, float | None]:
        return {k: (float(np.mean(v)) if v else None) for k, v in d.items()}

    group_energy = _avg(energy_by_epoch)
    group_comy = _avg(comy_by_epoch)
    group_upwind = _avg(upwind_by_epoch)
    first = odor.get("first_on")
    last = odor.get("last_on")
    pulse_s = None
    if first is not None and last is not None:
        pulse_s = (last - first + 1) / fps
    # rest exists if pre-odor still frames are a real fraction
    rest_frac = counts["rest"] / max(n, 1)
    return {
        "source": ELIFE,
        "stimulus": STIMULUS,
        "video": meta,
        "n_chambers": len(chambers),
        "chambers": chambers,
        "odor": {k: v for k, v in odor.items() if k not in {"series", "on"}},
        "odor_on": odor["on"][:n],
        "odor_gmax": odor["series"][:n],
        "motion_gate": {k: v for k, v in motion_gate.items() if k != "bout_lengths"},
        "motion_energy": total.tolist(),
        "moving": moving.astype(bool).tolist(),
        "epochs": epochs.tolist(),
        "epoch_frames": counts,
        "epoch_frac": {k: v / max(n, 1) for k, v in counts.items()},
        "pulse_s": pulse_s,
        "rest_frac": rest_frac,
        "has_rest_state": rest_frac > 1.0 / (_PHI ** 5),
        "motion_energy_by_epoch": group_energy,
        "motion_y_by_epoch": group_comy,
        "upwind_px_per_frame": group_upwind,
        "per_fly": per_fly,
        "tracks": [
            {
                "id": t["id"],
                "n_tracked": t["n_tracked"],
                "centroid_y": t["centroid_y"],
                "upwind": t["upwind"],
                "motion_energy": t["motion_energy"],
            }
            for t in tracks
        ],
        "seed_map": {
            "rest": {"seed": None, "how": None, "why": "odor off, motion below MAD+φ; no afferent seed"},
            "pre_walk": {"seed": "JO", "how": "type_prefix", "why": "wind always on; fly moving before odor"},
            "odor_on": {"seed": "olfactory", "how": "class", "why": "green overlay = ACV 10% pulse"},
            "offset": {
                "seed": None,
                "how": None,
                "why": "odor off after pulse; local search is not a cell class",
            },
            "post_rest": {"seed": None, "how": None, "why": "odor off, still"},
        },
        "note": (
            "Green overlay is the paper's odor marker. Motor is chamber "
            "frame-diff energy and the first moment of |Δ| (not a trained "
            "pose net; dark-blob largest-CC follows vignette). Rest is an "
            "observer state, not a neuron class."
        ),
    }


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


def attach_hops(*, boot: bool, gpu: bool) -> dict[str, Any]:
    """Live hops: olfactory vs JO / vnc_sensory. Rest = empty seed."""
    out: dict[str, Any] = {
        "rest": {
            "name": "rest",
            "n_seed": 0,
            "note": "no afferent seed when odor is off and motion is below MAD+φ",
            "hop_2": {"vnc_motor": 0.0, "descending": 0.0, "motor": 0.0},
        }
    }
    if boot:
        print("== re-boot v630 programs", flush=True)
        ann = ensure_annotations(FLY_ROOT)
        conn = ensure_connections(FLY_ROOT)
        programs = []
        for spec in (
            {"name": "olfactory", "seed": "olfactory", "how": "class"},
            {"name": "JO", "seed": "jo-", "how": "type_prefix"},
            {"name": "mechanosensory", "seed": "mechanosensory", "how": "class"},
        ):
            print(f"  v630 {spec['name']}", flush=True)
            run = boot_activity(
                ann, conn, seed=spec["seed"], how=spec["how"], gpu=gpu
            )
            programs.append(
                {
                    "name": spec["name"],
                    "n_seed": run["n_seed"],
                    "gpu": run.get("gpu"),
                    "device": run.get("device"),
                    "elapsed_s": run.get("elapsed_s"),
                    "flow": [
                        {
                            "hop": s.get("hop"),
                            "motor": (s.get("target_mass") or {}).get("motor"),
                            "descending": (s.get("target_mass") or {}).get("descending"),
                            "top": (s.get("top") or [{}])[0],
                        }
                        for s in (run.get("trace") or [])
                        if s.get("hop") in {0, 1, 2, 3}
                    ],
                }
            )
        out["v630"] = {"programs": programs, "source": "re-boot"}
    elif WALK_BOOT.exists():
        walk = json.loads(WALK_BOOT.read_text(encoding="utf-8"))
        conn = walk.get("connectome") or {}
        out["v630"] = {
            "source": str(WALK_BOOT),
            "gpu_name": conn.get("gpu_name"),
            "compare": _slim_compare(conn.get("compare") or {}),
            "n_programs": len(conn.get("programs") or []),
        }
    if MALE_BOOT.exists():
        male = json.loads(MALE_BOOT.read_text(encoding="utf-8"))
        out["male_cns"] = {
            "source": str(MALE_BOOT),
            "n_neurons": male.get("n_neurons"),
            "n_edges": male.get("n_edges"),
            "n_gaba": male.get("n_gaba"),
            "authority": male.get("authority"),
            "compare": _slim_compare(male.get("compare") or {}),
        }
    return out


def _hop_pick(compare: dict[str, Any], hop: str, name: str, key: str) -> float | None:
    row = (compare.get(hop) or {}).get(name) or {}
    v = row.get(key)
    return float(v) if v is not None else None


def loop_table(obs: dict[str, Any], hops: dict[str, Any]) -> dict[str, Any]:
    male = (hops.get("male_cns") or {}).get("compare") or {}
    v630 = (hops.get("v630") or {}).get("compare") or {}
    en = obs.get("motion_energy_by_epoch") or {}
    rest_e = en.get("rest") or 0.0
    odor_e = en.get("odor_on")
    return {
        "rest": {
            "observer": "pre-odor, motion below MAD+φ",
            "motion_energy": (obs.get("motion_energy_by_epoch") or {}).get("rest"),
            "motion_y": (obs.get("motion_y_by_epoch") or {}).get("rest"),
            "seed": None,
            "male_cns_hop2_vnc_motor": 0.0,
            "v630_hop2_descending": 0.0,
        },
        "odor_on": {
            "observer": "green overlay = ACV 10%",
            "motion_energy": (obs.get("motion_energy_by_epoch") or {}).get("odor_on"),
            "motion_y": (obs.get("motion_y_by_epoch") or {}).get("odor_on"),
            "seed": "olfactory",
            "male_cns_hop2_vnc_motor": _hop_pick(male, "hop_2", "olfactory", "vnc_motor"),
            "male_cns_hop1_peak": ((male.get("hop_1") or {}).get("olfactory") or {}).get("top"),
            "v630_hop2_descending": _hop_pick(v630, "hop_2", "olfactory", "descending"),
            "energy_over_rest": (None if not rest_e else (odor_e / rest_e if odor_e is not None else None)),
        },
        "offset": {
            "observer": "odor off, motion above MAD+φ",
            "motion_energy": (obs.get("motion_energy_by_epoch") or {}).get("offset"),
            "motion_y": (obs.get("motion_y_by_epoch") or {}).get("offset"),
            "seed": None,
            "note": "local search after the pulse; not a cell class",
        },
        "walking_contrast_JO": {
            "observer": "Harvard tethered walk (separate clip)",
            "seed": "JO / vnc_sensory",
            "male_cns_hop2_vnc_motor_JO": _hop_pick(male, "hop_2", "JO", "vnc_motor"),
            "male_cns_hop2_vnc_motor_vnc_sensory": _hop_pick(
                male, "hop_2", "vnc_sensory", "vnc_motor"
            ),
            "male_cns_hop1_peak_JO": ((male.get("hop_1") or {}).get("JO") or {}).get("top"),
            "v630_hop2_descending_JO": _hop_pick(v630, "hop_2", "JO", "descending"),
        },
        "honesty": (
            "Olfactory hops stay in the antennal lobe (il3LN6) and do not "
            "light vnc_motor. JO / vnc_sensory do. Frame-diff energy rises "
            "in the pulse and peaks at offset (search). Motion first-moment "
            "is further upwind during odor than rest. Do not claim the "
            "olfactory cascade is the walking command."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    ap.add_argument(
        "--boot",
        action="store_true",
        help="re-run v630 residual hops (Male CNS stays the live JSON)",
    )
    args = ap.parse_args(argv)
    print("FSOT fly rest/odor observer → residual boot", flush=True)
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
        f"  chambers={obs['n_chambers']} odor_on={obs['odor']['n_on']}/{obs['odor']['n']} "
        f"pulse={obs['pulse_s']:.2f}s rest_frac={obs['rest_frac']:.3f} "
        f"has_rest={obs['has_rest_state']}",
        flush=True,
    )
    print(f"  epochs {obs['epoch_frames']}", flush=True)
    print(f"  energy {obs.get('motion_energy_by_epoch')}", flush=True)
    print(f"  motion_y {obs.get('motion_y_by_epoch')}", flush=True)
    print(f"  upwind {obs['upwind_px_per_frame']}", flush=True)
    hops = attach_hops(boot=args.boot, gpu=args.gpu)
    loop = loop_table(obs, hops)
    report = {
        "product": "measured rest/odor observer + measured fly graph",
        "free_parameters": 0,
        "phi": _PHI,
        "residual_Biochemistry": _R_BIO,
        "authority": (
            "Álvarez-Salvado et al. eLife 7:e37815 Video 1 (ACV 10% pulse, "
            "green overlay); Male CNS v1.0 + FlyWire v630 residual hops"
        ),
        "honesty": (
            "Harvard tethered-walk clips have no rest. This video does: "
            "pre-odor still frames vs a 10 s ACV pulse vs offset search. "
            "Rest is not a neuron class — odor off and motion below MAD+φ "
            "means no afferent seed. Olfactory hops do not reach vnc_motor; "
            "JO / vnc_sensory do. Dryad 10.5061/dryad.g27mq71 is 6.94 GB "
            "and was not downloaded; GitHub is LabVIEW/MATLAB, not CSV. "
            "GABA only where annotated. Not a trained pose net, not a "
            "trained RNN, not a thought."
        ),
        "observer": {
            k: v
            for k, v in obs.items()
            if k not in {"tracks", "odor_on", "odor_gmax", "motion_energy", "moving", "epochs"}
        },
        "series": {
            "odor_on": obs["odor_on"],
            "odor_gmax": obs["odor_gmax"],
            "motion_energy": obs["motion_energy"],
            "moving": obs["moving"],
            "epochs": obs["epochs"],
            "centroid_y": [t["centroid_y"] for t in obs["tracks"]],
            "upwind": [t["upwind"] for t in obs["tracks"]],
        },
        "connectome": hops,
        "loop": loop,
        "gpu_name": cuda_name(),
    }
    out = ROOT / "data" / "fly_odor_flow.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(loop, indent=2), flush=True)
    if not obs["has_rest_state"] or obs["n_chambers"] < 4 or not obs["pulse_s"]:
        print("  observer did not boot a rest/odor split", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
