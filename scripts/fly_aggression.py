#!/usr/bin/env python3
"""Measured fly aggression observer → FSOT residual boot.

Gao et al. eLife 13:RP104212. Files on
D:\\FlyWire_Connectome\\behavior\\aggression (not git).

  fig4 video 1 — pC1SS2>CsChrimson, red-light overlay is the published
                 stimulus marker (like the odor green dot). Four wells.
  fig1 video 1 — wild-type Canton-S males tussling (motor corroboration).

Observer (not a trained pose net, not JAABA):
  title = near-black frames
  red-excess ≥ observer-max / φ in the center overlay = light ON
  MAD+φ on chamber frame-diff = tussling motor
  leftover dark CCs (drop ≥ max/φ) = the two males; distance = proximity

Seed is measured Male CNS genetics, not invented aggression cells:
  pC1_ types (paper: pC1SS2 promotes tussling), dsx, male-specific dimorphism.

pC1 also seeds courtship — same cells, different observer. 0 free parameters.
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

AGGR = FLY_ROOT / "behavior" / "aggression"
ELIFE = "https://doi.org/10.7554/eLife.104212"
GLENCOE = (
    "https://static-movie-usa.glencoesoftware.com/mp4/10.7554/81/"
    "4a24ece707af0d83caf6ee6b4a67fb9febf19a36/"
)
VIDEOS = [
    {
        "id": "pC1_optogenetic",
        "name": "elife-104212-fig4-video1.mp4",
        "what": "pC1SS2>UAS-CsChrimson; red-light overlay; 4 wells; 3x",
    },
    {
        "id": "wt_tussling",
        "name": "elife-104212-fig1-video1.mp4",
        "what": "Canton-S G14 males tussling; 2x",
    },
]

MALE_BOOT = ROOT / "data" / "male_cns_boot.json"
COURT_BOOT = ROOT / "data" / "fly_courtship_flow.json"

PROGRAMS = [
    {"name": "pC1", "seed": "pc1_", "how": "type_prefix"},
    {"name": "dsx", "seed": "dsx", "how": "fru_dsx"},
    {"name": "male_specific", "seed": "male-specific", "how": "dimorphism"},
]


def ensure_videos() -> list[Path]:
    AGGR.mkdir(parents=True, exist_ok=True)
    out = []
    for spec in VIDEOS:
        path = AGGR / spec["name"]
        if not (path.exists() and path.stat().st_size > 100_000):
            url = GLENCOE + spec["name"]
            print(f"  downloading {url}", flush=True)
            urllib.request.urlretrieve(url, path)
            print(f"  wrote {path} ({path.stat().st_size} bytes)", flush=True)
        out.append(path)
    return out


def _live_mask(vol: np.ndarray) -> np.ndarray:
    """Drop title cards: near-black frames."""
    g = vol.mean(axis=(1, 2, 3))
    peak = float(g.max()) + 1e-12
    return g >= (peak / _PHI)


def red_overlay(vol: np.ndarray) -> dict[str, Any]:
    """Published 'Red Light ON' overlay is red text in the center."""
    h, w = vol.shape[1], vol.shape[2]
    y0, y1 = int(h * 0.40), int(h * 0.60)
    x0, x1 = int(w * 0.15), int(w * 0.85)
    cen = vol[:, y0:y1, x0:x1, :]
    rex = cen[:, :, :, 0].astype(np.float64) - np.maximum(
        cen[:, :, :, 1], cen[:, :, :, 2]
    ).astype(np.float64)
    rmax = rex.max(axis=(1, 2))
    peak = float(rmax.max()) + 1e-12
    thr = peak / _PHI
    on = rmax >= thr
    w8 = on.astype(np.int8)
    edges = np.diff(np.concatenate(([0], w8, [0])))
    starts = np.where(edges == 1)[0]
    ends = np.where(edges == -1)[0]
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
        "bout_lengths": (ends - starts).tolist(),
        "series": rmax.tolist(),
        "on": on.astype(bool).tolist(),
    }


def _quadrants(h: int, w: int) -> list[dict[str, int]]:
    my, mx = h // 2, w // 2
    return [
        {"id": 0, "y0": 0, "y1": my, "x0": 0, "x1": mx},
        {"id": 1, "y0": 0, "y1": my, "x0": mx, "x1": w},
        {"id": 2, "y0": my, "y1": h, "x0": 0, "x1": mx},
        {"id": 3, "y0": my, "y1": h, "x0": mx, "x1": w},
    ]


def _chamber_energy(gray: np.ndarray, ch: dict[str, int]) -> np.ndarray:
    roi = gray[:, ch["y0"] : ch["y1"], ch["x0"] : ch["x1"]]
    peak = float(roi.max()) + 1e-12
    arena = roi[0] >= (peak / _PHI)
    d = np.abs(np.diff(roi, axis=0))
    return d[:, arena].mean(axis=1) if arena.any() else d.mean(axis=(1, 2))


def _two_fly_dist(frame: np.ndarray, ch: dict[str, int]) -> float | None:
    """Leftover dark CCs in one well. Drop components ≥ max/φ (the rim)."""
    from scipy import ndimage

    roi = frame[ch["y0"] : ch["y1"], ch["x0"] : ch["x1"]]
    peak = float(roi.max()) + 1e-12
    arena = roi >= (peak / _PHI)
    if int(arena.sum()) < 50:
        return None
    med = float(np.median(roi[arena]))
    mad = float(np.median(np.abs(roi[arena] - med))) + 1e-12
    dark = np.zeros(roi.shape, dtype=bool)
    dark[arena] = roi[arena] <= (med - _PHI * mad)
    if float(dark[arena].mean()) > 1.0 / _PHI:
        dark[arena] = roi[arena] <= (med - (_PHI * _PHI) * mad)
    labeled, nlab = ndimage.label(dark)
    if nlab < 2:
        return 0.0 if nlab == 1 else None
    sizes = ndimage.sum(dark, labeled, index=np.arange(1, nlab + 1))
    mx = float(np.max(sizes))
    leftover = [i + 1 for i, s in enumerate(sizes) if 2.0 <= float(s) < mx / _PHI]
    if len(leftover) == 1:
        return 0.0
    if len(leftover) < 2:
        return None
    leftover = sorted(leftover, key=lambda i: -float(sizes[i - 1]))[:2]
    coms = [ndimage.center_of_mass(dark, labeled, int(i)) for i in leftover]
    return float(np.hypot(coms[0][0] - coms[1][0], coms[0][1] - coms[1][1]))


def _mean_where(x: np.ndarray, mask: np.ndarray) -> float | None:
    hit = x[mask[: x.size]]
    hit = hit[np.isfinite(hit)]
    if hit.size == 0:
        return None
    return float(hit.mean())


def observe_pc1(vol: np.ndarray, meta: dict[str, Any]) -> dict[str, Any]:
    gray = vol.mean(axis=3)
    n = int(gray.shape[0])
    live = _live_mask(vol)
    odor = red_overlay(vol)  # name kept; it is the red-light overlay
    light = np.asarray(odor["on"], dtype=bool)[:n]
    chs = _quadrants(gray.shape[1], gray.shape[2])
    energies = [_chamber_energy(gray, ch) for ch in chs]
    n_e = min(e.size for e in energies)
    total = np.zeros(n, dtype=np.float64)
    for e in energies:
        total[1 : 1 + n_e] += e[:n_e]
    live_e = total[1:][live[1 : 1 + n_e]]
    gate = mad_phi_gate(live_e) if live_e.size else mad_phi_gate(np.zeros(1))
    moving = np.zeros(n, dtype=bool)
    moving[1 : 1 + n_e] = total[1 : 1 + n_e] >= gate["threshold"]
    fps = float(meta.get("fps") or 0.0) or 30.0
    on = light & live
    off = (~light) & live
    pulse_s = None
    if odor["first_on"] is not None and odor["last_on"] is not None:
        pulse_s = (odor["last_on"] - odor["first_on"] + 1) / fps
    dists = np.full((n, 4), np.nan)
    for t in np.flatnonzero(live):
        for i, ch in enumerate(chs):
            d = _two_fly_dist(gray[t], ch)
            dists[t, i] = np.nan if d is None else d
    per = []
    on_e: list[float] = []
    off_e: list[float] = []
    on_d: list[float] = []
    off_d: list[float] = []
    for i, e in enumerate(energies):
        ee = np.zeros(n, dtype=np.float64)
        ee[1 : 1 + e.size] += e
        dd = dists[:, i]
        row = {
            "id": i,
            "energy_light_on": _mean_where(ee, on),
            "energy_light_off": _mean_where(ee, off),
            "dist_light_on": _mean_where(dd, on),
            "dist_light_off": _mean_where(dd, off),
        }
        if row["energy_light_on"] is not None:
            on_e.append(row["energy_light_on"])
        if row["energy_light_off"] is not None:
            off_e.append(row["energy_light_off"])
        if row["dist_light_on"] is not None:
            on_d.append(row["dist_light_on"])
        if row["dist_light_off"] is not None:
            off_d.append(row["dist_light_off"])
        per.append(row)
    off_mean = float(np.mean(off_e)) if off_e else None
    on_mean = float(np.mean(on_e)) if on_e else None
    ratio = None
    if off_mean and on_mean is not None:
        ratio = on_mean / off_mean
    dist_on = float(np.mean(on_d)) if on_d else None
    dist_off = float(np.mean(off_d)) if off_d else None
    return {
        "id": "pC1_optogenetic",
        "source": ELIFE,
        "stimulus": "red light 0.02 mW/mm^2 on pC1SS2>CsChrimson",
        "video": meta,
        "n_live": int(live.sum()),
        "n_title": int((~live).sum()),
        "light": {k: v for k, v in odor.items() if k not in {"series", "on"}},
        "light_on": odor["on"][:n],
        "light_rex": odor["series"][:n],
        "pulse_s": pulse_s,
        "n_chambers": 4,
        "motion_gate": {k: v for k, v in gate.items() if k != "bout_lengths"},
        "motion_energy": total.tolist(),
        "moving": moving.astype(bool).tolist(),
        "energy_light_on": on_mean,
        "energy_light_off": off_mean,
        "energy_on_over_off": ratio,
        "dist_light_on": dist_on,
        "dist_light_off": dist_off,
        "dist_on_over_off": (
            None if not dist_off or dist_on is None else dist_on / dist_off
        ),
        "per_chamber": per,
        "seed_map": {
            "light_on": {
                "seed": "pC1_ / dsx / male-specific",
                "why": "published red overlay = pC1SS2 activation; paper: induces tussling",
            },
            "light_off": {"seed": None, "why": "pC1 opto off"},
        },
        "note": (
            "Red-excess overlay is the paper's light marker (same cut as the "
            "odor green dot: observer-max / φ). Motor is chamber frame-diff. "
            "Not a trained pose net and not JAABA."
        ),
    }


def observe_wt(vol: np.ndarray, meta: dict[str, Any]) -> dict[str, Any]:
    gray = vol.mean(axis=3)
    n = int(gray.shape[0])
    live = _live_mask(vol)
    d = np.abs(np.diff(gray, axis=0))
    energy = d.mean(axis=(1, 2))
    total = np.zeros(n, dtype=np.float64)
    total[1:] = energy
    live_e = energy[live[1 : 1 + energy.size]]
    gate = mad_phi_gate(live_e) if live_e.size else mad_phi_gate(np.zeros(1))
    moving = np.zeros(n, dtype=bool)
    moving[1 : 1 + energy.size] = energy >= gate["threshold"]
    return {
        "id": "wt_tussling",
        "source": ELIFE,
        "stimulus": "Canton-S G14 male-male tussling",
        "video": meta,
        "n_live": int(live.sum()),
        "n_title": int((~live).sum()),
        "motion_gate": {k: v for k, v in gate.items() if k != "bout_lengths"},
        "n_moving": int((moving & live).sum()),
        "paint_moving_live": float((moving & live).sum() / max(int(live.sum()), 1)),
        "energy_live_mean": _mean_where(total, live),
        "motion_energy": total.tolist(),
        "moving": moving.astype(bool).tolist(),
        "note": "Wild-type tussling clip after the title card. MAD+φ on frame-diff.",
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


def run_aggression_hops(*, gpu: bool) -> dict[str, Any]:
    print("== Male CNS aggression seeds", flush=True)
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
    contrast: dict[str, Any] = {}
    if MALE_BOOT.exists():
        male = json.loads(MALE_BOOT.read_text(encoding="utf-8"))
        contrast["walking_odor"] = _slim_compare(male.get("compare") or {})
    if COURT_BOOT.exists():
        court = json.loads(COURT_BOOT.read_text(encoding="utf-8"))
        contrast["courtship"] = _slim_compare(
            (court.get("connectome") or {}).get("compare") or {}
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


def loop_table(pc1: dict[str, Any], hops: dict[str, Any]) -> dict[str, Any]:
    cmp_ = hops.get("compare") or {}
    walking = (hops.get("contrast") or {}).get("walking_odor") or {}
    court = (hops.get("contrast") or {}).get("courtship") or {}
    return {
        "aggression_on": {
            "observer": "red overlay = pC1SS2 CsChrimson ON",
            "pulse_s": pc1.get("pulse_s"),
            "energy_on": pc1.get("energy_light_on"),
            "energy_off": pc1.get("energy_light_off"),
            "energy_on_over_off": pc1.get("energy_on_over_off"),
            "note": (
                "Tussling is in-place grappling: chamber translation energy "
                "falls during light-on. Leftover two-fly distance did not "
                "beat the well rim — overlay is the authority marker."
            ),
            "seed": "pC1_ / dsx / male-specific",
            "male_cns_hop2_descending_pC1": _hop_pick(cmp_, "hop_2", "pC1", "descending"),
            "male_cns_hop2_vnc_motor_pC1": _hop_pick(cmp_, "hop_2", "pC1", "vnc_motor"),
            "male_cns_hop1_peak_pC1": ((cmp_.get("hop_1") or {}).get("pC1") or {}).get("top"),
            "male_cns_hop2_descending_dsx": _hop_pick(cmp_, "hop_2", "dsx", "descending"),
            "male_cns_hop1_peak_dsx": ((cmp_.get("hop_1") or {}).get("dsx") or {}).get("top"),
            "male_cns_hop2_descending_male_specific": _hop_pick(
                cmp_, "hop_2", "male_specific", "descending"
            ),
            "male_cns_hop1_peak_male_specific": (
                (cmp_.get("hop_1") or {}).get("male_specific") or {}
            ).get("top"),
        },
        "courtship_contrast_TN1": {
            "seed": "TN1",
            "male_cns_hop2_vnc_motor_TN1": _hop_pick(court, "hop_2", "TN1", "vnc_motor"),
            "male_cns_hop1_peak_TN1": ((court.get("hop_1") or {}).get("TN1") or {}).get("top"),
            "note": "same pC1 cells also seed courtship; TN1 is the song MN contrast",
        },
        "walking_contrast_JO": {
            "seed": "JO",
            "male_cns_hop2_vnc_motor_JO": _hop_pick(walking, "hop_2", "JO", "vnc_motor"),
        },
        "odor_contrast_olfactory": {
            "seed": "olfactory",
            "male_cns_hop2_vnc_motor": _hop_pick(walking, "hop_2", "olfactory", "vnc_motor"),
        },
        "honesty": (
            "pC1_ is the measured type prefix (not LLPC1). The same cells "
            "promote tussling (this observer) and courtship (Movie S1). "
            "dsx and male-specific are measured annotations. Do not invent "
            "an aggression neuron class. Red overlay is the paper's light "
            "marker. GABA only where annotated."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gpu", dest="gpu", action="store_true", default=True)
    ap.add_argument("--cpu", dest="gpu", action="store_false")
    ap.add_argument("--skip-hops", action="store_true")
    args = ap.parse_args(argv)
    print("FSOT fly aggression observer → residual boot", flush=True)
    print(f"  phi={_PHI} R_BIO={_R_BIO} gpu={cuda_name()}", flush=True)
    paths = {spec["id"]: p for spec, p in zip(VIDEOS, ensure_videos())}

    print(f"  video pC1 {paths['pC1_optogenetic']}", flush=True)
    vol4, meta4 = decode_rgb(paths["pC1_optogenetic"])
    print(
        f"  frames={meta4['n_frames']} {meta4['width']}x{meta4['height']} "
        f"{meta4['duration_s']:.2f}s",
        flush=True,
    )
    pc1 = observe_pc1(vol4, meta4)
    del vol4
    print(
        f"  light_on={pc1['light']['n_on']}/{pc1['light']['n']} "
        f"pulse={pc1['pulse_s']:.2f}s energy_on/off={pc1['energy_on_over_off']} "
        f"dist_on/off={pc1.get('dist_on_over_off')}",
        flush=True,
    )

    print(f"  video WT {paths['wt_tussling']}", flush=True)
    vol1, meta1 = decode_rgb(paths["wt_tussling"])
    print(
        f"  frames={meta1['n_frames']} {meta1['width']}x{meta1['height']} "
        f"{meta1['duration_s']:.2f}s",
        flush=True,
    )
    wt = observe_wt(vol1, meta1)
    del vol1
    print(
        f"  live={wt['n_live']} moving={wt['n_moving']} "
        f"paint={wt['paint_moving_live']:.3f}",
        flush=True,
    )

    hops: dict[str, Any] = {}
    if not args.skip_hops:
        hops = run_aggression_hops(gpu=args.gpu)
    loop = loop_table(pc1, hops)
    report = {
        "product": "measured aggression observer + measured Male CNS genetics",
        "free_parameters": 0,
        "phi": _PHI,
        "residual_Biochemistry": _R_BIO,
        "authority": (
            "Gao et al. eLife 13:RP104212 fig4 (pC1SS2 red-light tussling) "
            "and fig1 (Canton-S tussling); Male CNS v1.0 pC1_ / dsx / dimorphism"
        ),
        "honesty": (
            "Red overlay is the published optogenetic marker. pC1_ types also "
            "seed courtship — same cells, different observer. male-specific "
            "and dsx are measured annotations, not invented fight neurons. "
            "GABA only where annotated. Not JAABA, not a trained pose net, "
            "not a thought."
        ),
        "observer": {
            "pC1_optogenetic": {
                k: v
                for k, v in pc1.items()
                if k not in {"light_on", "light_rex", "motion_energy", "moving"}
            },
            "wt_tussling": {
                k: v for k, v in wt.items() if k not in {"motion_energy", "moving"}
            },
        },
        "series": {
            "pC1_light_on": pc1["light_on"],
            "pC1_light_rex": pc1["light_rex"],
            "pC1_motion_energy": pc1["motion_energy"],
            "pC1_moving": pc1["moving"],
            "wt_motion_energy": wt["motion_energy"],
            "wt_moving": wt["moving"],
        },
        "connectome": hops,
        "loop": loop,
        "gpu_name": cuda_name(),
    }
    out = ROOT / "data" / "fly_aggression_flow.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    print(json.dumps(loop, indent=2), flush=True)
    if not pc1.get("pulse_s") or pc1["light"]["n_on"] < 10:
        print("  red overlay did not boot", flush=True)
        return 1
    if not hops.get("programs"):
        print("  hops did not boot", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
