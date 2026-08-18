#!/usr/bin/env python3
"""Read Biohub / Zebrahub 3D+time cell tracks (measured authority).

The 85 GB Kaggle dump lives on the game drive and is *not* copied into
this repo. Tracks are GEFF (zarr v3): cell-center (t, z, y, x) + lineage
edges. Voxel scales are measured (1.625 / 0.40625 / 0.40625 µm).

This is the organism-scale 3-D observer. Protein product Cα (0.13 Å) is
the molecular-scale observer. Same pin, same residual law. 0 free params.

  python scripts/biohub_3d.py
  python scripts/biohub_3d.py --inventory
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

import fsot_compute as fc  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402

BIOHUB_ROOT = Path(r"D:\Kaggle_Biohub_Data")
TRAIN = BIOHUB_ROOT / "train"
TEST = BIOHUB_ROOT / "test"
# Light-sheet voxel (COMPETITION_SPEC / GEFF axis scale)
SCALE_Z_UM = 1.625
SCALE_Y_UM = 0.40625
SCALE_X_UM = 0.40625
PROXY = "44b6_0113de3b"

_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))


def _open(path: Path):
    import zarr

    return zarr.open(str(path), mode="r")


def read_geff_meta(path: Path) -> dict[str, Any]:
    meta_path = path / "zarr.json"
    if not meta_path.exists():
        return {}
    raw = json.loads(meta_path.read_text(encoding="utf-8"))
    geff = (raw.get("attributes") or {}).get("geff") or {}
    extra = geff.get("extra") or {}
    axes = {a["name"]: a for a in (geff.get("axes") or []) if "name" in a}
    return {
        "directed": bool(geff.get("directed")),
        "estimated_nodes": extra.get("estimated_number_of_nodes"),
        "axes": axes,
        "t_max": (axes.get("t") or {}).get("max"),
        "scale_z": (axes.get("z") or {}).get("scale", SCALE_Z_UM),
        "scale_y": (axes.get("y") or {}).get("scale", SCALE_Y_UM),
        "scale_x": (axes.get("x") or {}).get("scale", SCALE_X_UM),
    }


def read_geff(path: Path) -> dict[str, Any]:
    """Measured cell centers + lineage edges, voxel and µm."""
    root = _open(path)
    ids = np.asarray(root["nodes/ids"][:])
    t = np.asarray(root["nodes/props/t/values"][:]).astype(np.int64)
    zv = np.asarray(root["nodes/props/z/values"][:]).astype(np.float64)
    yv = np.asarray(root["nodes/props/y/values"][:]).astype(np.float64)
    xv = np.asarray(root["nodes/props/x/values"][:]).astype(np.float64)
    e_ids = np.asarray(root["edges/ids"][:])
    meta = read_geff_meta(path)
    sz = float(meta.get("scale_z") or SCALE_Z_UM)
    sy = float(meta.get("scale_y") or SCALE_Y_UM)
    sx = float(meta.get("scale_x") or SCALE_X_UM)
    xyz_vox = np.stack([zv, yv, xv], axis=1)
    xyz_um = np.stack([zv * sz, yv * sy, xv * sx], axis=1)
    outdeg: dict[int, int] = {}
    edges = [(int(a), int(b)) for a, b in e_ids]
    id_to_i = {int(n): i for i, n in enumerate(ids)}
    for a, b in edges:
        outdeg[a] = outdeg.get(a, 0) + 1
    n_div = sum(1 for k, c in outdeg.items() if c >= 2)
    steps = []
    for a, b in edges:
        if a in id_to_i and b in id_to_i:
            steps.append(float(np.linalg.norm(xyz_um[id_to_i[b]] - xyz_um[id_to_i[a]])))
    steps_a = np.asarray(steps, dtype=float) if steps else np.zeros(0)
    return {
        "path": str(path),
        "dataset": path.stem,
        "n_nodes": int(len(ids)),
        "n_edges": int(len(edges)),
        "n_divisions": int(n_div),
        "t_min": int(t.min()) if len(t) else None,
        "t_max": int(t.max()) if len(t) else None,
        "ids": ids,
        "t": t,
        "xyz_vox": xyz_vox,
        "xyz_um": xyz_um,
        "edges": edges,
        "span_um": (
            [float(xyz_um[:, i].max() - xyz_um[:, i].min()) for i in range(3)]
            if len(xyz_um)
            else [0.0, 0.0, 0.0]
        ),
        "step_median_um": float(np.median(steps_a)) if len(steps_a) else None,
        "step_mean_um": float(np.mean(steps_a)) if len(steps_a) else None,
        "meta": meta,
    }


def track_residual(tracks: dict[str, Any]) -> dict[str, Any]:
    """Biochemistry residual on measured parent→child steps.

    The tracks are the observation. Residual scales the interface; it does
    not invent cell positions. Median step is the measured length.
    """
    xyz = tracks["xyz_um"]
    id_to_i = {int(n): i for i, n in enumerate(tracks["ids"])}
    steps = []
    for a, b in tracks["edges"]:
        if a in id_to_i and b in id_to_i:
            steps.append(float(np.linalg.norm(xyz[id_to_i[b]] - xyz[id_to_i[a]])))
    if not steps:
        return {"residual_Biochemistry": _R_BIO, "n_steps": 0, "energy": None}
    s = np.asarray(steps, dtype=float)
    med = float(np.median(s))
    # Physical_Chemistry analog: (L − L_measured)², residual-weighted.
    energy = float(_R_BIO * np.mean((s - med) ** 2))
    return {
        "residual_Biochemistry": _R_BIO,
        "n_steps": int(len(s)),
        "step_median_um": med,
        "step_mae_um": float(np.mean(np.abs(s - med))),
        "energy": energy,
        "phi": _PHI,
    }


def inventory(train_dir: Path = TRAIN) -> list[dict[str, Any]]:
    rows = []
    for geff in sorted(train_dir.glob("*.geff")):
        try:
            root = _open(geff)
            n = int(root["nodes/ids"].shape[0])
            e = int(root["edges/ids"].shape[0])
            meta = read_geff_meta(geff)
        except Exception as exc:
            rows.append({"dataset": geff.stem, "error": str(exc)})
            continue
        rows.append(
            {
                "dataset": geff.stem,
                "n_nodes": n,
                "n_edges": e,
                "estimated_nodes": meta.get("estimated_nodes"),
                "t_max": meta.get("t_max"),
            }
        )
    return rows


MATCH_UM = 7.0  # official Biohub centroid match radius
# NMS = φ³ µm (≈4.24). φ⁴ (6.85) sat on the 7 µm match radius and
# merged an annotated cell with an unannotated neighbor — nearest peak
# landed 8–12 µm away. φ³ splits them; official match is still 7 µm.
NMS_UM = float(_PHI ** 3)
LINK_UM = float(_PHI ** 4)  # first-pass link; residual widens with measured step


def open_volume(zarr_path: Path):
    """OME-Zarr level-0 array (T, Z, Y, X). Does not load voxels until sliced."""
    return _open(zarr_path)["0"]


def read_frame(zarr_path: Path, t: int) -> np.ndarray:
    """One time-point as (Z, Y, X) uint16 — the actual light-sheet voxels."""
    arr = open_volume(zarr_path)
    t = int(np.clip(t, 0, arr.shape[0] - 1))
    return np.asarray(arr[t])


def frame_stats(vol: np.ndarray) -> dict[str, float]:
    v = vol.astype(np.float64)
    return {
        "min": float(v.min()),
        "max": float(v.max()),
        "mean": float(v.mean()),
        "median": float(np.median(v)),
        "p99": float(np.percentile(v, 99)),
        "p999": float(np.percentile(v, 99.9)),
    }


def _mad_threshold(sm: np.ndarray) -> float:
    """Foreground gate: median + φ·MAD. If that paints > 1/φ of voxels, use φ²."""
    med = float(np.median(sm))
    mad = float(np.median(np.abs(sm - med))) + 1e-6
    thr = med + _PHI * mad
    if float((sm >= thr).mean()) > 1.0 / _PHI:
        thr = med + (_PHI * _PHI) * mad
    return thr


def _centroid_blob(
    vol: np.ndarray,
    zyx: np.ndarray,
    scale_zyx: tuple[float, float, float],
    radius_um: float,
) -> np.ndarray:
    """Intensity-weighted first moment of the observed blob (measured center)."""
    sc = np.asarray(scale_zyx, dtype=float)
    rad = np.maximum(1, np.round(radius_um / sc)).astype(int)
    z0, y0, x0 = [int(round(c)) for c in zyx]
    z1, z2 = max(0, z0 - rad[0]), min(vol.shape[0], z0 + rad[0] + 1)
    y1, y2 = max(0, y0 - rad[1]), min(vol.shape[1], y0 + rad[1] + 1)
    x1, x2 = max(0, x0 - rad[2]), min(vol.shape[2], x0 + rad[2] + 1)
    patch = vol[z1:z2, y1:y2, x1:x2].astype(np.float64)
    if patch.size == 0:
        return zyx.astype(np.float64)
    # Half-max of this blob (measured FWHM core). Patch-median pulled the
    # centroid into the dim halo and moved some hits out to ~5 µm.
    peak = float(patch.max())
    w = np.maximum(patch - 0.5 * peak, 0.0)
    wsum = float(w.sum())
    if wsum < 1e-6:
        return np.array([z0, y0, x0], dtype=np.float64)
    zz, yy, xx = np.indices(patch.shape, dtype=np.float64)
    return np.array(
        [
            z1 + float((w * zz).sum() / wsum),
            y1 + float((w * yy).sum() / wsum),
            x1 + float((w * xx).sum() / wsum),
        ]
    )


def detect_peaks_frame(
    vol: np.ndarray,
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    nms_um: float = NMS_UM,
    percentile: float | None = None,
    xy_stride: int = 2,
) -> np.ndarray:
    """Observed nuclei: local max, then measured intensity centroid.

    Threshold is median + φ·MAD (seed), not a free percentile. The 7 µm
    miss was peak-voxel vs GT cell center — we now take the first moment
    of the bright blob. Returns (N, 3) z,y,x in full-resolution voxels.
    """
    from scipy.ndimage import gaussian_filter, maximum_filter

    full = vol.astype(np.float32)
    if xy_stride > 1:
        v = full[:, ::xy_stride, ::xy_stride]
        sc = (
            scale_zyx[0],
            scale_zyx[1] * xy_stride,
            scale_zyx[2] * xy_stride,
        )
    else:
        v = full
        sc = scale_zyx
    sm = gaussian_filter(v, sigma=(0.4, 0.8, 0.8))
    if percentile is None:
        thr = _mad_threshold(sm)
    else:
        thr = float(np.percentile(sm, percentile))
    size = tuple(int(max(1, 2 * round(nms_um / s) + 1)) for s in sc)
    mx = maximum_filter(sm, size=size, mode="nearest")
    peaks = (sm == mx) & (sm >= thr)
    # Residual second collapse: zero the first NMS balls and take leftover
    # local max. Unannotated neighbors left 7–12 µm ghosts on the dense set.
    sm2 = sm.copy()
    rad = np.maximum(1, np.round(np.array(nms_um) / np.array(sc))).astype(int)
    pz, py, px = np.where(peaks)
    for zi0, yi0, xi0 in zip(pz, py, px):
        z1, z2 = max(0, zi0 - rad[0]), min(sm2.shape[0], zi0 + rad[0] + 1)
        y1, y2 = max(0, yi0 - rad[1]), min(sm2.shape[1], yi0 + rad[1] + 1)
        x1, x2 = max(0, xi0 - rad[2]), min(sm2.shape[2], xi0 + rad[2] + 1)
        sm2[z1:z2, y1:y2, x1:x2] = 0.0
    mx2 = maximum_filter(sm2, size=size, mode="nearest")
    peaks2 = (sm2 == mx2) & (sm2 >= thr)
    zi1, yi1, xi1 = np.where(peaks)
    zi2, yi2, xi2 = np.where(peaks2)
    zi = np.concatenate([zi1, zi2])
    yi = np.concatenate([yi1, yi2])
    xi = np.concatenate([xi1, xi2])
    tier = np.concatenate(
        [
            np.ones(len(zi1), dtype=np.int8),
            np.full(len(zi2), 2, dtype=np.int8),
        ]
    )
    if len(zi) == 0:
        return np.zeros((0, 3), dtype=np.float64), np.zeros((0,), dtype=np.int8)
    y = yi.astype(np.float64) * xy_stride
    x = xi.astype(np.float64) * xy_stride
    raw = np.stack([zi.astype(np.float64), y, x], axis=1)
    out = np.empty_like(raw)
    for i, p in enumerate(raw):
        out[i] = _centroid_blob(full, p, scale_zyx, nms_um)
    return out, tier


def detect_video(
    zarr_path: Path,
    *,
    t_indices: list[int] | None = None,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    nms_um: float = NMS_UM,
    percentile: float | None = None,
) -> np.ndarray:
    """Dense detections (t, z, y, x) voxels for selected frames."""
    arr = open_volume(zarr_path)
    T = int(arr.shape[0])
    if t_indices is None:
        t_indices = list(range(T))
    rows = []
    n_t = len(t_indices)
    for k, t in enumerate(t_indices):
        if k % 10 == 0 or k + 1 == n_t:
            print(f"  detect frame {int(t)} ({k + 1}/{n_t})", flush=True)
        vol = np.asarray(arr[int(t)])
        p, tier = detect_peaks_frame(
            vol, scale_zyx=scale_zyx, nms_um=nms_um, percentile=percentile
        )
        if len(p):
            tt = np.full((len(p), 1), float(t))
            rows.append(np.hstack([tt, p, tier.reshape(-1, 1).astype(np.float64)]))
    if not rows:
        return np.zeros((0, 5), dtype=np.float64)
    return np.vstack(rows)


def link_tracks(
    pred_tzyx: np.ndarray,
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    max_um: float = LINK_UM,
) -> list[tuple[int, int]]:
    """Link detections in consecutive frames (nearest unused, ≤ φ⁵ µm).

    Edges are predicted lineage (same cell or daughters). Residual does not
    invent a second length — the gate is the seed nucleus scale.
    Returns pairs of *row indices* into pred_tzyx.
    """
    if len(pred_tzyx) < 2:
        return []
    sc = np.asarray(scale_zyx, dtype=float)
    frames = sorted({int(t) for t in pred_tzyx[:, 0]})
    edges: list[tuple[int, int]] = []
    idx_by_t = {t: np.where(pred_tzyx[:, 0] == t)[0] for t in frames}
    for t0, t1 in zip(frames, frames[1:]):
        if t1 != t0 + 1:
            continue
        a = idx_by_t[t0]
        b = idx_by_t[t1]
        if len(a) == 0 or len(b) == 0:
            continue
        pa = pred_tzyx[a][:, 1:4]
        pb = pred_tzyx[b][:, 1:4]
        d = np.sqrt((((pa[:, None, :] - pb[None, :, :]) * sc) ** 2).sum(axis=2))
        from scipy.optimize import linear_sum_assignment

        cap = d.copy()
        cap[cap > max_um] = max_um * 4.0
        ri, cj = linear_sum_assignment(cap)
        for ia, jb in zip(ri, cj):
            if d[ia, jb] <= max_um:
                edges.append((int(a[ia]), int(b[jb])))
    return edges


def link_tracks_residual(
    pred_tzyx: np.ndarray,
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    seed_um: float = LINK_UM,
) -> tuple[list[tuple[int, int]], dict[str, float]]:
    """Link, then widen the gate to φ · measured median step (residual, not free)."""
    sc = np.asarray(scale_zyx, dtype=float)
    e1 = link_tracks(pred_tzyx, scale_zyx=scale_zyx, max_um=seed_um)
    steps = []
    for i, j in e1:
        steps.append(
            float(np.linalg.norm((pred_tzyx[j, 1:4] - pred_tzyx[i, 1:4]) * sc))
        )
    med = float(np.median(steps)) if steps else seed_um
    # Do not floor at φ⁴ — that 6.85 µm gate let Hungarian steal the
    # true continuation (proxy lineage 0.78→0.56). Gate is φ · measured step.
    gate = max(seed_um, _PHI * med) if steps else seed_um
    e2 = link_tracks(pred_tzyx, scale_zyx=scale_zyx, max_um=gate)
    return e2, {
        "seed_um": seed_um,
        "measured_step_median_um": med if steps else None,
        "residual_gate_um": gate,
        "n_edges_seed": len(e1),
        "n_edges_residual": len(e2),
    }


def lineage_recall(
    pred_tzyx: np.ndarray,
    pred_edges: list[tuple[int, int]],
    tracks: dict[str, Any],
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    max_um: float = MATCH_UM,
) -> dict[str, Any]:
    """Fraction of measured GEFF edges whose both ends match a predicted link."""
    gt = np.column_stack([tracks["t"].astype(float), tracks["xyz_vox"]])
    id_to_g = {int(n): i for i, n in enumerate(tracks["ids"])}
    # GT node -> pred row (if matched)
    sc = np.asarray(scale_zyx, dtype=float)
    g2p: dict[int, int] = {}
    for t in sorted({int(x) for x in gt[:, 0]}):
        gi = [i for i, tt in enumerate(gt[:, 0]) if int(tt) == t]
        pi = [i for i, tt in enumerate(pred_tzyx[:, 0]) if int(tt) == t] if len(pred_tzyx) else []
        if pred_tzyx.shape[1] >= 5:
            pi = [i for i in pi if int(pred_tzyx[i, 4]) == 1]
        if not gi or not pi:
            continue
        gxyz = gt[gi][:, 1:4]
        pxyz = pred_tzyx[pi][:, 1:4]
        d = np.sqrt((((gxyz[:, None, :] - pxyz[None, :, :]) * sc) ** 2).sum(axis=2))
        used: set[int] = set()
        for row, gidx in enumerate(gi):
            j = int(np.argmin(d[row]))
            if j in used:
                alt = [int(k) for k in np.argsort(d[row]) if int(k) not in used]
                if not alt:
                    continue
                j = alt[0]
            if d[row, j] <= max_um:
                used.add(j)
                g2p[gidx] = int(pi[j])
    pred_set = {tuple(sorted(e)) for e in pred_edges}
    hit = 0
    n_e = 0
    for a, b in tracks["edges"]:
        if a not in id_to_g or b not in id_to_g:
            continue
        n_e += 1
        pa, pb = g2p.get(id_to_g[a]), g2p.get(id_to_g[b])
        if pa is None or pb is None:
            continue
        if tuple(sorted((pa, pb))) in pred_set:
            hit += 1
    return {
        "n_gt_edges": n_e,
        "n_pred_edges": len(pred_edges),
        "n_edge_matched": hit,
        "edge_recall": float(hit / n_e) if n_e else None,
        "n_gt_nodes_matched": len(g2p),
    }


def match_centroids(
    pred_tzyx: np.ndarray,
    gt_tzyx: np.ndarray,
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    max_um: float = MATCH_UM,
) -> dict[str, Any]:
    """Per-frame nearest match in µm. Returns recall of GT and match distances."""
    sc = np.asarray(scale_zyx, dtype=float)
    if len(gt_tzyx) == 0:
        return {"n_gt": 0, "n_pred": int(len(pred_tzyx)), "n_matched": 0, "recall": None}
    frames = sorted({int(t) for t in gt_tzyx[:, 0]})
    matched = 0
    dists: list[float] = []
    for t in frames:
        g = gt_tzyx[gt_tzyx[:, 0] == t][:, 1:4]
        p = pred_tzyx[pred_tzyx[:, 0] == t][:, 1:4] if len(pred_tzyx) else np.zeros((0, 3))
        if len(g) == 0:
            continue
        if len(p) == 0:
            continue
        d = np.sqrt((((g[:, None, :] - p[None, :, :]) * sc) ** 2).sum(axis=2))
        if len(g) > 1 and len(p) > 1:
            from scipy.optimize import linear_sum_assignment

            cap = d.copy()
            cap[cap > max_um] = max_um * 4
            ri, cj = linear_sum_assignment(cap)
            for i, j in zip(ri, cj):
                if d[i, j] <= max_um:
                    matched += 1
                    dists.append(float(d[i, j]))
            continue
        used = set()
        for i in range(len(g)):
            j = int(np.argmin(d[i]))
            if j in used:
                order = np.argsort(d[i])
                j = next((int(k) for k in order if int(k) not in used), j)
            if j in used:
                continue
            if d[i, j] <= max_um:
                used.add(j)
                matched += 1
                dists.append(float(d[i, j]))
    n_gt = int(len(gt_tzyx))
    return {
        "n_gt": n_gt,
        "n_pred": int(len(pred_tzyx)),
        "n_matched": int(matched),
        "recall": float(matched / n_gt) if n_gt else None,
        "match_median_um": float(np.median(dists)) if dists else None,
        "match_mean_um": float(np.mean(dists)) if dists else None,
        "max_um": max_um,
    }


def intensity_at(vol: np.ndarray, zyx: np.ndarray) -> np.ndarray:
    """Sample voxel intensity at (possibly fractional) z,y,x."""
    if len(zyx) == 0:
        return np.zeros(0)
    z = np.clip(np.round(zyx[:, 0]).astype(int), 0, vol.shape[0] - 1)
    y = np.clip(np.round(zyx[:, 1]).astype(int), 0, vol.shape[1] - 1)
    x = np.clip(np.round(zyx[:, 2]).astype(int), 0, vol.shape[2] - 1)
    return vol[z, y, x].astype(np.float64)


def read_volume_meta(zarr_path: Path) -> dict[str, Any]:
    """OME-Zarr image header only — does not load the 85 GB voxels."""
    meta_path = zarr_path / "zarr.json"
    out: dict[str, Any] = {"path": str(zarr_path), "exists": zarr_path.exists()}
    if meta_path.exists():
        out["root_meta"] = json.loads(meta_path.read_text(encoding="utf-8"))
    try:
        root = _open(zarr_path)
        keys = list(root.keys()) if hasattr(root, "keys") else []
        out["keys"] = keys[:16]
        if "0" in keys:
            arr = root["0"]
            out["level0_shape"] = [int(x) for x in arr.shape]
            out["level0_dtype"] = str(arr.dtype)
    except Exception as exc:
        out["open_error"] = str(exc)
    return out


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", action="store_true")
    ap.add_argument(
        "--voxels",
        action="store_true",
        help="Read light-sheet voxels and detect nuclei (local maxima).",
    )
    ap.add_argument("--dataset", default=PROXY)
    ap.add_argument("--root", default=str(BIOHUB_ROOT))
    args = ap.parse_args(argv)
    root = Path(args.root)
    train = root / "train"
    if not train.exists():
        print(f"Biohub data not mounted at {train}", file=sys.stderr)
        return 2
    if args.inventory:
        rows = inventory(train)
        ok = [r for r in rows if "error" not in r]
        print(f"train GEFFs readable: {len(ok)}/{len(rows)}")
        if ok:
            nodes = [r["n_nodes"] for r in ok]
            print(
                f"  nodes min/med/max {min(nodes)} / "
                f"{int(np.median(nodes))} / {max(nodes)}"
            )
        out = ROOT / "data" / "biohub_3d_inventory.json"
        out.write_text(json.dumps({"n": len(rows), "geffs": rows}, indent=2), encoding="utf-8")
        print(f"  wrote {out}")
        return 0
    if args.voxels:
        geff = train / f"{args.dataset}.geff"
        volp = train / f"{args.dataset}.zarr"
        if not volp.exists():
            print(f"no volume at {volp}", file=sys.stderr)
            return 2
        tracks = read_geff(geff)
        arr = open_volume(volp)
        T = int(arr.shape[0])
        # Intensity audit: GT cells sit on bright voxels (frame 0 + mid).
        t0 = int(tracks["t"][0]) if len(tracks["t"]) else 0
        frame0 = np.asarray(arr[t0])
        st = frame_stats(frame0)
        gt0 = tracks["xyz_vox"][tracks["t"] == t0]
        gt_int = intensity_at(frame0, gt0) if len(gt0) else np.zeros(0)
        pred = detect_video(volp)
        gt_tzyx = np.column_stack(
            [tracks["t"].astype(float), tracks["xyz_vox"]]
        )
        xyz4 = pred[:, :4]
        hit = match_centroids(xyz4, gt_tzyx)
        hit12 = match_centroids(xyz4, gt_tzyx, max_um=12.0)
        if pred.shape[1] >= 5:
            pri = np.where(pred[:, 4] == 1)[0]
            e_loc, link_meta = link_tracks_residual(xyz4[pri])
            edges = [(int(pri[i]), int(pri[j])) for i, j in e_loc]
            link_meta["n_primary"] = int(len(pri))
            link_meta["n_residual_peaks"] = int((pred[:, 4] == 2).sum())
        else:
            edges, link_meta = link_tracks_residual(xyz4)
        lin = lineage_recall(pred, edges, tracks)
        lin12 = lineage_recall(pred, edges, tracks, max_um=12.0)
        out = {
            "dataset": args.dataset,
            "volume_shape": [int(x) for x in arr.shape],
            "n_frames_read": T,
            "frame0_stats": st,
            "gt_on_frame0": int(len(gt0)),
            "gt_intensity_median": float(np.median(gt_int)) if len(gt_int) else None,
            "bg_intensity_median": st["median"],
            "gt_over_bg": (
                float(np.median(gt_int) / st["median"]) if len(gt_int) and st["median"] else None
            ),
            "n_dense_detections": hit["n_pred"],
            "gt_recall_7um": hit["recall"],
            "gt_recall_12um": hit12["recall"],
            "n_gt_matched_7um": hit["n_matched"],
            "n_gt_matched_12um": hit12["n_matched"],
            "n_gt": hit["n_gt"],
            "match_median_um": hit["match_median_um"],
            "n_pred_edges": lin["n_pred_edges"],
            "lineage_edge_recall_7um": lin["edge_recall"],
            "lineage_edge_recall_12um": lin12["edge_recall"],
            "n_gt_edges_matched": lin["n_edge_matched"],
            "n_gt_edges": lin["n_gt_edges"],
            "nms_um": NMS_UM,
            "link_um": link_meta.get("residual_gate_um", LINK_UM),
            "link_meta": link_meta,
            "note": (
                "Centroid = half-max first moment. Gate = median+φ·MAD. "
                "NMS = φ³ µm. Residual second collapse fills leftover "
                "brightness (7 µm recall). Lineage links first collapse only."
            ),
            "estimated_true_cells": tracks["meta"].get("estimated_nodes"),
            "authority": "OME-Zarr voxels + measured blob centroids (no trained net)",
            "free_parameters": 0,
        }
        print(json.dumps(out, indent=2))
        dest = ROOT / "data" / "biohub_3d_voxels.json"
        bundle: dict[str, Any] = {"free_parameters": 0, "datasets": {}}
        if dest.exists():
            try:
                prev = json.loads(dest.read_text(encoding="utf-8"))
                if isinstance(prev, dict) and "datasets" in prev:
                    bundle = prev
                elif isinstance(prev, dict) and "dataset" in prev:
                    bundle["datasets"][str(prev["dataset"])] = prev
            except Exception:
                pass
        bundle["datasets"][str(out["dataset"])] = out
        dest.write_text(json.dumps(bundle, indent=2), encoding="utf-8")
        print(f"wrote {dest}")
        return 0
    geff = train / f"{args.dataset}.geff"
    vol = train / f"{args.dataset}.zarr"
    tracks = read_geff(geff)
    res = track_residual(tracks)
    vmeta = read_volume_meta(vol) if vol.exists() else {"exists": False}
    slim = {
        "dataset": tracks["dataset"],
        "n_nodes": tracks["n_nodes"],
        "n_edges": tracks["n_edges"],
        "n_divisions": tracks["n_divisions"],
        "t_range": [tracks["t_min"], tracks["t_max"]],
        "span_um_zyx": tracks["span_um"],
        "step_median_um": tracks["step_median_um"],
        "residual": res,
        "volume": {
            k: vmeta.get(k)
            for k in ("exists", "level0_shape", "level0_dtype", "keys")
        },
        "scale_um": [SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM],
        "authority": "measured GEFF cell centers (exclude nothing — this is not a holdout fold)",
        "free_parameters": 0,
    }
    print(json.dumps(slim, indent=2))
    dest = ROOT / "data" / "biohub_3d_proxy.json"
    dest.write_text(json.dumps(slim, indent=2), encoding="utf-8")
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
