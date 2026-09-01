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
# Detections stay on the game drive (not git). Re-link without re-detecting.
DETECT_CACHE = BIOHUB_ROOT / "_fsot_detect_cache"


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
        return (
            np.zeros((0, 3), dtype=np.float64),
            np.zeros((0,), dtype=np.int8),
            np.zeros((0,), dtype=np.float64),
        )
    y = yi.astype(np.float64) * xy_stride
    x = xi.astype(np.float64) * xy_stride
    raw = np.stack([zi.astype(np.float64), y, x], axis=1)
    out = np.empty_like(raw)
    inten = np.empty(len(raw), dtype=np.float64)
    for i, p in enumerate(raw):
        out[i] = _centroid_blob(full, p, scale_zyx, nms_um)
        zc, yc, xc = [int(round(c)) for c in out[i]]
        zc = int(np.clip(zc, 0, full.shape[0] - 1))
        yc = int(np.clip(yc, 0, full.shape[1] - 1))
        xc = int(np.clip(xc, 0, full.shape[2] - 1))
        inten[i] = float(full[zc, yc, xc])
    return out, tier, inten


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
        p, tier, inten = detect_peaks_frame(
            vol, scale_zyx=scale_zyx, nms_um=nms_um, percentile=percentile
        )
        if len(p):
            tt = np.full((len(p), 1), float(t))
            rows.append(
                np.hstack(
                    [
                        tt,
                        p,
                        tier.reshape(-1, 1).astype(np.float64),
                        inten.reshape(-1, 1),
                    ]
                )
            )
    if not rows:
        return np.zeros((0, 6), dtype=np.float64)
    return np.vstack(rows)


def detect_cache_path(dataset: str) -> Path:
    return DETECT_CACHE / f"{dataset}.npy"


def load_or_detect(dataset: str, volp: Path) -> np.ndarray:
    """Reuse measured peaks. Detection is the expensive observer, not the link."""
    cache = detect_cache_path(dataset)
    if cache.exists():
        pred = np.load(cache)
        print(f"  detect cache {cache} n={len(pred)}", flush=True)
        return pred
    pred = detect_video(volp)
    DETECT_CACHE.mkdir(parents=True, exist_ok=True)
    np.save(cache, pred)
    print(f"  wrote detect cache {cache}", flush=True)
    return pred


def _indices_at(pred: np.ndarray, t: int, tier: int | None = None) -> np.ndarray:
    m = pred[:, 0].astype(int) == int(t)
    if tier is not None and pred.shape[1] >= 5:
        m &= pred[:, 4].astype(int) == int(tier)
    return np.where(m)[0]


def _pair_cost(
    pred: np.ndarray,
    src: np.ndarray | list[int],
    dst: np.ndarray | list[int],
    sc: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """µm distance and intensity-identity cost between two index sets."""
    src = np.asarray(src, dtype=int)
    dst = np.asarray(dst, dtype=int)
    pa = pred[src][:, 1:4]
    pb = pred[dst][:, 1:4]
    d = np.sqrt((((pa[:, None, :] - pb[None, :, :]) * sc) ** 2).sum(axis=2))
    if pred.shape[1] >= 6:
        ia_ = pred[src, 5]
        ib_ = pred[dst, 5]
        rel = np.abs(ia_[:, None] - ib_[None, :]) / (
            np.maximum(ia_[:, None], ib_[None, :]) + 1e-6
        )
        cost = d * (1.0 + rel)
    else:
        cost = d
    return d, cost


def link_tracks(
    pred_tzyx: np.ndarray,
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    max_um: float = LINK_UM,
    use_velocity: bool = True,
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
    parent: dict[int, int] = {}
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
        # T2 linear: x̂_{t+1} = x_t + (x_t − x_{t−1}) when a parent exists.
        # 223 dense Jaccard FN were within φ⁴ of the true child but assigned
        # elsewhere — velocity can point at a neighbor. Position-only is the AB.
        pred_pos = pa.copy()
        if use_velocity:
            for k, i in enumerate(a):
                par = parent.get(int(i))
                if par is not None:
                    pred_pos[k] = pa[k] + (pa[k] - pred_tzyx[par, 1:4])
        d_act = np.sqrt((((pa[:, None, :] - pb[None, :, :]) * sc) ** 2).sum(axis=2))
        d_inn = np.sqrt((((pred_pos[:, None, :] - pb[None, :, :]) * sc) ** 2).sum(axis=2))
        from scipy.optimize import linear_sum_assignment

        if pred_tzyx.shape[1] >= 6:
            ia_ = pred_tzyx[a, 5]
            ib_ = pred_tzyx[b, 5]
            rel = np.abs(ia_[:, None] - ib_[None, :]) / (np.maximum(ia_[:, None], ib_[None, :]) + 1e-6)
            cost = d_inn * (1.0 + rel)
        else:
            cost = d_inn
        cap = cost.copy()
        far = (d_inn > max_um) | (d_act > float(_PHI ** 5))
        cap[far] = max_um * 4.0
        ri, cj = linear_sum_assignment(cap)
        for ia, jb in zip(ri, cj):
            if d_inn[ia, jb] <= max_um and d_act[ia, jb] <= float(_PHI ** 5):
                i, j = int(a[ia]), int(b[jb])
                edges.append((i, j))
                parent[j] = i
    return edges


def link_tracks_residual(
    pred_tzyx: np.ndarray,
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    seed_um: float = LINK_UM,
    use_velocity: bool = True,
) -> tuple[list[tuple[int, int]], dict[str, float]]:
    """Link, then widen the gate to φ · measured median step (residual, not free)."""
    sc = np.asarray(scale_zyx, dtype=float)
    e1 = link_tracks(
        pred_tzyx, scale_zyx=scale_zyx, max_um=seed_um, use_velocity=use_velocity
    )
    steps = []
    for i, j in e1:
        steps.append(
            float(np.linalg.norm((pred_tzyx[j, 1:4] - pred_tzyx[i, 1:4]) * sc))
        )
    med = float(np.median(steps)) if steps else seed_um
    # Do not floor at φ⁴ — that 6.85 µm gate let Hungarian steal the
    # true continuation (proxy lineage 0.78→0.56). Gate is φ · measured step.
    gate = max(seed_um, _PHI * med) if steps else seed_um
    e2 = link_tracks(
        pred_tzyx, scale_zyx=scale_zyx, max_um=gate, use_velocity=use_velocity
    )
    return e2, {
        "seed_um": seed_um,
        "measured_step_median_um": med if steps else None,
        "residual_gate_um": gate,
        "n_edges_seed": len(e1),
        "n_edges_residual": len(e2),
        "use_velocity": use_velocity,
    }


def product_detections(pred: np.ndarray, *, iso_um: float | None = None) -> np.ndarray:
    """Primary peaks plus isolated residual (farther than iso_um from every primary).

    Halo residual next to a primary steals the 7 µm GT match if it is a
    second node. Isolated leftover nuclei stay. Halo inside the NMS shell
    is folded into the primary centroid (intensity-weighted) so 7–12 µm
    leftover brightness can move the reported center without extra nodes.
    Default iso is NMS (φ³ µm).
    """
    if len(pred) == 0 or pred.shape[1] < 5:
        return pred
    sc = np.array([SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM])
    if iso_um is None:
        iso_um = float(NMS_UM)
    keep: list[int] = []
    for t in sorted({int(x) for x in pred[:, 0]}):
        pri = _indices_at(pred, t, 1)
        res = _indices_at(pred, t, 2)
        keep.extend(int(i) for i in pri)
        if len(res) == 0:
            continue
        if len(pri) == 0:
            keep.extend(int(i) for i in res)
            continue
        dmin = np.sqrt(
            (
                ((pred[res][:, 1:4][:, None, :] - pred[pri][None, :, 1:4]) * sc)
                ** 2
            ).sum(axis=2)
        ).min(axis=1)
        keep.extend(int(res[k]) for k in range(len(res)) if float(dmin[k]) > iso_um)
    if not keep:
        return pred[:0]
    idx = np.array(sorted(set(keep)), dtype=int)
    out = pred[idx].copy()
    # Fold in-shell residual into the primary centroid (not a second node).
    for t in sorted({int(x) for x in out[:, 0]}):
        po = np.where((out[:, 0].astype(int) == t) & (out[:, 4] == 1))[0]
        res = pred[(pred[:, 0].astype(int) == t) & (pred[:, 4] == 2)]
        if len(po) == 0 or len(res) == 0:
            continue
        d = np.sqrt(
            (
                ((out[po][:, 1:4][:, None, :] - res[None, :, 1:4]) * sc)
                ** 2
            ).sum(axis=2)
        )
        for k, row in enumerate(po):
            near = d[k] <= iso_um
            if not np.any(near):
                continue
            xyz = np.vstack([out[row, 1:4], res[near][:, 1:4]])
            if pred.shape[1] >= 6:
                w = np.concatenate([[out[row, 5]], res[near, 5]])
            else:
                w = np.ones(len(xyz))
            wsum = float(w.sum())
            if wsum <= 0:
                continue
            out[row, 1:4] = (xyz * w[:, None]).sum(axis=0) / wsum
    return out


def _hungarian_gate(
    pred: np.ndarray,
    src: list[int] | np.ndarray,
    dst: list[int] | np.ndarray,
    sc: np.ndarray,
    max_um: float,
) -> list[tuple[int, int]]:
    """1-1 links with intensity identity, reject if d > max_um."""
    from scipy.optimize import linear_sum_assignment

    src = [int(i) for i in src]
    dst = [int(i) for i in dst]
    if not src or not dst:
        return []
    d, cost = _pair_cost(pred, src, dst, sc)
    cap = cost.copy()
    cap[d > max_um] = max_um * 4.0
    ri, cj = linear_sum_assignment(cap)
    out: list[tuple[int, int]] = []
    for ia, jb in zip(ri, cj):
        if d[ia, jb] <= max_um:
            out.append((src[ia], dst[jb]))
    return out


def _pairwise_swap_edges(
    pred: np.ndarray,
    edges: list[tuple[int, int]],
    sc: np.ndarray,
    max_um: float,
) -> tuple[list[tuple[int, int]], int]:
    """2-opt on consecutive-frame pairs. Hungarian per stage can leave a
    cheaper swap across stages (223 dense FN were within φ⁴ of the true child).
    """
    from collections import defaultdict

    def dist(i, j) -> float:
        return float(np.linalg.norm((pred[j, 1:4] - pred[i, 1:4]) * sc))

    by_t: dict[int, list[tuple[int, int]]] = defaultdict(list)
    other: list[tuple[int, int]] = []
    for i, j in edges:
        if int(pred[j, 0]) == int(pred[i, 0]) + 1:
            by_t[int(pred[i, 0])].append((i, j))
        else:
            other.append((i, j))
    nswap = 0
    out: list[tuple[int, int]] = list(other)
    for t in sorted(by_t):
        pairs = list(by_t[t])
        changed = True
        while changed:
            changed = False
            for a in range(len(pairs)):
                i1, j1 = pairs[a]
                for b in range(a + 1, len(pairs)):
                    i2, j2 = pairs[b]
                    d11, d22 = dist(i1, j1), dist(i2, j2)
                    d12, d21 = dist(i1, j2), dist(i2, j1)
                    if (
                        d12 <= max_um
                        and d21 <= max_um
                        and d12 + d21 + 1e-9 < d11 + d22
                    ):
                        pairs[a] = (i1, j2)
                        pairs[b] = (i2, j1)
                        nswap += 1
                        changed = True
                        i1, j1 = pairs[a]
        out.extend(pairs)
    return out, nswap


def link_tracks_staged(
    pred: np.ndarray,
    *,
    use_velocity: bool = True,
) -> tuple[list[tuple[int, int]], dict[str, float]]:
    """First-collapse links, leftover primaries, isolated mix, residual–residual.

    Halo residual↔primary stole tracks. Isolated residual (farther than
    NMS from every primary) may meet an unmatched primary; leftover
    residual→residual recovers second-collapse nuclei that never mix.
    """
    if len(pred) == 0:
        return [], {}
    pri = np.where(pred[:, 4] == 1)[0] if pred.shape[1] >= 5 else np.arange(len(pred))
    e_loc, meta = link_tracks_residual(pred[pri], use_velocity=use_velocity)
    edges = [(int(pri[i]), int(pri[j])) for i, j in e_loc]
    has_n = {i for i, _j in edges}
    has_p = {j for _i, j in edges}
    leftover_um = float(_PHI ** 5)
    frames = sorted({int(t) for t in pred[:, 0]})
    sc = np.array([SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM])

    extra = 0
    for t0, t1 in zip(frames, frames[1:]):
        if t1 != t0 + 1:
            continue
        src = [int(i) for i in _indices_at(pred, t0, 1) if int(i) not in has_n]
        dst = [int(i) for i in _indices_at(pred, t1, 1) if int(i) not in has_p]
        for i, j in _hungarian_gate(pred, src, dst, sc, leftover_um):
            edges.append((i, j))
            has_n.add(i)
            has_p.add(j)
            extra += 1

    iso_um = float(NMS_UM)

    def _iso(res, pri):
        if not res:
            return []
        if len(pri) == 0:
            return list(res)
        dmin = np.sqrt(
            (
                ((pred[res][:, 1:4][:, None, :] - pred[pri][None, :, 1:4]) * sc)
                ** 2
            ).sum(axis=2)
        ).min(axis=1)
        return [res[k] for k in range(len(res)) if float(dmin[k]) > iso_um]

    # Isolated residual ↔ unmatched primary before residual–residual
    # consumes those leftovers. Halo residual↔primary stole tracks.
    mix_n = 0
    for t0, t1 in zip(frames, frames[1:]):
        if t1 != t0 + 1:
            continue
        pri0 = [int(i) for i in _indices_at(pred, t0, 1) if int(i) not in has_n]
        pri1 = [int(i) for i in _indices_at(pred, t1, 1) if int(i) not in has_p]
        res0 = [int(i) for i in _indices_at(pred, t0, 2) if int(i) not in has_n]
        res1 = [int(i) for i in _indices_at(pred, t1, 2) if int(i) not in has_p]
        iso0 = _iso(res0, _indices_at(pred, t0, 1))
        iso1 = _iso(res1, _indices_at(pred, t1, 1))
        for i, j in _hungarian_gate(pred, iso0, pri1, sc, leftover_um):
            edges.append((i, j))
            has_n.add(i)
            has_p.add(j)
            mix_n += 1
        pri0 = [i for i in pri0 if i not in has_n]
        for i, j in _hungarian_gate(pred, pri0, iso1, sc, leftover_um):
            edges.append((i, j))
            has_n.add(i)
            has_p.add(j)
            mix_n += 1

    iso_n = 0
    res_n = 0
    for t0, t1 in zip(frames, frames[1:]):
        if t1 != t0 + 1:
            continue
        res0 = [int(i) for i in _indices_at(pred, t0, 2) if int(i) not in has_n]
        res1 = [int(i) for i in _indices_at(pred, t1, 2) if int(i) not in has_p]
        src_iso = _iso(res0, _indices_at(pred, t0, 1))
        dst_iso = _iso(res1, _indices_at(pred, t1, 1))
        for i, j in _hungarian_gate(pred, src_iso, dst_iso, sc, leftover_um):
            edges.append((i, j))
            has_n.add(i)
            has_p.add(j)
            iso_n += 1
        src_r = [i for i in res0 if i not in has_n]
        dst_r = [i for i in res1 if i not in has_p]
        for i, j in _hungarian_gate(pred, src_r, dst_r, sc, leftover_um):
            edges.append((i, j))
            has_n.add(i)
            has_p.add(j)
            res_n += 1

    # Unused destinations still inside φ⁴ of some source: 223 dense Jaccard
    # FN were 1-1 steals (parent already had an out-edge). Fill dest-only.
    # Gate is first-pass φ⁴ so this is not a long-shotgun.
    fill_n = 0
    fill_um = float(LINK_UM)
    for t0, t1 in zip(frames, frames[1:]):
        if t1 != t0 + 1:
            continue
        src = [int(i) for i in _indices_at(pred, t0)]
        dst = [int(i) for i in _indices_at(pred, t1) if int(i) not in has_p]
        for i, j in _hungarian_gate(pred, src, dst, sc, fill_um):
            if j in has_p:
                continue
            edges.append((i, j))
            has_n.add(i)
            has_p.add(j)
            fill_n += 1

    edges, nswap = _pairwise_swap_edges(pred, edges, sc, float(LINK_UM))

    meta = dict(meta)
    meta["n_leftover_edges"] = extra
    meta["leftover_um"] = leftover_um
    meta["n_isolated_residual_edges"] = iso_n
    meta["n_residual_residual_edges"] = res_n
    meta["n_isolated_mix_edges"] = mix_n
    meta["n_unmatched_dest_fill"] = fill_n
    meta["dest_fill_um"] = fill_um
    meta["n_pairwise_swaps"] = nswap
    return edges, meta


def _greedy_match_frame(
    gxyz: np.ndarray,
    pxyz: np.ndarray,
    sc: np.ndarray,
    max_um: float,
    used_p: set[int] | None = None,
) -> dict[int, int]:
    """Greedy nearest unused pred for each GT row. Local indices."""
    if len(gxyz) == 0 or len(pxyz) == 0:
        return {}
    d = np.sqrt((((gxyz[:, None, :] - pxyz[None, :, :]) * sc) ** 2).sum(axis=2))
    used: set[int] = set() if used_p is None else set(used_p)
    out: dict[int, int] = {}
    for row in range(len(gxyz)):
        order = [int(k) for k in np.argsort(d[row]) if int(k) not in used]
        if not order:
            continue
        j = order[0]
        if d[row, j] <= max_um:
            used.add(j)
            out[row] = j
    return out


def _hungarian_match_frame(
    gxyz: np.ndarray,
    pxyz: np.ndarray,
    sc: np.ndarray,
    max_um: float,
) -> dict[int, int]:
    """1-1 nearest unused pred (same as node recall). Local indices."""
    if len(gxyz) == 0 or len(pxyz) == 0:
        return {}
    from scipy.optimize import linear_sum_assignment

    d = np.sqrt((((gxyz[:, None, :] - pxyz[None, :, :]) * sc) ** 2).sum(axis=2))
    cap = d.copy()
    cap[d > max_um] = max_um * 4.0
    ri, cj = linear_sum_assignment(cap)
    out: dict[int, int] = {}
    for r, c in zip(ri, cj):
        if d[r, c] <= max_um:
            out[int(r)] = int(c)
    return out


def _assign_frame(
    pred_tzyx: np.ndarray,
    gi: list[int],
    gxyz: np.ndarray,
    t: int,
    sc: np.ndarray,
    max_um: float,
    g2p: dict[int, int],
    *,
    primary_only: bool,
    fill_residual: bool,
) -> int:
    """Match unmatched GT in one frame. Returns residual-fill count."""
    n_fill = 0
    taken = set(g2p[g] for g in gi if g in g2p)
    leftover_rows = [r for r, g in enumerate(gi) if g not in g2p]
    if not leftover_rows:
        return 0
    g_left = gxyz[leftover_rows]
    if primary_only and pred_tzyx.shape[1] >= 5:
        pri = [
            i
            for i, tt in enumerate(pred_tzyx[:, 0])
            if int(tt) == t and int(pred_tzyx[i, 4]) == 1 and i not in taken
        ]
        hit = _hungarian_match_frame(
            g_left, pred_tzyx[pri][:, 1:4] if pri else np.zeros((0, 3)), sc, max_um
        )
        for row, j in hit.items():
            g2p[gi[leftover_rows[row]]] = int(pri[j])
            taken.add(int(pri[j]))
        if fill_residual:
            still = [r for r, g in enumerate(gi) if g not in g2p]
            if still:
                res = [
                    i
                    for i, tt in enumerate(pred_tzyx[:, 0])
                    if int(tt) == t and int(pred_tzyx[i, 4]) == 2 and i not in taken
                ]
                if res:
                    fill = _hungarian_match_frame(
                        gxyz[still], pred_tzyx[res][:, 1:4], sc, max_um
                    )
                    for row, j in fill.items():
                        g2p[gi[still[row]]] = int(res[j])
                        n_fill += 1
    else:
        pi = [
            i
            for i, tt in enumerate(pred_tzyx[:, 0])
            if int(tt) == t and i not in taken
        ]
        if pi:
            hit = _hungarian_match_frame(g_left, pred_tzyx[pi][:, 1:4], sc, max_um)
            for row, j in hit.items():
                g2p[gi[leftover_rows[row]]] = int(pi[j])
    return n_fill


def lineage_recall(
    pred_tzyx: np.ndarray,
    pred_edges: list[tuple[int, int]],
    tracks: dict[str, Any],
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    max_um: float = MATCH_UM,
    primary_only: bool = False,
    fill_residual: bool = False,
    seed_um: float | None = None,
    follow: bool = False,
) -> dict[str, Any]:
    """Fraction of measured GEFF edges recovered by predicted lineage.

    primary_only maps GT onto first-collapse peaks. fill_residual then gives
    leftover GT the residual second-collapse peaks — residual never competes
    with a primary match (that mix stole tracks). seed_um locks a tighter
    match so a wider radius cannot remap a cell onto a farther primary.

    follow: the parent's predicted child (outgoing link) lands within
    max_um of the measured next cell. Independent greedy pairing of both
    ends was matching a closer ghost that was not the continuation.
    """
    gt = np.column_stack([tracks["t"].astype(float), tracks["xyz_vox"]])
    id_to_g = {int(n): i for i, n in enumerate(tracks["ids"])}
    sc = np.asarray(scale_zyx, dtype=float)
    g2p: dict[int, int] = {}
    n_fill = 0
    radii = [float(max_um)] if seed_um is None else [float(seed_um), float(max_um)]
    # unique, increasing
    seen_r: list[float] = []
    for r in radii:
        if not seen_r or r > seen_r[-1] + 1e-9:
            seen_r.append(r)
    for t in sorted({int(x) for x in gt[:, 0]}):
        gi = [i for i, tt in enumerate(gt[:, 0]) if int(tt) == t]
        if not gi:
            continue
        gxyz = gt[gi][:, 1:4]
        for rad in seen_r:
            n_fill += _assign_frame(
                pred_tzyx,
                gi,
                gxyz,
                t,
                sc,
                rad,
                g2p,
                primary_only=primary_only,
                fill_residual=fill_residual,
            )
    pred_set = {tuple(sorted(e)) for e in pred_edges}
    out_of = {int(i): int(j) for i, j in pred_edges}
    sc_gt = sc
    xyz = tracks["xyz_vox"]
    hit = 0
    n_e = 0
    used_dest: set[int] = set()
    for a, b in tracks["edges"]:
        if a not in id_to_g or b not in id_to_g:
            continue
        n_e += 1
        ga, gb = id_to_g[a], id_to_g[b]
        if follow:
            pa = g2p.get(ga)
            if pa is None or pa not in out_of:
                continue
            dest = out_of[pa]
            if dest in used_dest:
                continue
            d = float(np.linalg.norm((pred_tzyx[dest, 1:4] - xyz[gb]) * sc_gt))
            if d <= max_um:
                used_dest.add(dest)
                hit += 1
            continue
        pa, pb = g2p.get(ga), g2p.get(gb)
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
        "n_gt_residual_fill": n_fill,
        "follow": bool(follow),
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


def edge_jaccard_official(
    pred_tzyx: np.ndarray,
    pred_edges: list[tuple[int, int]],
    tracks: dict[str, Any],
    *,
    scale_zyx: tuple[float, float, float] = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM),
    max_um: float = MATCH_UM,
    node_penalty_a: float = 0.1,
) -> dict[str, Any]:
    """Kaggle edge Jaccard (TP/FP/FN) + node-count penalty.

    Recall of GEFF edges is not this number. Unmatched predicted edges are
    ignored; a predicted edge that steals an annotated in/out is FP.
    Adjusted Jaccard taxes extra nodes vs estimated_number_of_nodes.
    """
    from scipy.optimize import linear_sum_assignment

    sc = np.asarray(scale_zyx, dtype=float)
    gt_t = tracks["t"].astype(int)
    gt_xyz = tracks["xyz_vox"]
    g2p: dict[int, int] = {}
    p2g: dict[int, int] = {}
    for t in sorted(set(gt_t.tolist())):
        gi = np.where(gt_t == t)[0]
        pi = np.where(pred_tzyx[:, 0].astype(int) == t)[0]
        if len(gi) == 0 or len(pi) == 0:
            continue
        d = np.sqrt((((gt_xyz[gi][:, None, :] - pred_tzyx[pi][:, 1:4][None, :, :]) * sc) ** 2).sum(axis=2))
        cap = d.copy()
        cap[d > max_um] = max_um * 4.0
        ri, cj = linear_sum_assignment(cap)
        for r, c in zip(ri, cj):
            if d[r, c] <= max_um:
                g2p[int(gi[r])] = int(pi[c])
                p2g[int(pi[c])] = int(gi[r])
    id_to_g = {int(n): i for i, n in enumerate(tracks["ids"])}
    gt_pair: set[tuple[int, int]] = set()
    gt_out: dict[int, set[int]] = {}
    gt_in: dict[int, set[int]] = {}
    for a, b in tracks["edges"]:
        if a not in id_to_g or b not in id_to_g:
            continue
        ga, gb = id_to_g[a], id_to_g[b]
        gt_pair.add((ga, gb))
        gt_out.setdefault(ga, set()).add(gb)
        gt_in.setdefault(gb, set()).add(ga)
    tp = 0
    fp = 0
    matched: set[tuple[int, int]] = set()
    for i, j in pred_edges:
        ga, gb = p2g.get(int(i)), p2g.get(int(j))
        if ga is None or gb is None:
            continue
        if (ga, gb) in gt_pair:
            tp += 1
            matched.add((ga, gb))
            continue
        if gb in gt_in and ga not in gt_in.get(gb, set()):
            fp += 1
        elif ga in gt_out and gb not in gt_out.get(ga, set()):
            fp += 1
    fn = len(gt_pair) - len(matched)
    den = tp + fp + fn
    jac = float(tp / den) if den else None
    t_pred = int(len(pred_tzyx))
    t_true = int(tracks["meta"].get("estimated_nodes") or 0)
    adj = (
        max(0.0, jac * (1.0 - node_penalty_a * (t_pred - t_true) / t_true))
        if jac is not None and t_true
        else jac
    )
    return {
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "n_gt_edges": int(len(gt_pair)),
        "edge_jaccard": jac,
        "adjusted_edge_jaccard": adj,
        "n_pred_nodes": t_pred,
        "n_est_true": t_true if t_true else None,
        "n_gt_nodes_matched": len(g2p),
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
        pred = load_or_detect(args.dataset, volp)
        gt_tzyx = np.column_stack(
            [tracks["t"].astype(float), tracks["xyz_vox"]]
        )
        xyz4 = pred[:, :4]
        hit = match_centroids(xyz4, gt_tzyx)
        hit12 = match_centroids(xyz4, gt_tzyx, max_um=12.0)
        if pred.shape[1] >= 5:
            prod = product_detections(pred)
            edges, link_meta = link_tracks_staged(prod)
            link_meta["n_primary"] = int((pred[:, 4] == 1).sum())
            link_meta["n_residual_peaks"] = int((pred[:, 4] == 2).sum())
            link_meta["n_product_nodes"] = int(len(prod))
            link_meta["n_isolated_residual_nodes"] = int(
                (prod[:, 4] == 2).sum() if prod.shape[1] >= 5 else 0
            )
        else:
            prod = pred
            edges, link_meta = link_tracks_residual(xyz4)
        # Observer find-recall stays on all peaks. Lineage / Jaccard are
        # the product graph (primary + isolated residual) so halo cannot
        # steal the 7 µm match off a linked first-collapse track.
        lin_pri = lineage_recall(prod, edges, tracks, primary_only=True)
        lin_pair = lineage_recall(
            prod, edges, tracks, primary_only=True, fill_residual=True
        )
        lin = lineage_recall(
            prod, edges, tracks, primary_only=True, fill_residual=True, follow=True
        )
        lin12 = lineage_recall(
            prod,
            edges,
            tracks,
            max_um=12.0,
            seed_um=MATCH_UM,
            primary_only=True,
            fill_residual=True,
            follow=True,
        )
        link_meta["lineage_primary_only_7um"] = lin_pri["edge_recall"]
        link_meta["lineage_pair_7um"] = lin_pair["edge_recall"]
        link_meta["n_gt_residual_fill"] = lin["n_gt_residual_fill"]
        jac = edge_jaccard_official(prod, edges, tracks)
        if pred.shape[1] >= 5:
            pri = pred[pred[:, 4] == 1]
            e_pri, _ = link_tracks_residual(pri)
            jac_pri = edge_jaccard_official(pri, e_pri, tracks)
        else:
            jac_pri = jac
        link_meta["edge_jaccard"] = jac["edge_jaccard"]
        link_meta["adjusted_edge_jaccard"] = jac["adjusted_edge_jaccard"]
        link_meta["edge_jaccard_primary"] = jac_pri["edge_jaccard"]
        link_meta["adjusted_edge_jaccard_primary"] = jac_pri["adjusted_edge_jaccard"]
        link_meta["jaccard_tp_fp_fn"] = [jac["tp"], jac["fp"], jac["fn"]]
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
            "n_pred_edges": len(edges),
            "n_product_nodes": int(len(prod)),
            "lineage_edge_recall_7um": lin["edge_recall"],
            "lineage_edge_recall_12um": lin12["edge_recall"],
            "n_gt_edges_matched": lin["n_edge_matched"],
            "n_gt_edges": lin["n_gt_edges"],
            "edge_jaccard": jac["edge_jaccard"],
            "adjusted_edge_jaccard": jac["adjusted_edge_jaccard"],
            "jaccard_tp_fp_fn": [jac["tp"], jac["fp"], jac["fn"]],
            "nms_um": NMS_UM,
            "link_um": link_meta.get("residual_gate_um", LINK_UM),
            "link_meta": link_meta,
            "note": (
                "Centroid = half-max first moment. Gate = median+φ·MAD. "
                "NMS = φ³ µm. Residual second collapse fills leftover "
                "brightness (7 µm find). Product graph = first collapse + "
                "isolated residual; in-shell residual folds into the primary "
                "centroid. Jaccard is the Kaggle reference metric."
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
