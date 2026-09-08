#!/usr/bin/env python3
"""Eye / relay: U-Net paints a field; FSOT reports the center. Not shipped.

The net is lab apparatus (cones). It must not emit annotator xyz.
Activation field × measured photons is the scene the relay reads:
MAD+φ gate, φ³ NMS, leftover second collapse, half-max centroid,
leftover-yield linker. Identity intensity is native photons.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import zarr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from biohub_3d import (  # noqa: E402
    BIOHUB_ROOT,
    MATCH_UM,
    NMS_UM,
    SCALE_X_UM,
    SCALE_Y_UM,
    SCALE_Z_UM,
    TRAIN,
    detect_peaks_frame,
    edge_jaccard_official,
    lineage_recall,
    link_tracks_staged,
    match_centroids,
    product_detections,
    read_geff,
)

UNET_ROOT = Path(r"C:\Users\damia\biohub-fsot-unet")
CELLMOT = UNET_ROOT / "vendor" / "kaggle-cell-tracking-competition"
FT_WEIGHTS = Path(
    r"D:\Kaggle_Biohub_Data\cellmot\cellmot-ft-detector-biohub\edge_predictor_best.pth"
)
EYE_CACHE = BIOHUB_ROOT / "_fsot_eye_cache"
SCALE = (SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM)


def _wire_cellmot() -> None:
    for p in (CELLMOT / "src", CELLMOT / "scripts", UNET_ROOT / "vendor"):
        s = str(p)
        if p.exists() and s not in sys.path:
            sys.path.insert(0, s)


def eye_field_path(dataset: str) -> Path:
    return EYE_CACHE / f"{dataset}_ft_sigmoid.npy"


def paint_eye(dataset: str, volp: Path) -> np.ndarray:
    """Downsampled sigmoid field (T, Z, Yd, Xd). Cached on the game drive."""
    cache = eye_field_path(dataset)
    if cache.exists():
        field = np.load(cache)
        print(f"  eye cache {cache} {field.shape}", flush=True)
        return field

    _wire_cellmot()
    import torch
    import torch.nn.functional as F
    from predict_unet_transformer import _load_frame, load_model

    zg = zarr.open_group(str(volp), mode="r")
    zarr_arr = zg["0"]
    qs = zg.attrs["image_statistics"]["quantiles"]
    q_low = float(qs["0.001"])
    q_high = float(qs["0.999"])
    T, Z, Y, X = [int(x) for x in zarr_arr.shape]
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, window_size, downsample = load_model(FT_WEIGHTS, device)
    dz, dy, dx = downsample
    target_shape = [Z // dz, Y // dy, X // dx]
    W = int(window_size)
    fields = np.zeros((T, target_shape[0], target_shape[1], target_shape[2]), dtype=np.float16)
    seen: set[int] = set()
    stride = max(W - 1, 1)
    starts = list(range(0, T - W + 1, stride))
    if not starts or starts[-1] + W < T:
        last = max(T - W, 0)
        if not starts or last != starts[-1]:
            starts.append(last)
    print(f"  eye encode device={device} windows={len(starts)} T={T}", flush=True)
    with torch.no_grad():
        for k, ws in enumerate(starts):
            if k % 10 == 0 or k + 1 == len(starts):
                print(f"  eye window {ws} ({k + 1}/{len(starts)})", flush=True)
            idxs = list(range(ws, ws + W))
            imgs = torch.stack(
                [_load_frame(zarr_arr, t, target_shape, downsample) for t in idxs]
            )
            imgs = ((imgs - q_low) / (q_high - q_low + 1e-6)).clamp(0.0)
            imgs = imgs.unsqueeze(0).to(device)
            _out, det_logits = model.encode(imgs)
            for f_idx, t in enumerate(idxs):
                if t in seen:
                    continue
                prob = torch.sigmoid(det_logits[f_idx][0, 0]).detach().cpu().numpy()
                fields[int(t)] = prob.astype(np.float16)
                seen.add(int(t))
            del imgs, _out, det_logits
    EYE_CACHE.mkdir(parents=True, exist_ok=True)
    np.save(cache, fields)
    print(f"  wrote {cache}", flush=True)
    return fields


def relay_pred_path(dataset: str) -> Path:
    return EYE_CACHE / f"{dataset}_ft_relay.npy"


def detect_eye_relay(dataset: str, volp: Path, fields: np.ndarray) -> np.ndarray:
    """FSOT observer on (photons × upsampled eye field)."""
    from scipy.ndimage import zoom

    zg = zarr.open_group(str(volp), mode="r")
    zarr_arr = zg["0"]
    T, Z, Y, X = [int(x) for x in zarr_arr.shape]
    _, Zd, Yd, Xd = fields.shape
    zy, zx = Y / Yd, X / Xd
    rows = []
    n_t = T
    for t in range(T):
        if t % 10 == 0 or t + 1 == n_t:
            print(f"  relay detect t={t} ({t + 1}/{n_t})", flush=True)
        photons = np.asarray(zarr_arr[t]).astype(np.float32)
        eye = np.asarray(fields[t], dtype=np.float32)
        if eye.shape != photons.shape:
            eye = zoom(eye, (Z / Zd, zy, zx), order=1).astype(np.float32)
            if eye.shape != photons.shape:
                # trilinear fallback to exact grid
                import torch
                import torch.nn.functional as F

                tsr = torch.from_numpy(eye)[None, None]
                eye = (
                    F.interpolate(
                        tsr, size=photons.shape, mode="trilinear", align_corners=False
                    )
                    .numpy()[0, 0]
                    .astype(np.float32)
                )
        scene = photons * eye
        p, tier, _inten = detect_peaks_frame(scene, xy_stride=1)
        if len(p) == 0:
            continue
        inten = np.empty(len(p), dtype=np.float64)
        for i, c in enumerate(p):
            zc, yc, xc = [int(round(v)) for v in c]
            zc = int(np.clip(zc, 0, photons.shape[0] - 1))
            yc = int(np.clip(yc, 0, photons.shape[1] - 1))
            xc = int(np.clip(xc, 0, photons.shape[2] - 1))
            inten[i] = float(photons[zc, yc, xc])
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


def score(ds: str, pred: np.ndarray, tracks: dict, tag: str) -> dict:
    gt = np.column_stack([tracks["t"].astype(float), tracks["xyz_vox"]])
    hit = match_centroids(pred[:, :4], gt, max_um=MATCH_UM)
    hit12 = match_centroids(pred[:, :4], gt, max_um=12.0)
    prod = product_detections(pred)
    hit_p = match_centroids(prod[:, :4], gt, max_um=MATCH_UM)
    edges, meta = link_tracks_staged(prod)
    jac = edge_jaccard_official(prod, edges, tracks)
    lin = lineage_recall(
        prod, edges, tracks, primary_only=True, fill_residual=True, follow=True
    )
    return {
        "tag": tag,
        "n_pred": int(len(pred)),
        "n_prod": int(len(prod)),
        "find7": hit["recall"],
        "find12": hit12["recall"],
        "prod_find7": hit_p["recall"],
        "match_median_um": hit["match_median_um"],
        "follow7": lin["edge_recall"],
        "jac": jac["edge_jaccard"],
        "adj": jac["adjusted_edge_jaccard"],
        "tp_fp_fn": [jac["tp"], jac["fp"], jac["fn"]],
        "n_yield": meta.get("n_leftover_yield"),
        "n_edges": len(edges),
        "nms_um": NMS_UM,
    }


def main() -> int:
    ds = sys.argv[1] if len(sys.argv) > 1 else "6bba_09961292"
    volp = TRAIN / f"{ds}.zarr"
    tracks = read_geff(BIOHUB_ROOT / "train" / f"{ds}.geff")
    pred_path = relay_pred_path(ds)
    if pred_path.exists():
        pred = np.load(pred_path)
        print(f"  relay pred cache {pred_path} n={len(pred)}", flush=True)
    else:
        fields = paint_eye(ds, volp)
        pred = detect_eye_relay(ds, volp, fields)
        EYE_CACHE.mkdir(parents=True, exist_ok=True)
        np.save(pred_path, pred)
        print(f"  wrote {pred_path} n={len(pred)}", flush=True)
    reports = [score(ds, pred, tracks, "eye_product_isolated")]
    if pred.shape[1] >= 5:
        pri = pred[pred[:, 4] == 1]
        reports.append(score(ds, pri, tracks, "eye_primary_only"))
    print(json.dumps({"dataset": ds, "reports": reports}, indent=2, default=str), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
