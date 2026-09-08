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
    _centroid_blob,
    detect_cache_path,
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


def fold_eye_into_native(
    native_prod: np.ndarray,
    eye: np.ndarray,
    *,
    iso_um: float | None = None,
) -> np.ndarray:
    """Move native centers with in-shell eye leftover. Do not add eye nodes.

    37k dense eye peaks sit inside native NMS; 26k isolated extra nodes
    paid the 0.1 tax and stole proxy tracks. Fold is residual brightness
    on an already-measured blob, same as in-shell photon residual.
    """
    sc = np.array([SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM])
    if iso_um is None:
        iso_um = float(NMS_UM)
    out = native_prod.copy()
    if len(out) == 0 or len(eye) == 0:
        return out
    for t in sorted({int(x) for x in out[:, 0]}):
        po = np.where(out[:, 0].astype(int) == t)[0]
        ey = eye[eye[:, 0].astype(int) == t]
        if len(po) == 0 or len(ey) == 0:
            continue
        d = np.sqrt(
            (((out[po][:, 1:4][:, None, :] - ey[None, :, 1:4]) * sc) ** 2).sum(axis=2)
        )
        for k, row in enumerate(po):
            near = d[k] <= iso_um
            if not np.any(near):
                continue
            xyz = np.vstack([out[row, 1:4], ey[near][:, 1:4]])
            if out.shape[1] >= 6 and ey.shape[1] >= 6:
                w = np.concatenate([[out[row, 5]], ey[near, 5]])
            else:
                w = np.ones(len(xyz))
            wsum = float(w.sum())
            if wsum <= 0:
                continue
            out[row, 1:4] = (xyz * w[:, None]).sum(axis=0) / wsum
    return out


def fill_isolated_eye(
    native_prod: np.ndarray,
    eye: np.ndarray,
    n_est: int,
    *,
    iso_um: float | None = None,
) -> np.ndarray:
    """Add isolated eye leftover up to estimated_number_of_nodes.

    Rank by measured photon intensity. Proxy is already over the estimate
    so it receives none; dense leftover (~3k) can take the 7 µm misses
    the eye painted without the 15k node tax.
    """
    sc = np.array([SCALE_Z_UM, SCALE_Y_UM, SCALE_X_UM])
    if iso_um is None:
        iso_um = float(NMS_UM)
    budget = int(n_est) - int(len(native_prod))
    if budget <= 0 or len(eye) == 0:
        return native_prod
    iso_rows: list[np.ndarray] = []
    for t in sorted({int(x) for x in eye[:, 0]}):
        ey = eye[eye[:, 0].astype(int) == t]
        nat = native_prod[native_prod[:, 0].astype(int) == t]
        if len(ey) == 0:
            continue
        if len(nat) == 0:
            iso_rows.append(ey)
            continue
        dmin = np.sqrt(
            (((ey[:, 1:4][:, None, :] - nat[None, :, 1:4]) * sc) ** 2).sum(axis=2)
        ).min(axis=1)
        keep = ey[dmin > iso_um]
        if len(keep):
            iso_rows.append(keep)
    if not iso_rows:
        return native_prod
    iso = np.vstack(iso_rows)
    if iso.shape[1] >= 6:
        order = np.argsort(-iso[:, 5])
        iso = iso[order]
    iso = iso[:budget].copy()
    if iso.shape[1] >= 5:
        iso[:, 4] = 2
    return np.vstack([native_prod, iso])


def reccentroid_native_on_field(
    native_prod: np.ndarray,
    fields: np.ndarray,
    volp: Path,
) -> np.ndarray:
    """Half-max first moment of photons×eye in the native NMS patch.

    Discrete in-shell eye peaks yanked centers toward neighbors (proxy
    1.00→0.94). The field modulates leftover brightness of the same blob.
    """
    from scipy.ndimage import zoom

    zg = zarr.open_group(str(volp), mode="r")
    zarr_arr = zg["0"]
    T, Z, Y, X = [int(x) for x in zarr_arr.shape]
    _, Zd, Yd, Xd = fields.shape
    zy, zx = Y / Yd, X / Xd
    out = native_prod.copy()
    for t in sorted({int(x) for x in out[:, 0]}):
        if t % 10 == 0:
            print(f"  field-centroid t={t}", flush=True)
        photons = np.asarray(zarr_arr[int(t)]).astype(np.float32)
        eye = np.asarray(fields[int(t)], dtype=np.float32)
        if eye.shape != photons.shape:
            eye = zoom(eye, (Z / Zd, zy, zx), order=1).astype(np.float32)
        scene = photons * eye
        rows = np.where(out[:, 0].astype(int) == int(t))[0]
        for i in rows:
            out[i, 1:4] = _centroid_blob(scene, out[i, 1:4], SCALE, float(NMS_UM))
            zc, yc, xc = [int(round(c)) for c in out[i, 1:4]]
            zc = int(np.clip(zc, 0, photons.shape[0] - 1))
            yc = int(np.clip(yc, 0, photons.shape[1] - 1))
            xc = int(np.clip(xc, 0, photons.shape[2] - 1))
            if out.shape[1] >= 6:
                out[i, 5] = float(photons[zc, yc, xc])
    return out


def score(
    ds: str,
    pred: np.ndarray,
    tracks: dict,
    tag: str,
    *,
    as_product: bool = True,
) -> dict:
    gt = np.column_stack([tracks["t"].astype(float), tracks["xyz_vox"]])
    hit = match_centroids(pred[:, :4], gt, max_um=MATCH_UM)
    hit12 = match_centroids(pred[:, :4], gt, max_um=12.0)
    prod = product_detections(pred) if as_product and pred.shape[1] >= 5 else pred
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
    mode = sys.argv[2] if len(sys.argv) > 2 else "fold"
    volp = TRAIN / f"{ds}.zarr"
    tracks = read_geff(BIOHUB_ROOT / "train" / f"{ds}.geff")
    if mode == "fill":
        native = np.load(detect_cache_path(ds))
        nat_prod = product_detections(native)
        eye = np.load(relay_pred_path(ds))
        n_est = int(tracks["meta"].get("estimated_nodes") or len(nat_prod))
        filled = fill_isolated_eye(nat_prod, eye, n_est)
        reports = [
            score(ds, nat_prod, tracks, "native_product", as_product=False),
            score(ds, filled, tracks, "native_plus_eye_fill", as_product=False),
        ]
        print(
            json.dumps(
                {
                    "dataset": ds,
                    "mode": mode,
                    "n_est": n_est,
                    "n_native": int(len(nat_prod)),
                    "n_filled": int(len(filled)),
                    "reports": reports,
                },
                indent=2,
                default=str,
            ),
            flush=True,
        )
        return 0
    if mode in ("fold", "field"):
        native = np.load(detect_cache_path(ds))
        nat_prod = product_detections(native)
        if mode == "fold":
            eye = np.load(relay_pred_path(ds))
            moved = fold_eye_into_native(nat_prod, eye)
            tag = "native_fold_eye_nms"
        else:
            fields = paint_eye(ds, volp)
            moved = reccentroid_native_on_field(nat_prod, fields, volp)
            tag = "native_field_centroid"
        reports = [
            score(ds, nat_prod, tracks, "native_product", as_product=False),
            score(ds, moved, tracks, tag, as_product=False),
        ]
        print(json.dumps({"dataset": ds, "mode": mode, "reports": reports}, indent=2, default=str), flush=True)
        return 0
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
    print(json.dumps({"dataset": ds, "mode": mode, "reports": reports}, indent=2, default=str), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
