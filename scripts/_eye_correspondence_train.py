#!/usr/bin/env python3
"""Fine-tune the eye so it agrees with the FSOT relay. Lab apparatus.

The FT detector already lights leftover GT (sigmoid ~1) and silences
ghosts (~0.002). We had ranked fill by photons and put the ghosts back.
This train step adds the steal-shell: around each labeled cell, do not
paint a second nucleus (NMS..leftover annulus). Hold out the bench
videos. Claim path stays 0 free params.
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
import zarr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from biohub_3d import (  # noqa: E402
    BIOHUB_ROOT,
    NMS_UM,
    TRAIN,
    _PHI,
    read_geff,
)
from _eye_relay import CELLMOT, FT_WEIGHTS, UNET_ROOT, _wire_cellmot  # noqa: E402

HOLD = {"6bba_09961292", "44b6_0113de3b"}
OUT = Path(r"D:\Kaggle_Biohub_Data\_fsot_eye_cache\correspondence")
DOWNSAMPLE = (1, 4, 4)
VOX_UM = 1.625  # after XY stride 4


def _stems() -> list[str]:
    return sorted(
        p.name[:-5]
        for p in TRAIN.glob("*.zarr")
        if (TRAIN / f"{p.name[:-5]}.geff").exists() and p.name[:-5] not in HOLD
    )


def _frame(zarr_arr, t: int, q_low: float, q_high: float) -> torch.Tensor:
    dz, dy, dx = DOWNSAMPLE
    raw = np.asarray(zarr_arr[int(t), ::dz, ::dy, ::dx], dtype=np.float32)
    x = torch.from_numpy(raw)
    return ((x - q_low) / (q_high - q_low + 1e-6)).clamp(0.0)


def _annulus_loss(
    logits: torch.Tensor,
    gt_zyx: np.ndarray,
) -> torch.Tensor:
    """Target 0 in the leftover shell around each GT (don't steal the track)."""
    # logits: (Z,Y,X)
    z, y, x = logits.shape
    zz, yy, xx = torch.meshgrid(
        torch.arange(z, device=logits.device),
        torch.arange(y, device=logits.device),
        torch.arange(x, device=logits.device),
        indexing="ij",
    )
    inner = float(NMS_UM) / VOX_UM
    outer = float(_PHI ** 5) / VOX_UM
    mask = torch.zeros_like(logits, dtype=torch.bool)
    for z0, y0, x0 in gt_zyx:
        d2 = (zz.float() - z0) ** 2 + (yy.float() - y0) ** 2 + (xx.float() - x0) ** 2
        mask |= (d2 > inner * inner) & (d2 <= outer * outer)
    if not bool(mask.any()):
        return logits.new_zeros(())
    return F.binary_cross_entropy_with_logits(logits[mask], torch.zeros_like(logits[mask]))


def main() -> int:
    _wire_cellmot()
    from predict_unet_transformer import load_model
    from train_unet_transformer import compute_detection_loss

    OUT.mkdir(parents=True, exist_ok=True)
    stems = _stems()
    rng = random.Random(1)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, window_size, _ds = load_model(FT_WEIGHTS, device)
    model.train()
    opt = torch.optim.Adam(model.parameters(), lr=5e-5)
    n_iters = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    log = []
    print(f"correspondence train device={device} stems={len(stems)} iters={n_iters}", flush=True)
    for it in range(n_iters):
        stem = rng.choice(stems)
        volp = TRAIN / f"{stem}.zarr"
        geff = BIOHUB_ROOT / "train" / f"{stem}.geff"
        try:
            tracks = read_geff(geff)
            zg = zarr.open_group(str(volp), mode="r")
            arr = zg["0"]
            qs = zg.attrs["image_statistics"]["quantiles"]
            q_low, q_high = float(qs["0.001"]), float(qs["0.999"])
        except Exception as exc:
            print(f"  skip {stem}: {exc}", flush=True)
            continue
        tcol = tracks["t"].astype(int)
        frames = sorted({int(t) for t in tcol.tolist()})
        frames = [t for t in frames if t + 1 in set(frames)]
        if not frames:
            continue
        t0 = int(rng.choice(frames))
        t1 = t0 + 1
        imgs = torch.stack([_frame(arr, t0, q_low, q_high), _frame(arr, t1, q_low, q_high)])
        imgs = imgs.unsqueeze(0).to(device)
        opt.zero_grad(set_to_none=True)
        _out, det_logits = model.encode(imgs)
        det_loss = imgs.new_zeros(())
        ann_loss = imgs.new_zeros(())
        spatial = det_logits[0].shape[2:]
        for f, t in enumerate((t0, t1)):
            xyz = tracks["xyz_vox"][tcol == t]
            if len(xyz) == 0:
                continue
            coords = np.stack(
                [xyz[:, 0], xyz[:, 1] / 4.0, xyz[:, 2] / 4.0], axis=1
            ).astype(np.float32)
            coords_t = torch.from_numpy(coords).unsqueeze(0).to(device)
            mask = torch.ones((1, len(coords)), dtype=torch.bool, device=device)
            det_loss = det_loss + compute_detection_loss(
                det_logits[f], coords_t, mask, neg_weight=0.1
            )
            ann_loss = ann_loss + _annulus_loss(det_logits[f][0, 0], coords)
        loss = det_loss + ann_loss
        if float(loss.detach()) == 0.0:
            continue
        loss.backward()
        opt.step()
        if it % 50 == 0 or it + 1 == n_iters:
            rec = {
                "iter": it,
                "loss": float(loss.detach()),
                "det": float(det_loss.detach()),
                "annulus": float(ann_loss.detach()),
                "stem": stem,
            }
            log.append(rec)
            print(json.dumps(rec), flush=True)
        del imgs, _out, det_logits, loss
    dest = OUT / "edge_predictor_correspondence.pth"
    torch.save(model.state_dict(), dest)
    (OUT / "train_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
    print(f"wrote {dest}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
