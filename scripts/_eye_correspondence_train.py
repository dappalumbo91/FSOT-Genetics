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


def _paint_frac(stem: str) -> float:
    """MAD-gate paint fraction on a mid frame (low = photon observer is blind)."""
    from scipy.ndimage import gaussian_filter
    from biohub_3d import _mad_threshold, open_volume

    volp = TRAIN / f"{stem}.zarr"
    tracks = read_geff(BIOHUB_ROOT / "train" / f"{stem}.geff")
    t = int(np.median(tracks["t"])) if len(tracks["t"]) else 0
    arr = open_volume(volp)
    vol = np.asarray(arr[t]).astype(np.float32)
    sm = gaussian_filter(vol[:, ::2, ::2], sigma=(0.4, 0.8, 0.8))
    thr = _mad_threshold(sm)
    return float((sm >= thr).mean())


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


def _gt_masks(logits: torch.Tensor, gt_zyx: np.ndarray) -> tuple[torch.Tensor, torch.Tensor]:
    """Nucleus ball (φ³) vs steal-shell (φ³..φ⁵) in downsampled voxels."""
    z, y, x = logits.shape
    zz, yy, xx = torch.meshgrid(
        torch.arange(z, device=logits.device),
        torch.arange(y, device=logits.device),
        torch.arange(x, device=logits.device),
        indexing="ij",
    )
    inner = float(NMS_UM) / VOX_UM
    outer = float(_PHI ** 5) / VOX_UM
    ball = torch.zeros_like(logits, dtype=torch.bool)
    shell = torch.zeros_like(logits, dtype=torch.bool)
    for z0, y0, x0 in gt_zyx:
        d2 = (zz.float() - z0) ** 2 + (yy.float() - y0) ** 2 + (xx.float() - x0) ** 2
        ball |= d2 <= inner * inner
        shell |= (d2 > inner * inner) & (d2 <= outer * outer)
    return ball, shell


def _ball_loss(logits: torch.Tensor, ball: torch.Tensor) -> torch.Tensor:
    if not bool(ball.any()):
        return logits.new_zeros(())
    return F.binary_cross_entropy_with_logits(logits[ball], torch.ones_like(logits[ball]))


def _annulus_loss(logits: torch.Tensor, shell: torch.Tensor) -> torch.Tensor:
    """Target 0 in the leftover shell around each GT (don't steal the track)."""
    if not bool(shell.any()):
        return logits.new_zeros(())
    return F.binary_cross_entropy_with_logits(logits[shell], torch.zeros_like(logits[shell]))


def main() -> int:
    _wire_cellmot()
    from predict_unet_transformer import load_model
    from train_unet_transformer import compute_detection_loss

    OUT.mkdir(parents=True, exist_ok=True)
    stems = _stems()
    rng = random.Random(1)
    lowcon = "lowcon" in sys.argv[1:]
    frac_path = OUT / "paint_frac.json"
    weights = {s: 1 for s in stems}
    if lowcon:
        if frac_path.exists():
            fr = json.loads(frac_path.read_text(encoding="utf-8"))
        else:
            fr = {}
            print("  scanning paint fractions", flush=True)
            for i, s in enumerate(stems):
                try:
                    fr[s] = _paint_frac(s)
                except Exception:
                    fr[s] = 1.0
                if i % 20 == 0:
                    print(f"  frac {i+1}/{len(stems)} {s} {fr[s]:.4f}", flush=True)
            OUT.mkdir(parents=True, exist_ok=True)
            frac_path.write_text(json.dumps(fr, indent=2), encoding="utf-8")
        cut = 1.0 / (_PHI ** 3)
        n_low = 0
        for s in stems:
            if float(fr.get(s, 1.0)) < cut:
                weights[s] = 4
                n_low += 1
        print(f"  low-contrast stems {n_low}/{len(stems)} (frac < 1/phi^3={cut:.3f})", flush=True)
    bag = []
    for s, w in weights.items():
        bag.extend([s] * int(w))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, window_size, _ds = load_model(FT_WEIGHTS, device)
    ckpt = OUT / "edge_predictor_correspondence.pth"
    if ckpt.exists():
        state = torch.load(ckpt, map_location=device, weights_only=True)
        model.load_state_dict(state)
        print(f"  resume {ckpt}", flush=True)
    model.train()
    opt = torch.optim.Adam(model.parameters(), lr=5e-5)
    n_iters = 800
    for a in sys.argv[1:]:
        if a.isdigit():
            n_iters = int(a)
            break
    log = []
    dest_name = (
        "edge_predictor_correspondence_lowcon.pth"
        if lowcon
        else "edge_predictor_correspondence_ball.pth"
    )
    print(
        f"correspondence train device={device} stems={len(stems)} "
        f"iters={n_iters} lowcon={lowcon}",
        flush=True,
    )
    for it in range(n_iters):
        stem = rng.choice(bag)
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
        ball_loss = imgs.new_zeros(())
        ann_loss = imgs.new_zeros(())
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
            ball, shell = _gt_masks(det_logits[f][0, 0], coords)
            if not lowcon:
                ball_loss = ball_loss + _ball_loss(det_logits[f][0, 0], ball)
            ann_loss = ann_loss + _annulus_loss(det_logits[f][0, 0], shell)
        loss = det_loss + ball_loss + ann_loss
        if float(loss.detach()) == 0.0:
            continue
        loss.backward()
        opt.step()
        if it % 50 == 0 or it + 1 == n_iters:
            rec = {
                "iter": it,
                "loss": float(loss.detach()),
                "det": float(det_loss.detach()),
                "ball": float(ball_loss.detach()),
                "annulus": float(ann_loss.detach()),
                "stem": stem,
            }
            log.append(rec)
            print(json.dumps(rec), flush=True)
        del imgs, _out, det_logits, loss
    dest = OUT / dest_name
    torch.save(model.state_dict(), dest)
    (OUT / "train_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
    print(f"wrote {dest}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
