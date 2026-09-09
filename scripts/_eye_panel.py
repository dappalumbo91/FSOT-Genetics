#!/usr/bin/env python3
"""Score correspondence fill on a diverse GT panel. Caches stay on D:\\.

Does not download. Train already has 199 labeled GEFFs. This is variety
of *known* information: native detect + steal-shell eye fill vs native.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from biohub_3d import (  # noqa: E402
    BIOHUB_ROOT,
    TRAIN,
    detect_cache_path,
    load_or_detect,
    product_detections,
    read_geff,
)
from _eye_relay import (  # noqa: E402
    CORR_WEIGHTS,
    EYE_CACHE,
    fill_isolated_eye,
    paint_eye,
    detect_eye_relay,
    score,
)

# Held-out bench already scored: 44b6_0113de3b, 6bba_09961292
PANEL = [
    "44b6_18ced818",  # hardest density (T_true ~78k, 100 GT)
    "44b6_0b24845f",  # sparse annot / hard5
    "6bba_05db0fb1",  # dense-ish + divisions (~1.2k GT, T_true ~70k)
    "6bba_48816121",  # mid T_true, divisions
]


def one(ds: str) -> dict:
    volp = TRAIN / f"{ds}.zarr"
    geff = BIOHUB_ROOT / "train" / f"{ds}.geff"
    print(f"\n=== {ds} ===", flush=True)
    tracks = read_geff(geff)
    native = load_or_detect(ds, volp)
    nat_prod = product_detections(native)
    fields = paint_eye(
        ds, volp, weights=CORR_WEIGHTS, cache_name=f"{ds}_corr_sigmoid.npy"
    )
    pred_c = EYE_CACHE / f"{ds}_corr_relay.npy"
    if pred_c.exists():
        eye = np.load(pred_c)
        print(f"  corr relay cache n={len(eye)}", flush=True)
    else:
        eye = detect_eye_relay(ds, volp, fields)
        EYE_CACHE.mkdir(parents=True, exist_ok=True)
        np.save(pred_c, eye)
        print(f"  wrote {pred_c} n={len(eye)}", flush=True)
    n_est = int(tracks["meta"].get("estimated_nodes") or len(nat_prod))
    filled = fill_isolated_eye(nat_prod, eye, n_est, fields=fields)
    native_s = score(ds, nat_prod, tracks, "native_product", as_product=False)
    fill_s = score(ds, filled, tracks, "corr_fill", as_product=False)
    row = {
        "dataset": ds,
        "n_gt": tracks["n_nodes"],
        "n_gt_edges": tracks["n_edges"],
        "n_divisions": tracks["n_divisions"],
        "n_est": n_est,
        "n_native": int(len(nat_prod)),
        "n_filled": int(len(filled)),
        "native": native_s,
        "fill": fill_s,
    }
    print(json.dumps(row, indent=2, default=str), flush=True)
    return row


def main() -> int:
    names = sys.argv[1:] if len(sys.argv) > 1 else PANEL
    rows = [one(ds) for ds in names]
    out = EYE_CACHE / "variety_panel.json"
    payload = {"n": len(rows), "holdout_note": "eye saw these except 6bba_09961292 / 44b6_0113de3b", "rows": rows}
    out.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(f"wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
