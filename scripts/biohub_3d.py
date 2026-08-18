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
