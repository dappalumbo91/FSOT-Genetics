#!/usr/bin/env python3
"""Read the FlyWire adult Drosophila connectome (measured authority).

Annotations live on the game drive, not in git:

  D:\\FlyWire_Connectome

  python scripts/fly_connectome.py --inventory

Same pin as protein product and Biohub: measured soma/anchor coordinates
and measured types/transmitters. Residual does not invent synapses.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402

FLY_ROOT = Path(r"D:\FlyWire_Connectome")
ANN_URL = (
    "https://raw.githubusercontent.com/flyconnectome/flywire_annotations/"
    "main/supplemental_files/Supplemental_file1_neuron_annotations.tsv"
)
ANN_NAME = "Supplemental_file1_neuron_annotations.tsv"
# FlyWire FAFB EM voxel (Schlegel et al. 2024)
SCALE_XY_NM = 4.0
SCALE_Z_NM = 40.0

_PHI = float(fc.PHI)
_R_BIO = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))


def ensure_annotations(dest_dir: Path = FLY_ROOT) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / ANN_NAME
    if path.exists() and path.stat().st_size > 1_000_000:
        return path
    print(f"  downloading {ANN_URL}", flush=True)
    urllib.request.urlretrieve(ANN_URL, path)
    print(f"  wrote {path} ({path.stat().st_size} bytes)", flush=True)
    return path


def _float(v: str) -> float | None:
    v = (v or "").strip()
    if not v or v.lower() in {"na", "nan", "none"}:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def read_annotations(path: Path) -> dict[str, Any]:
    """Measured neuron atlas: types, transmitters, soma xyz (voxels)."""
    import csv

    with path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    n = len(rows)
    super_c = Counter((r.get("super_class") or "").strip() or "unlabeled" for r in rows)
    flow = Counter((r.get("flow") or "").strip() or "unlabeled" for r in rows)
    nt = Counter((r.get("top_nt") or "").strip() or "unlabeled" for r in rows)
    side = Counter((r.get("side") or "").strip() or "unlabeled" for r in rows)
    dimorph = Counter((r.get("dimorphism") or "").strip() or "unlabeled" for r in rows)
    fru = Counter((r.get("fru_dsx") or "").strip() or "unlabeled" for r in rows)
    n_type = sum(1 for r in rows if (r.get("cell_type") or "").strip())
    n_soma = 0
    soma = []
    for r in rows:
        x, y, z = _float(r.get("soma_x")), _float(r.get("soma_y")), _float(r.get("soma_z"))
        if x is None or y is None or z is None:
            continue
        n_soma += 1
        soma.append((x * SCALE_XY_NM, y * SCALE_XY_NM, z * SCALE_Z_NM))
    span = None
    if soma:
        a = np.asarray(soma, dtype=np.float64)
        span = (a.max(axis=0) - a.min(axis=0)).tolist()
    return {
        "path": str(path),
        "n_neurons": n,
        "n_with_cell_type": n_type,
        "n_with_soma_xyz": n_soma,
        "super_class": dict(super_c.most_common(24)),
        "flow": dict(flow.most_common()),
        "top_nt": dict(nt.most_common()),
        "side": dict(side.most_common()),
        "dimorphism": dict(dimorph.most_common(12)),
        "fru_dsx": dict(fru.most_common(12)),
        "soma_span_nm": span,
        "voxel_nm": [SCALE_XY_NM, SCALE_XY_NM, SCALE_Z_NM],
        "authority": "FlyWire FAFB v783 + Schlegel et al. 2024 annotations",
        "free_parameters": 0,
        "phi": _PHI,
        "residual_Biochemistry": _R_BIO,
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", action="store_true")
    ap.add_argument("--root", default=str(FLY_ROOT))
    args = ap.parse_args(argv)
    root = Path(args.root)
    path = ensure_annotations(root)
    inv = read_annotations(path)
    print(json.dumps(inv, indent=2))
    out = ROOT / "data" / "fly_connectome_inventory.json"
    out.write_text(json.dumps(inv, indent=2), encoding="utf-8")
    print(f"  wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
