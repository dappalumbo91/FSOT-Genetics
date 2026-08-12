#!/usr/bin/env python3
"""FSOT-Genetics medical/research front door — one call, regime-selected structure.

Usage
-----
  python scripts/fsot_predict.py --seq MQIFV... --pdb-out model.pdb
  python scripts/fsot_predict.py --uniprot P0CG47 --pfam PF00240
  python scripts/fsot_predict.py --id 1UBQ

Regime policy (zero free params; templates/MSA = data):
  1. If a real homolog structure is found → template + packing fuse (MSA if deep)
  2. Else if deep MSA → bulk MSA-augmented F15
  3. Else → pure single-sequence bulk

Always reports free_parameters=0, structure_mode, confidence when MSA exists.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import predict_ca_coords, write_ca_pdb, clean_sequence  # noqa: E402
from msa_pipeline import build_msa_features, conservation_confidence  # noqa: E402
from msa_template_fuse import fuse_predict, select_regime  # noqa: E402

# optional template search (network)
try:
    from run_rcsb_template_holdout import best_template  # noqa: E402
except Exception:  # noqa: BLE001
    best_template = None  # type: ignore

BUILTINS = {
    "1UBQ": {
        "seq": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "pfam": "PF00240",
        "exclude_pdb": "1UBQ",
    },
    "1CRN": {
        "seq": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
        "pfam": "PF00304",
        "exclude_pdb": "1CRN",
    },
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seq", default=None)
    ap.add_argument("--id", default=None, help="Builtin short id")
    ap.add_argument("--uniprot", default=None)
    ap.add_argument("--pfam", default=None)
    ap.add_argument("--exclude-pdb", default=None, help="Self-PDB to exclude from templates")
    ap.add_argument("--no-template", action="store_true")
    ap.add_argument("--no-msa", action="store_true")
    ap.add_argument("--rounds", type=int, default=16)
    ap.add_argument("--pdb-out", default=None)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args(argv)

    builtin = BUILTINS.get((args.id or "").upper())
    if args.seq:
        seq = clean_sequence(args.seq)
        default_pfam = args.pfam
        default_excl = args.exclude_pdb or "XXXX"
    elif builtin:
        seq = clean_sequence(builtin["seq"])
        default_pfam = args.pfam or builtin.get("pfam")
        default_excl = args.exclude_pdb or builtin.get("exclude_pdb") or "XXXX"
    else:
        print("Provide --seq or --id (1UBQ/1CRN)", file=sys.stderr)
        return 2

    feat = None
    if not args.no_msa:
        try:
            feat = build_msa_features(seq, pfam=default_pfam, uniprot=args.uniprot)
        except Exception as exc:  # noqa: BLE001
            feat = None
            msa_err = str(exc)
        else:
            msa_err = None
    else:
        msa_err = "disabled"

    tmpl = None
    if not args.no_template and best_template is not None:
        try:
            tmpl = best_template(seq, default_excl)
        except Exception:
            tmpl = None

    regime = select_regime(tmpl is not None, feat if feat and feat.n_seqs else None)
    report: dict = {
        "sequence": seq,
        "length": len(seq),
        "free_parameters": 0,
        "formula": "S=K(T1+T2+T3)",
        "authority_pin": "D1D38A",
        "deploy_regime": regime,
        "msa": feat.summary() if feat and feat.n_seqs else None,
        "msa_error": msa_err,
    }

    if tmpl is not None:
        fused = fuse_predict(
            seq,
            tmpl["model"],
            feat if feat and feat.depth_ok else None,
            tertiary_contacts=tmpl.get("tertiary_contacts"),
        )
        X = fused["ca_coords"]
        report.update(
            {
                "structure_mode": fused.get("regime", "template"),
                "engine": fused.get("engine"),
                "template_pdb": tmpl.get("pdb_id"),
                "template_identity": tmpl.get("identity"),
                "template_coverage": tmpl.get("coverage"),
                "mean_confidence": float(np.mean(fused["confidence"]))
                if fused.get("confidence") is not None and len(fused["confidence"])
                else None,
            }
        )
        conf = fused.get("confidence")
    else:
        mode = "msa" if regime == "bulk_msa" else "single"
        pred = predict_ca_coords(
            seq,
            rounds=args.rounds,
            mode=mode,
            msa_features=feat if mode == "msa" else None,
            canonicalize_chirality=True,
        )
        X = pred["ca_coords"]
        conf = conservation_confidence(feat) if feat and feat.n_seqs else None
        report.update(
            {
                "structure_mode": pred.get("structure_mode", mode),
                "engine": pred.get("engine"),
                "rg_A": pred.get("rg_A"),
                "predict_ms": pred.get("predict_ms"),
                "mean_confidence": float(np.mean(conf)) if conf is not None else None,
            }
        )

    report["ca_coords"] = X.tolist()
    if conf is not None:
        report["per_residue_confidence"] = [float(c) for c in conf]

    print(
        f"FSOT predict  n={len(seq)}  regime={report.get('structure_mode')}  "
        f"deploy={regime}  free_parameters=0"
    )
    if report.get("template_pdb"):
        print(
            f"  template={report['template_pdb']}  "
            f"id={report.get('template_identity'):.2f}  "
            f"cov={report.get('template_coverage'):.2f}"
        )
    if report.get("msa"):
        print(
            f"  msa backend={report['msa'].get('backend')}  "
            f"n={report['msa'].get('n_seqs')}  neff={report['msa'].get('neff'):.1f}"
        )
    if report.get("mean_confidence") is not None:
        print(f"  mean confidence={report['mean_confidence']:.3f}")

    if args.pdb_out:
        write_ca_pdb(Path(args.pdb_out), seq, X, name="FSOT")
        print(f"  wrote {args.pdb_out}")
    if args.json_out:
        # drop huge coords optionally? keep for medical audit
        Path(args.json_out).write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
