#!/usr/bin/env python3
"""Dual-mode fold: pure single-sequence vs MSA-augmented F15.

Reports both modes side-by-side so claims stay honest:
  - single_sequence  = published de-novo path (default engine behaviour)
  - msa_augmented    = same arithmetic + optional coevolution channel

MSA sources (first hit):
  --msa FILE.sto|.a3m|.fasta
  --pfam PF#####          (InterPro public full alignment)
  --uniprot ACCESSION     (resolve Pfam then fetch)
  env FSOT_JACKHMMER_DB / FSOT_HHBLITS_DB for local search tools

Zero free parameters: evolutionary signal is data input; amplitude is F09-family.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import predict_ca_coords, clean_sequence  # noqa: E402
from msa_pipeline import (  # noqa: E402
    build_msa_features,
    conservation_confidence,
    tools_available,
)

# Classic short benchmarks (sequence only; optional RMSD if --pdb)
DEFAULTS = {
    "1UBQ": {
        "seq": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "pfam": "PF00240",
        "name": "ubiquitin",
    },
    "1CRN": {
        "seq": "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN",
        "pfam": "PF00304",
        "name": "crambin",
    },
}


def _kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(((Pc @ R.T - Qc) ** 2).sum(axis=1).mean()))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seq", default=None, help="Amino-acid sequence")
    ap.add_argument("--id", default="1UBQ", help="Builtin id (1UBQ, 1CRN) if --seq omitted")
    ap.add_argument("--pfam", default=None, help="Pfam accession for MSA")
    ap.add_argument("--uniprot", default=None, help="UniProt accession → Pfam resolve")
    ap.add_argument("--msa", default=None, help="Precomputed MSA path (sto/a3m/fasta)")
    ap.add_argument("--rounds", type=int, default=16)
    ap.add_argument("--no-network", action="store_true", help="Skip Pfam fetch; file/local only")
    ap.add_argument("--json-out", default=None, help="Write dual-mode report JSON")
    args = ap.parse_args(argv)

    if args.seq:
        seq = clean_sequence(args.seq)
        pfam = args.pfam
        name = "custom"
    else:
        hit = DEFAULTS.get(args.id.upper())
        if not hit:
            print(f"unknown --id {args.id}; pass --seq", file=sys.stderr)
            return 2
        seq = hit["seq"]
        pfam = args.pfam or hit["pfam"]
        name = hit["name"]

    print("FSOT dual-mode fold")
    print(f"  target={name}  n={len(seq)}  rounds={args.rounds}")
    print(f"  tools={tools_available()}")

    # ── single-sequence (default claim path) ──────────────────────────────
    r_single = predict_ca_coords(seq, rounds=args.rounds, mode="single")
    print(
        f"  single_sequence  rg={r_single['rg_A']:.2f} A  "
        f"mode={r_single.get('structure_mode')}  ms={r_single['predict_ms']:.1f}"
    )

    # ── MSA features ──────────────────────────────────────────────────────
    feat = None
    if args.msa or (not args.no_network and (pfam or args.uniprot)):
        try:
            feat = build_msa_features(
                seq,
                msa_path=args.msa,
                pfam=None if args.no_network else pfam,
                uniprot=None if args.no_network else args.uniprot,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"  MSA obtain failed: {exc}")
            feat = None

    report: dict = {
        "name": name,
        "length": len(seq),
        "single": {
            "structure_mode": r_single.get("structure_mode"),
            "rg_A": r_single["rg_A"],
            "rg_target_fsot_A": r_single["rg_target_fsot_A"],
            "predict_ms": r_single["predict_ms"],
            "free_parameters": r_single["free_parameters"],
            "engine": r_single["engine"],
        },
        "msa": None,
        "msa_augmented": None,
        "tools": tools_available(),
    }

    if feat is not None and feat.n_seqs > 0:
        conf = conservation_confidence(feat)
        print(
            f"  MSA  backend={feat.backend}  detail={feat.detail}  "
            f"n_seqs={feat.n_seqs}  neff={feat.neff:.1f}  "
            f"mean_cons={feat.conservation.mean():.3f}  "
            f"mean_conf={conf.mean():.3f}  depth_ok={feat.depth_ok}"
        )
        report["msa"] = feat.summary()
        report["msa"]["mean_evo_confidence"] = float(conf.mean())

        r_msa = predict_ca_coords(
            seq, rounds=args.rounds, mode="msa", msa_features=feat
        )
        print(
            f"  msa_augmented    rg={r_msa['rg_A']:.2f} A  "
            f"mode={r_msa.get('structure_mode')}  "
            f"msa={r_msa.get('msa')}  ms={r_msa['predict_ms']:.1f}"
        )
        report["msa_augmented"] = {
            "structure_mode": r_msa.get("structure_mode"),
            "rg_A": r_msa["rg_A"],
            "rg_target_fsot_A": r_msa["rg_target_fsot_A"],
            "predict_ms": r_msa["predict_ms"],
            "free_parameters": r_msa["free_parameters"],
            "msa": r_msa.get("msa"),
            "engine": r_msa["engine"],
        }
        # Coord drift between modes (not RMSD-to-native; topology shift probe)
        drift = float(
            np.linalg.norm(
                (r_msa["ca_coords"] - r_msa["ca_coords"].mean(0))
                - (r_single["ca_coords"] - r_single["ca_coords"].mean(0)),
                axis=1,
            ).mean()
        )
        report["mean_ca_drift_single_vs_msa_A"] = drift
        print(f"  mean |ΔCα| (centered) single↔msa = {drift:.3f} A")
    else:
        print("  MSA  not available — msa_augmented skipped (single-sequence only)")

    if args.json_out:
        out = Path(args.json_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"  wrote {out}")

    print("  free_parameters=0  (MSA is data input, not trained weights)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
