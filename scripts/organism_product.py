#!/usr/bin/env python3
"""Cell → gene → FSOT product Cα on the live organism graphs.

Named cells that already carry residual mass get the proteins that
sit on them (FlyBase / WormBase / UniProt). Fold only with a measured
homolog. Same leftover / elongated-homolog gates as the walking set.

  python scripts/organism_product.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_predict import main as predict_main  # noqa: E402

FASTA = Path(r"D:\FlyWire_Connectome\male_cns\next_proteins.fasta")
FASTA2 = Path(r"D:\FlyWire_Connectome\male_cns\open_proteins.fasta")
OUT_D = Path(r"D:\FlyWire_Connectome\male_cns\product")
OUT_GIT = ROOT / "data" / "organism_product_join.json"
WALKING = ROOT / "data" / "fly_walking_product.json"

# Measured sit-on: not invented. JO/GABA/mechanosensory from the live boots.
NEW = [
    {
        "symbol": "ChAT",
        "uniprot": "P07668",
        "organism": "Drosophila melanogaster",
        "sits_on": "cholinergic neurons (most DNs, many sensory); ACh is the default excitatory edge",
        "graph": "male-cns / FlyWire",
        "function": "acetylcholine synthesis",
    },
    {
        "symbol": "unc-25",
        "uniprot": "G5EDB7",
        "organism": "Caenorhabditis elegans",
        "sits_on": "GABAergic motor neurons DD/VD/RME (inhibitory residual on the worm graph)",
        "graph": "C. elegans SI5",
        "function": "GABA synthesis (Gad1 homolog)",
    },
    {
        "symbol": "mec-4",
        "uniprot": "P24612",
        "organism": "Caenorhabditis elegans",
        "sits_on": "ALML/ALMR touch neurons (worm hop-0 peak)",
        "graph": "C. elegans SI5",
        "function": "DEG/ENaC mechanotransduction",
    },
    {
        "symbol": "myo-3",
        "uniprot": "P12844",
        "organism": "Caenorhabditis elegans",
        "sits_on": "body-wall muscles (worm hop-4 effector mass)",
        "graph": "C. elegans SI5",
        "function": "muscle myosin heavy chain",
    },
    {
        "symbol": "VGlut",
        "uniprot": "Q9VQC0",
        "organism": "Drosophila melanogaster",
        "sits_on": "glutamatergic motor neurons (vnc_motor NMJ)",
        "graph": "male-cns",
        "function": "vesicular glutamate transporter (FlyBase FBgn0031424; TrEMBL)",
    },
    {
        "symbol": "Mhc",
        "uniprot": "P05661",
        "organism": "Drosophila melanogaster",
        "sits_on": "muscle (effector of vnc_motor; not a CNS cell)",
        "graph": "fly body, not in Male CNS dump",
        "function": "muscle myosin heavy chain",
    },
]


def _parse_fasta(path: Path) -> dict[str, str]:
    acc: dict[str, str] = {}
    cur = None
    seq: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur is not None:
                acc[cur] = "".join(seq)
            parts = line[1:].split("|")
            cur = parts[1] if len(parts) > 1 else line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if cur is not None:
        acc[cur] = "".join(seq)
    return acc


def _fold_one(g: dict, seq: str) -> dict:
    pdb_out = OUT_D / f"{g['symbol']}_{g['uniprot']}.pdb"
    json_out = OUT_D / f"{g['symbol']}_{g['uniprot']}.json"
    print(f"== {g['symbol']} {g['uniprot']} n={len(seq)} {g['organism']}", flush=True)
    rc = predict_main(
        [
            "--seq",
            seq,
            "--uniprot",
            g["uniprot"],
            "--pdb-out",
            str(pdb_out),
            "--json-out",
            str(json_out),
        ]
    )
    rec = {
        **g,
        "length": len(seq),
        "predict_rc": rc,
        "free_parameters": 0,
    }
    if json_out.exists():
        full = json.loads(json_out.read_text(encoding="utf-8"))
        rec.update(
            {
                "structure_mode": full.get("structure_mode"),
                "deploy_regime": full.get("deploy_regime"),
                "template_pdb": full.get("template_pdb"),
                "template_identity": full.get("template_identity"),
                "template_coverage": full.get("template_coverage"),
                "mean_confidence": full.get("mean_confidence"),
                "rg_target_A": full.get("rg_target_A"),
                "engine": full.get("engine"),
                "pdb": str(pdb_out) if pdb_out.exists() else None,
            }
        )
    print(
        f"  {g['symbol']} mode={rec.get('structure_mode')} "
        f"tmpl={rec.get('template_pdb')} id={rec.get('template_identity')}",
        flush=True,
    )
    return rec


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args(argv)
    if not FASTA.exists():
        raise SystemExit(f"missing {FASTA}")
    seqs = _parse_fasta(FASTA)
    if FASTA2.exists():
        seqs.update(_parse_fasta(FASTA2))
    OUT_D.mkdir(parents=True, exist_ok=True)
    want = {s.lower() for s in args.only} if args.only else None
    prior = []
    if OUT_GIT.exists() and want:
        try:
            prior = json.loads(OUT_GIT.read_text(encoding="utf-8")).get("new_folds", [])
        except Exception:
            prior = []
    prior_map = {r["symbol"]: r for r in prior}
    new_folds = []
    for g in NEW:
        if want and g["symbol"].lower() not in want:
            if g["symbol"] in prior_map:
                new_folds.append(prior_map[g["symbol"]])
            continue
        seq = seqs.get(g["uniprot"])
        if not seq:
            print(f"  skip {g['symbol']}: no sequence", flush=True)
            continue
        new_folds.append(_fold_one(g, seq))
    walking = []
    if WALKING.exists():
        walking = json.loads(WALKING.read_text(encoding="utf-8")).get("genes", [])
    join = {
        "product": "named cell on measured graph → UniProt sequence → FSOT product Cα",
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3); residual Biochemistry on synapses; product Cα when a homolog exists",
        "cells": [
            {
                "cell": "JO-* / Johnston organ",
                "graph": "male-cns / FlyWire",
                "boot": "JO seed → DNg29 → vnc_motor",
                "proteins": ["iav", "nan"],
            },
            {
                "cell": "mechanosensory / vnc_sensory",
                "graph": "male-cns",
                "boot": "vnc_sensory seed → vnc_motor",
                "proteins": ["nompC"],
            },
            {
                "cell": "GABAergic (fly) / DD VD RME (worm)",
                "graph": "both",
                "boot": "GABA edges inhibitory",
                "proteins": ["Gad1", "unc-25"],
            },
            {
                "cell": "cholinergic DN / motor",
                "graph": "male-cns",
                "boot": "ACh default excitatory",
                "proteins": ["ChAT"],
            },
            {
                "cell": "ALML/ALMR",
                "graph": "C. elegans SI5",
                "boot": "hop 0 peak on sensory seed",
                "proteins": ["mec-4"],
            },
            {
                "cell": "body-wall muscle",
                "graph": "C. elegans SI5",
                "boot": "hop 4 effector mass",
                "proteins": ["myo-3"],
            },
            {
                "cell": "vnc_motor (glutamatergic NMJ)",
                "graph": "male-cns",
                "boot": "vnc_sensory / JO → vnc_motor",
                "proteins": ["VGlut"],
            },
            {
                "cell": "leg/body muscle (not in CNS EM)",
                "graph": "fly body",
                "boot": "effector of vnc_motor",
                "proteins": ["Mhc"],
            },
        ],
        "walking_folds": walking,
        "new_folds": new_folds,
    }
    OUT_GIT.write_text(json.dumps(join, indent=2), encoding="utf-8")
    print(f"  wrote {OUT_GIT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
