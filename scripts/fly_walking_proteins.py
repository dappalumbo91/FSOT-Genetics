#!/usr/bin/env python3
"""Fold walking-program fly proteins on the FSOT product path.

iav / nan / nompC sit on JO and mechanosensory cells that seeded the
Male CNS boot. Gad1 is the GABA enzyme already used as inhibitory sign.

Measured homolog → product Cα. No homolog → no_measured_map (Rg +
secondary only). 0 free parameters. Sequences from UniProt on D:.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_predict import main as predict_main  # noqa: E402

FASTA = Path(r"D:\FlyWire_Connectome\male_cns\fly_walking_proteins.fasta")
OUT_D = Path(r"D:\FlyWire_Connectome\male_cns\product")
OUT_GIT = ROOT / "data" / "fly_walking_product.json"

GENES = [
    {"symbol": "Gad1", "uniprot": "P20228", "sits_on": "GABAergic neurons"},
    {"symbol": "nan", "uniprot": "Q9VUD5", "sits_on": "Johnston organ"},
    {"symbol": "iav", "uniprot": "Q9W3W0", "sits_on": "Johnston organ"},
    {"symbol": "nompC", "uniprot": "Q7KIQ2", "sits_on": "mechanosensory transduction"},
]


def _parse_fasta(path: Path) -> dict[str, str]:
    acc: dict[str, str] = {}
    cur = None
    seq: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur is not None:
                acc[cur] = "".join(seq)
            # >sp|P20228|DCE_DROME ...
            parts = line[1:].split("|")
            cur = parts[1] if len(parts) > 1 else line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if cur is not None:
        acc[cur] = "".join(seq)
    return acc


def main() -> int:
    if not FASTA.exists():
        raise SystemExit(f"missing {FASTA}")
    seqs = _parse_fasta(FASTA)
    OUT_D.mkdir(parents=True, exist_ok=True)
    rows = []
    for g in GENES:
        acc = g["uniprot"]
        seq = seqs.get(acc)
        if not seq:
            print(f"  skip {g['symbol']}: no sequence for {acc}", flush=True)
            continue
        pdb_out = OUT_D / f"{g['symbol']}_{acc}.pdb"
        json_out = OUT_D / f"{g['symbol']}_{acc}.json"
        print(f"== {g['symbol']} {acc} n={len(seq)}", flush=True)
        rc = predict_main(
            [
                "--seq",
                seq,
                "--uniprot",
                acc,
                "--pdb-out",
                str(pdb_out),
                "--json-out",
                str(json_out),
            ]
        )
        rec: dict = {
            "symbol": g["symbol"],
            "uniprot": acc,
            "sits_on": g["sits_on"],
            "length": len(seq),
            "predict_rc": rc,
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
                    "free_parameters": full.get("free_parameters", 0),
                }
            )
        rows.append(rec)
        print(
            f"  {g['symbol']} mode={rec.get('structure_mode')} "
            f"tmpl={rec.get('template_pdb')} id={rec.get('template_identity')}",
            flush=True,
        )
    report = {
        "product": "walking-program fly proteins on FSOT product Cα",
        "pin": "D1D38A",
        "free_parameters": 0,
        "authority": "UniProt sequences; measured homologs via RCSB; no MDS fold",
        "genes": rows,
    }
    OUT_GIT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {OUT_GIT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
