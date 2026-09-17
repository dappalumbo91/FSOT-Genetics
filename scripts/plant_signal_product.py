#!/usr/bin/env python3
"""Product Cα for proteins that sit on the Arabidopsis signaling graph.

Same product rule. These are the receptors/TFs/transporters the residual
hops actually land on. Not a connectome.

  python scripts/plant_signal_product.py
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_predict import main as predict_main  # noqa: E402
from plant_product import _fetch_fasta, _parse_fasta  # noqa: E402

OUT_D = Path(r"D:\FlyWire_Connectome\plants\product")
FASTA = Path(r"D:\FlyWire_Connectome\plants\arabidopsis_signal.fasta")
OUT_GIT = ROOT / "data" / "plant_signal_product.json"

PANEL = [
    {"symbol": "PHOT1", "uniprot": "Q2V2M9", "sits_on": "blue-light phototropin", "system": "light receptor"},
    {"symbol": "PHYB", "uniprot": "P14713", "sits_on": "red-light phytochrome", "system": "light receptor"},
    {"symbol": "CRY1", "uniprot": "Q43125", "sits_on": "blue-light cryptochrome", "system": "light receptor"},
    {"symbol": "PIF3", "uniprot": "Q495N3", "sits_on": "phytochrome-interacting bHLH", "system": "light transcription"},
    {"symbol": "OST1", "uniprot": "Q940H6", "sits_on": "ABA SnRK2 kinase", "system": "ABA / stomata"},
    {"symbol": "ABI1", "uniprot": "P49597", "sits_on": "PP2C ABA phosphatase", "system": "ABA / stomata"},
    {"symbol": "PIN1", "uniprot": "Q9C6B8", "sits_on": "auxin efflux carrier", "system": "auxin transport"},
    {"symbol": "TIR1", "uniprot": "Q570C0", "sits_on": "auxin F-box receptor", "system": "auxin receptor"},
]


def main() -> int:
    FASTA.parent.mkdir(parents=True, exist_ok=True)
    OUT_D.mkdir(parents=True, exist_ok=True)
    seqs: dict[str, str] = {}
    if FASTA.exists():
        seqs.update(_parse_fasta(FASTA))
    chunks: list[str] = []
    rows = []
    for g in PANEL:
        acc = g["uniprot"]
        seq = seqs.get(acc)
        if not seq:
            print(f"  fetch {g['symbol']} {acc}", flush=True)
            seq = _fetch_fasta(acc)
            seqs[acc] = seq
        chunks.append(f">{g['symbol']}|{acc}|Arabidopsis thaliana\n")
        for i in range(0, len(seq), 60):
            chunks.append(seq[i : i + 60] + "\n")
        pdb_out = OUT_D / f"{g['symbol']}_{acc}.pdb"
        json_out = OUT_D / f"{g['symbol']}_{acc}.json"
        rec = {
            **g,
            "organism": "Arabidopsis thaliana",
            "taxid": 3702,
            "length": len(seq),
            "free_parameters": 0,
        }
        if json_out.exists():
            full = json.loads(json_out.read_text(encoding="utf-8"))
            rec.update(
                {
                    "predict_rc": 0,
                    "structure_mode": full.get("structure_mode"),
                    "template_pdb": full.get("template_pdb"),
                    "template_identity": full.get("template_identity"),
                    "template_coverage": full.get("template_coverage"),
                    "pdb": str(pdb_out) if pdb_out.exists() else None,
                }
            )
            print(f"  skip existing {json_out.name}", flush=True)
        else:
            print(f"== {g['symbol']} {acc} n={len(seq)}", flush=True)
            rc = predict_main(
                ["--seq", seq, "--uniprot", acc, "--pdb-out", str(pdb_out), "--json-out", str(json_out)]
            )
            rec["predict_rc"] = rc
            if json_out.exists():
                full = json.loads(json_out.read_text(encoding="utf-8"))
                rec.update(
                    {
                        "structure_mode": full.get("structure_mode"),
                        "template_pdb": full.get("template_pdb"),
                        "template_identity": full.get("template_identity"),
                        "template_coverage": full.get("template_coverage"),
                        "pdb": str(pdb_out) if pdb_out.exists() else None,
                    }
                )
            print(
                f"  {g['symbol']} mode={rec.get('structure_mode')} "
                f"tmpl={rec.get('template_pdb')} id={rec.get('template_identity')}",
                flush=True,
            )
        rows.append(rec)
    FASTA.write_text("".join(chunks), encoding="utf-8")
    join = {
        "product": "Arabidopsis signaling proteins on the measured PPI graph → FSOT product Cα",
        "pin": "D1D38A",
        "free_parameters": 0,
        "not": "Not a plant connectome. These proteins sit on IntAct physical edges.",
        "organism": "Arabidopsis thaliana",
        "taxid": 3702,
        "genes": rows,
        "n_folded": sum(1 for r in rows if r.get("predict_rc") == 0),
    }
    OUT_GIT.write_text(json.dumps(join, indent=2), encoding="utf-8")
    print(f"  wrote {OUT_GIT} n_folded={join['n_folded']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
