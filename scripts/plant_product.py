#!/usr/bin/env python3
"""Arabidopsis measured-homolog product Cα — not a plant connectome.

Plants have genomes, proteomes, and crystals. They do not have a fly-class
EM wiring diagram. Same product rule as animals: fold only with a measured
homolog. Same pin, 0 free parameters.

  python scripts/plant_product.py

Sequences on D:\\FlyWire_Connectome\\plants (not git).
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_predict import main as predict_main  # noqa: E402

OUT_D = Path(r"D:\FlyWire_Connectome\plants\product")
FASTA = Path(r"D:\FlyWire_Connectome\plants\arabidopsis_panel.fasta")
OUT_GIT = ROOT / "data" / "plant_product.json"

_UA = "FSOT-Genetics plant product (mailto:local)"

# Reviewed Arabidopsis thaliana. PDB xrefs on the UniProt page or a
# universal measured class (actin). Not invented.
PANEL = [
    {
        "symbol": "rbcL",
        "uniprot": "O03042",
        "sits_on": "chloroplast stroma; carbon fixation (Calvin cycle)",
        "system": "photosynthesis / carbon fixation",
        "pdb_on_entry": ["5IU0", "9MUR", "9N37"],
    },
    {
        "symbol": "psbA",
        "uniprot": "P83755",
        "sits_on": "photosystem II reaction center D1",
        "system": "photosynthesis (light reactions)",
        "pdb_on_entry": ["5MDX", "7OUI", "9LE7"],
    },
    {
        "symbol": "LHCB1.3",
        "uniprot": "P04778",
        "sits_on": "PSII light-harvesting complex II",
        "system": "photosynthesis (antenna)",
        "pdb_on_entry": ["5MDX", "7OUI", "9LE7"],
    },
    {
        "symbol": "GAPA1",
        "uniprot": "P25856",
        "sits_on": "chloroplast GAPDH (Calvin cycle)",
        "system": "photosynthesis / carbon fixation",
        "pdb_on_entry": ["3K2B", "3QV1", "3RVD"],
    },
    {
        "symbol": "ACT2",
        "uniprot": "Q96292",
        "sits_on": "cytoskeleton (universal actin)",
        "system": "cell shape / trafficking",
        "pdb_on_entry": [],
    },
    {
        "symbol": "CESA3",
        "uniprot": "Q941L0",
        "sits_on": "plasma membrane cellulose synthase",
        "system": "cell wall",
        "pdb_on_entry": ["7CK1", "7CK2", "7CK3"],
    },
]


def _fetch_fasta(acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=90) as fh:
        raw = fh.read().decode("utf-8")
    return "".join(ln.strip() for ln in raw.splitlines() if ln and not ln.startswith(">"))


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


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args(argv)
    FASTA.parent.mkdir(parents=True, exist_ok=True)
    OUT_D.mkdir(parents=True, exist_ok=True)
    want = {s.lower() for s in args.only} if args.only else None
    chunks: list[str] = []
    seqs: dict[str, str] = {}
    if FASTA.exists():
        seqs.update(_parse_fasta(FASTA))
    rows = []
    for g in PANEL:
        if want and g["symbol"].lower() not in want:
            continue
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
        print(f"== {g['symbol']} {acc} n={len(seq)} Arabidopsis thaliana", flush=True)
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
        rec = {
            **g,
            "organism": "Arabidopsis thaliana",
            "taxid": 3702,
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
        rows.append(rec)
    FASTA.write_text("".join(chunks), encoding="utf-8")
    join = {
        "product": "Arabidopsis measured homolog → FSOT product Cα",
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3); product Cα when a homolog exists",
        "not": (
            "Not a plant connectome. Arabidopsis has a genome and crystals. "
            "No fly-class EM wiring diagram. Predict proteins, not invented synapses."
        ),
        "organism": "Arabidopsis thaliana",
        "taxid": 3702,
        "proteome": "UP000006548",
        "genes": rows,
        "n_folded": sum(1 for r in rows if r.get("predict_rc") == 0),
        "fasta": str(FASTA),
    }
    OUT_GIT.write_text(json.dumps(join, indent=2), encoding="utf-8")
    print(f"  wrote {OUT_GIT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
