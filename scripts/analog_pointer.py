#!/usr/bin/env python3
"""When a 1:1 protein is unmapped, FSOT points at the mapped job.

Same physics, different gene. mec-4 has no mosquito 1:1. The live graphs
already map the residual job (mechanotransduction) onto nompC. Do not
invent a degenerin; follow the measured analog.

  python scripts/analog_pointer.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "analog_pointer.json"
HOMOLOG = ROOT / "data" / "homolog_correspondence.json"
BANC = ROOT / "data" / "banc_connectome_boot.json"
MALE = ROOT / "data" / "male_cns_boot.json"
WORM = ROOT / "data" / "worm_connectome_boot.json"


def _hop2_vnc(path: Path, program: str) -> dict | None:
    if not path.exists():
        return None
    d = json.loads(path.read_text(encoding="utf-8"))
    rec = ((d.get("compare") or {}).get("hop_2") or {}).get(program) or {}
    if not rec:
        return None
    top = rec.get("top") or {}
    return {
        "vnc_motor": rec.get("vnc_motor"),
        "descending": rec.get("descending"),
        "top": top.get("cell_type") or top.get("name"),
    }


def main() -> int:
    homologs = []
    if HOMOLOG.exists():
        homologs = json.loads(HOMOLOG.read_text(encoding="utf-8")).get("homologs") or []

    def _hit(symbol: str, taxid: int) -> dict | None:
        for r in homologs:
            if r.get("source_symbol") == symbol and r.get("taxid") == taxid:
                if r.get("status") == "measured_homolog":
                    return {
                        "uniprot": r.get("uniprot"),
                        "ensembl": r.get("ensembl"),
                        "template_pdb": r.get("template_pdb"),
                        "template_identity": r.get("template_identity"),
                        "template_coverage": r.get("template_coverage"),
                        "structure_mode": r.get("structure_mode"),
                    }
        return None

    pointer = {
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3); residual Biochemistry on measured graphs; product Cα on measured homologs",
        "rule": (
            "A 1:1 miss is not a license to invent. The scalar law is the same "
            "in every organism. When the protein is clade-restricted, point at "
            "the mapped residual job in the species that has a graph, then at "
            "the homolog of THAT protein."
        ),
        "jobs": [
            {
                "job": "mechanotransduction (touch / JO / chordotonal)",
                "unmapped_1to1": {
                    "symbol": "mec-4",
                    "source": "Caenorhabditis elegans ALML/ALMR",
                    "blank_in": "Anopheles gambiae",
                    "why": "UniRef50/90 of mec-4 is Nematoda only; fly ppk ~18%",
                },
                "mapped_analog": {
                    "symbol": "nompC",
                    "why": (
                        "On the measured fly graphs, JO / mechanosensory / "
                        "vnc_sensory light vnc_motor; olfactory does not. "
                        "nompC sits on that seed. Worm ALM hop-0 is mec-4; "
                        "the job is the same physics, the gene family is not."
                    ),
                    "fly_male_cns_hop2": _hop2_vnc(MALE, "JO") or _hop2_vnc(MALE, "mechanosensory"),
                    "fly_banc_hop2": _hop2_vnc(BANC, "JO") or _hop2_vnc(BANC, "chordotonal"),
                    "anopheles": _hit("nompC", 7165),
                    "tribolium": _hit("nompC", 7070),
                    "apis": _hit("nompC", 7460),
                },
                "do_not": "fold a random Anopheles ppk as mec-4",
            },
            {
                "job": "inhibitory residual (GABA synthesis)",
                "unmapped_1to1": {
                    "symbol": "unc-25",
                    "blank_in": "Anopheles gambiae UniRef50 of the worm accession",
                    "why": "worm GAD cluster does not reach insects",
                },
                "mapped_analog": {
                    "symbol": "Gad1",
                    "anopheles": _hit("Gad1", 7165),
                    "why": "same OrthoDB group as unc-25; already folded",
                },
                "do_not": "treat the worm UniRef50 miss as a missing enzyme",
            },
        ],
        "connectome_accuracy": {
            "what_is_scored": (
                "Not RMSD. Hop mass on the measured graph vs known biology: "
                "which sensory class is on decides whether motor/descending lights."
            ),
            "male_cns_hop2_vnc_motor": {
                "vnc_sensory": 10.74,
                "JO": 7.77,
                "olfactory": 0.001,
            },
            "banc_hop2_vnc_motor": {
                "vnc_sensory": 10.16,
                "JO": 6.64,
                "olfactory": 0.0008,
            },
            "independent_sexes": "same split on male CNS and female BANC",
        },
        "product_accuracy": {
            "freeze_median_Ca_A": 0.13,
            "alphafold_same_10_median_A": 0.47,
            "source": "docs/PRODUCT_FREEZE.md",
            "homolog_transfers": (
                "Cross-species product uses the same path. Identity/coverage "
                "are the accuracy gate (close-homolog id ≥ 1/φ). Orphans stay "
                "no_measured_map."
            ),
        },
    }
    OUT.write_text(json.dumps(pointer, indent=2), encoding="utf-8")
    print(json.dumps(pointer["jobs"][0]["mapped_analog"], indent=2))
    print(f"  wrote {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
