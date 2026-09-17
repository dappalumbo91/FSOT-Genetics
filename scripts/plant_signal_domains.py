#!/usr/bin/env python3
"""Domain-split product for plant signaling proteins that lack a full-chain map.

PHOT1 Q2V2M9 is 1422 aa, full-chain no_measured_map. Phototropin LOV/kinase
domains have measured crystals. Same product rule per domain; inter-domain
pose is not a claimed fold. UVR8 was missing from the IntAct hop graph —
fold it only if a measured homolog exists. Not a connectome.

  python scripts/plant_signal_domains.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from domain_split_assemble import (  # noqa: E402
    assemble_domains,
    fetch_interpro_domains,
)
from fsot_predict import main as predict_main  # noqa: E402
from plant_product import _fetch_fasta, _parse_fasta  # noqa: E402

OUT_D = Path(r"D:\FlyWire_Connectome\plants\product")
FASTA = Path(r"D:\FlyWire_Connectome\plants\arabidopsis_signal.fasta")
OUT_GIT = ROOT / "data" / "plant_signal_domains.json"
PHI = (1.0 + 5.0 ** 0.5) / 2.0
CLOSE = 1.0 / PHI


def _seq(acc: str, fasta: dict[str, str]) -> str:
    if acc in fasta and fasta[acc]:
        return fasta[acc]
    return _fetch_fasta(acc)


def main() -> int:
    OUT_D.mkdir(parents=True, exist_ok=True)
    seqs = _parse_fasta(FASTA) if FASTA.exists() else {}

    phot_acc = "O48963"
    phot_seq = _seq(phot_acc, seqs)
    print(f"PHOT1 {phot_acc} n={len(phot_seq)} InterPro domains", flush=True)
    domains = fetch_interpro_domains(phot_acc)
    print(f"  n_domains={len(domains)}", flush=True)
    for d in domains:
        print(f"    {d.pfam} {d.name} {d.start}-{d.end} L={d.length}", flush=True)

    assembly = assemble_domains(phot_seq, domains, exclude_pdb="XXXX", identity_cap=1.0)
    phot_domains = []
    n_product = 0
    for d in assembly.get("domains") or []:
        ident = d.get("template_identity")
        cov = d.get("template_coverage")
        mode = d.get("source") or "none"
        product = bool(
            d.get("template_pdb")
            and ident is not None
            and float(ident) >= CLOSE
        )
        if product:
            n_product += 1
        phot_domains.append(
            {
                "pfam": d.get("pfam"),
                "name": d.get("name"),
                "start": d.get("start"),
                "end": d.get("end"),
                "length": d.get("length"),
                "source": mode,
                "template_pdb": d.get("template_pdb"),
                "template_identity": ident,
                "template_coverage": cov,
                "close_homolog": product,
            }
        )
        print(
            f"  {d.get('name')} {d.get('start')}-{d.get('end')} "
            f"tmpl={d.get('template_pdb')} id={ident} cov={cov} "
            f"{'product' if product else 'not close-homolog'}",
            flush=True,
        )

    # UVR8: missing from IntAct hops; still a light receptor with crystals.
    uvr_acc = "Q9FN03"
    uvr_seq = _seq(uvr_acc, seqs)
    uvr_pdb = OUT_D / f"UVR8_{uvr_acc}.pdb"
    uvr_json = OUT_D / f"UVR8_{uvr_acc}.json"
    print(f"UVR8 {uvr_acc} n={len(uvr_seq)}", flush=True)
    if not uvr_json.exists():
        rc = predict_main(
            [
                "--seq",
                uvr_seq,
                "--uniprot",
                uvr_acc,
                "--pdb-out",
                str(uvr_pdb),
                "--json-out",
                str(uvr_json),
            ]
        )
    else:
        rc = 0
        print(f"  skip existing {uvr_json.name}", flush=True)
    uvr: dict = {
        "symbol": "UVR8",
        "uniprot": uvr_acc,
        "length": len(uvr_seq),
        "predict_rc": rc,
        "on_intact_graph": False,
        "note": "Missing from IntAct photoreceptor hop (seed_missing). Product only if homolog exists.",
    }
    if uvr_json.exists():
        full = json.loads(uvr_json.read_text(encoding="utf-8"))
        uvr.update(
            {
                "structure_mode": full.get("structure_mode"),
                "template_pdb": full.get("template_pdb"),
                "template_identity": full.get("template_identity"),
                "template_coverage": full.get("template_coverage"),
                "pdb": str(uvr_pdb) if uvr_pdb.exists() else None,
            }
        )
    print(
        f"  UVR8 mode={uvr.get('structure_mode')} tmpl={uvr.get('template_pdb')} "
        f"id={uvr.get('template_identity')}",
        flush=True,
    )

    report = {
        "product": "Plant signaling domain product (not a connectome)",
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3); per-domain measured homolog; inter-domain pose not claimed",
        "close_homolog_floor": CLOSE,
        "PHOT1": {
            "uniprot": phot_acc,
            "length": len(phot_seq),
            "full_chain": "leftover_template_5HZI_id_0.61_cov_0.46",
            "n_domains": len(domains),
            "n_close_homolog_domains": n_product,
            "covered_fraction": (assembly.get("covered_fraction") if assembly else None),
            "joint_template": assembly.get("joint_template"),
            "domains": phot_domains,
            "note": "Full-chain leftover map 5HZI id 0.61 (<1/φ) cov 0.46. Domain crystals: LOV2 4HHD id 1.00. Kinase 4L3J id 0.55 is not close-homolog. Inter-domain pose not claimed.",
        },
        "UVR8": uvr,
        "not": "Not plant synapses. Not a 13 Å MDS PHOT1 brain.",
    }
    OUT_GIT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"  wrote {OUT_GIT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
