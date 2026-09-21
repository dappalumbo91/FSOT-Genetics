#!/usr/bin/env python3
"""FSOT map for chains that still have no full measured structure.

Coordinates only where a domain crystal clears the close-homolog floor
(identity ≥ 1/φ). Everything else stays an observable: Rg target,
secondary, and the measured interaction it sits on.

Plant edges are IntAct physical PPIs (Biochemistry residual).
DNA-binding domains are Electromagnetism observers when a crystal exists.
Insect mec-4 sequence hits with no close map stay on the analog pointer
(nompC), not a transferred DEG/ENaC fold.

  python scripts/orphan_system_map.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from domain_interface import domain_slice  # noqa: E402
from fsot_structure_engine import sequence_observables  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402
from domain_split_assemble import assemble_domains, fetch_interpro_domains  # noqa: E402
from plant_product import _fetch_fasta  # noqa: E402
from plant_signal import CACHE, CLASSES, load_graph  # noqa: E402

PHI = float(fc.PHI)
CLOSE = 1.0 / PHI
DATA = ROOT / "data"
OUT = DATA / "orphan_system_map.json"

# DNA-binding names from InterPro. The nucleic acid is the observer.
_DNA = ("dna", "bzip", "helix-loop-helix", "b3 ")


def _chem(name: str, close: bool) -> dict[str, Any]:
    dna = any(k in name.lower() for k in _DNA)
    domain = "Electromagnetism" if dna else "Biochemistry"
    sl = domain_slice(domain)
    if not close:
        role = "no_measured_map"
    elif dna:
        role = "dna_observer"
    else:
        role = "measured_domain"
    return {
        "domain": domain,
        "D_eff": int(sl.D_eff),
        "role": role,
        "residual": residual_scale(abs(float(fc.domain_scalar(domain)))),
    }


def _neighbors(graph: dict[str, Any], acc: str, k: int = 8) -> list[dict[str, Any]]:
    idx = graph["idx"]
    if acc not in idx:
        return []
    i = idx[acc]
    row = graph["W"].getrow(i)
    nodes = graph["nodes"]
    genes = graph["genes"]
    pairs = []
    for j, w in zip(row.indices, row.data):
        b = nodes[int(j)]
        info = CLASSES.get(b) or {}
        pairs.append(
            (
                float(w),
                {
                    "uniprot": b,
                    "gene": info.get("symbol") or genes.get(b) or "",
                    "class": info.get("class") or "unlabeled",
                    "n_reports": int(w),
                    "chem_link": "tertiary_biochem",
                    "domain": "Biochemistry",
                    "D_eff": 13,
                },
            )
        )
    pairs.sort(key=lambda t: -t[0])
    return [p[1] for p in pairs[:k]]


def _span_card(symbol: str, acc: str, domains: list[dict[str, Any]], graph: dict[str, Any]) -> dict[str, Any]:
    seq = _fetch_fasta(acc)
    obs = sequence_observables(seq)
    close = [d for d in domains if d.get("close_homolog")]
    covered = sum(int(d.get("length") or 0) for d in close)
    return {
        "symbol": symbol,
        "uniprot": acc,
        "length": len(seq),
        "full_chain": "no_measured_map",
        "rg_target_A": obs["rg_target_A"],
        "secondary_head": (obs.get("secondary") or "")[:48],
        "n_close_domains": len(close),
        "close_residues": covered,
        "product_fraction": covered / len(seq) if seq else 0.0,
        "domains": [
            {
                **d,
                "interface": _chem(str(d.get("name") or ""), bool(d.get("close_homolog"))),
            }
            for d in domains
        ],
        "measured_partners": _neighbors(graph, acc),
        "on_intact_graph": acc in graph["idx"],
    }


def main() -> int:
    if not CACHE.is_file():
        raise SystemExit(f"missing IntAct cache {CACHE}")
    graph = load_graph(CACHE)
    bio = domain_slice("Biochemistry")
    r_bio = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))
    domains_doc = json.loads((DATA / "plant_signal_domains.json").read_text(encoding="utf-8"))
    orphans = domains_doc.get("orphans") or {}
    cards = []
    for symbol in ("PIF3", "HY5", "ARF5"):
        block = orphans.get(symbol) or {}
        card = _span_card(symbol, block["uniprot"], block.get("domains") or [], graph)
        print(
            f"{symbol} on_graph={card['on_intact_graph']} "
            f"close_domains={card['n_close_domains']} "
            f"product_fraction={card['product_fraction']:.3f} "
            f"top={card['measured_partners'][:1]}",
            flush=True,
        )
        cards.append(card)

    # Insect sequence hits that never received a coordinate map.
    homo = json.loads((DATA / "homolog_correspondence.json").read_text(encoding="utf-8"))
    insects = []
    for row in homo.get("homologs") or []:
        if row.get("structure_mode") != "no_measured_map":
            continue
        if row.get("source_symbol") != "mec-4":
            continue
        acc = str(row.get("uniprot") or "")
        print(f"mec-4 {acc} domain search", flush=True)
        seq = _fetch_fasta(acc) if acc else ""
        doms = fetch_interpro_domains(acc) if acc else []
        asm = assemble_domains(seq, doms, exclude_pdb="XXXX", identity_cap=1.0) if seq and doms else {}
        drows = []
        n_close = 0
        for d in asm.get("domains") or []:
            ident = d.get("template_identity")
            close = bool(d.get("template_pdb") and ident is not None and float(ident) >= CLOSE)
            if close:
                n_close += 1
            drows.append(
                {
                    "name": d.get("name"),
                    "start": d.get("start"),
                    "end": d.get("end"),
                    "template_pdb": d.get("template_pdb"),
                    "template_identity": ident,
                    "template_coverage": d.get("template_coverage"),
                    "close_homolog": close,
                    "source": d.get("source"),
                }
            )
            print(
                f"  {d.get('name')} tmpl={d.get('template_pdb')} id={ident} "
                f"{'close' if close else 'not close'}",
                flush=True,
            )
        insects.append(
            {
                "source": "mec-4",
                "uniprot": acc,
                "organism": row.get("target_organism"),
                "length": len(seq) or row.get("length"),
                "full_chain": "no_measured_map",
                "n_close_domains": n_close,
                "domains": drows,
                "why": (
                    "No close measured map on the full chain. "
                    "A domain crystal is product only at identity ≥ 1/φ. "
                    "The live gentle-touch protein on the fly graph is nompC "
                    "(5VKQ). Do not transfer a distant DEG/ENaC crystal onto mec-4."
                ),
            }
        )

    out = {
        "product": "FSOT system map for chains without a full measured structure",
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3)",
        "close_homolog": CLOSE,
        "residual_Biochemistry": r_bio,
        "biochem_D_eff": int(bio.D_eff),
        "not": (
            "Not a docked complex. Not a 13 Å MDS fill. "
            "Partner edges are IntAct counts. Domain Cα only if identity ≥ 1/φ."
        ),
        "proteins": cards,
        "insect_sequence_hits": insects,
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
