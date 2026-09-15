#!/usr/bin/env python3
"""Arabidopsis panel → rice / maize measured homologs → product Cα.

Same product rule as animals: fold only with a measured homolog. Plants have
genomes and crystals, not a fly-class EM connectome. Do not invent synapses.

  python scripts/plant_homolog.py
  python scripts/plant_homolog.py --resolve-only

Sequences on D:\\FlyWire_Connectome\\plants\\homologs (not git).
"""
from __future__ import annotations

import json
import sys
import urllib.parse
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_predict import main as predict_main  # noqa: E402
from homolog_correspondence import (  # noqa: E402
    _MIN_LEN_FRAC,
    _ensembl_id,
    _gene_name,
    _get_json,
    _isoform_accessions,
    _orthodb_id,
    _parse_fasta,
    _seq_of,
    _uniref50_cluster_ids,
)
from plant_product import PANEL  # noqa: E402

OUT_D = Path(r"D:\FlyWire_Connectome\plants\homologs\product")
FASTA = Path(r"D:\FlyWire_Connectome\plants\homologs\plant_homologs.fasta")
OUT_GIT = ROOT / "data" / "plant_homolog.json"

# Japonica (39947) is the reviewed rice proteome. Species-level 4530 is fallback
# when a chloroplast entry is tagged only at Oryza sativa.
TAXA = [
    {
        "organism": "Oryza sativa Japonica Group",
        "short": "rice",
        "taxid": 39947,
        "taxid_fallback": 4530,
        "proteome": "UP000059680",
    },
    {
        "organism": "Zea mays",
        "short": "maize",
        "taxid": 4577,
        "taxid_fallback": None,
        "proteome": "UP000007305",
    },
]


def _ensembl_plants(rec: dict[str, Any]) -> str:
    for x in rec.get("uniProtKBCrossReferences") or []:
        if x.get("database") in ("EnsemblPlants", "Ensembl"):
            return str(x.get("id") or "")
    return _ensembl_id(rec)


def _gene_stem(symbol: str) -> str:
    # LHCB1.3 is subtype LHCB1 (not Lhcb5). CESA3/ACT2/GAPA1 drop the paralog digit.
    if "." in symbol:
        return "".join(ch for ch in symbol.split(".")[0].upper() if ch.isalnum())
    s = "".join(ch for ch in symbol.upper() if ch.isalnum())
    i = len(s)
    while i and s[i - 1].isdigit():
        i -= 1
    return s[:i] or s


def _norm_gene(name: str) -> str:
    return "".join(ch for ch in (name or "").upper() if ch.isalnum())


def _protein_name(rec: dict[str, Any]) -> str:
    rec_name = ((rec.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName") or {}
    if isinstance(rec_name, dict):
        return str(rec_name.get("value") or "")
    return str(rec_name or "")


def _ordered_locus(rec: dict[str, Any]) -> str:
    for g in rec.get("genes") or []:
        for n in g.get("orderedLocusNames") or []:
            val = n.get("value") if isinstance(n, dict) else n
            if val:
                return str(val)
    return ""


def _is_reviewed(rec: dict[str, Any]) -> bool:
    return rec.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"


def _len_of(rec: dict[str, Any]) -> int:
    return int(((rec.get("sequence") or {}).get("length")) or 0)


def _rank_hit(rec: dict[str, Any], src_len: int) -> tuple:
    """Reviewed first, then closest length to the source, then longer."""
    return (
        0 if _is_reviewed(rec) else 1,
        abs(_len_of(rec) - src_len),
        -_len_of(rec),
    )


def _pick_plant_hit(
    results: list[dict[str, Any]],
    *,
    symbol: str,
    min_length: int,
    src_len: int,
) -> dict[str, Any] | None:
    """Do not let OrthoDB superfamilies pick CSLD as CESA or GAPC as GAPA."""
    pool = [r for r in results if _len_of(r) >= min_length]
    if not pool:
        return None
    stem = _gene_stem(symbol)
    want = _norm_gene(symbol)

    def names(rec: dict[str, Any]) -> str:
        return _norm_gene(_gene_name(rec))

    exact = [r for r in pool if names(r) == want]
    if exact:
        return sorted(exact, key=lambda r: _rank_hit(r, src_len))[0]
    same = [r for r in pool if names(r).startswith(stem) and stem]
    # CESA stem must not accept CSLD/CSLA/CSLF; GAPA must not accept GAPC.
    if stem == "CESA":
        same = [r for r in same if names(r).startswith("CESA") and not names(r).startswith("CSL")]
    if stem == "GAPA":
        same = [r for r in same if names(r).startswith("GAPA") and not names(r).startswith("GAPC")]
    if same:
        return sorted(same, key=lambda r: _rank_hit(r, src_len))[0]
    if stem == "GAPA":
        pool = [r for r in pool if not names(r).startswith("GAPC")]
        if not pool:
            return None
        return sorted(pool, key=lambda r: _rank_hit(r, src_len))[0]
    if stem == "CESA":
        cesa = [r for r in pool if names(r).startswith("CESA") and not names(r).startswith("CSL")]
        if cesa:
            return sorted(cesa, key=lambda r: _rank_hit(r, src_len))[0]
    return sorted(pool, key=lambda r: _rank_hit(r, src_len))[0]


def _search_taxon(query_core: str, taxid: int, size: int = 50) -> list[dict[str, Any]]:
    q = f"({query_core}) AND (organism_id:{taxid})"
    url = (
        "https://rest.uniprot.org/uniprotkb/search?query="
        + urllib.parse.quote(q)
        + "&fields=accession,gene_names,organism_name,protein_name,sequence,reviewed,xref_ensemblplants"
        + f"&size={size}&format=json"
    )
    return _get_json(url).get("results") or []


def _search_with_fallback(
    query_core: str,
    tx: dict[str, Any],
    *,
    symbol: str,
    min_length: int,
    src_len: int,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]], int]:
    used = int(tx["taxid"])
    hits = _search_taxon(query_core, used)
    hit = _pick_plant_hit(hits, symbol=symbol, min_length=min_length, src_len=src_len)
    if hit or not tx.get("taxid_fallback"):
        return hit, hits, used
    used = int(tx["taxid_fallback"])
    hits = _search_taxon(query_core, used)
    hit = _pick_plant_hit(hits, symbol=symbol, min_length=min_length, src_len=src_len)
    return hit, hits, used


def _named_gene_search(
    symbol: str, tx: dict[str, Any], min_length: int, src_len: int
) -> tuple[dict[str, Any] | None, list[dict[str, Any]], int, str]:
    """UniProt gene:SYMBOL in taxon. Same measured bar as NCBI named nompC."""
    stem = _gene_stem(symbol)
    # Exact gene symbol only. Stem wildcards pick the wrong CESA/LHCB paralog.
    queries = [f"gene:{symbol}"]
    if symbol.upper() == "ACT2":
        queries.append("gene:ACT1")
    if _norm_gene(symbol).startswith("LHCB1"):
        # CAB1 / CAB1R is the UniProt gene for LHCII type I in grasses.
        queries.extend(["gene:CAB1", "gene:CAB1R"])
    for q in queries:
        hit, hits, used = _search_with_fallback(
            q, tx, symbol=symbol, min_length=min_length, src_len=src_len
        )
        if hit:
            gn = _norm_gene(_gene_name(hit))
            want = _norm_gene(symbol)
            via = "gene_name"
            if gn != want:
                if want.startswith("LHCB1") and gn in ("CAB1", "CAB1R"):
                    via = "gene_name"
                elif gn.startswith(stem):
                    via = "gene_name_paralog"
                else:
                    via = "gene_name_family"
            return hit, hits, used, via
    return None, [], int(tx["taxid"]), "gene_name"


def resolve() -> dict[str, Any]:
    FASTA.parent.mkdir(parents=True, exist_ok=True)
    fasta_chunks: list[str] = []
    sources_out = []
    rows: list[dict[str, Any]] = []
    for src in PANEL:
        print(f"== {src['symbol']} {src['uniprot']}", flush=True)
        entry = _get_json(f"https://rest.uniprot.org/uniprotkb/{src['uniprot']}.json")
        odb = _orthodb_id(entry)
        src_len = int(((entry.get("sequence") or {}).get("length")) or 0)
        min_len = int(round(src_len * _MIN_LEN_FRAC))
        clusters: list[str] = []
        for iso_acc in _isoform_accessions(entry, src["uniprot"])[:12]:
            for cid in _uniref50_cluster_ids(iso_acc):
                if cid not in clusters:
                    clusters.append(cid)
        sources_out.append(
            {
                **src,
                "source_organism": "Arabidopsis thaliana",
                "source_taxid": 3702,
                "orthodb": odb,
                "length": src_len,
                "uniref50_clusters": clusters,
            }
        )
        for tx in TAXA:
            hit = None
            hits: list[dict[str, Any]] = []
            via = "gene_name"
            used_tax = int(tx["taxid"])
            try:
                hit, hits, used_tax, via = _named_gene_search(
                    src["symbol"], tx, min_len, src_len
                )
            except Exception as exc:  # noqa: BLE001
                print(f"  {tx['organism']} gene-name fail ({exc})", flush=True)
                hit = None
            if not hit:
                via = "uniref50"
                for cid in clusters or [f"UniRef50_{src['uniprot']}"]:
                    try:
                        hit, hits, used_tax = _search_with_fallback(
                            f"uniref_cluster_50:{cid}",
                            tx,
                            symbol=src["symbol"],
                            min_length=min_len,
                            src_len=src_len,
                        )
                    except Exception as exc:  # noqa: BLE001
                        print(f"  {tx['organism']} UniRef50 {cid} fail ({exc})", flush=True)
                        hit, hits = None, []
                    if hit:
                        via = f"uniref50:{cid}"
                        break
            if not hit and odb:
                via = "orthodb"
                try:
                    hit, hits, used_tax = _search_with_fallback(
                        f"xref:orthodb-{odb}",
                        tx,
                        symbol=src["symbol"],
                        min_length=min_len,
                        src_len=src_len,
                    )
                except Exception as exc:  # noqa: BLE001
                    print(f"  {tx['organism']} OrthoDB fail ({exc})", flush=True)
                    rows.append(
                        {
                            "source_symbol": src["symbol"],
                            "source_uniprot": src["uniprot"],
                            "source_organism": "Arabidopsis thaliana",
                            "sits_on": src["sits_on"],
                            "system": src["system"],
                            "target_organism": tx["organism"],
                            "taxid": tx["taxid"],
                            "orthodb": odb,
                            "status": "uniprot_error",
                            "reason": str(exc),
                        }
                    )
                    continue
            if not hit:
                print(f"  {tx['organism']}: no named-gene/UniRef50/OrthoDB member", flush=True)
                rows.append(
                    {
                        "source_symbol": src["symbol"],
                        "source_uniprot": src["uniprot"],
                        "source_organism": "Arabidopsis thaliana",
                        "sits_on": src["sits_on"],
                        "system": src["system"],
                        "target_organism": tx["organism"],
                        "taxid": tx["taxid"],
                        "orthodb": odb,
                        "status": "no_measured_homolog",
                        "reason": (
                            "no UniProt named gene, UniRef50, or OrthoDB member in taxon. "
                            "Chloroplast UniRef50 of Arabidopsis rbcL/psbA is a split cluster; "
                            "named gene is the measured recover."
                        ),
                        "n_hits": 0,
                    }
                )
                continue
            acc = hit["primaryAccession"]
            seq = _seq_of(hit)
            gn = _gene_name(hit)
            oln = _ordered_locus(hit)
            reviewed = hit.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"
            print(
                f"  {tx['organism']}: {acc} gene={gn or oln or '-'} n={len(seq)} "
                f"reviewed={reviewed} via={via} tax={used_tax} hits={len(hits)}",
                flush=True,
            )
            fasta_chunks.append(f">{src['symbol']}|{acc}|{tx['organism']}\n")
            for i in range(0, len(seq), 60):
                fasta_chunks.append(seq[i : i + 60] + "\n")
            rows.append(
                {
                    "source_symbol": src["symbol"],
                    "source_uniprot": src["uniprot"],
                    "source_organism": "Arabidopsis thaliana",
                    "sits_on": src["sits_on"],
                    "system": src["system"],
                    "orthodb": odb,
                    "correspondence": via,
                    "target_organism": tx["organism"],
                    "taxid": used_tax,
                    "query_taxid": tx["taxid"],
                    "uniprot": acc,
                    "gene": gn or oln,
                    "ordered_locus": oln or None,
                    "length": len(seq),
                    "reviewed": reviewed,
                    "n_hits": len(hits),
                    "status": "measured_homolog",
                    "ensembl": _ensembl_plants(hit),
                    "uniref50_cluster": via.split(":", 1)[1] if via.startswith("uniref50:") else None,
                    "free_parameters": 0,
                }
            )
    FASTA.write_text("".join(fasta_chunks), encoding="utf-8")
    print(f"  wrote {FASTA}", flush=True)
    return {
        "product": "Arabidopsis panel → rice/maize measured homolog → FSOT product Cα",
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3); product Cα only where a measured homolog exists",
        "not": (
            "Not a plant connectome. Rice and maize have genomes and crystals. "
            "No fly-class EM wiring diagram. Predict proteins, not invented synapses."
        ),
        "source_organism": "Arabidopsis thaliana",
        "source_taxid": 3702,
        "source_proteome": "UP000006548",
        "sources": sources_out,
        "taxa": TAXA,
        "homologs": rows,
        "fasta": str(FASTA),
    }


def fold_rows(join: dict[str, Any]) -> dict[str, Any]:
    seqs = _parse_fasta(FASTA) if FASTA.exists() else {}
    OUT_D.mkdir(parents=True, exist_ok=True)
    folded: list[dict[str, Any]] = []
    cache: dict[str, dict[str, Any]] = {}
    for row in join.get("homologs") or []:
        if row.get("status") != "measured_homolog":
            folded.append(row)
            continue
        acc = row["uniprot"]
        seq = seqs.get(acc)
        if not seq:
            folded.append({**row, "predict_rc": None, "reason": "missing fasta"})
            continue
        pdb_out = OUT_D / f"{row['source_symbol']}_{row['taxid']}_{acc}.pdb"
        json_out = OUT_D / f"{row['source_symbol']}_{row['taxid']}_{acc}.json"
        if json_out.exists() and acc not in cache:
            full = json.loads(json_out.read_text(encoding="utf-8"))
            cache[acc] = {
                "predict_rc": 0,
                "pdb": str(pdb_out) if pdb_out.exists() else None,
                "structure_mode": full.get("structure_mode"),
                "deploy_regime": full.get("deploy_regime"),
                "template_pdb": full.get("template_pdb"),
                "template_identity": full.get("template_identity"),
                "template_coverage": full.get("template_coverage"),
                "mean_confidence": full.get("mean_confidence"),
                "rg_target_A": full.get("rg_target_A"),
                "engine": full.get("engine"),
            }
            print(f"  skip existing {json_out.name}", flush=True)
        if acc in cache:
            rec = {**row, **cache[acc], "pdb": str(pdb_out) if pdb_out.exists() else cache[acc].get("pdb")}
            folded.append(rec)
            print(f"  reuse {acc} for {row['source_symbol']} {row['target_organism']}", flush=True)
            continue
        print(
            f"== fold {row['source_symbol']} {row['target_organism']} {acc} n={len(seq)}",
            flush=True,
        )
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
            **row,
            "predict_rc": rc,
            "pdb": str(pdb_out) if pdb_out.exists() else None,
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
                }
            )
        cache[acc] = {
            k: rec[k]
            for k in (
                "predict_rc",
                "pdb",
                "structure_mode",
                "deploy_regime",
                "template_pdb",
                "template_identity",
                "template_coverage",
                "mean_confidence",
                "rg_target_A",
                "engine",
            )
            if k in rec
        }
        print(
            f"  {row['source_symbol']} {acc} mode={rec.get('structure_mode')} "
            f"tmpl={rec.get('template_pdb')} id={rec.get('template_identity')}",
            flush=True,
        )
        folded.append(rec)
    join["homologs"] = folded
    join["n_measured_homologs"] = sum(1 for r in folded if r.get("status") == "measured_homolog")
    join["n_folded"] = sum(1 for r in folded if r.get("predict_rc") == 0)
    join["n_miss"] = sum(1 for r in folded if r.get("status") == "no_measured_homolog")
    return join


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--resolve-only", action="store_true")
    args = ap.parse_args(argv)
    join = resolve()
    if not args.resolve_only:
        join = fold_rows(join)
    OUT_GIT.write_text(json.dumps(join, indent=2), encoding="utf-8")
    print(f"  wrote {OUT_GIT}", flush=True)
    n_ok = sum(1 for r in join["homologs"] if r.get("status") == "measured_homolog")
    n_miss = sum(1 for r in join["homologs"] if r.get("status") == "no_measured_homolog")
    print(f"  measured={n_ok} missing={n_miss}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
