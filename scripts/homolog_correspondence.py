#!/usr/bin/env python3
"""Measured homologs of residual-mass proteins in other insects.

Predict *up* only where UniProt/OrthoDB lists a homolog. Fold with the
existing product Cα path. A genome is not a connectome — this does not
invent a bee, mosquito, or beetle wiring diagram.

  python scripts/homolog_correspondence.py
  python scripts/homolog_correspondence.py --resolve-only

Sequences stay on D:\\FlyWire_Connectome\\homologs (not git).
"""
from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from fsot_predict import main as predict_main  # noqa: E402

_PHI = float(fc.PHI)
# Leftover coverage analog: hit shorter than source/φ² is a fragment, not a homolog.
_MIN_LEN_FRAC = 1.0 / (_PHI ** 2)

OUT_D = Path(r"D:\FlyWire_Connectome\homologs\product")
FASTA = Path(r"D:\FlyWire_Connectome\homologs\homologs.fasta")
OUT_GIT = ROOT / "data" / "homolog_correspondence.json"

# Residual-mass proteins already on the live fly/worm graphs.
SOURCES = [
    {
        "symbol": "Gad1",
        "uniprot": "P20228",
        "organism": "Drosophila melanogaster",
        "sits_on": "GABAergic neurons (inhibitory residual)",
    },
    {
        "symbol": "nompC",
        "uniprot": "Q7KIQ2",
        "organism": "Drosophila melanogaster",
        "sits_on": "mechanosensory transduction",
    },
    {
        "symbol": "iav",
        "uniprot": "Q9W3W0",
        "organism": "Drosophila melanogaster",
        "sits_on": "Johnston organ TRPV",
    },
    {
        "symbol": "nan",
        "uniprot": "Q9VUD5",
        "organism": "Drosophila melanogaster",
        "sits_on": "Johnston organ TRPV",
    },
    {
        "symbol": "ChAT",
        "uniprot": "P07668",
        "organism": "Drosophila melanogaster",
        "sits_on": "cholinergic neurons (default excitatory edge)",
    },
    {
        "symbol": "VGlut",
        "uniprot": "Q9VQC0",
        "organism": "Drosophila melanogaster",
        "sits_on": "glutamatergic motor neurons (vnc_motor NMJ)",
    },
    {
        "symbol": "Mhc",
        "uniprot": "P05661",
        "organism": "Drosophila melanogaster",
        "sits_on": "muscle (effector of vnc_motor; not a CNS cell)",
    },
    {
        "symbol": "unc-25",
        "uniprot": "G5EDB7",
        "organism": "Caenorhabditis elegans",
        "sits_on": "GABAergic motor neurons DD/VD/RME",
    },
    {
        "symbol": "mec-4",
        "uniprot": "P24612",
        "organism": "Caenorhabditis elegans",
        "sits_on": "ALML/ALMR touch neurons",
    },
    {
        "symbol": "myo-3",
        "uniprot": "P12844",
        "organism": "Caenorhabditis elegans",
        "sits_on": "body-wall muscle myosin",
    },
]

TAXA = [
    {"organism": "Apis mellifera", "taxid": 7460},
    {"organism": "Anopheles gambiae", "taxid": 7165},
    {"organism": "Tribolium castaneum", "taxid": 7070},
]

_UA = "FSOT-Genetics homolog correspondence (mailto:local)"


def _get_json(url: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=90) as fh:
        return json.loads(fh.read().decode("utf-8"))


def _orthodb_id(entry: dict[str, Any]) -> str | None:
    for x in entry.get("uniProtKBCrossReferences") or []:
        if x.get("database") == "OrthoDB":
            return str(x.get("id") or "") or None
    return None


def _gene_name(rec: dict[str, Any]) -> str:
    genes = rec.get("genes") or []
    if not genes:
        return ""
    gn = genes[0].get("geneName") or {}
    return str(gn.get("value") or "")


def _seq_of(rec: dict[str, Any]) -> str:
    return str(((rec.get("sequence") or {}).get("value")) or "")


def _pick_hit(
    results: list[dict[str, Any]],
    *,
    min_length: int = 0,
) -> dict[str, Any] | None:
    if not results:
        return None
    long_enough = [
        r
        for r in results
        if int(((r.get("sequence") or {}).get("length")) or 0) >= min_length
    ]
    pool = long_enough or []
    if not pool:
        return None
    reviewed = [r for r in pool if r.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"]
    pool = reviewed or pool
    pool = sorted(pool, key=lambda r: int(((r.get("sequence") or {}).get("length")) or 0), reverse=True)
    return pool[0]


def _isoform_accessions(entry: dict[str, Any], acc: str) -> list[str]:
    out = [acc]
    for c in entry.get("comments") or []:
        if c.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for iso in c.get("isoforms") or []:
            for iid in iso.get("isoformIds") or []:
                if iid and iid not in out:
                    out.append(str(iid))
    return out


def _uniref50_cluster_ids(acc: str) -> list[str]:
    """Actual UniRef50 cluster IDs, including isoform clusters (P05661-2).

    Querying UniRef50_{accession} misses when the representative is an
    isoform or a different member (nompC Q7KIQ2 lives in UniRef50_Q9VMR4).
    """
    url = (
        "https://rest.uniprot.org/uniref/search?query="
        + urllib.parse.quote(f"uniprot_id:{acc} AND identity:0.5")
        + "&size=10&format=json"
    )
    recs = _get_json(url).get("results") or []
    out: list[str] = []
    for rec in recs:
        cid = rec.get("id")
        if cid and cid not in out:
            out.append(str(cid))
    return out


def _ensembl_id(rec: dict[str, Any]) -> str:
    for x in rec.get("uniProtKBCrossReferences") or []:
        if x.get("database") == "EnsemblMetazoa":
            return str(x.get("id") or "")
    return ""


def resolve() -> dict[str, Any]:
    FASTA.parent.mkdir(parents=True, exist_ok=True)
    fasta_chunks: list[str] = []
    sources_out = []
    rows = []
    for src in SOURCES:
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
        src_row = {
            **src,
            "orthodb": odb,
            "length": src_len,
            "uniref50_clusters": clusters,
        }
        sources_out.append(src_row)
        if not odb:
            print("  no OrthoDB xref — no_measured_homolog", flush=True)
            for tx in TAXA:
                rows.append(
                    {
                        **src,
                        "target_organism": tx["organism"],
                        "taxid": tx["taxid"],
                        "status": "no_measured_homolog",
                        "reason": "source has no OrthoDB xref",
                    }
                )
            continue
        for tx in TAXA:
            q = f"(xref:orthodb-{odb}) AND (organism_id:{tx['taxid']})"
            url = (
                "https://rest.uniprot.org/uniprotkb/search?query="
                + urllib.parse.quote(q)
                + "&fields=accession,gene_names,organism_name,protein_name,sequence,reviewed,xref_ensemblmetazoa"
                + "&size=50&format=json"
            )
            via = "orthodb"
            try:
                hits = _get_json(url).get("results") or []
            except Exception as exc:
                print(f"  {tx['organism']} UniProt fail ({exc})", flush=True)
                rows.append(
                    {
                        **src,
                        "target_organism": tx["organism"],
                        "taxid": tx["taxid"],
                        "orthodb": odb,
                        "status": "uniprot_error",
                        "reason": str(exc),
                    }
                )
                continue
            hit = _pick_hit(hits, min_length=min_len)
            if not hit:
                hits = []
                via = "uniref50"
                for cid in clusters or [f"UniRef50_{src['uniprot']}"]:
                    q50 = f"uniref_cluster_50:{cid} AND organism_id:{tx['taxid']}"
                    url50 = (
                        "https://rest.uniprot.org/uniprotkb/search?query="
                        + urllib.parse.quote(q50)
                        + "&fields=accession,gene_names,organism_name,protein_name,sequence,reviewed,xref_ensemblmetazoa"
                        + "&size=50&format=json"
                    )
                    try:
                        hits = _get_json(url50).get("results") or []
                    except Exception as exc:
                        print(f"  {tx['organism']} UniRef50 {cid} fail ({exc})", flush=True)
                        hits = []
                    hit = _pick_hit(hits, min_length=min_len)
                    if hit:
                        via = f"uniref50:{cid}"
                        break
            if not hit:
                print(f"  {tx['organism']}: no OrthoDB/UniRef50 member", flush=True)
                rows.append(
                    {
                        **src,
                        "target_organism": tx["organism"],
                        "taxid": tx["taxid"],
                        "orthodb": odb,
                        "status": "no_measured_homolog",
                        "reason": "no UniProt OrthoDB or UniRef50 member in taxon",
                        "n_hits": 0,
                    }
                )
                continue
            acc = hit["primaryAccession"]
            seq = _seq_of(hit)
            gn = _gene_name(hit)
            reviewed = hit.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"
            print(
                f"  {tx['organism']}: {acc} gene={gn or '-'} n={len(seq)} "
                f"reviewed={reviewed} via={via} hits={len(hits)}",
                flush=True,
            )
            fasta_chunks.append(f">{src['symbol']}|{acc}|{tx['organism']}\n")
            for i in range(0, len(seq), 60):
                fasta_chunks.append(seq[i : i + 60] + "\n")
            rows.append(
                {
                    "source_symbol": src["symbol"],
                    "source_uniprot": src["uniprot"],
                    "source_organism": src["organism"],
                    "sits_on": src["sits_on"],
                    "orthodb": odb,
                    "correspondence": via,
                    "target_organism": tx["organism"],
                    "taxid": tx["taxid"],
                    "uniprot": acc,
                    "gene": gn,
                    "length": len(seq),
                    "reviewed": reviewed,
                    "n_hits": len(hits),
                    "status": "measured_homolog",
                    "ensembl": _ensembl_id(hit),
                    "uniref50_cluster": via.split(":", 1)[1] if via.startswith("uniref50:") else None,
                    "free_parameters": 0,
                }
            )
    # Same OrthoDB group already hit in this taxon under another source
    # (nan vs iav, unc-25 vs Gad1, myo-3 vs Mhc).
    measured = [
        r
        for r in rows
        if r.get("status") == "measured_homolog"
    ]
    for r in rows:
        if r.get("status") != "no_measured_homolog":
            continue
        odb = r.get("orthodb")
        tax = r.get("taxid")
        cover = next(
            (
                m
                for m in measured
                if m.get("orthodb") == odb
                and m.get("taxid") == tax
                and m.get("source_uniprot") != r.get("source_uniprot")
                and m.get("source_uniprot") != r.get("uniprot")
            ),
            None,
        )
        if not cover:
            continue
        r["status"] = "covered_by_orthodb_paralog"
        r["covered_by"] = {
            "source_symbol": cover.get("source_symbol"),
            "uniprot": cover.get("uniprot"),
            "ensembl": cover.get("ensembl"),
        }
        r["reason"] = (
            f"same OrthoDB {odb} already measured as "
            f"{cover.get('source_symbol')} {cover.get('uniprot')}"
        )
        print(
            f"  cover {r.get('source_symbol') or r.get('symbol')} "
            f"{r.get('target_organism')} via {cover.get('source_symbol')} "
            f"{cover.get('uniprot')}",
            flush=True,
        )
    FASTA.write_text("".join(fasta_chunks), encoding="utf-8")
    print(f"  wrote {FASTA}", flush=True)
    return {
        "product": "OrthoDB/UniProt homolog of residual-mass protein → FSOT product Cα",
        "pin": "D1D38A",
        "free_parameters": 0,
        "law": "S=K(T1+T2+T3); product Cα only where a measured homolog exists",
        "not": (
            "Not a bee/mosquito/beetle connectome. Genomes without synapses "
            "are not wiring diagrams. Predict up only on the protein."
        ),
        "sources": sources_out,
        "taxa": TAXA,
        "homologs": rows,
        "fasta": str(FASTA),
        "miss_audit": {
            "no_measured_homolog": [
                {
                    "source": r.get("source_symbol") or r.get("symbol"),
                    "target": r.get("target_organism"),
                    "reason": r.get("reason"),
                }
                for r in rows
                if r.get("status") == "no_measured_homolog"
            ],
            "covered_by_orthodb_paralog": [
                {
                    "source": r.get("source_symbol") or r.get("symbol"),
                    "target": r.get("target_organism"),
                    "covered_by": r.get("covered_by"),
                }
                for r in rows
                if r.get("status") == "covered_by_orthodb_paralog"
            ],
            "note": (
                "UniRef50_{accession} misses isoform clusters "
                "(P05661-2, Q9VMR4). Fragments shorter than source/φ² dropped. "
                "DEG/ENaC family in Anopheles is not a mec-4 1:1. "
                "Tribolium nompC has no full-length UniRef50 member."
            ),
        },
    }


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


def fold_rows(join: dict[str, Any]) -> dict[str, Any]:
    seqs = _parse_fasta(FASTA) if FASTA.exists() else {}
    OUT_D.mkdir(parents=True, exist_ok=True)
    folded = []
    cache: dict[str, dict[str, Any]] = {}
    for row in join.get("homologs") or []:
        if row.get("status") != "measured_homolog":
            folded.append(row)
            continue
        acc = row["uniprot"]
        seq = seqs.get(acc)
        if not seq:
            row = {**row, "predict_rc": None, "reason": "missing fasta"}
            folded.append(row)
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
    n_ok = sum(1 for r in folded if r.get("status") == "measured_homolog")
    n_fold = sum(1 for r in folded if r.get("predict_rc") == 0)
    join["n_measured_homologs"] = n_ok
    join["n_folded"] = n_fold
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
