#!/usr/bin/env python3
"""Probe: does real-homolog template transfer generalize from protein to RNA?

Selects fresh RNA chains from RCSB, finds an RNA homolog structure (excluding the
target), transfers its C1' coordinates onto the aligned query, and measures RMSD.
Same real-data, zero-trained-weight idea as the protein template pipeline, adapted
to nucleotides (C1' backbone, RNA sequence search).
"""

from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_fsot_vs_alphafold_structure import kabsch_rmsd  # noqa: E402
from run_rcsb_template_holdout import nw_align  # noqa: E402

TCACHE = Path.home() / ".cache" / "fsot-genetics" / "rna_template"
TCACHE.mkdir(parents=True, exist_ok=True)
NT3 = {"A": "A", "U": "U", "G": "G", "C": "C", "RA": "A", "RU": "U", "RG": "G",
       "RC": "C", "ADE": "A", "URA": "U", "GUA": "G", "CYT": "C"}


def _post(url, body, timeout=60):
    req = urllib.request.Request(url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "fsot"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read())


def _pdb(pid):
    fp = TCACHE / f"{pid}.pdb"
    if fp.exists():
        return fp.read_text(encoding="utf-8", errors="replace")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60) as r:
        txt = r.read().decode("utf-8", "replace")
    fp.write_text(txt, encoding="utf-8")
    return txt


def parse_rna_c1(text, chain):
    seq, coords = [], []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "C1'":
            one = NT3.get(line[17:20].strip())
            if one is None:
                continue
            seq.append(one)
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return "".join(seq), np.array(coords, dtype=np.float64)


def rna_chains(text):
    order, seen = [], set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[12:16].strip() == "C1'":
            c = line[21]
            if c not in seen and line[17:20].strip() in NT3:
                seen.add(c)
                order.append(c)
    return order


def select_rna(n=25):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "entity_poly.rcsb_entity_polymer_type", "operator": "exact_match", "value": "RNA"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined", "operator": "less", "value": 3.0}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.deposited_polymer_monomer_count", "operator": "range",
                "value": {"from": 30, "to": 120}}}]},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": 120},
                             "sort": [{"sort_by": "rcsb_id", "direction": "asc"}],
                             "results_content_type": ["experimental"]}}
    return [h["identifier"] for h in _post("https://search.rcsb.org/rcsbsearch/v2/query", q).get("result_set", [])]


def rna_homologs(seq):
    q = {"query": {"type": "terminal", "service": "sequence", "parameters": {
            "evalue_cutoff": 1.0, "identity_cutoff": 0.4, "sequence_type": "rna", "value": seq}},
         "return_type": "polymer_entity",
         "request_options": {"paginate": {"start": 0, "rows": 30}, "results_content_type": ["experimental"]}}
    try:
        d = _post("https://search.rcsb.org/rcsbsearch/v2/query", q)
    except Exception:
        return []
    ids, seen = [], set()
    for h in d.get("result_set", []):
        p = h["identifier"].split("_")[0].upper()
        if p not in seen:
            seen.add(p)
            ids.append(p)
    return ids


def main():
    rows = []
    for pid in select_rna():
        if len(rows) >= 12:
            break
        try:
            txt = _pdb(pid)
        except Exception:
            continue
        chs = rna_chains(txt)
        if not chs:
            continue
        seq, X = parse_rna_c1(txt, chs[0])
        if not (30 <= len(seq) <= 120) or len(X) != len(seq):
            continue
        best = None
        for hp in rna_homologs(seq):
            if hp == pid.upper():
                continue
            try:
                htxt = _pdb(hp)
            except Exception:
                continue
            for hc in rna_chains(htxt):
                hseq, hX = parse_rna_c1(htxt, hc)
                if len(hseq) < 15:
                    continue
                pairs = nw_align(seq, hseq)
                if len(pairs) < 10:
                    continue
                ident = sum(1 for qi, ti in pairs if seq[qi] == hseq[ti]) / len(pairs)
                cov = len(pairs) / len(seq)
                if ident > 0.95 or cov < 0.6:
                    continue
                if best is None or cov * ident > best[0]:
                    best = (cov * ident, hp, hc, ident, cov, pairs, hX)
            if best and best[0] > 0.7:
                break
        if not best:
            print(f"{pid} n={len(seq):3d}  no RNA homolog")
            continue
        _s, hp, hc, ident, cov, pairs, hX = best
        qi = [a for a, b in pairs]
        ti = [b for a, b in pairs]
        rmsd = kabsch_rmsd(hX[ti], X[qi])
        rows.append(rmsd)
        print(f"{pid} n={len(seq):3d}  RNA template {hp} id={ident:.2f} cov={cov:.2f} -> {rmsd:.2f} A")
    if rows:
        print("--- RNA template-transfer median C1' RMSD: %.2f A  (n=%d) ---" % (
            float(np.median(rows)), len(rows)))


if __name__ == "__main__":
    raise SystemExit(main())
