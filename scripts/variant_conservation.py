#!/usr/bin/env python3
"""Variant-effect engine: evolutionary conservation from real homolog structures.

Conservation is the field-standard signal for pathogenicity (SIFT/PolyPhen): a
position that is invariant across evolution is intolerant to mutation. We build a
per-position conservation profile for p53 from real homolog structures (the same
homolog pipeline used for template modeling), then score variants as
conservation x trinary-opcode change and check that the known cancer drivers
surface at the top. Zero trained weights; real evolutionary data as input.
"""

from __future__ import annotations

import gzip
import sys
import urllib.request
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_rcsb_template_holdout import _get, nw_align  # noqa: E402
from run_fsot_vs_alphafold_structure import fetch_pdb  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_headtohead"
MSA_CACHE = Path.home() / ".cache" / "fsot-genetics" / "pfam_msa"
MSA_CACHE.mkdir(parents=True, exist_ok=True)
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}
AA20 = "ARNDCQEGHILKMFPSTWYV"
HOTSPOTS = [(175, "R", "H"), (245, "G", "S"), (248, "R", "Q"),
            (249, "R", "S"), (273, "R", "H"), (282, "R", "W")]


def resnums(text, chain):
    seq, nums, seen = [], [], set()
    for line in text.splitlines():
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
            key = line[22:26].strip()
            if (chain, key) in seen:
                continue
            seen.add((chain, key))
            aa = AA3.get(line[17:20].strip())
            if aa:
                seq.append(aa)
                nums.append(int(key))
    return "".join(seq), nums


def pfam_accession(pdb_id):
    for ent in ("1", "2", "3", "4"):
        try:
            u = _get(f"https://data.rcsb.org/rest/v1/core/uniprot/{pdb_id}/{ent}")
            acc = u[0]["rcsb_uniprot_container_identifiers"]["uniprot_id"]
            fd = _get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/")
            return fd["results"][0]["metadata"]["accession"]
        except Exception:
            continue
    return None


def fetch_msa(pfam, kind="full"):
    fp = MSA_CACHE / f"{pfam}.{kind}.sto"
    if fp.exists():
        return fp.read_text(encoding="utf-8", errors="replace")
    url = f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam}/?annotation=alignment:{kind}"
    req = urllib.request.Request(url, headers={"User-Agent": "fsot"})
    with urllib.request.urlopen(req, timeout=180) as r:
        raw = r.read()
    try:
        txt = gzip.decompress(raw).decode("utf-8", "replace")
    except Exception:
        txt = raw.decode("utf-8", "replace")
    fp.write_text(txt, encoding="utf-8")
    return txt


def parse_stockholm(txt):
    rows: dict[str, str] = {}
    for line in txt.splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) == 2:
            rows[parts[0]] = rows.get(parts[0], "") + parts[1]
    return list(rows.values())


def conservation_profile(seq, self_pdb, pfam=None):
    """Per-position identity conservation from the diverse Pfam family MSA."""
    pfam = pfam or pfam_accession(self_pdb)
    if not pfam:
        return np.zeros(len(seq)), [], 0, None
    rows = parse_stockholm(fetch_msa(pfam))
    # pick the alignment row whose ungapped sequence best matches our query
    best = None
    for row in rows:
        ung = row.replace(".", "").replace("-", "").upper()
        if len(ung) < 20:
            continue
        pairs = nw_align(seq, ung)
        score = sum(1 for a, b in pairs if seq[a] == ung[b])
        if best is None or score > best[0]:
            best = (score, row)
    if best is None:
        return np.zeros(len(seq)), [], 0, pfam
    ref = best[1]
    # map alignment columns -> query residue index via the ref row
    ref_ung = ref.replace(".", "").replace("-", "").upper()
    qpairs = dict(nw_align(seq, ref_ung))  # query_idx -> ref_ungapped_idx
    ung_to_col = [c for c, ch in enumerate(ref) if ch not in ".-"]
    match = np.zeros(len(seq))
    cover = np.zeros(len(seq))
    from collections import Counter
    colcount = [Counter() for _ in range(len(seq))]
    cols = [list(r) for r in rows]
    for qi, ui in qpairs.items():
        if ui >= len(ung_to_col):
            continue
        col = ung_to_col[ui]
        for r in cols:
            if col >= len(r):
                continue
            ch = r[col].upper()
            if ch in AA20:
                cover[qi] += 1
                colcount[qi][ch] += 1
                if ch == seq[qi]:
                    match[qi] += 1
    cons = np.where(cover > 0, match / np.maximum(cover, 1), 0.0)
    freq = [{a: c[a] / max(sum(c.values()), 1) for a in c} for c in colcount]
    return cons, freq, len(rows), pfam


def main() -> int:
    fetch_pdb("1TUP", "A", CACHE)
    raw = (CACHE / "1TUP.pdb").read_text(encoding="utf-8", errors="replace")
    seq, nums = resnums(raw, "A")
    idx = {num: i for i, num in enumerate(nums)}
    cons, freq, nrows, pfam = conservation_profile(seq, "1TUP", pfam="PF00870")
    print(f"p53 1TUP/A  n={len(seq)}  Pfam {pfam}  MSA rows={nrows}  "
          f"mean conservation={cons.mean():.2f}\n")

    cons_scores = cons.copy()
    # SIFT-style substitution-specific impact: conserved position AND mutant residue disallowed
    var_scores = []
    for num, i in idx.items():
        for mut in AA20:
            if mut == seq[i]:
                continue
            fmut = freq[i].get(mut, 0.0) if i < len(freq) else 0.0
            var_scores.append(cons[i] * (1.0 - fmut))
    var_scores = np.array(var_scores)

    print(f"{'hotspot':<9}{'conservation':>13}{'cons-pctile':>12}{'variant-pctile':>15}")
    print("-" * 49)
    cpct, vpct = [], []
    for num, w, mut in HOTSPOTS:
        if num not in idx:
            continue
        i = idx[num]
        cp = float((cons_scores < cons[i]).mean()) * 100
        fmut = freq[i].get(mut, 0.0) if i < len(freq) else 0.0
        vimp = cons[i] * (1.0 - fmut)
        vp = float((var_scores < vimp).mean()) * 100
        cpct.append(cp)
        vpct.append(vp)
        print(f"{w}{num}{mut:<5}{cons[i]:>13.2f}{cp:>11.0f}%{vp:>14.0f}%")
    print("-" * 49)
    print(f"mean hotspot percentile:  conservation {np.mean(cpct):.1f}%   "
          f"variant {np.mean(vpct):.1f}%   (random 50%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
