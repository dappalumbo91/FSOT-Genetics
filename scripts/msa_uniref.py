#!/usr/bin/env python3
"""Protein-specific deep MSA from UniRef clusters (real observables).

Why this exists
---------------
Broad Pfam MSAs mix distant family members (e.g. all small GTPases for KRAS),
which dilutes site conservation at classic driver positions (G12, etc.) and
leaves N-terminal / domain-edge columns unmapped.

UniRef50/90 clusters for the *query UniProt accession* are protein-specific
homolog sets published by UniProt — real sequence observables, not fitted
weights. We fetch cluster members, align each to the query, and compute
per-position frequencies + Shannon/identity conservation.

API: UniProt REST (credential-free)
  GET /uniref/search?query={acc}+AND+identity:0.5
  GET /uniref/{clusterId}/members?size=500  (+ Link pagination)
"""

from __future__ import annotations

import json
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT / "scripts"))

from run_rcsb_template_holdout import nw_align  # noqa: E402

CACHE = Path.home() / ".cache" / "fsot-genetics" / "uniref_msa"
CACHE.mkdir(parents=True, exist_ok=True)

AA20 = "ARNDCQEGHILKMFPSTWYV"
UA = "fsot-genetics-uniref/1.0"


def _http_json(url: str, timeout: int = 120) -> tuple[Any, dict]:
    req = urllib.request.Request(
        url, headers={"User-Agent": UA, "Accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        headers = {k.lower(): v for k, v in r.headers.items()}
        return json.loads(r.read().decode("utf-8", "replace")), headers


def _next_link(headers: dict) -> str | None:
    link = headers.get("link") or headers.get("Link") or ""
    # <url>; rel="next"
    m = re.search(r'<([^>]+)>;\s*rel="next"', link)
    return m.group(1) if m else None


def resolve_uniref_cluster(uniprot: str, identity: float = 0.5) -> str | None:
    """Return UniRef cluster id for accession (0.5 → UniRef50, 0.9 → UniRef90)."""
    # identity filter: 0.5 or 0.9
    id_tag = "0.5" if identity <= 0.55 else ("0.9" if identity <= 0.95 else "1.0")
    url = (
        f"https://rest.uniprot.org/uniref/search?"
        f"query={uniprot}+AND+identity:{id_tag}&format=json&size=5"
    )
    cache = CACHE / f"cluster_{uniprot}_{id_tag}.json"
    if cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
    else:
        data, _ = _http_json(url)
        cache.write_text(json.dumps(data), encoding="utf-8")
    results = data.get("results") or []
    if not results:
        return None
    # Prefer cluster whose seed/name matches accession when possible
    for r in results:
        rid = r.get("id") or ""
        if uniprot.upper() in rid.upper():
            return rid
    return results[0].get("id")


def _parse_fasta(txt: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    mid, chunks = None, []
    for line in txt.splitlines():
        if line.startswith(">"):
            if mid and chunks:
                seq = "".join(chunks)
                seq = "".join(c for c in seq.upper() if c in AA20)
                if len(seq) >= 20:
                    out.append((mid, seq))
            # >sp|P01116|RASK_HUMAN ...
            parts = line[1:].split("|")
            mid = parts[1] if len(parts) > 1 else line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if mid and chunks:
        seq = "".join(c for c in "".join(chunks).upper() if c in AA20)
        if len(seq) >= 20:
            out.append((mid, seq))
    return out


def fetch_uniref_sequences(
    cluster_id: str, *, max_members: int = 2000
) -> list[tuple[str, str]]:
    """Return list of (accession, sequence) from a UniRef cluster.

    Uses UniProtKB search ``uniref_cluster_{50|90|100}:{id}`` in FASTA form —
    the members JSON endpoint often omits sequences for non-seed entries.
    """
    cache = CACHE / f"fasta_{cluster_id.replace('/', '_')}_{max_members}.json"
    if cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
        return [(a, s) for a, s in data]

    # map UniRef50_P01116 → uniref_cluster_50:UniRef50_P01116
    if cluster_id.startswith("UniRef50_"):
        qfield = f"uniref_cluster_50:{cluster_id}"
    elif cluster_id.startswith("UniRef90_"):
        qfield = f"uniref_cluster_90:{cluster_id}"
    elif cluster_id.startswith("UniRef100_"):
        qfield = f"uniref_cluster_100:{cluster_id}"
    else:
        qfield = f"uniref_cluster_50:{cluster_id}"

    out: list[tuple[str, str]] = []
    # page with size=500; UniProt uses Link: rel="next"
    url = (
        f"https://rest.uniprot.org/uniprotkb/search?"
        f"query={qfield}&format=fasta&size=500"
    )
    while url and len(out) < max_members:
        req = urllib.request.Request(
            url, headers={"User-Agent": UA, "Accept": "text/plain"}
        )
        with urllib.request.urlopen(req, timeout=180) as r:
            headers = {k.lower(): v for k, v in r.headers.items()}
            txt = r.read().decode("utf-8", "replace")
        batch = _parse_fasta(txt)
        out.extend(batch)
        if len(out) >= max_members:
            out = out[:max_members]
            break
        url = _next_link(headers)
        if not batch:
            break
    cache.write_text(json.dumps(out), encoding="utf-8")
    return out


def _align_freq(
    query: str, homologs: list[str], *, max_rows: int = 1500
) -> tuple[np.ndarray, list[dict], int]:
    """Map homologs onto query via NW; return identity-to-query cons, freqs, n."""
    n = len(query)
    match = np.zeros(n)
    cover = np.zeros(n)
    colcount = [Counter() for _ in range(n)]
    if len(homologs) > max_rows:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(homologs), max_rows, replace=False)
        homologs = [homologs[i] for i in idx]
    used = 0
    for h in homologs:
        # length gate — protein-specific clusters are close; skip outliers
        if not (0.5 * n <= len(h) <= 1.8 * n):
            continue
        pairs = nw_align(query, h)
        if len(pairs) < max(10, int(0.3 * n)):
            continue
        used += 1
        for qi, hi in pairs:
            aa = h[hi]
            if aa not in AA20:
                continue
            cover[qi] += 1
            colcount[qi][aa] += 1
            if aa == query[qi]:
                match[qi] += 1
    cons = np.where(cover > 0, match / np.maximum(cover, 1), 0.0)
    freq = [{a: c[a] / max(sum(c.values()), 1) for a in c} for c in colcount]
    return cons, freq, used


def shannon_from_freq(freq_i: dict) -> float:
    if not freq_i:
        return 0.0
    vals = np.array([float(v) for v in freq_i.values() if v > 0], dtype=float)
    if vals.size == 0:
        return 0.0
    vals = vals / vals.sum()
    H = float(-(vals * np.log(vals)).sum())
    return max(0.0, 1.0 - H / math.log(20.0))


def conservation_profile_uniref(
    sequence: str,
    uniprot: str,
    *,
    identity: float = 0.5,
    max_members: int = 2000,
    max_align: int = 1200,
) -> tuple[np.ndarray, list[dict], int, str | None, dict]:
    """Per-position conservation from UniRef protein-specific homologs.

    Returns (identity_cons, freq, n_aligned, cluster_id, meta).
    """
    seq = "".join(c for c in sequence.upper() if c in AA20)
    meta: dict[str, Any] = {
        "backend": "uniref",
        "uniprot": uniprot,
        "identity": identity,
        "free_parameters": 0,
    }
    cluster = resolve_uniref_cluster(uniprot, identity=identity)
    if not cluster:
        meta["error"] = "no_cluster"
        return np.zeros(len(seq)), [{} for _ in range(len(seq))], 0, None, meta
    meta["cluster_id"] = cluster
    members = fetch_uniref_sequences(cluster, max_members=max_members)
    meta["n_members_fetched"] = len(members)
    # Always include the query itself once
    homologs = [seq] + [s for _id, s in members]
    cons, freq, n_used = _align_freq(seq, homologs, max_rows=max_align)
    meta["n_aligned"] = n_used
    meta["mean_identity_cons"] = float(cons.mean()) if len(cons) else 0.0
    meta["mean_shannon"] = float(
        np.mean([shannon_from_freq(f) for f in freq]) if freq else 0.0
    )
    meta["n_covered_positions"] = int((cons > 0).sum())
    return cons, freq, n_used, cluster, meta


def build_uniref_msa_features(
    sequence: str,
    uniprot: str,
    *,
    pfam: str | None = None,
    max_members: int = 600,
) -> "Any":
    """Build MsaFeatures from UniRef cluster for tertiary chem-link bridge."""
    from msa_pipeline import MsaFeatures, map_msa_to_query, mutual_information_apc  # noqa: WPS433

    seq = "".join(c for c in sequence.upper() if c in AA20)
    cons, freq, n, meta = best_conservation_profile(seq, uniprot=uniprot, pfam=pfam)
    cluster = meta.get("cluster_id") or resolve_uniref_cluster(uniprot, 0.5)
    rows: list[str] = [seq]
    if cluster:
        members = fetch_uniref_sequences(cluster, max_members=max_members)
        rows.extend(s for _i, s in members)
    coev = np.zeros((len(seq), len(seq)))
    covered = list(range(len(seq)))
    if len(rows) >= 20:
        mapped = map_msa_to_query(seq, rows)
        if mapped is not None:
            mat, covered, _ref = mapped
            if len(covered) >= 5:
                coev = mutual_information_apc(mat, covered)
    cons_a = np.asarray(cons, dtype=float)
    if len(cons_a) != len(seq):
        cons_a = np.zeros(len(seq))
    return MsaFeatures(
        sequence=seq,
        n_seqs=max(n, len(rows)),
        neff=float(max(n, 1)),
        backend="uniref",
        detail=str(cluster),
        conservation=cons_a,
        gap_frac=np.zeros(len(seq)),
        entropy=np.full(len(seq), math.log(20.0)),
        aa_freq=np.zeros((len(seq), 20)),
        coevolution=coev,
        covered=covered,
    )


def best_conservation_profile(
    sequence: str,
    uniprot: str | None = None,
    pfam: str | None = None,
    self_pdb: str | None = None,
) -> tuple[np.ndarray, list[dict], int, dict]:
    """Prefer UniRef protein-specific MSA; fall back to Pfam.

    Returns (cons, freq, n, meta) where cons is max(identity, shannon) blend
    stored separately in meta arrays when useful.
    """
    from variant_conservation import conservation_profile  # noqa: WPS433

    meta: dict[str, Any] = {"sources_tried": []}
    seq = "".join(c for c in sequence.upper() if c in AA20)

    if uniprot:
        for ident in (0.5, 0.9):
            cons, freq, n, cluster, m = conservation_profile_uniref(
                seq, uniprot, identity=ident
            )
            meta["sources_tried"].append(m)
            # Accept if enough depth and coverage of sequence
            cov = float((cons > 0).mean()) if len(cons) else 0.0
            if n >= 20 and cov >= 0.5:
                meta["chosen"] = "uniref"
                meta["cluster_id"] = cluster
                meta["uniref"] = m
                # blended conservation for callers that use cons as strength
                blend = np.array(
                    [
                        max(float(cons[i]), shannon_from_freq(freq[i]))
                        for i in range(len(seq))
                    ]
                )
                meta["identity_cons_mean"] = float(cons.mean())
                return blend, freq, n, meta

    if pfam or self_pdb:
        cons, freq, n, pf = conservation_profile(
            seq, self_pdb or "XXXX", pfam=pfam, max_rows=3000
        )
        meta["sources_tried"].append({"backend": "pfam", "pfam": pf, "n": n})
        meta["chosen"] = "pfam"
        meta["pfam"] = pf
        blend = np.array(
            [
                max(float(cons[i]), shannon_from_freq(freq[i] if i < len(freq) else {}))
                for i in range(len(seq))
            ]
        )
        return blend, freq, n, meta

    meta["chosen"] = "empty"
    return np.zeros(len(seq)), [{} for _ in range(len(seq))], 0, meta


if __name__ == "__main__":
    # Smoke: KRAS G12 should be highly conserved in UniRef50
    seq = (
        "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY"
        "SAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTK"
        "QAQDLARSYGIPFIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGCVKIKKCIIM"
    )
    cons, freq, n, cluster, meta = conservation_profile_uniref(seq, "P01116")
    print("cluster", cluster, "n", n, "meta", meta)
    for pos in (12, 13, 61):  # 1-based
        i = pos - 1
        print(
            f"  p.{pos}{seq[i]}  id={cons[i]:.2f}  shan={shannon_from_freq(freq[i]):.2f}  "
            f"top={sorted(freq[i].items(), key=lambda x: -x[1])[:4]}"
        )
