#!/usr/bin/env python3
"""Template-based FSOT modeling: real homolog structures + coordinate transfer.

For each query chain, search RCSB for homolog structures, EXCLUDE the target PDB
(and, by default, near-identical redeposits above an identity cap), transfer the
best homolog's real Cα coordinates onto aligned positions, and model gap loops.
When no template exists, fall back to the zero-parameter bulk fold. This uses
real observed data as input to a zero-parameter map; no trained weights.
"""

from __future__ import annotations

import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from fsot_structure_engine import CA_CA, predict_ca_coords, target_rg_fsot  # noqa: E402
from run_fsot_vs_alphafold_structure import kabsch_rmsd, parse_pdb_ca  # noqa: E402
from run_rcsb_holdout import bootstrap_median, fetch_pdb, git_commit  # noqa: E402

MANIFEST_EVAL = ROOT / "data" / "rcsb_live_api_holdout_eval.json"
OUTPUT = ROOT / "data" / "rcsb_template_holdout_eval.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_live_api_holdout"
TCACHE = Path.home() / ".cache" / "fsot-genetics" / "template_pdb"
TCACHE.mkdir(parents=True, exist_ok=True)
IDENTITY_CAP = 0.95  # exclude near-identical redeposits of the same protein
MIN_IDENTITY = 0.45  # below this a homolog alignment is unreliable
MIN_COVERAGE = 0.65  # need the template to span most of the query
MAX_CANDIDATES = 8


def _post(url: str, body: dict, timeout: int = 60) -> dict:
    request = urllib.request.Request(
        url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "fsot"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read())


def _get(url: str, timeout: int = 40) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": "fsot", "Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read())


def pfam_family_pdbs(query_pdb: str) -> list[str]:
    """PDB structures in the query's Pfam family (remote structural homologs)."""
    try:
        u = _get(f"https://data.rcsb.org/rest/v1/core/uniprot/{query_pdb}/1")
        acc = u[0]["rcsb_uniprot_container_identifiers"]["uniprot_id"]
        fd = _get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/")
        fam = fd["results"][0]["metadata"]["accession"]
        sd = _get(f"https://www.ebi.ac.uk/interpro/api/structure/PDB/entry/pfam/{fam}/?page_size=100")
        out = []
        for row in sd.get("results", []):
            acc_pdb = row.get("metadata", {}).get("accession")
            if acc_pdb:
                out.append(acc_pdb.upper())
        return out
    except Exception:
        return []


def fetch_template_pdb(pdb_id: str) -> str:
    path = TCACHE / f"{pdb_id}.pdb"
    if path.exists():
        return path.read_text(encoding="utf-8", errors="replace")
    with urllib.request.urlopen(
        f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60
    ) as response:
        text = response.read().decode("utf-8", "replace")
    path.write_text(text, encoding="utf-8")
    return text


def chains_of(text: str) -> list[str]:
    order, seen = [], set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            chain = line[21]
            if chain not in seen:
                seen.add(chain)
                order.append(chain)
    return order


def nw_align(a: str, b: str) -> list[tuple[int, int]]:
    n, m = len(a), len(b)
    dp = np.zeros((n + 1, m + 1))
    dp[:, 0] = -np.arange(n + 1)
    dp[0, :] = -np.arange(m + 1)
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            match = dp[i - 1, j - 1] + (1.0 if ai == b[j - 1] else -1.0)
            dp[i, j] = max(match, dp[i - 1, j] - 1.0, dp[i, j - 1] - 1.0)
    i, j, pairs = n, m, []
    while i > 0 and j > 0:
        if dp[i, j] == dp[i - 1, j - 1] + (1.0 if a[i - 1] == b[j - 1] else -1.0):
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif dp[i, j] == dp[i - 1, j] - 1.0:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def homolog_ids(sequence: str) -> list[str]:
    query = {
        "query": {"type": "terminal", "service": "sequence", "parameters": {
            "evalue_cutoff": 0.1, "identity_cutoff": 0.25,
            "sequence_type": "protein", "value": sequence}},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": 30},
                            "results_content_type": ["experimental"]},
    }
    try:
        data = _post("https://search.rcsb.org/rcsbsearch/v2/query", query)
    except Exception:
        return []
    ids, seen = [], set()
    for hit in data.get("result_set", []):
        pdb = hit["identifier"].split("_")[0].upper()
        if pdb not in seen:
            seen.add(pdb)
            ids.append(pdb)
    return ids


def best_template(sequence: str, exclude_pdb: str) -> dict | None:
    best = None
    # Candidates from sequence search AND the query's Pfam family (remote homologs).
    candidates, seen = [], set()
    for pdb in homolog_ids(sequence) + pfam_family_pdbs(exclude_pdb):
        if pdb != exclude_pdb.upper() and pdb not in seen:
            seen.add(pdb)
            candidates.append(pdb)
    for pdb in candidates[:90]:
        try:
            text = fetch_template_pdb(pdb)
        except Exception:
            continue
        for chain in chains_of(text):
            try:
                tseq, tcoords = parse_pdb_ca(text, chain)
            except Exception:
                continue
            if len(tseq) < 20:
                continue
            pairs = nw_align(sequence, tseq)
            if len(pairs) < 10:
                continue
            identity = sum(1 for qi, ti in pairs if sequence[qi] == tseq[ti]) / len(pairs)
            coverage = len(pairs) / len(sequence)
            if identity > IDENTITY_CAP or identity < MIN_IDENTITY or coverage < MIN_COVERAGE:
                continue
            score = coverage * identity
            if best is not None and score <= best["score"]:
                continue
            model = build_from_template(len(sequence), tcoords, pairs)
            if not model_is_sane(model, len(sequence)):
                continue  # only keep templates whose transferred geometry is physical
            best = {"score": score, "pdb_id": pdb, "chain": chain, "model": model,
                    "identity": identity, "coverage": coverage}
        if best is not None and best["score"] > 0.85:
            break
    return best


def build_from_template(n: int, tcoords: np.ndarray, pairs: list[tuple[int, int]]) -> np.ndarray:
    coord = np.full((n, 3), np.nan)
    for qi, ti in pairs:
        coord[qi] = tcoords[ti]
    aligned = sorted(qi for qi, _ in pairs)
    for a, b in zip(aligned, aligned[1:]):
        if b > a + 1:
            pa, pb = coord[a], coord[b]
            for k, qi in enumerate(range(a + 1, b), 1):
                coord[qi] = pa + (pb - pa) * (k / (b - a))
    first, last = aligned[0], aligned[-1]
    if first > 0:
        step = (coord[aligned[1]] - coord[first]) if len(aligned) > 1 else np.array([CA_CA, 0.0, 0.0])
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(first - 1, -1, -1):
            coord[qi] = coord[qi + 1] - CA_CA * step
    if last < n - 1:
        step = (coord[last] - coord[aligned[-2]]) if len(aligned) > 1 else np.array([CA_CA, 0.0, 0.0])
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(last + 1, n):
            coord[qi] = coord[qi - 1] + CA_CA * step
    return coord - coord.mean(axis=0)


def model_is_sane(coord: np.ndarray, n: int) -> bool:
    """Intrinsic reject gate: catches assemblies, wrong chains, gap spikes."""
    rg = float(np.sqrt(((coord - coord.mean(axis=0)) ** 2).sum(axis=1).mean()))
    target = target_rg_fsot(n)
    bonds = np.linalg.norm(coord[1:] - coord[:-1], axis=1)
    frac_broken = float(np.mean(bonds > 5.0))
    return (0.5 * target < rg < 2.0 * target) and frac_broken < 0.1


def main() -> int:
    manifest = json.loads(MANIFEST_EVAL.read_text(encoding="utf-8"))
    rows = []
    for entry in manifest["results"]:
        pdb_id, chain = entry["pdb_id"], entry["chain"]
        text, _ = fetch_pdb(pdb_id, CACHE)
        sequence, native = parse_pdb_ca(text, chain)
        bulk = predict_ca_coords(sequence, canonicalize_chirality=True, observer_bulk_dim=25)
        bulk_rmsd = kabsch_rmsd(bulk["ca_coords"], native)
        template = best_template(sequence, pdb_id)
        if template is not None:
            template_rmsd = kabsch_rmsd(template["model"], native)
            meta = {"template_pdb": template["pdb_id"], "template_chain": template["chain"],
                    "identity": template["identity"], "coverage": template["coverage"]}
        else:
            template_rmsd = None
            meta = None
        rows.append({
            "pdb_id": pdb_id, "chain": chain, "length": len(sequence),
            "baseline_rmsd_A": entry["rmsd_A"]["baseline"],
            "bulk_rmsd_A": bulk_rmsd,
            "template_rmsd_A": template_rmsd,
            "template": meta,
            "best_rmsd_A": template_rmsd if template_rmsd is not None else bulk_rmsd,
        })
        tag = (f"template {meta['template_pdb']} id={meta['identity']:.2f} -> {template_rmsd:.2f}"
               if template_rmsd is not None else "no template -> bulk %.2f" % bulk_rmsd)
        print(f"{pdb_id}:{chain} n={len(sequence):3d} bulk={bulk_rmsd:5.2f}  {tag}")

    covered = [r for r in rows if r["template_rmsd_A"] is not None]
    uncovered = [r for r in rows if r["template_rmsd_A"] is None]
    aggregate = {
        "n_chains": len(rows),
        "n_template_covered": len(covered),
        "identity_cap": IDENTITY_CAP,
        "baseline_rmsd_A": bootstrap_median([r["baseline_rmsd_A"] for r in rows]),
        "bulk_rmsd_A": bootstrap_median([r["bulk_rmsd_A"] for r in rows]),
        "best_rmsd_A": bootstrap_median([r["best_rmsd_A"] for r in rows]),
        "template_covered_rmsd_A": bootstrap_median([r["template_rmsd_A"] for r in covered]) if covered else None,
        "template_free_bulk_rmsd_A": bootstrap_median([r["bulk_rmsd_A"] for r in uncovered]) if uncovered else None,
    }
    output = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "purpose": "Template-based FSOT modeling; real homolog data, zero trained weights",
        "candidate_commit": git_commit(),
        "controls": {"exclude_self_pdb": True, "identity_cap": IDENTITY_CAP,
                     "template_source": "RCSB sequence search v2"},
        "aggregate": aggregate,
        "results": rows,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(aggregate, indent=2))
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
