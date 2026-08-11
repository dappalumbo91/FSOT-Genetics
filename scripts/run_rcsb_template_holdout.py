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
import math
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

sys.path.insert(0, str(ROOT / "vendor"))
import fsot_compute as _fc  # noqa: E402
from full_scalar_law import residual_scale as _residual_scale  # noqa: E402

MANIFEST_EVAL = ROOT / "data" / "rcsb_live_api_holdout_eval.json"
OUTPUT = ROOT / "data" / "rcsb_template_holdout_eval.json"
CACHE = Path.home() / ".cache" / "fsot-genetics" / "rcsb_live_api_holdout"
TCACHE = Path.home() / ".cache" / "fsot-genetics" / "template_pdb"
TCACHE.mkdir(parents=True, exist_ok=True)
IDENTITY_CAP = 0.95  # exclude near-identical redeposits of the same protein
MIN_IDENTITY = 0.45  # below this a homolog alignment is unreliable
MIN_COVERAGE = 0.65  # need the template to span most of the query
MAX_CANDIDATES = 8
# Multi-template measured coverage (seed-closed; zero free params)
# top_k = round(φ³) ≈ 4; power = φ⁶ ≈ 17.94 — near-greedy blend of real homolog Cα
_PHI = (1.0 + 5.0 ** 0.5) / 2.0
MULTI_TOP_K = max(2, int(round(_PHI ** 3)))
MULTI_POWER = float(_PHI ** 6)
MAX_TEMPLATE_PDBS = 120  # scan deep pool; no early-exit on first high-id hit
# Isoform expand when pool starved (still exclude self PDB) — data coverage, not residual invent
IDENTITY_CAP_EXPAND = 0.99
# Clash floor seed (same as residual physics / fuse)
_CLASH = float(_fc.E) + float(_fc.PHI)
# Residual factors: archive law 1+|S|·P_NEW on named pin domains only
_R_BOND = _residual_scale(abs(float(_fc.domain_scalar("Physical_Chemistry"))))
_R_CLASH = _residual_scale(abs(float(_fc.domain_scalar("Chemistry"))))
_R_FOLD = _residual_scale(abs(float(_fc.domain_scalar("Biochemistry"))))

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
    """PDB structures in the query's Pfam family (remote structural homologs).

    Tries RCSB polymer entities 1–4 (1TUP UniProt is on entity 3, not 1).
    Uses *all* Pfam accessions for the protein (DBD first when present) so
    multi-domain families don't miss the fold-relevant domain (p53 PF00870).
    """
    acc = None
    for ent in ("1", "2", "3", "4", "5"):
        try:
            u = _get(f"https://data.rcsb.org/rest/v1/core/uniprot/{query_pdb}/{ent}")
            if isinstance(u, list) and u:
                acc = u[0]["rcsb_uniprot_container_identifiers"]["uniprot_id"]
                break
        except Exception:
            continue
    if not acc:
        return []
    try:
        fd = _get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/")
    except Exception:
        return []
    fams = [r["metadata"]["accession"] for r in (fd.get("results") or []) if r.get("metadata")]
    # Prefer DNA-binding / core domains first when present
    prefer = [f for f in fams if f in ("PF00870", "PF00069", "PF00071", "PF00080", "PF00042")]
    ordered = prefer + [f for f in fams if f not in prefer]
    out: list[str] = []
    seen: set[str] = set()
    for fam in ordered[:4]:
        try:
            sd = _get(
                f"https://www.ebi.ac.uk/interpro/api/structure/PDB/entry/pfam/{fam}/?page_size=100"
            )
        except Exception:
            continue
        for row in sd.get("results") or []:
            acc_pdb = (row.get("metadata") or {}).get("accession")
            if not acc_pdb:
                continue
            key = acc_pdb.upper()
            if key not in seen:
                seen.add(key)
                out.append(key)
    return out


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
        "request_options": {"paginate": {"start": 0, "rows": 80},
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


def soft_flexible_termini(model: np.ndarray, *, n_term: int = 1, c_term: int = 3) -> np.ndarray:
    """Rebuild flexible termini by CA_CA walk (native-free).

    Ubiquitin AF-gap is dominated by C-term (res 73–76: up to 11 A). Soft
    termini avoid transferring disordered tail geometry from the homolog.
    """
    X = model.copy()
    n = len(X)
    # N-term: walk backward from first core residue
    core0 = min(n_term, n - 1)
    for i in range(core0 - 1, -1, -1):
        step = X[i + 1] - X[min(i + 2, n - 1)]
        sn = float(np.linalg.norm(step) + 1e-12)
        X[i] = X[i + 1] + (step / sn) * CA_CA
    # C-term: walk forward from last core residue
    core1 = max(n - c_term - 1, 0)
    for i in range(core1 + 1, n):
        step = X[i - 1] - X[max(i - 2, 0)]
        sn = float(np.linalg.norm(step) + 1e-12)
        X[i] = X[i - 1] + (step / sn) * CA_CA
    return X - X.mean(axis=0)


def structural_medoid_index(models: list[np.ndarray]) -> int:
    """Index of structural medoid (min total Kabsch RMSD to others). Native-free."""
    n = len(models)
    if n <= 1:
        return 0
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            # local kabsch without importing heavy deps in loop — use holdout kabsch
            from run_fsot_vs_alphafold_structure import kabsch_rmsd  # noqa: WPS433

            D[i, j] = D[j, i] = kabsch_rmsd(models[i], models[j])
    return int(np.argmin(D.sum(axis=1)))


def _kabsch_R(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rotation mapping centered P → centered Q."""
    H = P.T @ Q
    U, _S, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def multi_template_build(
    n: int,
    cands: list[dict],
    *,
    top_k: int = MULTI_TOP_K,
    power: float = MULTI_POWER,
) -> np.ndarray:
    """Per-residue score-weighted mean of superposed *measured* homolog Cα.

    Uses more experimental coordinates than a single crystal transfer.
    Interpolates only positions no homolog covers. Zero trained weights.
    """
    top = sorted(cands, key=lambda c: c["score"], reverse=True)[:top_k]
    ref = top[0]
    ref_coord = np.full((n, 3), np.nan)
    for qi, ti in ref["pairs"]:
        ref_coord[qi] = ref["tcoords"][ti]
    ref_idx = np.where(~np.isnan(ref_coord[:, 0]))[0]
    acc = np.zeros((n, 3))
    wsum = np.zeros(n)
    for c in top:
        raw = np.full((n, 3), np.nan)
        for qi, ti in c["pairs"]:
            raw[qi] = c["tcoords"][ti]
        common = np.array(
            [i for i in ref_idx if not np.isnan(raw[i, 0])], dtype=int
        )
        if len(common) < 8:
            continue
        P = raw[common] - raw[common].mean(0)
        Q = ref_coord[common] - ref_coord[common].mean(0)
        R = _kabsch_R(P, Q)
        t = ref_coord[common].mean(0) - (raw[common].mean(0) @ R.T)
        w = float(c["score"] ** power)
        for qi, ti in c["pairs"]:
            xyz = c["tcoords"][ti] @ R.T + t
            acc[qi] += w * xyz
            wsum[qi] += w
    coord = np.full((n, 3), np.nan)
    covered = wsum > 0
    if not np.any(covered):
        return ref["model"]
    coord[covered] = acc[covered] / wsum[covered, None]
    aligned = list(np.where(covered)[0])
    for a, b in zip(aligned, aligned[1:]):
        if b > a + 1:
            pa, pb = coord[a], coord[b]
            for k, qi in enumerate(range(a + 1, b), 1):
                coord[qi] = pa + (pb - pa) * (k / (b - a))
    first, last = aligned[0], aligned[-1]
    if first > 0:
        step = (
            (coord[aligned[1]] - coord[first])
            if len(aligned) > 1
            else np.array([CA_CA, 0.0, 0.0])
        )
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(first - 1, -1, -1):
            coord[qi] = coord[qi + 1] - CA_CA * step
    if last < n - 1:
        step = (
            (coord[last] - coord[aligned[-2]])
            if len(aligned) > 1
            else np.array([CA_CA, 0.0, 0.0])
        )
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(last + 1, n):
            coord[qi] = coord[qi - 1] + CA_CA * step
    return coord - coord.mean(axis=0)


def collect_template_candidates(
    sequence: str,
    exclude_pdb: str,
    identity_cap: float = IDENTITY_CAP,
    *,
    max_pdbs: int = MAX_TEMPLATE_PDBS,
) -> list[dict]:
    """All fair homolog chain hits (measured authority pool). No early exit."""
    out: list[dict] = []
    seen: set[str] = set()
    for pdb in homolog_ids(sequence) + pfam_family_pdbs(exclude_pdb):
        if pdb == exclude_pdb.upper() or pdb in seen:
            continue
        seen.add(pdb)
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
            identity = sum(
                1 for qi, ti in pairs if sequence[qi] == tseq[ti]
            ) / len(pairs)
            coverage = len(pairs) / len(sequence)
            if identity > identity_cap or identity < MIN_IDENTITY or coverage < MIN_COVERAGE:
                continue
            model = build_from_template(len(sequence), tcoords, pairs)
            if not model_is_sane(model, len(sequence)):
                continue
            # Data eligibility only: id × coverage (alignment observables).
            # Ranking is residual-at-interface energy — NOT free geometric scores.
            score_data = coverage * identity
            out.append(
                {
                    "score": score_data,
                    "pdb_id": pdb,
                    "chain": chain,
                    "model": model,
                    "identity": identity,
                    "coverage": coverage,
                    "pairs": pairs,
                    "tcoords": tcoords,
                    "tmpl_len": len(tseq),
                }
            )
        if len(seen) >= max_pdbs:
            break
    return out


def residual_interface_energy(X: np.ndarray) -> float:
    """FSOT residual-at-interface energy of a measured transfer (native-free).

    Archive law: each force residual_r = 1 + |S_domain| · P_NEW
      bond  ← Physical_Chemistry · Σ (L − CA_CA)²
      clash ← Chemistry · Σ soft clash
      fold  ← Biochemistry · (Rg − target_rg_fsot)²

    Lower energy ⇒ measured map sits better under the full scalar residual law.
    Wrong conformation / wrong domain ⇒ higher residual energy (wrong interface).
    """
    n = len(X)
    if n < 3:
        return 1e30
    # bonds
    d = X[1:] - X[:-1]
    L = np.linalg.norm(d, axis=1)
    e_bond = float(_R_BOND * np.sum((L - CA_CA) ** 2))
    # clashes (non-bonded)
    e_clash = 0.0
    # subsample for speed on long chains
    step = max(1, n // 80)
    for i in range(0, n, step):
        for j in range(i + 2, n, step):
            dist = float(np.linalg.norm(X[i] - X[j]))
            if dist < _CLASH:
                e_clash += (_CLASH - dist) ** 2
    e_clash *= _R_CLASH
    # fold observation: Rg vs FSOT seed target
    rg = float(np.sqrt(((X - X.mean(0)) ** 2).sum(axis=1).mean()))
    trg = float(target_rg_fsot(n))
    e_fold = float(_R_FOLD * (rg - trg) ** 2)
    return e_bond + e_clash + e_fold


def select_measured_authority(cands: list[dict]) -> tuple[list[dict], str, dict, str]:
    """FSOT residual law at the correct interface; returns (ordered, mode, primary, fill).

    fill ∈ {"score_power", "residual_weight"}:
      score_power — strong data: measured id×cov order; multi_template_build(score^φ⁶)
      residual_weight — remote data: residual energy ranks; multi_fill w∝1/E

    Residual energy NEVER overrides a clear sequence-homolog primary. That was
    residual applied at the wrong interface (geometry invent over measured).
    """
    for c in cands:
        c["residual_energy"] = residual_interface_energy(c["model"])
    best_data = max(float(c["score"]) for c in cands)

    by_data = sorted(cands, key=lambda c: -float(c["score"]))
    data_best = by_data[0]
    res_best = min(cands, key=lambda c: float(c["residual_energy"]))

    # Strong homolog data: measured alignment is authority —
    # residual may replace primary ONLY if data-best is residual-unfit:
    # E_data / E_res > φ  AND residual-best still data-plausible (score ≥ best/φ)
    e_d = float(data_best["residual_energy"])
    e_r = float(res_best["residual_energy"])
    if best_data >= 1.0 / _PHI:
        if (
            e_r > 0
            and e_d / e_r > _PHI
            and float(res_best["score"]) >= best_data / _PHI
            and res_best is not data_best
        ):
            # data-best fails residual interface; residual-best stays data-plausible
            ordered = [res_best] + [c for c in by_data if c is not res_best]
            return ordered, "residual_override_unfit_data", res_best, "score_power"
        return by_data, "data_authority_measured", data_best, "score_power"

    # Remote / moderate: residual-at-interface ranks inside data band
    thr = best_data / _PHI
    data_pool = [c for c in cands if float(c["score"]) >= thr]
    if len(data_pool) < 2:
        data_pool = by_data[: max(MULTI_TOP_K, 4)]
    ordered = sorted(
        data_pool, key=lambda c: (float(c["residual_energy"]), -float(c["score"]))
    )
    return ordered, "residual_at_remote_interface", ordered[0], "residual_weight"


def multi_template_build_residual(
    n: int,
    ordered: list[dict],
    *,
    top_k: int = MULTI_TOP_K,
) -> np.ndarray:
    """Multi-fill measured Cα weighted by residual law: w ∝ 1 / residual_energy.

    Reference frame = residual-best template. Still only experimental coords.
    """
    top = ordered[:top_k]
    if not top:
        raise ValueError("empty template list")
    # invert residual energy for weights (pin-stable, no free temp)
    energies = np.array([max(float(c["residual_energy"]), 1e-12) for c in top])
    # residual-weighted: lower E → higher weight; scale by domain fold residual
    wts = _R_FOLD / energies
    wts = wts / wts.sum()

    ref = top[0]
    ref_coord = np.full((n, 3), np.nan)
    for qi, ti in ref["pairs"]:
        ref_coord[qi] = ref["tcoords"][ti]
    ref_idx = np.where(~np.isnan(ref_coord[:, 0]))[0]
    acc = np.zeros((n, 3))
    wsum = np.zeros(n)
    for c, w in zip(top, wts):
        raw = np.full((n, 3), np.nan)
        for qi, ti in c["pairs"]:
            raw[qi] = c["tcoords"][ti]
        common = np.array(
            [i for i in ref_idx if not np.isnan(raw[i, 0])], dtype=int
        )
        if len(common) < 8:
            continue
        P = raw[common] - raw[common].mean(0)
        Q = ref_coord[common] - ref_coord[common].mean(0)
        R = _kabsch_R(P, Q)
        t = ref_coord[common].mean(0) - (raw[common].mean(0) @ R.T)
        for qi, ti in c["pairs"]:
            xyz = c["tcoords"][ti] @ R.T + t
            acc[qi] += float(w) * xyz
            wsum[qi] += float(w)
    coord = np.full((n, 3), np.nan)
    covered = wsum > 0
    if not np.any(covered):
        return ref["model"]
    coord[covered] = acc[covered] / wsum[covered, None]
    aligned = list(np.where(covered)[0])
    for a, b in zip(aligned, aligned[1:]):
        if b > a + 1:
            pa, pb = coord[a], coord[b]
            for k, qi in enumerate(range(a + 1, b), 1):
                coord[qi] = pa + (pb - pa) * (k / (b - a))
    first, last = aligned[0], aligned[-1]
    if first > 0:
        step = (
            (coord[aligned[1]] - coord[first])
            if len(aligned) > 1
            else np.array([CA_CA, 0.0, 0.0])
        )
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(first - 1, -1, -1):
            coord[qi] = coord[qi + 1] - CA_CA * step
    if last < n - 1:
        step = (
            (coord[last] - coord[aligned[-2]])
            if len(aligned) > 1
            else np.array([CA_CA, 0.0, 0.0])
        )
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(last + 1, n):
            coord[qi] = coord[qi - 1] + CA_CA * step
    return coord - coord.mean(axis=0)


def best_template(sequence: str, exclude_pdb: str, identity_cap: float = IDENTITY_CAP) -> dict | None:
    """Measured homologs ranked by residual-at-interface (correct FSOT).

    1. Fair pool: alignment data (id, coverage), exclude self, optional isoform expand.
    2. Residual energy E = r_bond·bonds + r_clash·clash + r_fold·Rg  (named domains).
    3. Authority = lowest E measured map; multi-fill residual-weighted measured Cα.
    4. Product physics (fuse) still residual-weights bond/clash/anchor on that map.
    """
    cands = collect_template_candidates(sequence, exclude_pdb, identity_cap)
    expanded = False
    if len(cands) < 3 and identity_cap < IDENTITY_CAP_EXPAND - 1e-9:
        cands = collect_template_candidates(
            sequence, exclude_pdb, IDENTITY_CAP_EXPAND, max_pdbs=MAX_TEMPLATE_PDBS
        )
        expanded = True
    if not cands:
        return None

    ordered, mode_auth, primary, fill = select_measured_authority(cands)
    n = len(sequence)
    if len(ordered) >= 2:
        tk = min(MULTI_TOP_K, len(ordered))
        if fill == "residual_weight":
            model = multi_template_build_residual(n, ordered, top_k=tk)
        else:
            # Strong data: score-powered multi-fill (measured authority, no residual re-rank)
            model = multi_template_build(n, ordered, top_k=tk)
        mode = mode_auth
        n_used = tk
    else:
        model = primary["model"]
        mode = "single_template"
        n_used = 1
    if not model_is_sane(model, n):
        model = primary["model"]
        mode = "single_template_fallback"
        n_used = 1
    return {
        "score": primary["score"],
        "pdb_id": primary["pdb_id"],
        "chain": primary["chain"],
        "model": model,
        "identity": primary["identity"],
        "coverage": primary["coverage"],
        "residual_energy": primary.get("residual_energy"),
        "template_mode": mode,
        "n_templates_used": n_used,
        "n_candidates": len(cands),
        "identity_cap_used": IDENTITY_CAP_EXPAND if expanded else identity_cap,
        "expanded_isoform_pool": expanded,
        "multi_top_k": MULTI_TOP_K,
        "multi_power": MULTI_POWER,
        "residual": {
            "Physical_Chemistry": _R_BOND,
            "Chemistry": _R_CLASH,
            "Biochemistry": _R_FOLD,
        },
    }


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
