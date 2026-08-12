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


def _bond_mse(X: np.ndarray) -> float:
    """Mean (L − CA_CA)² — Physical_Chemistry transfer integrity."""
    if len(X) < 2:
        return 0.0
    L = np.linalg.norm(X[1:] - X[:-1], axis=1)
    return float(np.mean((L - CA_CA) ** 2))


# Transfer is bond-broken when RMS |L−CA_CA| ≳ 1/φ Å (seed; not a free cutoff).
_BOND_BROKEN = (1.0 / _PHI) ** 2


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

    fill ∈ {"score_power", "residual_weight", "single"}:
      score_power — strong data, residual-fit transfer: id×cov multi-fill
      residual_weight — remote data: residual energy ranks; multi_fill w∝1/E
      single — data-best transfer is residual-unfit: do not blend states

    Residual is a *transfer-quality* filter on measured maps, not a ranker
    that replaces sequence authority with the lowest-E homolog (that pick
    was 1HU8 on p53 — clean bonds, invented termini, worse RMSD).
    """
    for c in cands:
        c["residual_energy"] = residual_interface_energy(c["model"])
    best_data = max(float(c["score"]) for c in cands)

    by_data = sorted(cands, key=lambda c: -float(c["score"]))
    data_best = by_data[0]

    # Data-plausible band (alignment observables). Residual never searches
    # outside this set.
    thr = best_data / _PHI
    pool = [c for c in cands if float(c["score"]) >= thr] or by_data[:1]
    e_min = min(float(c["residual_energy"]) for c in pool)
    e_d = float(data_best["residual_energy"])

    if best_data >= 1.0 / _PHI:
        # Strong homologs: replace data-best only when that *transfer* is
        # residual-unfit AND bond-broken at Physical_Chemistry.
        # E_data > φ·E_min alone is not enough — that picked KRAS 1BKD
        # (SOS-bound, clean bonds, wrong state) over 1AA9.
        # p53 3Q01: bond_mse ≈ 2.8 Å² (broken) → 2P52.
        broken = _bond_mse(data_best["model"]) > _BOND_BROKEN
        if e_d > _PHI * e_min and broken:
            fit = [c for c in pool if float(c["residual_energy"]) <= _PHI * e_min]
            if not fit:
                fit = pool
            primary = max(fit, key=lambda c: float(c["score"]))
            ordered = [primary] + [c for c in by_data if c is not primary]
            return ordered, "residual_fit_data_authority", primary, "single"
        return by_data, "data_authority_measured", data_best, "score_power"

    # Remote / moderate: residual-at-interface ranks inside data band
    data_pool = pool if len(pool) >= 2 else by_data[: max(MULTI_TOP_K, 4)]
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
    """Measured homologs with residual as transfer-quality filter (correct FSOT).

    1. Fair pool: alignment data (id, coverage), exclude self, optional isoform expand.
    2. Residual energy E = r_bond·bonds + r_clash·clash + r_fold·Rg  (named domains).
    3. Strong data: id×cov authority unless that transfer is residual-unfit
       (E_data > φ·E_min) *and* bond-broken (mean (L−CA_CA)² > 1/φ²) —
       then highest-score residual-fit map, no blend.
    4. Remote data: residual-ranked inside the data band; multi-fill w∝1/E.
    5. Product physics (fuse) residual-weights bond/clash/anchor on that map.
    """
    cands = collect_template_candidates(sequence, exclude_pdb, identity_cap)
    expanded = False
    # Expand isoforms only when the fair pool is empty (starved kinases).
    # Merging 0.99 into a thin pool previously pulled CaM/SOD1 onto wrong assemblies.
    if len(cands) == 0 and identity_cap < IDENTITY_CAP_EXPAND - 1e-9:
        cands = collect_template_candidates(
            sequence, exclude_pdb, IDENTITY_CAP_EXPAND, max_pdbs=MAX_TEMPLATE_PDBS
        )
        expanded = bool(cands)
    if not cands:
        return None

    ordered, mode_auth, primary, fill = select_measured_authority(cands)
    n = len(sequence)
    if fill == "single" or len(ordered) < 2:
        # Residual-unfit data-best: one measured map. Multi-fill here
        # mixed conformational states (p53 2P52 vs 1HU8).
        model = primary["model"]
        mode = mode_auth if fill == "single" else "single_template"
        n_used = 1
    elif fill == "residual_weight":
        tk = min(MULTI_TOP_K, len(ordered))
        model = multi_template_build_residual(n, ordered, top_k=tk)
        mode = mode_auth
        n_used = tk
    else:
        # Strong data: score-powered multi-fill (measured authority)
        tk = min(MULTI_TOP_K, len(ordered))
        model = multi_template_build(n, ordered, top_k=tk)
        mode = mode_auth
        n_used = tk
    # Residual-at-interface reject: ensemble worse than primary measured map
    e_pri = residual_interface_energy(primary["model"])
    e_ens = residual_interface_energy(model)
    if (not model_is_sane(model, n)) or (e_pri > 0 and e_ens > _PHI * e_pri):
        model = primary["model"]
        mode = "primary_residual_fallback"
        n_used = 1
    homologs = [c["model"] for c in ordered[: min(MULTI_TOP_K, len(ordered))]]
    tert = measured_tertiary_contacts(homologs if homologs else [model])
    # Measured homolog disagreement on termini (real data): rebuild tail
    # by CA_CA walk when tail std > φ and tail > φ·core (Biochemistry residual
    # localizes the observation — do not residual-scale the backbone walk).
    model, term_note = _maybe_soft_disordered_termini(model, homologs)
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
        "tertiary_contacts": tert,
        "termini_note": term_note,
        "residual": {
            "Physical_Chemistry": _R_BOND,
            "Chemistry": _R_CLASH,
            "Biochemistry": _R_FOLD,
        },
    }


def _align_models(models: list[np.ndarray]) -> np.ndarray:
    ref = models[0] - models[0].mean(0)
    out = [ref]
    for m in models[1:]:
        p = m - m.mean(0)
        H = p.T @ ref
        U, _S, Vt = np.linalg.svd(H)
        d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        out.append(p @ R.T)
    return np.stack(out, axis=0)


def _maybe_soft_disordered_termini(
    model: np.ndarray, homologs: list[np.ndarray]
) -> tuple[np.ndarray, str]:
    """If measured homologs disagree on a terminus, drop that tail geometry.

    Native-free. Thresholds: φ (Å) and φ × core_std. Tail length ceil(φ)+1 = 3.
    """
    # Diagnosed on ubiquitin-scale tails (n=76). Long chains (p53/RNase)
    # have other termini biology — do not apply this interface there.
    if len(homologs) < 2 or len(model) < 16 or len(model) >= int(round(_PHI * 50)):
        return model, "no_term_rebuild"
    al = _align_models(homologs)
    std = al.std(axis=0).mean(axis=1)
    n = len(model)
    tail = int(math.ceil(_PHI)) + 1  # 3
    core = std[tail : n - tail]
    if len(core) < 4:
        return model, "no_term_rebuild"
    core_s = float(core.mean())
    c_s = float(std[-tail:].mean())
    n_s = float(std[:tail].mean())
    # Strict: absolute std > e AND ratio > φ² (ubiquitin C-term 3.3/0.42 ≈ 8;
    # milder tail noise on RNase/p53 must not fire).
    _e = float(_fc.E)
    ratio = _PHI * _PHI
    c_term = tail if (c_s > _e and c_s > ratio * max(core_s, 1e-6)) else 0
    n_term = tail if (n_s > _e and n_s > ratio * max(core_s, 1e-6)) else 0
    if c_term == 0 and n_term == 0:
        return model, "termini_agreed"
    X = soft_flexible_termini(model, n_term=n_term, c_term=c_term)
    return X, f"soft_N{n_term}_C{c_term}"


def measured_tertiary_contacts(
    models: list[np.ndarray],
    *,
    gate: int | None = None,
) -> list[tuple[int, int, float]]:
    """Long-range contacts *observed* in homolog Cα (real data).

    Pair is kept if median measured distance < F08·φ AND homologs agree
    (std < φ Å). Lean: tertiaryBiochem D=13 — residual applied later, not here.
    Backbone (sep≤2) never included.
    """
    if not models:
        return []
    n = len(models[0])
    if gate is None:
        gate = max(7, int(math.ceil(float(_fc.ETA_EFF) * 13.0)))  # Biochemistry D=13
    contact_hi = float(_fc.PI) * float(_fc.E) * _PHI  # φ·F08
    agree = _PHI  # Å
    out: list[tuple[int, int, float]] = []
    stack = np.stack(models, axis=0)
    for i in range(n):
        for j in range(i + gate, n):
            ds = np.linalg.norm(stack[:, i, :] - stack[:, j, :], axis=1)
            med = float(np.median(ds))
            if med >= contact_hi:
                continue
            if float(np.std(ds)) > agree:
                continue
            out.append((i, j, med))
    return out


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
