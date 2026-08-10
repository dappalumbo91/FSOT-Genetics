#!/usr/bin/env python3
"""MSA front-end for FSOT-Genetics: generation path + evolutionary features.

Design (zero free parameters):
  - Core F01–F15 engine stays single-sequence by default.
  - MSA is *optional input data* (like templates), never trained weights.
  - Features are closed-form statistics of the alignment:
      conservation (1 - H/log20), gap fraction, MI+APC coevolution,
      per-position AA frequency, effective depth (Neff).
  - Injection amplitudes reuse domain scalars already in F09
      (evo_amp = |S_biochem| · P_NEW · C_EFF / φ), same family as region_amp.

Generation backends (first available wins):
  1. Precomputed MSA file (Stockholm / A3M / FASTA)
  2. Local jackhmmer / hhblits against a user database (if on PATH)
  3. Public Pfam full alignment via InterPro (network; proven in this repo)
  4. Empty → features are zeros; mode falls back to pure single-sequence

Authority: FSOT seeds only for scales; evolutionary signal is data, not fit.
"""

from __future__ import annotations

import gzip
import math
import os
import shutil
import subprocess
import tempfile
import urllib.request
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
P_NEW = float(fc.P_NEW)
C_EFF = float(fc.C_EFF)

try:
    S_BIOCHEM = abs(float(fc.domain_scalar("Biochemistry")))
except Exception:
    S_BIOCHEM = 0.3

# F09-family amplitude for evolutionary channel (seed-closed, not free).
EVO_AMP = S_BIOCHEM * P_NEW * C_EFF / PHI

AA20 = "ARNDCQEGHILKMFPSTWYV"
A2I = {a: i for i, a in enumerate(AA20)}
GAP = 20  # index for gap/other in 21-state alphabet
ALPH21 = AA20 + "-"

MSA_CACHE = Path.home() / ".cache" / "fsot-genetics" / "msa"
MSA_CACHE.mkdir(parents=True, exist_ok=True)

# ── parsers ───────────────────────────────────────────────────────────────


def parse_stockholm(txt: str) -> list[str]:
    rows: dict[str, str] = {}
    for line in txt.splitlines():
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) == 2:
            rows[parts[0]] = rows.get(parts[0], "") + parts[1]
    return list(rows.values())


def parse_a3m_or_fasta(txt: str) -> list[str]:
    """Parse A3M or FASTA. A3M lowercase = insertions relative to query; drop them."""
    seqs: list[str] = []
    cur: list[str] = []
    for line in txt.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
            continue
        # A3M: uppercase + gaps only for match states
        cleaned = "".join(ch for ch in line if not ch.islower())
        cur.append(cleaned)
    if cur:
        seqs.append("".join(cur))
    return seqs


def parse_msa_text(txt: str, hint: str | None = None) -> list[str]:
    head = (txt[:200] if txt else "").lstrip()
    kind = (hint or "").lower()
    if kind in ("sto", "stockholm") or head.startswith("# STOCKHOLM"):
        return parse_stockholm(txt)
    if kind in ("a3m", "fasta", "fa", "fas") or head.startswith(">"):
        return parse_a3m_or_fasta(txt)
    # heuristic
    if head.startswith("# STOCKHOLM"):
        return parse_stockholm(txt)
    if ">" in head[:80]:
        return parse_a3m_or_fasta(txt)
    return parse_stockholm(txt) or parse_a3m_or_fasta(txt)


def load_msa_file(path: Path | str) -> list[str]:
    path = Path(path)
    raw = path.read_bytes()
    try:
        txt = gzip.decompress(raw).decode("utf-8", "replace")
    except Exception:
        txt = raw.decode("utf-8", "replace")
    suf = path.suffix.lower().lstrip(".")
    if suf == "gz":
        suf = path.stem.split(".")[-1].lower()
    return parse_msa_text(txt, hint=suf)


# ── generation backends ───────────────────────────────────────────────────


def _have(cmd: str) -> str | None:
    return shutil.which(cmd)


def fetch_pfam_msa(pfam: str, kind: str = "full") -> list[str]:
    """Public Pfam full/seed alignment via InterPro (same path as variant_conservation)."""
    fp = MSA_CACHE / f"{pfam}.{kind}.sto"
    if not fp.exists():
        url = (
            f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam}/"
            f"?annotation=alignment:{kind}"
        )
        req = urllib.request.Request(url, headers={"User-Agent": "fsot-genetics-msa"})
        with urllib.request.urlopen(req, timeout=180) as r:
            raw = r.read()
        try:
            txt = gzip.decompress(raw).decode("utf-8", "replace")
        except Exception:
            txt = raw.decode("utf-8", "replace")
        fp.write_text(txt, encoding="utf-8")
    return parse_stockholm(fp.read_text(encoding="utf-8", errors="replace"))


def pfam_for_uniprot(uniprot: str) -> str | None:
    try:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{uniprot}/"
        req = urllib.request.Request(url, headers={"User-Agent": "fsot-genetics-msa"})
        with urllib.request.urlopen(req, timeout=60) as r:
            import json

            data = json.loads(r.read().decode("utf-8", "replace"))
        results = data.get("results") or []
        if results:
            return results[0]["metadata"]["accession"]
    except Exception:
        return None
    return None


def run_jackhmmer(
    query_fasta: Path,
    db: Path,
    *,
    n_iter: int = 2,
    evalue: float = 1e-3,
    cpu: int = 2,
) -> list[str] | None:
    """Optional local JackHMMER → Stockholm. Returns None if tool/db missing."""
    bin_ = _have("jackhmmer")
    if not bin_ or not db.exists():
        return None
    with tempfile.TemporaryDirectory(prefix="fsot_jh_") as td:
        out = Path(td) / "out.sto"
        cmd = [
            bin_,
            "-N",
            str(n_iter),
            "-E",
            str(evalue),
            "--cpu",
            str(cpu),
            "-A",
            str(out),
            "--noali",
            str(query_fasta),
            str(db),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=600)
        except (subprocess.SubprocessError, OSError, TimeoutError):
            return None
        if not out.exists():
            return None
        return parse_stockholm(out.read_text(encoding="utf-8", errors="replace"))


def run_hhblits(
    query_fasta: Path,
    db_prefix: Path,
    *,
    n_iter: int = 2,
    evalue: float = 1e-3,
    cpu: int = 2,
) -> list[str] | None:
    """Optional local HHblits → A3M. db_prefix is the HH-suite database stem."""
    bin_ = _have("hhblits")
    if not bin_:
        return None
    # HH-suite DBs are multi-file; require at least one sidecar
    if not any(Path(str(db_prefix) + suf).exists() for suf in ("", "_a3m.ffdata", ".cs219")):
        return None
    with tempfile.TemporaryDirectory(prefix="fsot_hh_") as td:
        out = Path(td) / "out.a3m"
        cmd = [
            bin_,
            "-i",
            str(query_fasta),
            "-d",
            str(db_prefix),
            "-oa3m",
            str(out),
            "-n",
            str(n_iter),
            "-e",
            str(evalue),
            "-cpu",
            str(cpu),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=900)
        except (subprocess.SubprocessError, OSError, TimeoutError):
            return None
        if not out.exists():
            return None
        return parse_a3m_or_fasta(out.read_text(encoding="utf-8", errors="replace"))


def write_fasta(path: Path, seq: str, name: str = "query") -> None:
    path.write_text(f">{name}\n{seq}\n", encoding="utf-8")


@dataclass
class MsaSource:
    rows: list[str]
    backend: str  # file | jackhmmer | hhblits | pfam | empty
    detail: str = ""


def obtain_msa(
    sequence: str,
    *,
    msa_path: Path | str | None = None,
    pfam: str | None = None,
    uniprot: str | None = None,
    jackhmmer_db: Path | str | None = None,
    hhblits_db: Path | str | None = None,
    prefer_local: bool = True,
) -> MsaSource:
    """Resolve an MSA for *sequence* using the first available backend."""
    seq = "".join(c for c in sequence.upper() if c in AA20)
    if msa_path is not None:
        rows = load_msa_file(msa_path)
        return MsaSource(rows=rows, backend="file", detail=str(msa_path))

    if prefer_local:
        with tempfile.TemporaryDirectory(prefix="fsot_msa_q_") as td:
            fa = Path(td) / "q.fasta"
            write_fasta(fa, seq)
            if hhblits_db:
                rows = run_hhblits(fa, Path(hhblits_db))
                if rows:
                    return MsaSource(rows=rows, backend="hhblits", detail=str(hhblits_db))
            if jackhmmer_db:
                rows = run_jackhmmer(fa, Path(jackhmmer_db))
                if rows:
                    return MsaSource(
                        rows=rows, backend="jackhmmer", detail=str(jackhmmer_db)
                    )

    pfam_id = pfam
    if not pfam_id and uniprot:
        pfam_id = pfam_for_uniprot(uniprot)
    if pfam_id:
        try:
            rows = fetch_pfam_msa(pfam_id)
            if rows:
                return MsaSource(rows=rows, backend="pfam", detail=pfam_id)
        except Exception as exc:
            return MsaSource(rows=[], backend="empty", detail=f"pfam_fail:{exc}")

    return MsaSource(rows=[], backend="empty", detail="no_msa_source")


# ── map MSA columns → query ───────────────────────────────────────────────


def _nw_align(a: str, b: str) -> list[tuple[int, int]]:
    """Lightweight Needleman–Wunsch (same spirit as run_rcsb_template_holdout)."""
    # Prefer shared implementation when present
    try:
        from run_rcsb_template_holdout import nw_align  # noqa: WPS433

        return nw_align(a, b)
    except Exception:
        pass
    na, nb = len(a), len(b)
    gap = -1
    match = 1
    mismatch = -1
    H = np.zeros((na + 1, nb + 1), dtype=np.int32)
    ptr = np.zeros((na + 1, nb + 1), dtype=np.int8)  # 1=diag 2=up 3=left
    for i in range(1, na + 1):
        H[i, 0] = i * gap
        ptr[i, 0] = 2
    for j in range(1, nb + 1):
        H[0, j] = j * gap
        ptr[0, j] = 3
    for i in range(1, na + 1):
        for j in range(1, nb + 1):
            s = match if a[i - 1] == b[j - 1] else mismatch
            diag = H[i - 1, j - 1] + s
            up = H[i - 1, j] + gap
            left = H[i, j - 1] + gap
            best = diag
            p = 1
            if up > best:
                best, p = up, 2
            if left > best:
                best, p = left, 3
            H[i, j] = best
            ptr[i, j] = p
    i, j = na, nb
    pairs: list[tuple[int, int]] = []
    while i > 0 or j > 0:
        p = ptr[i, j]
        if p == 1:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif p == 2:
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def map_msa_to_query(
    sequence: str, rows: Sequence[str], *, max_ref_scan: int = 400
) -> tuple[np.ndarray, list[int], str] | None:
    """Return (n_seq × L_query int8 matrix, covered query indices, ref_row).

    Matrix alphabet: 0..19 = AA20, 20 = gap/other.
    """
    seq = "".join(c for c in sequence.upper() if c in AA20)
    if not rows or not seq:
        return None
    best = None
    for r in list(rows)[:max_ref_scan]:
        ung = r.replace(".", "").replace("-", "").upper()
        if not (0.5 * len(seq) <= len(ung) <= 2.0 * len(seq)):
            continue
        pairs = _nw_align(seq, ung)
        score = sum(1 for a, b in pairs if seq[a] == ung[b])
        if best is None or score > best[0]:
            best = (score, r, pairs)
        if best[0] > 0.9 * len(seq):
            break
    if best is None:
        return None
    _score, ref, pairs = best
    ung_to_col = [c for c, ch in enumerate(ref) if ch not in ".-"]
    ref_ung = ref.replace(".", "").replace("-", "").upper()
    qpairs = dict(_nw_align(seq, ref_ung))
    n_rows = min(len(rows), 12000)  # hard cap for memory
    mat = np.full((n_rows, len(seq)), GAP, dtype=np.int8)
    listrows = [list(r) for r in rows[:n_rows]]
    for qi, ui in qpairs.items():
        if ui >= len(ung_to_col):
            continue
        col = ung_to_col[ui]
        for ri, r in enumerate(listrows):
            if col < len(r):
                mat[ri, qi] = A2I.get(r[col].upper(), GAP)
    covered = sorted(qpairs.keys())
    return mat, covered, ref


# ── feature extraction ────────────────────────────────────────────────────


@dataclass
class MsaFeatures:
    """Per-residue and pairwise evolutionary features aligned to the query."""

    sequence: str
    n_seqs: int
    neff: float
    backend: str
    detail: str
    conservation: np.ndarray  # (L,) in [0,1]  1=invariant
    gap_frac: np.ndarray  # (L,)
    entropy: np.ndarray  # (L,) bits
    aa_freq: np.ndarray  # (L, 20)
    coevolution: np.ndarray  # (L, L) MI-APC, non-negative after clip
    covered: list[int] = field(default_factory=list)

    @property
    def depth_ok(self) -> bool:
        return self.n_seqs >= int(math.ceil(PI * E)) and self.neff >= PHI

    def summary(self) -> dict[str, Any]:
        return {
            "backend": self.backend,
            "detail": self.detail,
            "n_seqs": self.n_seqs,
            "neff": float(self.neff),
            "length": len(self.sequence),
            "mean_conservation": float(self.conservation.mean()) if len(self.conservation) else 0.0,
            "mean_gap_frac": float(self.gap_frac.mean()) if len(self.gap_frac) else 0.0,
            "n_covered": len(self.covered),
            "depth_ok": self.depth_ok,
            "evo_amp": float(EVO_AMP),
            "free_parameters": 0,
        }


def empty_features(sequence: str, *, backend: str = "empty", detail: str = "") -> MsaFeatures:
    seq = "".join(c for c in sequence.upper() if c in AA20)
    n = len(seq)
    return MsaFeatures(
        sequence=seq,
        n_seqs=0,
        neff=0.0,
        backend=backend,
        detail=detail,
        conservation=np.zeros(n),
        gap_frac=np.ones(n),
        entropy=np.full(n, math.log(20.0)),
        aa_freq=np.zeros((n, 20)),
        coevolution=np.zeros((n, n)),
        covered=[],
    )


def _column_stats(mat: np.ndarray, covered: list[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """conservation, gap_frac, entropy, aa_freq for length = mat.shape[1]."""
    L = mat.shape[1]
    conservation = np.zeros(L)
    gap_frac = np.ones(L)
    entropy = np.full(L, math.log(20.0))
    aa_freq = np.zeros((L, 20))
    log20 = math.log(20.0)
    for j in covered:
        col = mat[:, j]
        nongap = col[col != GAP]
        if nongap.size == 0:
            continue
        gap_frac[j] = 1.0 - (nongap.size / max(mat.shape[0], 1))
        counts = np.bincount(nongap.astype(np.int64), minlength=20).astype(float)
        tot = counts.sum()
        if tot <= 0:
            continue
        p = counts / tot
        aa_freq[j] = p
        nz = p > 0
        H = float(-(p[nz] * np.log(p[nz])).sum())
        entropy[j] = H
        conservation[j] = max(0.0, 1.0 - H / log20)
    return conservation, gap_frac, entropy, aa_freq


def _neff(mat: np.ndarray, id_thresh: float | None = None) -> float:
    """Effective sequence count (diversity-aware).

    Two seed-closed pieces, then take the max (honest depth, not free fit):
      1) Greedy identity clusters at threshold 1/e ≈ 0.37 (classic MSA cut).
      2) Analytic N_eff ≈ N / (1 + (N-1)·mean_id) from a random pair sample.
    """
    if id_thresh is None:
        id_thresh = 1.0 / E
    n = mat.shape[0]
    if n == 0:
        return 0.0
    if n == 1:
        return 1.0
    rng = np.random.default_rng(0)
    work = mat
    n_work = n
    if n > 2000:
        idx = rng.choice(n, 2000, replace=False)
        work = mat[idx]
        n_work = work.shape[0]

    # (1) greedy clusters
    used = np.zeros(n_work, dtype=bool)
    n_clusters = 0
    for i in range(n_work):
        if used[i]:
            continue
        n_clusters += 1
        used[i] = True
        ai = work[i]
        for j in range(i + 1, n_work):
            if used[j]:
                continue
            bj = work[j]
            valid = (ai != GAP) & (bj != GAP)
            if int(valid.sum()) < 10:
                continue
            if float((ai[valid] == bj[valid]).mean()) >= id_thresh:
                used[j] = True
    # scale clusters back to full N if subsampled
    cluster_neff = float(n_clusters) * (n / n_work)

    # (2) mean pairwise identity → analytic Neff
    n_pairs = min(2000, n_work * (n_work - 1) // 2)
    ids: list[float] = []
    for _ in range(n_pairs):
        i = int(rng.integers(0, n_work))
        j = int(rng.integers(0, n_work))
        if i == j:
            continue
        ai, bj = work[i], work[j]
        valid = (ai != GAP) & (bj != GAP)
        if int(valid.sum()) < 10:
            continue
        ids.append(float((ai[valid] == bj[valid]).mean()))
    if ids:
        mean_id = float(np.mean(ids))
        mean_id = min(max(mean_id, 0.0), 1.0 - 1e-9)
        analytic = float(n) / (1.0 + (n - 1) * mean_id)
    else:
        analytic = float(n_clusters)

    return float(max(cluster_neff, analytic, 1.0))


def mutual_information_apc(
    mat: np.ndarray, covered: list[int], *, max_rows: int = 4000
) -> np.ndarray:
    """Pairwise MI with average-product correction; zeros on diagonal/uncovered."""
    n = mat.shape[1]
    mi = np.zeros((n, n), dtype=np.float64)
    if len(covered) < 2:
        return mi
    rng = np.random.default_rng(0)
    if mat.shape[0] > max_rows:
        mat = mat[rng.choice(mat.shape[0], max_rows, replace=False)]
    cols = {a: mat[:, a].astype(np.int64) for a in covered}
    gap = GAP
    for ai, a in enumerate(covered):
        ca = cols[a]
        va = ca != gap
        for b in covered[ai + 1 :]:
            cb = cols[b]
            valid = va & (cb != gap)
            if int(valid.sum()) < 20:
                continue
            xa, xb = ca[valid], cb[valid]
            # 20x20 joint (AA only)
            xa = np.clip(xa, 0, 19)
            xb = np.clip(xb, 0, 19)
            joint = np.bincount(xa * 20 + xb, minlength=400).reshape(20, 20).astype(float)
            tot = joint.sum()
            if tot < 20:
                continue
            joint /= tot
            pa = joint.sum(1)
            pb = joint.sum(0)
            nz = joint > 0
            outer = np.outer(pa, pb)
            m = joint[nz] * np.log((joint[nz]) / (outer[nz] + 1e-12))
            mi[a, b] = mi[b, a] = float(m.sum())
    # APC
    mi_mean = float(mi.mean())
    if mi_mean <= 0:
        return np.maximum(mi, 0.0)
    row_mean = mi.mean(1, keepdims=True)
    apc = (row_mean @ row_mean.T) / (mi_mean + 1e-12)
    corrected = mi - apc
    np.fill_diagonal(corrected, 0.0)
    return np.maximum(corrected, 0.0)


def extract_features(
    sequence: str,
    source: MsaSource,
    *,
    compute_coevolution: bool = True,
) -> MsaFeatures:
    seq = "".join(c for c in sequence.upper() if c in AA20)
    if not source.rows:
        return empty_features(seq, backend=source.backend, detail=source.detail)
    mapped = map_msa_to_query(seq, source.rows)
    if mapped is None:
        return empty_features(seq, backend=source.backend, detail=source.detail + "|map_fail")
    mat, covered, _ref = mapped
    conservation, gap_frac, entropy, aa_freq = _column_stats(mat, covered)
    neff = _neff(mat)
    if compute_coevolution and len(covered) >= int(math.ceil(PI + E)):
        coev = mutual_information_apc(mat, covered)
    else:
        coev = np.zeros((len(seq), len(seq)))
    return MsaFeatures(
        sequence=seq,
        n_seqs=int(mat.shape[0]),
        neff=float(neff),
        backend=source.backend,
        detail=source.detail,
        conservation=conservation,
        gap_frac=gap_frac,
        entropy=entropy,
        aa_freq=aa_freq,
        coevolution=coev,
        covered=covered,
    )


def build_msa_features(
    sequence: str,
    *,
    msa_path: Path | str | None = None,
    pfam: str | None = None,
    uniprot: str | None = None,
    jackhmmer_db: Path | str | None = None,
    hhblits_db: Path | str | None = None,
    compute_coevolution: bool = True,
) -> MsaFeatures:
    """One-call convenience: obtain MSA + extract features."""
    src = obtain_msa(
        sequence,
        msa_path=msa_path,
        pfam=pfam,
        uniprot=uniprot,
        jackhmmer_db=jackhmmer_db,
        hhblits_db=hhblits_db,
    )
    return extract_features(sequence, src, compute_coevolution=compute_coevolution)


# ── F-series injection (seed-closed) ──────────────────────────────────────


def evo_proximity_boost(features: MsaFeatures, gate: int = 7) -> np.ndarray:
    """Additive long-range proximity term from coevolution × conservation.

    M_evo[i,j] = evo_amp · Ĉ_ij · √(c_i c_j) · (1 - ḡ)
    where Ĉ is coevolution normalized by its top-L mean (data scale, not free param),
    c = conservation, ḡ = mean gap of the pair.
    Only |i-j| ≥ gate. Amplitude evo_amp is the F09-family scalar above.
    """
    n = len(features.sequence)
    M = np.zeros((n, n), dtype=np.float64)
    if not features.depth_ok or features.coevolution.max() <= 0:
        return M
    C = features.coevolution.copy()
    # Normalize by mean of top-L off-diagonal scores (L = n) — data-driven scale
    tri = C[np.triu_indices(n, k=gate)]
    if tri.size == 0 or float(tri.max()) <= 0:
        return M
    top_l = np.partition(tri, max(len(tri) - n, 0))[-min(n, len(tri)) :]
    scale = float(top_l.mean()) + 1e-12
    C = C / scale
    c = features.conservation
    g = features.gap_frac
    for i in range(n):
        for j in range(i + gate, n):
            w = math.sqrt(max(c[i], 0.0) * max(c[j], 0.0))
            gap_pen = 1.0 - 0.5 * (g[i] + g[j])
            if gap_pen <= 0:
                continue
            val = EVO_AMP * float(C[i, j]) * w * gap_pen
            if val > 0:
                M[i, j] = M[j, i] = val
    return M


def conservation_confidence(features: MsaFeatures) -> np.ndarray:
    """Per-residue evolutionary confidence in [0,1].

    conf = c · (1 - gap) · depth_factor
    depth_factor = 1 - e^{-Neff/φ}  (saturates with effective depth)
    """
    if features.n_seqs == 0:
        return np.zeros(len(features.sequence))
    depth = 1.0 - math.exp(-features.neff / PHI)
    return np.clip(
        features.conservation * (1.0 - features.gap_frac) * depth, 0.0, 1.0
    )


def tools_available() -> dict[str, Any]:
    return {
        "jackhmmer": _have("jackhmmer"),
        "hhblits": _have("hhblits"),
        "jackhmmer_db": os.environ.get("FSOT_JACKHMMER_DB"),
        "hhblits_db": os.environ.get("FSOT_HHBLITS_DB"),
        "pfam_network": True,
        "evo_amp": EVO_AMP,
    }


if __name__ == "__main__":
    # Self-check with a tiny synthetic MSA (no network)
    q = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    # Build a shallow fake MSA: query + a few mutated copies
    rows = [q]
    rng = np.random.default_rng(1)
    for _ in range(30):
        s = list(q)
        for i in rng.choice(len(s), size=5, replace=False):
            s[i] = AA20[int(rng.integers(0, 20))]
        rows.append("".join(s))
    # Align as ungapped equal-length (stockholm-like rows)
    feat = extract_features(q, MsaSource(rows=rows, backend="synthetic", detail="selfcheck"))
    print("selfcheck", feat.summary())
    print("conf mean", float(conservation_confidence(feat).mean()))
    print("evo boost nnz", int((evo_proximity_boost(feat) > 0).sum() // 2))
    print("tools", tools_available())
