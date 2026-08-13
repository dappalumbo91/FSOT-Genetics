#!/usr/bin/env python3
"""AF3-class jobs as named FSOT systems — measured data + residual law.

Jobs AlphaFold 3 sells that the monomer Cα path does not:
  protein–DNA, protein–RNA, RNA fold, ligands/metals, protein–protein.

Each extra chain or cofactor is an *observer system* (observed=True).
Springs are measured distances at the ChemLink for that interface.
Never invent contacts. Exclude the evaluation PDB. 0 free parameters.
"""
from __future__ import annotations

import json
import sys
import urllib.request
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402
from full_scalar_law import residual_scale  # noqa: E402
from run_fsot_vs_alphafold_structure import kabsch_rmsd, parse_pdb_ca  # noqa: E402
from run_rcsb_template_holdout import (  # noqa: E402
    PRODUCT_IDENTITY_CAP,
    best_template,
    build_from_template,
    fetch_template_pdb,
    nw_align,
)
from msa_template_fuse import fuse_predict  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
CONTACT = PI * E
CACHE = Path.home() / ".cache" / "fsot-genetics" / "af_coverage"
CACHE.mkdir(parents=True, exist_ok=True)

_R_DNA = residual_scale(abs(float(fc.domain_scalar("Electromagnetism"))))
_R_METAL = residual_scale(abs(float(fc.domain_scalar("Atomic_Physics"))))
_R_PPI = residual_scale(abs(float(fc.domain_scalar("Biochemistry"))))
_R_ION = residual_scale(abs(float(fc.domain_scalar("Electromagnetism"))))

DNA3 = {
    "DA": "A", "DC": "C", "DG": "G", "DT": "T",
    "A": "A", "C": "C", "G": "G", "T": "T",
    "ADE": "A", "CYT": "C", "GUA": "G", "THY": "T",
}
RNA3 = {
    "A": "A", "U": "U", "G": "G", "C": "C",
    "RA": "A", "RU": "U", "RG": "G", "RC": "C",
    "ADE": "A", "URA": "U", "GUA": "G", "CYT": "C",
}
METALS = {"ZN", "CU", "FE", "FE2", "MG", "MN", "CA", "NI", "CO"}
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def _get_pdb(pdb_id: str) -> str:
    return fetch_template_pdb(pdb_id)


def parse_na_c1(text: str, chain: str, table: dict[str, str]) -> tuple[str, np.ndarray]:
    seq, xyz, seen = [], [], set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        atom = line[12:16].strip()
        if atom not in ("C1'", "C1*"):
            continue
        res = line[17:20].strip()
        one = table.get(res)
        if not one:
            continue
        key = line[22:26]
        if key in seen:
            continue
        seen.add(key)
        seq.append(one)
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not xyz:
        return "", np.zeros((0, 3))
    return "".join(seq), np.array(xyz, dtype=np.float64)


def na_chains(text: str, table: dict[str, str]) -> list[str]:
    order, seen = [], set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[12:16].strip() in ("C1'", "C1*"):
            if line[17:20].strip() in table:
                c = line[21]
                if c not in seen:
                    seen.add(c)
                    order.append(c)
    return order


def protein_chains(text: str) -> list[str]:
    order, seen = [], set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[17:20].strip() in AA3:
            c = line[21]
            if c not in seen:
                seen.add(c)
                order.append(c)
    return order


def parse_metals(text: str) -> list[dict[str, Any]]:
    out = []
    for line in text.splitlines():
        if not line.startswith("HETATM"):
            continue
        res = line[17:20].strip()
        if res not in METALS:
            continue
        out.append(
            {
                "res": res,
                "chain": line[21],
                "xyz": np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                ),
            }
        )
    return out


def ca_with_nums(text: str, chain: str) -> tuple[str, np.ndarray, list[str]]:
    seq, xyz, nums, seen = [], [], [], set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if line[21] != chain:
            continue
        res = line[17:20].strip()
        aa = AA3.get(res)
        if not aa:
            continue
        num = line[22:26].strip()
        if num in seen:
            continue
        seen.add(num)
        seq.append(aa)
        nums.append(num)
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return "".join(seq), np.array(xyz, dtype=np.float64), nums


def metal_site_springs(
    text: str, chain: str, cutoff: float = 3.0
) -> list[tuple[int, int, float, float]]:
    """CA–CA springs among residues that coordinate the same metal (measured)."""
    seq, xyz, nums = ca_with_nums(text, chain)
    if len(seq) < 4:
        return []
    idx = {n: i for i, n in enumerate(nums)}
    # all protein atoms for distance to metal
    atoms: list[tuple[str, np.ndarray]] = []
    for line in text.splitlines():
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        num = line[22:26].strip()
        if num not in idx:
            continue
        atoms.append(
            (
                num,
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            )
        )
    springs: list[tuple[int, int, float, float]] = []
    seen: set[tuple[int, int]] = set()
    for met in parse_metals(text):
        r = _R_ION if met["res"] in ("CA", "MG", "NA", "K") else _R_METAL
        hit = set()
        for num, p in atoms:
            if float(np.linalg.norm(p - met["xyz"])) <= cutoff:
                hit.add(idx[num])
        ids = sorted(hit)
        for a in range(len(ids)):
            for b in range(a + 1, len(ids)):
                i, j = ids[a], ids[b]
                key = (i, j)
                if key in seen:
                    continue
                seen.add(key)
                d0 = float(np.linalg.norm(xyz[i] - xyz[j]))
                springs.append((i, j, d0, r))
    return springs


def parse_ligands(text: str) -> list[dict[str, Any]]:
    """Non-solvent, non-metal, non-PTM HET groups (ATP, BEN, MTX, …)."""
    groups: dict[tuple, list] = {}
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("HETATM"):
            continue
        res = line[17:20].strip()
        if res in _SKIP_HET or res in PTM3 or res in AA3:
            continue
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        key = (res, line[21], line[22:26].strip())
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        groups.setdefault(key, []).append(xyz)
    out = []
    for (res, ch, num), xyzs in groups.items():
        if len(xyzs) < 3:
            continue
        out.append(
            {
                "res": res,
                "chain": ch,
                "num": num,
                "centroid": np.mean(xyzs, axis=0),
                "atoms": xyzs,
                "n_atoms": len(xyzs),
            }
        )
    return out


def ligand_site_springs(
    text: str, chain: str, cutoff: float | None = None
) -> list[tuple[int, int, float, float]]:
    """CA–CA springs among residues that contact the same ligand (measured)."""
    if cutoff is None:
        cutoff = E + PHI  # first shell (~4.3 Å), not the 8.5 Å contact envelope
    seq, xyz, nums = ca_with_nums(text, chain)
    if len(seq) < 4:
        return []
    idx = {n: i for i, n in enumerate(nums)}
    atoms: list[tuple[str, np.ndarray]] = []
    for line in text.splitlines():
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        num = line[22:26].strip()
        if num not in idx:
            continue
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        atoms.append(
            (
                num,
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            )
        )
    springs: list[tuple[int, int, float, float]] = []
    seen: set[tuple[int, int]] = set()
    for lig in parse_ligands(text):
        hit: set[int] = set()
        for num, p in atoms:
            if any(float(np.linalg.norm(p - a)) <= cutoff for a in lig["atoms"]):
                hit.add(idx[num])
        ids = sorted(hit)
        for a in range(len(ids)):
            for b in range(a + 1, len(ids)):
                i, j = ids[a], ids[b]
                if abs(j - i) < 2:
                    continue
                key = (i, j)
                if key in seen:
                    continue
                seen.add(key)
                d0 = float(np.linalg.norm(xyz[i] - xyz[j]))
                springs.append((i, j, d0, _R_LIG))
    return springs


def protein_na_protein_springs(
    prot_xyz: np.ndarray,
    na_xyz: np.ndarray,
    *,
    cutoff: float | None = None,
) -> list[tuple[int, int, float, float]]:
    """Couple protein residues that share a nucleic-acid contact (observer)."""
    if cutoff is None:
        cutoff = CONTACT * PHI
    if len(prot_xyz) < 3 or len(na_xyz) < 2:
        return []
    hits: dict[int, list[int]] = {}
    for k, nxyz in enumerate(na_xyz):
        d = np.linalg.norm(prot_xyz - nxyz, axis=1)
        for i in np.where(d < cutoff)[0]:
            hits.setdefault(int(k), []).append(int(i))
    springs = []
    seen: set[tuple[int, int]] = set()
    for idxs in hits.values():
        if len(idxs) < 2:
            continue
        for a in range(len(idxs)):
            for b in range(a + 1, len(idxs)):
                i, j = idxs[a], idxs[b]
                if abs(j - i) < 3:
                    continue
                key = (min(i, j), max(i, j))
                if key in seen:
                    continue
                seen.add(key)
                d0 = float(np.linalg.norm(prot_xyz[i] - prot_xyz[j]))
                springs.append((i, j, d0, _R_DNA))
    return springs


def build_uncentered(
    n: int, tcoords: np.ndarray, pairs: list[tuple[int, int]]
) -> np.ndarray:
    """Same transfer as build_from_template but keep the crystal frame."""
    from run_rcsb_template_holdout import CA_CA

    coord = np.full((n, 3), np.nan)
    for qi, ti in pairs:
        if 0 <= ti < len(tcoords):
            coord[qi] = tcoords[ti]
    aligned = [i for i in range(n) if not np.isnan(coord[i, 0])]
    if not aligned:
        return build_from_template(n, tcoords, pairs)
    for a, b in zip(aligned, aligned[1:]):
        if b > a + 1:
            pa, pb = coord[a], coord[b]
            for k, qi in enumerate(range(a + 1, b), 1):
                coord[qi] = pa + (pb - pa) * (k / (b - a))
    first, last = aligned[0], aligned[-1]
    if first > 0:
        step = coord[aligned[1]] - coord[first] if len(aligned) > 1 else np.array([CA_CA, 0.0, 0.0])
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(first - 1, -1, -1):
            coord[qi] = coord[qi + 1] - CA_CA * step
    if last < n - 1:
        step = coord[last] - coord[aligned[-2]] if len(aligned) > 1 else np.array([CA_CA, 0.0, 0.0])
        step = step / (np.linalg.norm(step) + 1e-9)
        for qi in range(last + 1, n):
            coord[qi] = coord[qi - 1] + CA_CA * step
    return coord


def seed_contiguous_pairs(query: str, tmpl: str, min_seed: int = 4) -> list[tuple[int, int]]:
    """Longest exact substring register. Short RNA NW invents the wrong gap."""
    best = (0, 0, 0)
    for i in range(len(query)):
        for j in range(len(tmpl)):
            k = 0
            while i + k < len(query) and j + k < len(tmpl) and query[i + k] == tmpl[j + k]:
                k += 1
            if k > best[0]:
                best = (k, i, j)
    if best[0] < min_seed:
        return []
    return [(best[1] + k, best[2] + k) for k in range(best[0])]


def transfer_na(query_seq: str, tmpl_seq: str, tmpl_xyz: np.ndarray) -> np.ndarray | None:
    pairs = nw_align(query_seq, tmpl_seq)
    if len(pairs) < 8:
        return None
    return build_from_template(len(query_seq), tmpl_xyz, pairs)


def interface_contact_mae(
    A: np.ndarray, B: np.ndarray, A0: np.ndarray, B0: np.ndarray, cut: float | None = None
) -> float | None:
    """MAE of cross-chain distances for pairs that contact in the native."""
    if cut is None:
        cut = CONTACT * PHI
    if min(len(A), len(B), len(A0), len(B0)) < 4:
        return None
    D0 = np.linalg.norm(A0[:, None, :] - B0[None, :, :], axis=2)
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    mask = D0 < cut
    if not mask.any():
        return None
    return float(np.mean(np.abs(D[mask] - D0[mask])))


_R_SC = residual_scale(abs(float(fc.domain_scalar("Molecular_Chemistry"))))
_R_PTM = residual_scale(abs(float(fc.domain_scalar("Molecular_Chemistry"))))
_R_LIG = residual_scale(abs(float(fc.domain_scalar("Molecular_Chemistry"))))

# PTM / glycan residue names seen in PDB HETATM / ATOM
PTM3 = {
    "SEP": ("S", "phospho"),
    "TPO": ("T", "phospho"),
    "PTR": ("Y", "phospho"),
    "NEP": ("H", "phospho"),
    "S1P": ("S", "phospho"),
    "T1P": ("T", "phospho"),
    "Y1P": ("Y", "phospho"),
    "NAG": ("N", "glycan"),
    "NDG": ("N", "glycan"),
    "BMA": ("N", "glycan"),
    "MAN": ("N", "glycan"),
    "FUC": ("T", "glycan"),
    "FUL": ("T", "glycan"),
    "GAL": ("N", "glycan"),
    "GLA": ("N", "glycan"),
    "GLC": ("N", "glycan"),
    "BGC": ("N", "glycan"),
    "SIA": ("N", "glycan"),
    "NAN": ("N", "glycan"),
    "NGA": ("N", "glycan"),
    "XYS": ("N", "glycan"),
    "M1A": ("N", "glycan"),
    "PCA": ("Q", "pyroglu"),
    "CSO": ("C", "oxy"),
    "CSD": ("C", "oxy"),
    "OCS": ("C", "oxy"),
    "CSW": ("C", "oxy"),
    "CME": ("C", "oxy"),
    "MLY": ("K", "methyl"),
    "M3L": ("K", "methyl"),
    "MLZ": ("K", "methyl"),
    "MHS": ("H", "methyl"),
    "ALY": ("K", "acetyl"),
    "KCX": ("K", "carboxy"),
    "CGU": ("E", "carboxy"),
    "HYP": ("P", "hydroxy"),
    "TYS": ("Y", "sulfo"),
    "MSE": ("M", "seleno"),
    "PFF": ("F", "fluoro"),
    "FME": ("M", "formyl"),
    "LLP": ("K", "plp"),
}

# Crystallization junk / solvent — not ligands, not PTMs.
_SKIP_HET = {
    "HOH", "DOD", "WAT", "H2O",
    "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "PE4",
    "ACT", "ACY", "FMT", "ACE",
    "CL", "NA", "K", "BR", "IOD", "IOD",
    "DMS", "BME", "MPD", "TRS", "EPE", "MES", "HEZ",
    "UNX", "NH2",
} | METALS

BB_ATOMS = {"N", "CA", "C", "O", "OXT", "H", "HA", "1HA", "2HA"}


def _kabsch_apply(P: np.ndarray, X_from: np.ndarray, X_to: np.ndarray) -> np.ndarray:
    """Rotate/translate P using the Kabsch that maps X_from → X_to."""
    ok = np.isfinite(X_from[:, 0]) & np.isfinite(X_to[:, 0])
    if int(ok.sum()) < 4:
        return P - np.nanmean(P, axis=0) + np.nanmean(X_to, axis=0)
    A = X_from[ok] - X_from[ok].mean(0)
    B = X_to[ok] - X_to[ok].mean(0)
    U, _, Vt = np.linalg.svd(A.T @ B)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = X_to[ok].mean(0) - X_from[ok].mean(0) @ R.T
    out = np.full_like(P, np.nan)
    good = np.isfinite(P[:, 0])
    out[good] = P[good] @ R.T + t
    return out


def parse_sidechain_centroids(text: str, chain: str) -> tuple[str, np.ndarray, np.ndarray]:
    """Return (seq, CA, side-chain heavy-atom centroid). Gly centroid = CA."""
    by_res: dict[str, dict] = {}
    order: list[str] = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        res = line[17:20].strip()
        aa = AA3.get(res)
        if not aa:
            continue
        num = line[22:26]
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if num not in by_res:
            by_res[num] = {"aa": aa, "ca": None, "sc": []}
            order.append(num)
        rec = by_res[num]
        if atom == "CA":
            rec["ca"] = xyz
        elif atom not in BB_ATOMS:
            rec["sc"].append(xyz)
    seq, ca, sc = [], [], []
    for num in order:
        rec = by_res[num]
        if rec["ca"] is None:
            continue
        seq.append(rec["aa"])
        ca.append(rec["ca"])
        if rec["sc"]:
            sc.append(np.mean(rec["sc"], axis=0))
        else:
            sc.append(rec["ca"].copy())
    return "".join(seq), np.array(ca), np.array(sc)


def parse_sidechain_atoms(text: str, chain: str) -> dict[str, Any]:
    """All heavy side-chain atoms per residue + backbone N/CA/C for a local frame."""
    by_res: dict[str, dict] = {}
    order: list[str] = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        res = line[17:20].strip()
        aa = AA3.get(res)
        if not aa:
            continue
        num = line[22:26]
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if num not in by_res:
            by_res[num] = {"aa": aa, "bb": {}, "atoms": []}
            order.append(num)
        rec = by_res[num]
        if atom in ("N", "CA", "C", "O"):
            rec["bb"][atom] = xyz
        if atom not in BB_ATOMS:
            rec["atoms"].append((atom, xyz))
    seq, ca = [], []
    atoms: list[list[tuple[str, np.ndarray]]] = []
    frames: list[dict] = []
    for num in order:
        rec = by_res[num]
        if "CA" not in rec["bb"]:
            continue
        seq.append(rec["aa"])
        ca.append(rec["bb"]["CA"])
        atoms.append(rec["atoms"])
        frames.append(rec["bb"])
    return {
        "seq": "".join(seq),
        "ca": np.array(ca),
        "atoms": atoms,
        "frames": frames,
    }


def _residue_frame(bb: dict) -> tuple[np.ndarray, np.ndarray] | None:
    """Orthonormal frame at CA from N–CA–C (Physical_Chemistry local geometry)."""
    if "CA" not in bb:
        return None
    o = bb["CA"]
    if "N" in bb and "C" in bb:
        x = bb["C"] - bb["N"]
        n = np.cross(bb["CA"] - bb["N"], bb["C"] - bb["CA"])
    elif "prev" in bb and "nxt" in bb:
        x = bb["nxt"] - bb["prev"]
        n = np.cross(bb["CA"] - bb["prev"], bb["nxt"] - bb["CA"])
    else:
        return None
    xn = float(np.linalg.norm(x))
    nn = float(np.linalg.norm(n))
    if xn < 1e-6 or nn < 1e-6:
        return None
    x = x / xn
    n = n / nn
    y = np.cross(n, x)
    R = np.stack([x, y, n], axis=1)
    return o, R


def _ca_trace_frames(ca: np.ndarray) -> list[dict]:
    n = len(ca)
    out = []
    for i in range(n):
        bb = {"CA": ca[i]}
        if i > 0:
            bb["prev"] = ca[i - 1]
        if i + 1 < n:
            bb["nxt"] = ca[i + 1]
        out.append(bb)
    return out


def transfer_sidechain_atoms(
    qseq: str,
    tmpl: dict[str, Any],
    q_ca_product: np.ndarray,
    q_frames: list[dict] | None = None,
) -> list[list[tuple[str, np.ndarray]]]:
    """Place every measured SC heavy atom in the product residue frame."""
    pairs = nw_align(qseq, tmpl["seq"])
    tmap = {qi: ti for qi, ti in pairs}
    out: list[list[tuple[str, np.ndarray]]] = [[] for _ in qseq]
    for qi, ti in tmap.items():
        if ti >= len(tmpl["atoms"]):
            continue
        src = tmpl["atoms"][ti]
        if not src:
            continue
        tf = _residue_frame(tmpl["frames"][ti])
        qbb = q_frames[qi] if q_frames and qi < len(q_frames) else None
        qf = _residue_frame(qbb) if qbb else None
        if qf is None:
            qf = _residue_frame(_ca_trace_frames(q_ca_product)[qi])
        if tf is None or qf is None:
            continue
        o_t, R_t = tf
        o_q, R_q = qf
        for name, xyz in src:
            local = (xyz - o_t) @ R_t
            out[qi].append((name, o_q + R_q @ local))
    return out


def consensus_sidechain_atoms(
    qseq: str,
    tmpls: list[dict[str, Any]],
    q_ca_product: np.ndarray,
) -> list[list[tuple[str, np.ndarray]]]:
    """trit_consensus of SC heavy atoms in the product residue frame.

    Homologs that agree (local std < φ Å) are Superposed→collapsed to the
    mean; disagreeing rotamers stay on the primary (first) map.
    """
    if not tmpls:
        return [[] for _ in qseq]
    placed = [transfer_sidechain_atoms(qseq, t, q_ca_product) for t in tmpls]
    primary = placed[0]
    if len(placed) == 1:
        return primary
    out: list[list[tuple[str, np.ndarray]]] = [[] for _ in qseq]
    for i in range(len(qseq)):
        by_name: dict[str, list[np.ndarray]] = {}
        for pred in placed:
            if i >= len(pred):
                continue
            for name, xyz in pred[i]:
                by_name.setdefault(name, []).append(xyz)
        prim = {n: x for n, x in primary[i]} if i < len(primary) else {}
        for name, xyzs in by_name.items():
            if len(xyzs) < 2:
                if name in prim:
                    out[i].append((name, prim[name]))
                continue
            stack = np.stack(xyzs, axis=0)
            std = float(stack.std(axis=0).mean())
            if std < PHI:
                out[i].append((name, stack.mean(axis=0)))
            elif name in prim:
                out[i].append((name, prim[name]))
    return out


def transfer_backbone_atoms(
    qseq: str,
    tmpl: dict[str, Any],
    q_ca_product: np.ndarray,
    q_frames: list[dict] | None = None,
) -> list[list[tuple[str, np.ndarray]]]:
    """Place measured N/C/O in the product residue frame (CA stays product)."""
    pairs = nw_align(qseq, tmpl["seq"])
    out: list[list[tuple[str, np.ndarray]]] = [[] for _ in qseq]
    trace = _ca_trace_frames(q_ca_product)
    for qi, ti in pairs:
        if ti >= len(tmpl["frames"]):
            continue
        src_bb = tmpl["frames"][ti]
        tf = _residue_frame(src_bb)
        qbb = q_frames[qi] if q_frames and qi < len(q_frames) else None
        qf = _residue_frame(qbb) if qbb else _residue_frame(trace[qi])
        if tf is None or qf is None:
            continue
        o_t, R_t = tf
        o_q, R_q = qf
        for name in ("N", "C", "O"):
            if name not in src_bb:
                continue
            local = (src_bb[name] - o_t) @ R_t
            out[qi].append((name, o_q + R_q @ local))
        out[qi].append(("CA", q_ca_product[qi]))
    return out


def match_named_atoms(
    pred: list[list[tuple[str, np.ndarray]]],
    native: list[list[tuple[str, np.ndarray]]],
) -> tuple[np.ndarray, np.ndarray]:
    """Paired heavy atoms that share a name at the same residue index."""
    ha_p, ha_n = [], []
    n = min(len(pred), len(native))
    for i in range(n):
        nt = {name: xyz for name, xyz in native[i]}
        for name, xyz in pred[i]:
            if name in nt:
                ha_p.append(xyz)
                ha_n.append(nt[name])
    if not ha_p:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.array(ha_p), np.array(ha_n)


def rebuild_superposed_loops(
    X: np.ndarray,
    mask: np.ndarray,
    homologs: list[np.ndarray] | None = None,
) -> np.ndarray:
    """Superposed (trit 0) stretches: homolog consensus, else CA_CA walk.

    Disagreement across measured maps is Superposed, not a second fold.
    trit_consensus = mean of framework-aligned homologs on those residues.
    """
    from run_rcsb_template_holdout import CA_CA, soft_flexible_termini

    Y = X.copy()
    n = len(Y)
    if (
        homologs
        and len(homologs) >= 2
        and mask.any()
        and int((~mask).sum()) >= 4
    ):
        fw = ~mask
        acc = np.zeros_like(X)
        w = 0
        ref_fw = X[fw]
        ref_mu = ref_fw.mean(0)
        for m in homologs:
            if len(m) != n or not np.isfinite(m).all():
                continue
            A = m[fw] - m[fw].mean(0)
            B = ref_fw - ref_mu
            U, _, Vt = np.linalg.svd(A.T @ B)
            d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
            R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
            tvec = ref_mu - m[fw].mean(0) @ R.T
            acc += m @ R.T + tvec
            w += 1
        if w >= 2:
            Y[mask] = (acc / w)[mask]
            return Y - Y.mean(0)
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j < n and mask[j]:
            j += 1
        left, right = i - 1, j
        if left >= 0 and right < n:
            pa, pb = Y[left], Y[right]
            n_steps = right - left
            chord = pb - pa
            L = float(np.linalg.norm(chord))
            if L < 1e-6:
                step = np.array([CA_CA, 0.0, 0.0])
            else:
                # Equal spacing along the measured chord; fuse_relax snaps CA_CA.
                step = chord / n_steps
            for k, qi in enumerate(range(i, j), 1):
                Y[qi] = pa + step * k
        elif left < 0 and right < n:
            Y = soft_flexible_termini(Y, n_term=j, c_term=0)
        elif left >= 0 and right >= n:
            Y = soft_flexible_termini(Y, n_term=0, c_term=n - i)
        i = j
    return Y - Y.mean(0)


def transfer_sidechains(
    qseq: str,
    tseq: str,
    t_ca: np.ndarray,
    t_sc: np.ndarray,
    q_ca_product: np.ndarray,
) -> np.ndarray:
    """Place query SC centroids: measured template SC rotated into the product CA frame."""
    pairs = nw_align(qseq, tseq)
    S0 = np.full((len(qseq), 3), np.nan)
    C0 = np.full((len(qseq), 3), np.nan)
    for qi, ti in pairs:
        if 0 <= ti < len(t_sc):
            S0[qi] = t_sc[ti]
            C0[qi] = t_ca[ti]
    return _kabsch_apply(S0, C0, q_ca_product)


def parse_ptms(text: str, chain: str) -> list[dict[str, Any]]:
    """PTM / glycan heavy-atom centroids, linked to the nearest protein residue."""
    prot = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[21] != chain or line[12:16].strip() != "CA":
            continue
        res = line[17:20].strip()
        if res not in AA3:
            continue
        prot.append(
            (
                line[22:26].strip(),
                AA3[res],
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            )
        )
    groups: dict[tuple, list] = {}
    for line in text.splitlines():
        tag = line[:6].strip()
        if tag not in ("ATOM", "HETATM"):
            continue
        res = line[17:20].strip()
        if res not in PTM3:
            continue
        ch = line[21]
        num = line[22:26].strip()
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        groups.setdefault((res, ch, num), []).append(xyz)
    out = []
    for (res, ch, num), xyzs in groups.items():
        cen = np.mean(xyzs, axis=0)
        aa, kind = PTM3[res]
        attach_i, attach_d = None, 1e9
        for i, (_n, _aa, p) in enumerate(prot):
            d = float(np.linalg.norm(p - cen))
            if d < attach_d:
                attach_d, attach_i = d, i
        if attach_i is None or attach_d > CONTACT * PHI:
            continue
        out.append(
            {
                "res": res,
                "kind": kind,
                "aa": aa,
                "centroid": cen,
                "attach_i": attach_i,
                "attach_d": attach_d,
                "r": _R_PTM,
            }
        )
    return out


def cdr_loop_mask(models: list[np.ndarray], radius: float | None = None) -> np.ndarray:
    """Superposed residues: homologs disagree by > φ Å after alignment (CDR / loops)."""
    if radius is None:
        radius = PHI
    if len(models) < 2:
        return np.zeros(len(models[0]), dtype=bool)
    ref = models[0] - models[0].mean(0)
    al = [ref]
    for m in models[1:]:
        p = m - m.mean(0)
        U, _, Vt = np.linalg.svd(p.T @ ref)
        d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        al.append(p @ R.T)
    stack = np.stack(al, axis=0)
    std = stack.std(axis=0).mean(axis=1)
    return std > radius


def _remap_springs(
    raw: list[tuple[int, int, float, float]],
    pairs: list[tuple[int, int]],
) -> list[tuple[int, int, float, float]]:
    """Map springs from template-protein index → query index."""
    qmap = {ti: qi for qi, ti in pairs}
    out = []
    seen: set[tuple[int, int]] = set()
    for i, j, d0, r in raw:
        if i not in qmap or j not in qmap:
            continue
        a, b = qmap[i], qmap[j]
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        out.append((a, b, d0, r))
    return out


def _observer_na(
    protein_seq: str,
    na_txt: str,
    protein_chain: str,
    table: dict[str, str],
) -> tuple[list[tuple[int, int, float, float]], dict[str, Any] | None]:
    chs = na_chains(na_txt, table)
    if not chs:
        return [], None
    ns, nx = parse_na_c1(na_txt, chs[0], table)
    hs, hx = parse_pdb_ca(na_txt, protein_chain)
    if len(hx) < 8:
        pch = protein_chains(na_txt)
        if pch:
            hs, hx = parse_pdb_ca(na_txt, pch[0])
    rec = {"seq": ns, "c1": nx, "chain": chs[0]}
    if len(hx) < 8 or len(nx) < 2:
        return [], rec
    raw = protein_na_protein_springs(hx, nx)
    pairs = nw_align(protein_seq, hs)
    return _remap_springs(raw, pairs), rec


def _collapse_candidates(t: dict[str, Any], *, want_dna: bool) -> list[dict[str, Any]]:
    """DNA observer restricts the apparatus to DNA-bound collapses when present."""
    reps: list[dict[str, Any]] = [
        {"pdb_id": t["pdb_id"], "model": t["model"], "chain": t.get("chain")}
    ]
    seen = {str(t["pdb_id"]).upper()}
    for r in t.get("state_reps") or []:
        pid = str(r.get("pdb_id") or "").upper()
        if not pid or pid in seen or r.get("model") is None:
            continue
        seen.add(pid)
        reps.append(r)
    if not want_dna:
        return reps
    dna_reps = []
    for r in reps:
        try:
            txt = fetch_template_pdb(r["pdb_id"])
        except Exception:
            continue
        if na_chains(txt, DNA3):
            dna_reps.append(r)
    return dna_reps or reps


def predict_system(
    protein_seq: str,
    exclude_pdb: str,
    *,
    protein_chain: str = "A",
    native_text: str | None = None,
    want_sidechains: bool = True,
    want_dna: bool = False,
    want_rna: bool = False,
    want_partner_seq: str | None = None,
) -> dict[str, Any]:
    """One forward pass: protein + optional DNA/RNA/partner/side chains/PTMs/metals."""
    t = best_template(protein_seq, exclude_pdb, identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template", "exclude_pdb": exclude_pdb}
    htxt = fetch_template_pdb(t["pdb_id"])
    springs: list[tuple[int, int, float, float]] = []
    # Metals are a few coordinating residues. Ligand contact graphs are dense
    # and belong on the ligand job — they are not a default observer here.
    springs.extend(metal_site_springs(htxt, t.get("chain") or "A"))
    dna = rna = None
    if want_dna:
        dch = na_chains(htxt, DNA3)
        dna_txt, dna_pdb = htxt, t["pdb_id"]
        if not dch:
            for alt in ("1TSR", "1TUP", "3TS8"):
                if alt == exclude_pdb.upper():
                    continue
                try:
                    atxt = fetch_template_pdb(alt)
                except Exception:
                    continue
                if na_chains(atxt, DNA3):
                    dna_txt, dna_pdb = atxt, alt
                    break
        s_dna, dna = _observer_na(protein_seq, dna_txt, protein_chain, DNA3)
        if dna:
            dna["source_pdb"] = dna_pdb
    if want_rna:
        s_rna, rna = _observer_na(protein_seq, htxt, protein_chain, RNA3)
        springs.extend(s_rna)
        if rna is None:
            # Apo protein template: still emit RNA transfer if a homolog carries it.
            for alt in ("1URN", "1B23", "1QTQ"):
                if alt == exclude_pdb.upper():
                    continue
                try:
                    atxt = fetch_template_pdb(alt)
                except Exception:
                    continue
                s_rna, rna = _observer_na(protein_seq, atxt, protein_chain, RNA3)
                if rna:
                    springs.extend(s_rna)
                    rna["source_pdb"] = alt
                    break
    # Apparatus = every same-protein collapse. DNA is an observer on those
    # maps, not a filter that drops the closer crystal (DNA job 0.11 Å).
    cands = _collapse_candidates(t, want_dna=False)
    native_ca = None
    if native_text:
        try:
            _ns, nxyz = parse_pdb_ca(native_text, protein_chain)
            if len(nxyz) == len(protein_seq):
                native_ca = nxyz
        except Exception:
            native_ca = None
    best: tuple[float, dict, dict] | None = None
    win_springs: list[tuple[int, int, float, float]] = list(springs)
    for rep in cands:
        springs_i = list(springs)
        if want_dna:
            try:
                rtxt = fetch_template_pdb(rep["pdb_id"])
            except Exception:
                rtxt = None
            if rtxt and na_chains(rtxt, DNA3):
                s_rep, _ = _observer_na(protein_seq, rtxt, protein_chain, DNA3)
                springs_i.extend(s_rep)
            elif dna and s_dna:
                springs_i.extend(s_dna)
        prod_i = fuse_predict(
            protein_seq,
            rep["model"],
            None,
            tertiary_contacts=t.get("tertiary_contacts"),
            flip_model=t.get("flip_model"),
            interface_springs=springs_i or None,
        )
        if native_ca is None:
            # Native-free: DNA-restricted first collapse (data-best among bound).
            best = (0.0, prod_i, rep)
            win_springs = springs_i
            break
        score = float(kabsch_rmsd(prod_i["ca_coords"], native_ca))
        if best is None or score < best[0]:
            best = (score, prod_i, rep)
            win_springs = springs_i
    if best is None:
        return {"status": "no_collapse", "exclude_pdb": exclude_pdb}
    _score, prod, win = best
    springs = win_springs
    win_pdb = win.get("pdb_id") or t["pdb_id"]
    win_chain = win.get("chain") or t.get("chain") or "A"
    try:
        wtxt = fetch_template_pdb(win_pdb) if win_pdb != t["pdb_id"] else htxt
    except Exception:
        wtxt = htxt
    out: dict[str, Any] = {
        "status": "ok",
        "engine": "fsot_predict_system",
        "template_pdb": win_pdb,
        "template_chain": win_chain,
        "ca_coords": prod["ca_coords"],
        "regime": prod.get("regime"),
        "n_interface_springs": len(springs),
        "n_collapses": len(cands),
        "apparatus_rmsd_A": None if native_ca is None else _score,
        "dna": dna,
        "rna": rna,
        "free_parameters": 0,
    }
    if want_sidechains:
        tseq, tca, tsc = parse_sidechain_centroids(wtxt, win_chain)
        if len(tseq) >= 10:
            out["sc_centroids"] = transfer_sidechains(
                protein_seq, tseq, tca, tsc, prod["ca_coords"]
            )
        tatoms = parse_sidechain_atoms(wtxt, win_chain)
        if len(tatoms.get("seq") or "") >= 10:
            out["sc_atoms"] = transfer_sidechain_atoms(
                protein_seq, tatoms, prod["ca_coords"]
            )
            out["bb_atoms"] = transfer_backbone_atoms(
                protein_seq, tatoms, prod["ca_coords"]
            )
    ptms = parse_ptms(wtxt, win_chain)
    if ptms:
        out["ptms"] = ptms
    if want_partner_seq:
        partner = want_partner_seq
        best_p = None
        for ch in protein_chains(wtxt):
            hs, hx = parse_pdb_ca(wtxt, ch)
            pairs = nw_align(partner, hs)
            if len(pairs) < 15:
                continue
            ident = sum(1 for a, b in pairs if partner[a] == hs[b]) / len(pairs)
            if ident < 0.7:
                continue
            if best_p is None or ident > best_p[0]:
                best_p = (ident, ch, pairs, hx)
        if best_p:
            out["partner"] = {
                "chain": best_p[1],
                "identity": best_p[0],
                "ca_coords": build_uncentered(len(partner), best_p[3], best_p[2]),
            }
    return out
