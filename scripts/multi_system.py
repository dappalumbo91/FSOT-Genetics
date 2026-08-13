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
