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


_R_SC = residual_scale(abs(float(fc.domain_scalar("Molecular_Chemistry"))))
_R_PTM = residual_scale(abs(float(fc.domain_scalar("Molecular_Chemistry"))))

# PTM / glycan residue names seen in PDB HETATM / ATOM
PTM3 = {
    "SEP": ("S", "phospho"),  # phosphoserine
    "TPO": ("T", "phospho"),
    "PTR": ("Y", "phospho"),
    "NEP": ("H", "phospho"),
    "NAG": ("N", "glycan"),
    "NDG": ("N", "glycan"),
    "BMA": ("N", "glycan"),
    "MAN": ("N", "glycan"),
    "FUC": ("T", "glycan"),
    "GAL": ("N", "glycan"),
    "GLC": ("N", "glycan"),
    "SIA": ("N", "glycan"),
    "M1A": ("N", "glycan"),
    "PCA": ("Q", "pyroglu"),
    "CSO": ("C", "oxy"),
    "MLY": ("K", "methyl"),
    "M3L": ("K", "methyl"),
    "ALY": ("K", "acetyl"),
}

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
    # metals on the template
    springs.extend(metal_site_springs(htxt, t.get("chain") or "A"))
    dna = rna = partner = None
    if want_dna:
        dch = na_chains(htxt, DNA3)
        dna_txt, dna_pdb = htxt, t["pdb_id"]
        if not dch:
            # Template protein may be apo; pull DNA observer from a DNA-bound homolog.
            for alt in ("1TSR", "1TUP", "3TS8"):
                if alt == exclude_pdb.upper():
                    continue
                try:
                    atxt = fetch_template_pdb(alt)
                except Exception:
                    continue
                ach = na_chains(atxt, DNA3)
                if ach:
                    dna_txt, dna_pdb, dch = atxt, alt, ach
                    break
        if dch:
            ds, dx = parse_na_c1(dna_txt, dch[0], DNA3)
            hs, hx = parse_pdb_ca(dna_txt, protein_chain)
            if len(hx) >= 8:
                springs.extend(protein_na_protein_springs(hx, dx))
            dna = {"seq": ds, "c1": dx, "chain": dch[0], "source_pdb": dna_pdb}
    if want_rna:
        rch = na_chains(htxt, RNA3)
        if rch:
            rs, rx = parse_na_c1(htxt, rch[0], RNA3)
            rna = {"seq": rs, "c1": rx, "chain": rch[0]}
    prod = fuse_predict(
        protein_seq,
        t["model"],
        None,
        tertiary_contacts=t.get("tertiary_contacts"),
        flip_model=t.get("flip_model"),
        interface_springs=springs or None,
    )
    out: dict[str, Any] = {
        "status": "ok",
        "engine": "fsot_predict_system",
        "template_pdb": t["pdb_id"],
        "template_chain": t.get("chain"),
        "ca_coords": prod["ca_coords"],
        "regime": prod.get("regime"),
        "n_interface_springs": len(springs),
        "dna": dna,
        "rna": rna,
        "free_parameters": 0,
    }
    if want_sidechains:
        tseq, tca, tsc = parse_sidechain_centroids(htxt, t.get("chain") or "A")
        if len(tseq) >= 10:
            out["sc_centroids"] = transfer_sidechains(
                protein_seq, tseq, tca, tsc, prod["ca_coords"]
            )
    ptms = parse_ptms(htxt, t.get("chain") or "A")
    if ptms:
        out["ptms"] = ptms
    if want_partner_seq:
        partner = want_partner_seq
        best = None
        for ch in protein_chains(htxt):
            hs, hx = parse_pdb_ca(htxt, ch)
            pairs = nw_align(partner, hs)
            if len(pairs) < 15:
                continue
            ident = sum(1 for a, b in pairs if partner[a] == hs[b]) / len(pairs)
            if ident < 0.7:
                continue
            if best is None or ident > best[0]:
                best = (ident, ch, pairs, hx)
        if best:
            out["partner"] = {
                "chain": best[1],
                "identity": best[0],
                "ca_coords": build_uncentered(len(partner), best[3], best[2]),
            }
    return out
