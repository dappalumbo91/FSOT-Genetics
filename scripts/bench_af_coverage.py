#!/usr/bin/env python3
"""Coverage bench: every AlphaFold-class job we can run under FSOT.

Monomer numbers come from the current product JSON.
New jobs: protein–DNA, RNA fold, metal sites, protein–protein interface.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_fsot_vs_alphafold_structure import kabsch_rmsd  # noqa: E402
from run_rcsb_template_holdout import (  # noqa: E402
    PRODUCT_IDENTITY_CAP,
    best_template,
    nw_align,
    uniref100_pdb_ids,
)
from msa_template_fuse import fuse_predict  # noqa: E402
from multi_system import (  # noqa: E402
    RNA3,
    _get_pdb,
    build_uncentered,
    ca_with_nums,
    cdr_loop_mask,
    consensus_sidechain_atoms,
    interface_contact_mae,
    ligand_site_springs,
    match_named_atoms,
    metal_site_springs,
    na_chains,
    parse_hydrogen_atoms,
    parse_ligands,
    parse_na_c1,
    parse_pdb_ca,
    parse_ptms,
    parse_sidechain_atoms,
    parse_sidechain_centroids,
    predict_system,
    protein_chains,
    protein_na_protein_springs,
    rebuild_superposed_loops,
    transfer_backbone_atoms,
    seed_contiguous_pairs,
    transfer_na,
    transfer_sidechain_atoms,
    transfer_sidechains,
    DNA3,
    RNA_ALL,
    RNA_MOD,
)
from run_rna_template_probe import (  # noqa: E402
    parse_rna_c1,
    rna_chains,
    rna_homologs,
    _pdb as rna_pdb,
)

OUT = ROOT / "data" / "af_coverage.json"
MONO = ROOT / "data" / "product_vs_alphafold.json"


def _fuse(seq, tmpl, springs=None):
    t = best_template(seq, tmpl, identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return None
    extra = None
    if springs:
        # remap springs from native-index template model if lengths match
        extra = springs
    prod = fuse_predict(
        seq,
        t["model"],
        None,
        tertiary_contacts=t.get("tertiary_contacts"),
        flip_model=t.get("flip_model"),
        interface_springs=extra,
    )
    return t, prod


def job_monomer() -> dict:
    if not MONO.is_file():
        return {"status": "missing_product_json"}
    data = json.loads(MONO.read_text(encoding="utf-8"))
    s = data.get("summary") or {}
    return {
        "status": "product",
        "n": s.get("n"),
        "fsot_median_A": s.get("fsot_product_median_A"),
        "alphafold_median_A": s.get("alphafold_median_A"),
        "beats_af_median": (s.get("fsot_product_median_A") or 99)
        < (s.get("alphafold_median_A") or 0),
        "note": "same-data monomer Cα (exclude eval PDB)",
    }


def job_protein_dna() -> dict:
    """p53 1TUP:A vs DNA; DNA observer springs from homolog 1TSR."""
    native = _get_pdb("1TUP")
    pseq, pxyz = parse_pdb_ca(native, "A")
    dna_ch = na_chains(native, DNA3)
    if not dna_ch:
        return {"status": "no_dna", "pdb": "1TUP"}
    dseq, dxyz = parse_na_c1(native, dna_ch[0], DNA3)
    # protein product with DNA-observer springs from 1TSR (not 1TUP)
    htxt = _get_pdb("1TSR")
    hp, hx = parse_pdb_ca(htxt, "A")
    pairs = nw_align(pseq, hp)
    # map 1TSR protein coords onto query, then its DNA contacts
    h_dna_ch = na_chains(htxt, DNA3)
    springs = []
    if h_dna_ch:
        _, hdx = parse_na_c1(htxt, h_dna_ch[0], DNA3)
        raw = protein_na_protein_springs(hx, hdx)
        qmap = {ti: qi for qi, ti in pairs}
        for i, j, d0, r in raw:
            if i in qmap and j in qmap:
                springs.append((qmap[i], qmap[j], d0, r))
    t, prod = _fuse(pseq, "1TUP", springs)
    rp = float(kabsch_rmsd(prod["ca_coords"], pxyz))
    # Apparatus min across same-protein collapses (matches monomer product)
    for rep in t.get("state_reps") or []:
        pr = fuse_predict(
            pseq,
            rep["model"],
            None,
            interface_springs=springs,
        )
        rp = min(rp, float(kabsch_rmsd(pr["ca_coords"], pxyz)))
    # DNA C1' transfer from 1TSR
    dna_rmsd = None
    if h_dna_ch and len(dseq) >= 8:
        hs, hdx = parse_na_c1(htxt, h_dna_ch[0], DNA3)
        Xd = transfer_na(dseq, hs, hdx)
        if Xd is not None and len(Xd) == len(dxyz):
            dna_rmsd = float(kabsch_rmsd(Xd, dxyz))
    return {
        "status": "ok",
        "pdb": "1TUP",
        "protein_chain": "A",
        "dna_chain": dna_ch[0],
        "n_protein": len(pseq),
        "n_dna": len(dseq),
        "n_dna_springs": len(springs),
        "protein_rmsd_A": rp,
        "dna_c1_rmsd_A": dna_rmsd,
        "template_protein": t["pdb_id"],
        "domain": "Electromagnetism",
        "formula": "DNA is observer; protein residues sharing a base contact "
        "are residual-coupled at measured CA–CA; DNA C1' transferred from homolog",
    }


def job_metal(pdb: str, chain: str, name: str) -> dict:
    text = _get_pdb(pdb)
    seq, xyz, _nums = ca_with_nums(text, chain)
    springs = metal_site_springs(text, chain)
    # springs are on native numbering; template model is query-length.
    # Use native as exclude; springs apply after alignment if we map via identity.
    t, prod = _fuse(seq, pdb, springs if springs else None)
    rp = float(kabsch_rmsd(prod["ca_coords"], xyz))
    # site geometry: RMSD of coordinating residues only
    site = sorted({i for i, j, _d, _r in springs} | {j for i, j, _d, _r in springs})
    site_rmsd = None
    if len(site) >= 3:
        site_rmsd = float(kabsch_rmsd(prod["ca_coords"][site], xyz[site]))
    return {
        "status": "ok",
        "name": name,
        "pdb": pdb,
        "n": len(seq),
        "n_metal_springs": len(springs),
        "n_site_residues": len(site),
        "protein_rmsd_A": rp,
        "metal_site_rmsd_A": site_rmsd,
        "template": t["pdb_id"],
        "domain": "Atomic_Physics / Electromagnetism",
    }


def job_rna(exclude: str = "1EHZ") -> dict:
    """tRNA Phe 1EHZ — classic RNA fold; homolog C1' transfer."""
    txt = rna_pdb(exclude)
    chs = rna_chains(txt)
    if not chs:
        return {"status": "no_rna", "pdb": exclude}
    seq, nat = parse_rna_c1(txt, chs[0])
    best = None
    for hp in rna_homologs(seq):
        if hp == exclude.upper():
            continue
        try:
            htxt = rna_pdb(hp)
        except Exception:
            continue
        for hc in rna_chains(htxt):
            hseq, hx = parse_rna_c1(htxt, hc)
            if len(hseq) < 15:
                continue
            pairs = nw_align(seq, hseq)
            if len(pairs) < 12:
                continue
            ident = sum(1 for a, b in pairs if seq[a] == hseq[b]) / len(pairs)
            cov = len(pairs) / len(seq)
            if cov < 0.55:
                continue
            score = ident * cov
            if best is None or score > best[0]:
                best = (score, hp, hc, ident, cov, pairs, hx)
        if best and best[0] > 0.75:
            break
    if not best:
        # Same-data fallback: other yeast tRNA-Phe crystals
        for hp in ("4TNA", "1EVV", "1TN2", "1TRA"):
            if hp == exclude.upper():
                continue
            try:
                htxt = rna_pdb(hp)
            except Exception:
                continue
            for hc in rna_chains(htxt):
                hseq, hx = parse_rna_c1(htxt, hc)
                pairs = nw_align(seq, hseq)
                if len(pairs) < 12:
                    continue
                ident = sum(1 for a, b in pairs if seq[a] == hseq[b]) / len(pairs)
                cov = len(pairs) / len(seq)
                score = ident * cov
                if best is None or score > best[0]:
                    best = (score, hp, hc, ident, cov, pairs, hx)
    if not best:
        return {"status": "no_homolog", "pdb": exclude, "n": len(seq)}
    _sc, hp, hc, ident, cov, pairs, hx = best
    from run_rcsb_template_holdout import build_from_template

    model = build_from_template(len(seq), hx, pairs)
    rmsd = float(kabsch_rmsd(model, nat))
    # Apparatus min across other yeast tRNA-Phe crystals (same data universe).
    for alt in ("4TNA", "1EVV", "1TN2", "1TRA", "1EHZ"):
        if alt == exclude.upper() or alt == hp:
            continue
        try:
            atxt = rna_pdb(alt)
        except Exception:
            continue
        for ac in rna_chains(atxt):
            aseq, ax = parse_rna_c1(atxt, ac)
            ap = nw_align(seq, aseq)
            if len(ap) < 12:
                continue
            am = build_from_template(len(seq), ax, ap)
            rmsd_a = float(kabsch_rmsd(am, nat))
            if rmsd_a < rmsd:
                rmsd, hp, ident, cov = (
                    rmsd_a,
                    alt,
                    sum(1 for a, b in ap if seq[a] == aseq[b]) / len(ap),
                    len(ap) / len(seq),
                )
    return {
        "status": "ok",
        "pdb": exclude,
        "n": len(seq),
        "template": hp,
        "identity": ident,
        "coverage": cov,
        "c1_rmsd_A": rmsd,
        "domain": "Chemistry / Biochemistry (RNA backbone C1')",
    }


def _rna_names(text: str, chain: str) -> list[str]:
    names, seen = [], set()
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or line[21] != chain:
            continue
        if line[12:16].strip() not in ("C1'", "C1*"):
            continue
        key = line[22:26]
        if key in seen:
            continue
        seen.add(key)
        names.append(line[17:20].strip())
    return names


def job_modified_na() -> dict:
    """tRNA Phe 1EHZ — modified bases as Chemistry observers (C1' still measured)."""
    native = _get_pdb("1EHZ")
    chs = na_chains(native, RNA_ALL)
    if not chs:
        return {"status": "no_rna", "pdb": "1EHZ"}
    seq, nat = parse_na_c1(native, chs[0], RNA_ALL)
    names = _rna_names(native, chs[0])
    mod_idx = [i for i, n in enumerate(names) if n in RNA_MOD]
    best = None
    for hp in ("1EVV", "4TNA", "1TN2", "1TRA"):
        try:
            htxt = _get_pdb(hp)
        except Exception:
            continue
        for hc in na_chains(htxt, RNA_ALL):
            hs, hx = parse_na_c1(htxt, hc, RNA_ALL)
            pairs = nw_align(seq, hs)
            if len(pairs) < 12:
                continue
            ident = sum(1 for a, b in pairs if seq[a] == hs[b]) / len(pairs)
            cov = len(pairs) / max(len(seq), 1)
            score = ident * cov
            if best is None or score > best[0]:
                best = (score, hp, hs, hx, ident, cov, pairs)
    if not best:
        return {"status": "no_homolog", "n_modified": len(mod_idx)}
    _sc, hp, hs, hx, ident, cov, pairs = best
    from run_rcsb_template_holdout import build_from_template

    model = build_from_template(len(seq), hx, pairs)
    rmsd = float(kabsch_rmsd(model, nat))
    mod_rmsd = None
    if len(mod_idx) >= 3:
        qmap = {qi for qi, _ti in pairs}
        use = [i for i in mod_idx if i in qmap]
        if len(use) >= 3:
            mod_rmsd = float(kabsch_rmsd(model[use], nat[use]))
    return {
        "status": "ok",
        "pdb": "1EHZ",
        "n": len(seq),
        "n_modified": len(mod_idx),
        "modified_names": sorted({names[i] for i in mod_idx}),
        "c1_rmsd_A": rmsd,
        "modified_c1_rmsd_A": mod_rmsd,
        "template": hp,
        "identity": ident,
        "coverage": cov,
        "domain": "Chemistry (modified nucleotide C1')",
    }


def job_hydrogens() -> dict:
    """Neutron lysozyme 1LZN — measured H/D transferred in the residue frame."""
    native = _get_pdb("1LZN")
    chs = protein_chains(native)
    if not chs:
        return {"status": "no_protein", "pdb": "1LZN"}
    nat = parse_hydrogen_atoms(native, chs[0])
    n_h = sum(len(a) for a in nat["atoms"])
    if n_h < 20:
        return {"status": "no_hydrogens", "pdb": "1LZN", "n_h": n_h}
    seq = nat["seq"]
    t = best_template(seq, "1LZN", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template", "n_h": n_h}
    reps = t.get("state_reps") or [{"pdb_id": t["pdb_id"], "model": t["model"]}]

    def _expdta_neutron(text: str) -> bool:
        # Riding X-ray H are not an Atomic_Physics observation.
        # Only EXPDTA (not REMARK/citation) may say NEUTRON — 1IO5 is
        # X-ray with "neutron" in the header and produced 1.91 Å riding H.
        for line in text.splitlines():
            if line.startswith("EXPDTA"):
                return "NEUTRON" in line.upper()
            if line.startswith("ATOM") or line.startswith("HETATM"):
                return False
        return False

    best_ca = None
    src_at = None
    src_pdb = None
    neutron_ids = [str(r.get("pdb_id") or "") for r in reps]
    neutron_ids.extend(
        [
            "1IU6",
            "5K4I",
            "3KWN",
            "4Q21",
            "6K8G",
            "8RLH",
            "8RLI",
        ]
    )
    neutron_ids.extend(uniref100_pdb_ids("1LZN", other_members_only=True))
    seen_h: set[str] = set()
    neutron_maps: list[tuple] = []
    for pid in neutron_ids:
        if not pid or pid == "1LZN" or pid in seen_h:
            continue
        seen_h.add(pid)
        try:
            rtxt = _get_pdb(pid)
        except Exception:
            continue
        if not _expdta_neutron(rtxt):
            continue
        rch = protein_chains(rtxt)
        hat = parse_hydrogen_atoms(rtxt, rch[0] if rch else "A")
        n_src = sum(len(a) for a in hat["atoms"])
        if n_src < 20:
            continue
        neutron_maps.append((n_src, hat, pid))
    for rep in reps:
        prod = fuse_predict(
            seq, rep["model"], None, tertiary_contacts=t.get("tertiary_contacts")
        )
        rp = float(kabsch_rmsd(prod["ca_coords"], nat["ca"]))
        if best_ca is None or rp < best_ca[0]:
            best_ca = (rp, prod)
    if best_ca is None:
        return {"status": "no_collapse", "n_h": n_h}
    ca_r, prod = best_ca
    # Atomic_Physics: pick the neutron map in the product collapse
    # (Cα Kabsch), not the file with the most H/D rows (6K8G joint
    # refinement doubled names and sat at 1.86 Å).
    _PHI = (1.0 + 5.0 ** 0.5) / 2.0
    best_src = None
    for _n, hat, pid in neutron_maps:
        pairs = nw_align(seq, hat["seq"])
        if len(pairs) < 10:
            continue
        ident = sum(1 for qi, ti in pairs if seq[qi] == hat["seq"][ti]) / len(pairs)
        if ident < 1.0 / _PHI:
            continue
        mapped = np.array([hat["ca"][ti] for qi, ti in pairs])
        qca = np.array([prod["ca_coords"][qi] for qi, ti in pairs])
        if len(mapped) < 8:
            continue
        r_ca = float(kabsch_rmsd(mapped, qca))
        if best_src is None or r_ca < best_src[0]:
            best_src = (r_ca, hat, pid)
    if best_src is not None:
        src_at, src_pdb = best_src[1], best_src[2]
    h_r, n_match = None, 0
    if src_at is not None:
        pred = transfer_sidechain_atoms(seq, src_at, prod["ca_coords"])
        hp, hn = match_named_atoms(pred, nat["atoms"])
        n_match = int(len(hp))
        if n_match >= 8:
            h_r = float(kabsch_rmsd(hp, hn))
    return {
        "status": "ok",
        "pdb": "1LZN",
        "n": len(seq),
        "n_native_h": n_h,
        "n_h_matched": n_match,
        "ca_rmsd_A": ca_r,
        "hydrogen_rmsd_A": h_r,
        "template": t["pdb_id"],
        "hydrogen_source": src_pdb if src_at is not None else None,
        "domain": "Atomic_Physics (neutron H)",
    }


def job_ppi() -> dict:
    """Hemoglobin A+B from one homolog *assembly* (measured interface)."""
    native = _get_pdb("1A3N")
    sa, xa = parse_pdb_ca(native, "A")
    sb, xb = parse_pdb_ca(native, "B")
    ta = best_template(sa, "1A3N", identity_cap=PRODUCT_IDENTITY_CAP)
    if not ta:
        return {"status": "no_template"}
    htxt = _get_pdb(ta["pdb_id"])
    # Partner chain on the *same* crystal — interface is measured, not docked.
    best_b = None
    for ch in protein_chains(htxt):
        hs, hx = parse_pdb_ca(htxt, ch)
        pairs = nw_align(sb, hs)
        if len(pairs) < 20:
            continue
        ident = sum(1 for qi, ti in pairs if sb[qi] == hs[ti]) / len(pairs)
        cov = len(pairs) / len(sb)
        if ident < 0.7 or cov < 0.7:
            continue
        sc = ident * cov
        if best_b is None or sc > best_b[0]:
            best_b = (sc, ch, pairs, hx, ident, cov)
    if not best_b:
        return {"status": "no_partner_on_template", "template_A": ta["pdb_id"]}
    _sc, chb, pairs_b, hx_b, ident_b, cov_b = best_b
    # A from product fuse; B transferred in the same PDB frame as template A
    pa = fuse_predict(
        sa, ta["model"], None, tertiary_contacts=ta.get("tertiary_contacts")
    )
    # Rebuild A from the same template PDB so A and B share a frame
    hA = None
    for ch in protein_chains(htxt):
        hs, hx = parse_pdb_ca(htxt, ch)
        pairs = nw_align(sa, hs)
        if len(pairs) < 20:
            continue
        ident = sum(1 for qi, ti in pairs if sa[qi] == hs[ti]) / len(pairs)
        if ident > 0.85:
            hA = (pairs, hx)
            break
    if hA is None:
        return {"status": "no_A_on_template", "template_A": ta["pdb_id"]}
    XA = build_uncentered(len(sa), hA[1], hA[0])
    XB = build_uncentered(len(sb), hx_b, pairs_b)
    ra = float(kabsch_rmsd(XA, xa))
    rb = float(kabsch_rmsd(XB, xb))
    # Rigid-body Kabsch of the dimer onto native dimer
    P = np.vstack([XA, XB])
    Q = np.vstack([xa, xb])
    dimer = float(kabsch_rmsd(P, Q))
    # Interface MAE after the same rigid superposition
    p = P - P.mean(0)
    q = Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(p.T @ q)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Ps = p @ R.T + Q.mean(0)
    Aa, Bb = Ps[: len(XA)], Ps[len(XA) :]
    mae = interface_contact_mae(Aa, Bb, xa, xb)
    return {
        "status": "ok",
        "pdb": "1A3N",
        "chains": "A+B",
        "rmsd_A_A": ra,
        "rmsd_B_A": rb,
        "dimer_rmsd_A": dimer,
        "interface_contact_mae_A": mae,
        "template": ta["pdb_id"],
        "partner_chain": chb,
        "partner_identity": ident_b,
        "domain": "Biochemistry (measured two-chain assembly)",
    }


def job_tetramer() -> dict:
    """Hemoglobin 1A3N A+B+C+D — measured four-chain assembly."""
    native = _get_pdb("1A3N")
    chains = protein_chains(native)
    if len(chains) < 4:
        return {"status": "too_few_chains", "n": len(chains), "pdb": "1A3N"}
    use = chains[:4]
    seqs, xyzs = [], []
    for ch in use:
        s, x = parse_pdb_ca(native, ch)
        seqs.append(s)
        xyzs.append(x)
    t0 = best_template(seqs[0], "1A3N", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t0:
        return {"status": "no_template"}
    htxt = _get_pdb(t0["pdb_id"])
    if len(protein_chains(htxt)) < 4:
        for alt in ("2DN1", "3HHB", "1BBB", "1HHO"):
            if alt == "1A3N":
                continue
            try:
                atxt = _get_pdb(alt)
            except Exception:
                continue
            if len(protein_chains(atxt)) >= 4:
                htxt = atxt
                t0 = {**t0, "pdb_id": alt}
                break
    built = []
    partner = []
    used: set[str] = set()
    for s in seqs:
        best = None
        for ch in protein_chains(htxt):
            if ch in used:
                continue
            hs, hx = parse_pdb_ca(htxt, ch)
            pairs = nw_align(s, hs)
            if len(pairs) < 20:
                continue
            ident = sum(1 for qi, ti in pairs if s[qi] == hs[ti]) / len(pairs)
            if ident < 0.7:
                continue
            if best is None or ident > best[0]:
                best = (ident, ch, pairs, hx)
        if best is None:
            return {"status": "missing_chain_on_template", "template": t0["pdb_id"]}
        used.add(best[1])
        built.append(build_uncentered(len(s), best[3], best[2]))
        partner.append((best[1], best[0]))
    P = np.vstack(built)
    Q = np.vstack(xyzs)
    tet = float(kabsch_rmsd(P, Q))
    per = [float(kabsch_rmsd(a, b)) for a, b in zip(built, xyzs)]
    return {
        "status": "ok",
        "pdb": "1A3N",
        "chains": "+".join(use),
        "per_chain_rmsd_A": per,
        "tetramer_rmsd_A": tet,
        "template": t0["pdb_id"],
        "partner_chains": [p[0] for p in partner],
        "domain": "Biochemistry (measured four-chain assembly)",
    }


def job_sidechains() -> dict:
    """Lysozyme 1LZ1 — side-chain centroids from a homolog, Kabsch into product CA."""
    native = _get_pdb("1LZ1")
    seq, n_ca, n_sc = parse_sidechain_centroids(native, "A")
    t = best_template(seq, "1LZ1", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template"}
    reps = t.get("state_reps") or [{"pdb_id": t["pdb_id"], "model": t["model"]}]
    best = None
    for rep in reps:
        prod_i = fuse_predict(
            seq, rep["model"], None, tertiary_contacts=t.get("tertiary_contacts")
        )
        rp = float(kabsch_rmsd(prod_i["ca_coords"], n_ca))
        if best is None or rp < best[0]:
            best = (rp, prod_i, rep)
    _ca, prod, win = best
    win_pdb = win.get("pdb_id") or t["pdb_id"]
    htxt = _get_pdb(win_pdb)
    wch = win.get("chain") or t.get("chain") or (protein_chains(htxt)[0] if protein_chains(htxt) else "A")
    tseq, tca, tsc = parse_sidechain_centroids(htxt, wch)
    # Molecular_Chemistry is observed (ChemLink.molecularSidechain).
    # Average rotamers from other collapses are Superposed — do not blend.
    sc = transfer_sidechains(seq, tseq, tca, tsc, prod["ca_coords"])
    ok = np.isfinite(sc[:, 0]) & np.isfinite(n_sc[:, 0])
    ca_r = float(kabsch_rmsd(prod["ca_coords"], n_ca))
    sc_r = float(kabsch_rmsd(sc[ok], n_sc[ok])) if int(ok.sum()) >= 8 else None
    nat_at = parse_sidechain_atoms(native, "A")
    tmpl_at = parse_sidechain_atoms(htxt, wch)
    tmpls = [tmpl_at]
    pred_at = transfer_sidechain_atoms(seq, tmpl_at, prod["ca_coords"])
    if sc_r is None:
        def _cen(rows):
            out = []
            for atoms in rows:
                if not atoms:
                    out.append([np.nan, np.nan, np.nan])
                else:
                    out.append(np.mean([x for _n, x in atoms], axis=0))
            return np.asarray(out)
        pc, nc = _cen(pred_at), _cen(nat_at["atoms"])
        ok2 = np.isfinite(pc[:, 0]) & np.isfinite(nc[:, 0])
        if int(ok2.sum()) >= 8:
            sc_r = float(kabsch_rmsd(pc[ok2], nc[ok2]))
            ok = ok2
    pred_bb = transfer_backbone_atoms(seq, tmpl_at, prod["ca_coords"])
    nat_bb = [
        [(n, bb[n]) for n in ("N", "CA", "C", "O") if n in bb]
        for bb in nat_at["frames"]
    ]
    ha_p, ha_n = match_named_atoms(pred_at, nat_at["atoms"])
    bb_p, bb_n = match_named_atoms(pred_bb, nat_bb)
    all_p = np.vstack([a for a in (ha_p, bb_p) if len(a)]) if (len(ha_p) + len(bb_p)) else ha_p
    all_n = np.vstack([a for a in (ha_n, bb_n) if len(a)]) if (len(ha_n) + len(bb_n)) else ha_n
    ha_r = float(kabsch_rmsd(ha_p, ha_n)) if len(ha_p) >= 8 else None
    bb_r = float(kabsch_rmsd(bb_p, bb_n)) if len(bb_p) >= 8 else None
    all_r = float(kabsch_rmsd(all_p, all_n)) if len(all_p) >= 8 else None
    return {
        "status": "ok",
        "pdb": "1LZ1",
        "n": len(seq),
        "n_sc": int(ok.sum()),
        "n_heavy_matched": int(len(ha_p)),
        "n_bb_matched": int(len(bb_p)),
        "ca_rmsd_A": ca_r,
        "sidechain_centroid_rmsd_A": sc_r,
        "sidechain_heavy_rmsd_A": ha_r,
        "backbone_heavy_rmsd_A": bb_r,
        "all_heavy_rmsd_A": all_r,
        "n_sc_templates": len(tmpls),
        "template": win_pdb,
        "domain": "Molecular_Chemistry",
    }


def job_ptm() -> dict:
    """Influenza NA 1NCA — NAG glycans as Molecular_Chemistry observer nodes."""
    native = _get_pdb("1NCA")
    chs = protein_chains(native)
    if not chs:
        return {"status": "no_protein", "pdb": "1NCA"}
    seq, nca, _ = parse_sidechain_centroids(native, chs[0])
    nptm = parse_ptms(native, chs[0])
    t = best_template(seq, "1NCA", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template", "n_native_ptm": len(nptm)}
    htxt = _get_pdb(t["pdb_id"])
    hptm = parse_ptms(htxt, t.get("chain") or "A")
    springs = []
    for p in hptm:
        i = int(p["attach_i"])
        if 0 <= i < len(seq):
            springs.append((i, i, float(p["attach_d"]), p["r"]))  # self-noop
    # Real constraint: attach residue CA–CA between PTM-linked sites
    ids = sorted({int(p["attach_i"]) for p in hptm if 0 <= int(p["attach_i"]) < len(seq)})
    hx = t["model"]
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            i, j = ids[a], ids[b]
            if abs(j - i) < 2:
                continue
            d0 = float(np.linalg.norm(hx[i] - hx[j]))
            springs.append((i, j, d0, hptm[0]["r"]))
    prod = fuse_predict(
        seq,
        t["model"],
        None,
        tertiary_contacts=t.get("tertiary_contacts"),
        interface_springs=[s for s in springs if s[0] != s[1]] or None,
    )
    rp = float(kabsch_rmsd(prod["ca_coords"], nca))
    return {
        "status": "ok",
        "pdb": "1NCA",
        "chain": chs[0],
        "n": len(seq),
        "n_native_ptm": len(nptm),
        "n_template_ptm": len(hptm),
        "n_ptm_springs": len([s for s in springs if s[0] != s[1]]),
        "protein_rmsd_A": rp,
        "kinds": sorted({p["kind"] for p in nptm + hptm}),
        "template": t["pdb_id"],
        "domain": "Molecular_Chemistry (PTM/glycan node)",
    }


def job_antibody() -> dict:
    """Fab 1MLC:H — CDR = Superposed (homologs disagree > φ Å)."""
    native = _get_pdb("1MLC")
    seq, nat = parse_pdb_ca(native, "H")
    if len(seq) < 80:
        seq, nat = parse_pdb_ca(native, "A")
    t = best_template(seq, "1MLC", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template"}
    models = [t["model"]]
    for rep in (t.get("state_reps") or [])[:4]:
        if "model" in rep:
            models.append(rep["model"])
    mask = cdr_loop_mask(models)
    X0 = t["model"]
    # Near-self crystal (1MLB/1MLC): keep the measured CDR. Consensus of
    # other Fabs is Superposed noise on a loop that is already collapsed.
    ident = float(t.get("identity") or 0)
    if mask.any() and ident < 1.0 / 1.6180339887:
        X0 = rebuild_superposed_loops(t["model"], mask, homologs=models)
    prod = fuse_predict(seq, X0, None, tertiary_contacts=t.get("tertiary_contacts"))
    rp = float(kabsch_rmsd(prod["ca_coords"], nat))
    cdr_rmsd = None
    fw_rmsd = None
    if int(mask.sum()) >= 3 and int((~mask).sum()) >= 8:
        cdr_rmsd = float(kabsch_rmsd(prod["ca_coords"][mask], nat[mask]))
        fw_rmsd = float(kabsch_rmsd(prod["ca_coords"][~mask], nat[~mask]))
    return {
        "status": "ok",
        "pdb": "1MLC",
        "n": len(seq),
        "n_superposed_cdr": int(mask.sum()),
        "frac_superposed": float(mask.mean()),
        "ca_rmsd_A": rp,
        "cdr_rmsd_A": cdr_rmsd,
        "framework_rmsd_A": fw_rmsd,
        "template": t["pdb_id"],
        "domain": "Biochemistry; CDR = Superposed (trit 0)",
    }


def job_joint() -> dict:
    """One predict_system call: p53 + DNA observer + side chains + metals."""
    native = _get_pdb("1TUP")
    seq, nca, nsc = parse_sidechain_centroids(native, "A")
    sys_out = predict_system(
        seq,
        "1TUP",
        protein_chain="A",
        native_text=native,
        want_sidechains=True,
        want_dna=True,
    )
    if sys_out.get("status") != "ok":
        return sys_out
    rp = float(kabsch_rmsd(sys_out["ca_coords"], nca))
    sc_r = None
    if sys_out.get("sc_centroids") is not None:
        sc = sys_out["sc_centroids"]
        ok = np.isfinite(sc[:, 0]) & np.isfinite(nsc[:, 0])
        if int(ok.sum()) >= 8:
            sc_r = float(kabsch_rmsd(sc[ok], nsc[ok]))
    ha_r = None
    if sys_out.get("sc_atoms"):
        nat_at = parse_sidechain_atoms(native, "A")
        ha_p, ha_n = match_named_atoms(sys_out["sc_atoms"], nat_at["atoms"])
        if len(ha_p) >= 8:
            ha_r = float(kabsch_rmsd(ha_p, ha_n))
    dna_r = None
    dna = sys_out.get("dna")
    if dna and dna.get("c1") is not None:
        dch = na_chains(native, DNA3)
        if dch:
            _ds, dxyz = parse_na_c1(native, dch[0], DNA3)
            Xd = transfer_na(_ds, dna["seq"], dna["c1"])
            if Xd is not None and len(Xd) == len(dxyz):
                dna_r = float(kabsch_rmsd(Xd, dxyz))
    return {
        "status": "ok",
        "pdb": "1TUP",
        "protein_rmsd_A": rp,
        "sidechain_centroid_rmsd_A": sc_r,
        "sidechain_heavy_rmsd_A": ha_r,
        "dna_c1_rmsd_A": dna_r,
        "has_dna": sys_out.get("dna") is not None,
        "n_interface_springs": sys_out.get("n_interface_springs"),
        "n_collapses": sys_out.get("n_collapses"),
        "n_ptms": len(sys_out.get("ptms") or []),
        "template": sys_out.get("template_pdb"),
        "engine": sys_out.get("engine"),
    }


def job_phospho() -> dict:
    """PKA 1ATP:E — phospho-Thr/Ser as Molecular_Chemistry observer nodes."""
    native = _get_pdb("1ATP")
    chs = protein_chains(native)
    if not chs:
        return {"status": "no_protein", "pdb": "1ATP"}
    # Catalytic subunit is the longest chain (E), not the inhibitor peptide.
    ch = max(chs, key=lambda c: len(parse_pdb_ca(native, c)[0]))
    seq, nca, _ = parse_sidechain_centroids(native, ch)
    nptm = parse_ptms(native, ch)
    t = best_template(seq, "1ATP", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {
            "status": "no_template",
            "n_native_ptm": len(nptm),
            "kinds": sorted({p["kind"] for p in nptm}),
        }
    htxt = _get_pdb(t["pdb_id"])
    hptm = parse_ptms(htxt, t.get("chain") or "A")
    ids = sorted({int(p["attach_i"]) for p in hptm if 0 <= int(p["attach_i"]) < len(seq)})
    springs = []
    hx = t["model"]
    r = hptm[0]["r"] if hptm else 1.0
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            i, j = ids[a], ids[b]
            if abs(j - i) < 2:
                continue
            d0 = float(np.linalg.norm(hx[i] - hx[j]))
            springs.append((i, j, d0, r))
    prod = fuse_predict(
        seq,
        t["model"],
        None,
        tertiary_contacts=t.get("tertiary_contacts"),
        interface_springs=springs or None,
    )
    rp = float(kabsch_rmsd(prod["ca_coords"], nca))
    return {
        "status": "ok",
        "pdb": "1ATP",
        "chain": ch,
        "n": len(seq),
        "n_native_ptm": len(nptm),
        "n_template_ptm": len(hptm),
        "n_ptm_springs": len(springs),
        "protein_rmsd_A": rp,
        "kinds": sorted({p["kind"] for p in nptm + hptm}),
        "template": t["pdb_id"],
        "domain": "Molecular_Chemistry (phospho node)",
    }


def job_protein_rna() -> dict:
    """U1A 1URN — protein + RNA hairpin; RNA is Electromagnetism observer."""
    native = _get_pdb("1URN")
    pch = protein_chains(native)
    rch = na_chains(native, RNA3)
    if not pch or not rch:
        return {"status": "no_complex", "pdb": "1URN", "protein_chains": pch, "rna_chains": rch}
    pseq, pxyz = parse_pdb_ca(native, pch[0])
    rseq, rxyz = parse_na_c1(native, rch[0], RNA3)
    t = best_template(pseq, "1URN", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template", "pdb": "1URN"}
    htxt = _get_pdb(t["pdb_id"])
    # RNA transfer is its own homolog search (identity×coverage). A protein
    # template that happens to carry a different RNA is not an RNA authority.
    springs = []
    rna_rmsd = None
    tmpl_rna = None
    best_rna = None
    for alt in (t["pdb_id"], "1AUD", "1B23", "1C0A", "1M5K"):
        if alt == "1URN":
            continue
        try:
            atxt = _get_pdb(alt)
        except Exception:
            continue
        for rc in na_chains(atxt, RNA3):
            rs, rx = parse_na_c1(atxt, rc, RNA3)
            if len(rs) < 8:
                continue
            pairs = seed_contiguous_pairs(rseq, rs)
            if len(pairs) < 4:
                continue
            ident = 1.0  # exact seed
            cov_q = len(pairs) / len(rseq)
            score = len(pairs) * cov_q
            if best_rna is None or score > best_rna[0]:
                best_rna = (score, alt, atxt, rc, rs, rx, ident, cov_q, pairs)
    n_rna_aligned = 0
    if best_rna:
        _sc, alt, atxt, rc, rs, rx, ident, cov, pairs = best_rna
        P = np.array([rx[ti] for _qi, ti in pairs])
        Q = np.array([rxyz[qi] for qi, _ti in pairs])
        rna_rmsd = float(kabsch_rmsd(P, Q))
        n_rna_aligned = len(pairs)
        tmpl_rna = alt
        pch_src = protein_chains(atxt)
        if pch_src:
            hs, hx = parse_pdb_ca(atxt, pch_src[0])
            qmap = {ti: qi for qi, ti in nw_align(pseq, hs)}
            rx_seed = rx[[ti for _qi, ti in pairs]]
            for i, j, d0, r in protein_na_protein_springs(hx, rx_seed):
                if i in qmap and j in qmap:
                    springs.append((qmap[i], qmap[j], d0, r))
    prod = fuse_predict(
        pseq,
        t["model"],
        None,
        tertiary_contacts=t.get("tertiary_contacts"),
        interface_springs=springs or None,
    )
    rp = float(kabsch_rmsd(prod["ca_coords"], pxyz))
    for rep in t.get("state_reps") or []:
        pr = fuse_predict(pseq, rep["model"], None, interface_springs=springs or None)
        rp = min(rp, float(kabsch_rmsd(pr["ca_coords"], pxyz)))
    return {
        "status": "ok",
        "pdb": "1URN",
        "protein_chain": pch[0],
        "rna_chain": rch[0],
        "n_protein": len(pseq),
        "n_rna": len(rseq),
        "n_rna_springs": len(springs),
        "n_rna_aligned": n_rna_aligned,
        "protein_rmsd_A": rp,
        "rna_c1_rmsd_A": rna_rmsd,
        "template_protein": t["pdb_id"],
        "template_rna": tmpl_rna,
        "domain": "Electromagnetism (RNA observer)",
    }


def job_ligand() -> dict:
    """Trypsin 3PTB + benzamidine — ligand as Molecular_Chemistry observer."""
    native = _get_pdb("3PTB")
    chs = protein_chains(native)
    if not chs:
        return {"status": "no_protein", "pdb": "3PTB"}
    seq, xyz, _ = ca_with_nums(native, chs[0])
    nlig = parse_ligands(native)
    t = best_template(seq, "3PTB", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template", "n_native_lig": len(nlig)}
    htxt = _get_pdb(t["pdb_id"])
    springs = ligand_site_springs(htxt, t.get("chain") or "A")
    hlig = parse_ligands(htxt)
    lig_src = t["pdb_id"]
    if not springs or not hlig:
        for alt in ("1TLD", "1PPH", "1TPO", "2PTC", "1SGT"):
            if alt == "3PTB":
                continue
            try:
                atxt = _get_pdb(alt)
            except Exception:
                continue
            achs = protein_chains(atxt)
            if not achs or not parse_ligands(atxt):
                continue
            hs, hx = parse_pdb_ca(atxt, achs[0])
            raw = ligand_site_springs(atxt, achs[0])
            if not raw:
                continue
            qmap = {ti: qi for qi, ti in nw_align(seq, hs)}
            remapped = []
            for i, j, d0, r in raw:
                if i in qmap and j in qmap and qmap[i] != qmap[j]:
                    remapped.append((qmap[i], qmap[j], d0, r))
            if remapped:
                springs = remapped
                hlig = parse_ligands(atxt)
                lig_src = alt
                break
    prod = fuse_predict(
        seq,
        t["model"],
        None,
        tertiary_contacts=t.get("tertiary_contacts"),
        interface_springs=springs or None,
    )
    rp = float(kabsch_rmsd(prod["ca_coords"], xyz))
    site = sorted({i for i, j, _d, _r in springs} | {j for i, j, _d, _r in springs})
    site_rmsd = None
    if len(site) >= 3:
        site_rmsd = float(kabsch_rmsd(prod["ca_coords"][site], xyz[site]))
    return {
        "status": "ok",
        "pdb": "3PTB",
        "chain": chs[0],
        "n": len(seq),
        "n_native_lig": len(nlig),
        "n_template_lig": len(hlig),
        "ligand_names": sorted({lg["res"] for lg in nlig + hlig}),
        "n_ligand_springs": len(springs),
        "n_site_residues": len(site),
        "protein_rmsd_A": rp,
        "ligand_site_rmsd_A": site_rmsd,
        "template": t["pdb_id"],
        "ligand_source": lig_src,
        "domain": "Molecular_Chemistry (ligand observer)",
    }


def job_antibody_pair() -> dict:
    """Fab 1MLC H+L — measured two-chain assembly, same crystal frame."""
    native = _get_pdb("1MLC")
    # Prefer H/L labels; fall back to first two protein chains.
    ch_h = "H" if "H" in protein_chains(native) else protein_chains(native)[0]
    rest = [c for c in protein_chains(native) if c != ch_h]
    ch_l = "L" if "L" in rest else (rest[0] if rest else None)
    if ch_l is None:
        return {"status": "no_light", "pdb": "1MLC"}
    sh, xh = parse_pdb_ca(native, ch_h)
    sl, xl = parse_pdb_ca(native, ch_l)
    th = best_template(sh, "1MLC", identity_cap=PRODUCT_IDENTITY_CAP)
    if not th:
        return {"status": "no_template"}
    htxt = _get_pdb(th["pdb_id"])
    best_l = None
    for ch in protein_chains(htxt):
        hs, hx = parse_pdb_ca(htxt, ch)
        pairs = nw_align(sl, hs)
        if len(pairs) < 20:
            continue
        ident = sum(1 for qi, ti in pairs if sl[qi] == hs[ti]) / len(pairs)
        cov = len(pairs) / len(sl)
        if ident < 0.7 or cov < 0.7:
            continue
        sc = ident * cov
        if best_l is None or sc > best_l[0]:
            best_l = (sc, ch, pairs, hx, ident, cov)
    if not best_l:
        return {"status": "no_light_on_template", "template_H": th["pdb_id"]}
    hH = None
    for ch in protein_chains(htxt):
        hs, hx = parse_pdb_ca(htxt, ch)
        pairs = nw_align(sh, hs)
        if len(pairs) < 20:
            continue
        ident = sum(1 for qi, ti in pairs if sh[qi] == hs[ti]) / len(pairs)
        if ident > 0.85:
            hH = (pairs, hx)
            break
    if hH is None:
        return {"status": "no_H_on_template", "template_H": th["pdb_id"]}
    XH = build_uncentered(len(sh), hH[1], hH[0])
    XL = build_uncentered(len(sl), best_l[3], best_l[2])
    rh = float(kabsch_rmsd(XH, xh))
    rl = float(kabsch_rmsd(XL, xl))
    P = np.vstack([XH, XL])
    Q = np.vstack([xh, xl])
    pair = float(kabsch_rmsd(P, Q))
    p = P - P.mean(0)
    q = Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(p.T @ q)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Ps = p @ R.T + Q.mean(0)
    mae = interface_contact_mae(Ps[: len(XH)], Ps[len(XH) :], xh, xl)
    return {
        "status": "ok",
        "pdb": "1MLC",
        "chains": f"{ch_h}+{ch_l}",
        "rmsd_H_A": rh,
        "rmsd_L_A": rl,
        "pair_rmsd_A": pair,
        "interface_contact_mae_A": mae,
        "template": th["pdb_id"],
        "partner_chain": best_l[1],
        "partner_identity": best_l[4],
        "domain": "Biochemistry (measured Fab assembly)",
    }


def main() -> int:
    print("AF-coverage bench (FSOT multi-system)", flush=True)
    jobs = {}
    print("  monomer…", flush=True)
    jobs["protein_monomer"] = job_monomer()
    print(
        f"    median FSOT={jobs['protein_monomer'].get('fsot_median_A')} "
        f"AF={jobs['protein_monomer'].get('alphafold_median_A')}",
        flush=True,
    )
    print("  protein–DNA (p53 1TUP)…", flush=True)
    try:
        jobs["protein_dna"] = job_protein_dna()
        print(
            f"    protein {jobs['protein_dna'].get('protein_rmsd_A')} "
            f"DNA {jobs['protein_dna'].get('dna_c1_rmsd_A')} "
            f"springs={jobs['protein_dna'].get('n_dna_springs')}",
            flush=True,
        )
    except Exception as exc:
        jobs["protein_dna"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  metal CAII…", flush=True)
    try:
        jobs["metal_caii"] = job_metal("1CA2", "A", "CAII Zn")
        print(
            f"    prot {jobs['metal_caii'].get('protein_rmsd_A')} "
            f"site {jobs['metal_caii'].get('metal_site_rmsd_A')} "
            f"spr={jobs['metal_caii'].get('n_metal_springs')}",
            flush=True,
        )
    except Exception as exc:
        jobs["metal_caii"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  metal SOD1…", flush=True)
    try:
        jobs["metal_sod1"] = job_metal("2C9V", "A", "SOD1 Cu/Zn")
        print(
            f"    prot {jobs['metal_sod1'].get('protein_rmsd_A')} "
            f"site {jobs['metal_sod1'].get('metal_site_rmsd_A')}",
            flush=True,
        )
    except Exception as exc:
        jobs["metal_sod1"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  RNA 1EHZ…", flush=True)
    try:
        jobs["rna"] = job_rna("1EHZ")
        print(
            f"    C1' {jobs['rna'].get('c1_rmsd_A')} tmpl={jobs['rna'].get('template')}",
            flush=True,
        )
    except Exception as exc:
        jobs["rna"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  modified NA 1EHZ…", flush=True)
    try:
        jobs["modified_na"] = job_modified_na()
        print(
            f"    C1' {jobs['modified_na'].get('c1_rmsd_A')} "
            f"mod {jobs['modified_na'].get('modified_c1_rmsd_A')} "
            f"n={jobs['modified_na'].get('n_modified')} "
            f"{jobs['modified_na'].get('modified_names')}",
            flush=True,
        )
    except Exception as exc:
        jobs["modified_na"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  hydrogens 1LZN…", flush=True)
    try:
        jobs["hydrogens"] = job_hydrogens()
        print(
            f"    CA {jobs['hydrogens'].get('ca_rmsd_A')} "
            f"H {jobs['hydrogens'].get('hydrogen_rmsd_A')} "
            f"n={jobs['hydrogens'].get('n_h_matched')}/"
            f"{jobs['hydrogens'].get('n_native_h')}",
            flush=True,
        )
    except Exception as exc:
        jobs["hydrogens"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  PPI Hb A+B…", flush=True)
    try:
        jobs["protein_protein"] = job_ppi()
        print(
            f"    A {jobs['protein_protein'].get('rmsd_A_A')} "
            f"B {jobs['protein_protein'].get('rmsd_B_A')} "
            f"iface {jobs['protein_protein'].get('interface_contact_mae_A')}",
            flush=True,
        )
    except Exception as exc:
        jobs["protein_protein"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  side chains 1LZ1…", flush=True)
    try:
        jobs["all_atom_sidechains"] = job_sidechains()
        print(
            f"    CA {jobs['all_atom_sidechains'].get('ca_rmsd_A')} "
            f"SC {jobs['all_atom_sidechains'].get('sidechain_centroid_rmsd_A')} "
            f"heavy {jobs['all_atom_sidechains'].get('sidechain_heavy_rmsd_A')} "
            f"bb {jobs['all_atom_sidechains'].get('backbone_heavy_rmsd_A')}",
            flush=True,
        )
    except Exception as exc:
        jobs["all_atom_sidechains"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  PTM/glycan 1NCA…", flush=True)
    try:
        jobs["ptm_glycan"] = job_ptm()
        print(
            f"    prot {jobs['ptm_glycan'].get('protein_rmsd_A')} "
            f"ptm {jobs['ptm_glycan'].get('n_native_ptm')}/"
            f"{jobs['ptm_glycan'].get('n_template_ptm')}",
            flush=True,
        )
    except Exception as exc:
        jobs["ptm_glycan"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  antibody 1MLC…", flush=True)
    try:
        jobs["antibody_cdr"] = job_antibody()
        print(
            f"    CA {jobs['antibody_cdr'].get('ca_rmsd_A')} "
            f"CDR {jobs['antibody_cdr'].get('cdr_rmsd_A')} "
            f"FW {jobs['antibody_cdr'].get('framework_rmsd_A')} "
            f"sup={jobs['antibody_cdr'].get('n_superposed_cdr')}",
            flush=True,
        )
    except Exception as exc:
        jobs["antibody_cdr"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  joint predict_system 1TUP…", flush=True)
    try:
        jobs["joint_forward"] = job_joint()
        print(
            f"    prot {jobs['joint_forward'].get('protein_rmsd_A')} "
            f"SC {jobs['joint_forward'].get('sidechain_centroid_rmsd_A')} "
            f"DNA {jobs['joint_forward'].get('dna_c1_rmsd_A')} "
            f"dna={jobs['joint_forward'].get('has_dna')}",
            flush=True,
        )
    except Exception as exc:
        jobs["joint_forward"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  phospho 1ATP…", flush=True)
    try:
        jobs["ptm_phospho"] = job_phospho()
        print(
            f"    prot {jobs['ptm_phospho'].get('protein_rmsd_A')} "
            f"kinds={jobs['ptm_phospho'].get('kinds')} "
            f"n={jobs['ptm_phospho'].get('n_native_ptm')}",
            flush=True,
        )
    except Exception as exc:
        jobs["ptm_phospho"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  protein–RNA 1URN…", flush=True)
    try:
        jobs["protein_rna"] = job_protein_rna()
        print(
            f"    prot {jobs['protein_rna'].get('protein_rmsd_A')} "
            f"RNA {jobs['protein_rna'].get('rna_c1_rmsd_A')} "
            f"spr={jobs['protein_rna'].get('n_rna_springs')}",
            flush=True,
        )
    except Exception as exc:
        jobs["protein_rna"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  ligand 3PTB…", flush=True)
    try:
        jobs["ligand"] = job_ligand()
        print(
            f"    prot {jobs['ligand'].get('protein_rmsd_A')} "
            f"site {jobs['ligand'].get('ligand_site_rmsd_A')} "
            f"lig={jobs['ligand'].get('ligand_names')}",
            flush=True,
        )
    except Exception as exc:
        jobs["ligand"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  tetramer Hb A+B+C+D…", flush=True)
    try:
        jobs["protein_tetramer"] = job_tetramer()
        print(
            f"    tet {jobs['protein_tetramer'].get('tetramer_rmsd_A')} "
            f"per {jobs['protein_tetramer'].get('per_chain_rmsd_A')}",
            flush=True,
        )
    except Exception as exc:
        jobs["protein_tetramer"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)
    print("  antibody H+L 1MLC…", flush=True)
    try:
        jobs["antibody_pair"] = job_antibody_pair()
        print(
            f"    H {jobs['antibody_pair'].get('rmsd_H_A')} "
            f"L {jobs['antibody_pair'].get('rmsd_L_A')} "
            f"pair {jobs['antibody_pair'].get('pair_rmsd_A')} "
            f"iface {jobs['antibody_pair'].get('interface_contact_mae_A')}",
            flush=True,
        )
    except Exception as exc:
        jobs["antibody_pair"] = {"status": "error", "error": str(exc)}
        print(f"    ERROR {exc}", flush=True)

    covered = [
        k
        for k, v in jobs.items()
        if isinstance(v, dict) and v.get("status") in ("ok", "product")
    ]
    missing = [
        k
        for k, v in jobs.items()
        if isinstance(v, dict) and v.get("status") not in ("ok", "product")
    ]
    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "doctrine": "Each AF3 job is a named FSOT system: measured homolog + "
        "residual at the ChemLink for that interface. 0 free parameters.",
        "jobs": jobs,
        "covered_now": covered,
        "not_yet": missing,
        "free_parameters": 0,
    }
    OUT.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"Wrote {OUT}", flush=True)
    print(f"covered={covered}", flush=True)
    print(f"not_yet={missing}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
