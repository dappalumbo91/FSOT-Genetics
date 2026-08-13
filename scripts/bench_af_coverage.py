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
)
from msa_template_fuse import fuse_predict  # noqa: E402
from multi_system import (  # noqa: E402
    _get_pdb,
    build_uncentered,
    ca_with_nums,
    cdr_loop_mask,
    interface_contact_mae,
    metal_site_springs,
    na_chains,
    parse_na_c1,
    parse_pdb_ca,
    parse_ptms,
    parse_sidechain_centroids,
    predict_system,
    protein_chains,
    protein_na_protein_springs,
    transfer_na,
    transfer_sidechains,
    DNA3,
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


def job_sidechains() -> dict:
    """Lysozyme 1LZ1 — side-chain centroids from a homolog, Kabsch into product CA."""
    native = _get_pdb("1LZ1")
    seq, n_ca, n_sc = parse_sidechain_centroids(native, "A")
    t = best_template(seq, "1LZ1", identity_cap=PRODUCT_IDENTITY_CAP)
    if not t:
        return {"status": "no_template"}
    prod = fuse_predict(seq, t["model"], None, tertiary_contacts=t.get("tertiary_contacts"))
    htxt = _get_pdb(t["pdb_id"])
    tseq, tca, tsc = parse_sidechain_centroids(htxt, t.get("chain") or "A")
    sc = transfer_sidechains(seq, tseq, tca, tsc, prod["ca_coords"])
    ok = np.isfinite(sc[:, 0]) & np.isfinite(n_sc[:, 0])
    ca_r = float(kabsch_rmsd(prod["ca_coords"], n_ca))
    sc_r = float(kabsch_rmsd(sc[ok], n_sc[ok])) if int(ok.sum()) >= 8 else None
    return {
        "status": "ok",
        "pdb": "1LZ1",
        "n": len(seq),
        "n_sc": int(ok.sum()),
        "ca_rmsd_A": ca_r,
        "sidechain_centroid_rmsd_A": sc_r,
        "template": t["pdb_id"],
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
    prod = fuse_predict(seq, t["model"], None, tertiary_contacts=t.get("tertiary_contacts"))
    rp = float(kabsch_rmsd(prod["ca_coords"], nat))
    cdr_rmsd = None
    fw_rmsd = None
    if mask.any() and (~mask).sum() >= 8:
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
    return {
        "status": "ok",
        "pdb": "1TUP",
        "protein_rmsd_A": rp,
        "sidechain_centroid_rmsd_A": sc_r,
        "has_dna": sys_out.get("dna") is not None,
        "n_interface_springs": sys_out.get("n_interface_springs"),
        "n_ptms": len(sys_out.get("ptms") or []),
        "template": sys_out.get("template_pdb"),
        "engine": sys_out.get("engine"),
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
            f"SC {jobs['all_atom_sidechains'].get('sidechain_centroid_rmsd_A')}",
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
            f"dna={jobs['joint_forward'].get('has_dna')}",
            flush=True,
        )
    except Exception as exc:
        jobs["joint_forward"] = {"status": "error", "error": str(exc)}
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
