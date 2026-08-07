#!/usr/bin/env python3
"""Multi-scale D_eff interface routing for protein fold (FSOT pin table only).

Doctrine:
  - D_eff comes from named DomainConfig rows — never free continuous D.
  - Different F-layers attach to different physical domains (multi-interface).
  - Residual-at-interface: wrong domain → worse residual; report honestly.
  - Zero free parameters.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "vendor"))

import fsot_compute as fc  # noqa: E402

PI = float(fc.PI)
E = float(fc.E)
PHI = float(fc.PHI)
P_NEW = float(fc.P_NEW)
C_EFF = float(fc.C_EFF)
ETA_EFF = float(fc.ETA_EFF)

# Protein-relevant named domains (subset of 35-domain table)
PROTEIN_SCALE_DOMAINS = (
    "Atomic_Physics",
    "Physical_Chemistry",
    "Chemistry",
    "Molecular_Chemistry",
    "Materials_Science",
    "Biology",
    "Biochemistry",
    "Condensed_Matter",
    "Neuroscience",
)


@dataclass(frozen=True)
class DomainSlice:
    name: str
    D_eff: int
    hits: int
    delta_psi: float
    observed: bool
    S: float
    abs_S: float


def domain_slice(name: str) -> DomainSlice:
    d = fc.DOMAINS[name]
    s = float(fc.domain_scalar(name))
    return DomainSlice(
        name=name,
        D_eff=int(d.D_eff),
        hits=int(d.hits),
        delta_psi=float(d.delta_psi),
        observed=bool(d.observed),
        S=s,
        abs_S=abs(s),
    )


def long_range_gate(D_eff: int) -> int:
    """F13 gate: ceil(η_eff · D) — same law as protein derivations, variable D."""
    return int(math.ceil(ETA_EFF * float(D_eff)))


@dataclass(frozen=True)
class InterfaceRouting:
    """Named multi-interface routing — theory labels, not free dials."""

    name: str
    chem_domain: str
    ss_domain: str
    region_domain: str
    packing_domain: str
    notes: str

    def resolve(self) -> dict[str, Any]:
        chem = domain_slice(self.chem_domain)
        ss = domain_slice(self.ss_domain)
        reg = domain_slice(self.region_domain)
        pack = domain_slice(self.packing_domain)
        chem_amp = chem.abs_S * P_NEW
        ss_amp = ss.abs_S * P_NEW
        region_amp = reg.abs_S * P_NEW * C_EFF
        packing_amp = pack.abs_S * P_NEW
        gate = long_range_gate(reg.D_eff)
        return {
            "routing": self.name,
            "notes": self.notes,
            "chem": asdict(chem),
            "ss": asdict(ss),
            "region": asdict(reg),
            "packing": asdict(pack),
            "chem_amp": chem_amp,
            "ss_amp": ss_amp,
            "region_amp": region_amp,
            "packing_amp": packing_amp,
            "long_range_gate": gate,
            "free_parameters": 0,
        }


# Theory-first routings (all domain names from pin table)
ROUTINGS: dict[str, InterfaceRouting] = {
    # Legacy v7 single-pair amplitudes
    "legacy_v7": InterfaceRouting(
        name="legacy_v7",
        chem_domain="Molecular_Chemistry",
        ss_domain="Molecular_Chemistry",
        region_domain="Biochemistry",
        packing_domain="Biochemistry",
        notes="Original protein derivations: chem+SS → molchem; region → biochem D=13",
    ),
    # Multi-scale ladder (default v9)
    "multi_scale_v9": InterfaceRouting(
        name="multi_scale_v9",
        chem_domain="Molecular_Chemistry",
        ss_domain="Chemistry",
        region_domain="Biochemistry",
        packing_domain="Condensed_Matter",
        notes="Residue chem D=9; H-bond SS D=8; tertiary D=13; packing D=14",
    ),
    # Polymer / physical chemistry heavy
    "polymer_physchem": InterfaceRouting(
        name="polymer_physchem",
        chem_domain="Physical_Chemistry",
        ss_domain="Chemistry",
        region_domain="Materials_Science",
        packing_domain="Condensed_Matter",
        notes="Backbone/polymer-leaning interfaces; packing condensed",
    ),
    # Cellular / biology context for regions
    "bio_context": InterfaceRouting(
        name="bio_context",
        chem_domain="Molecular_Chemistry",
        ss_domain="Chemistry",
        region_domain="Biology",
        packing_domain="Biochemistry",
        notes="Region gate at Biology D=12 (cellular context); packing biochem",
    ),
    # Dense packing emphasis
    "packing_dense": InterfaceRouting(
        name="packing_dense",
        chem_domain="Molecular_Chemistry",
        ss_domain="Chemistry",
        region_domain="Condensed_Matter",
        packing_domain="Condensed_Matter",
        notes="Tertiary+packing both Condensed_Matter D=14",
    ),
    # Atomic local + biochem global
    "atomic_to_biochem": InterfaceRouting(
        name="atomic_to_biochem",
        chem_domain="Atomic_Physics",
        ss_domain="Physical_Chemistry",
        region_domain="Biochemistry",
        packing_domain="Condensed_Matter",
        notes="Local bond-scale chem; global biochem regions",
    ),
}

# Claim default = protein derivation authority (legacy_v7).
# multi_scale_v9 and others remain for residual-at-interface diagnosis
# (see run_deff_interface_probe.py) — not auto-switched by lowest RMSD.
DEFAULT_ROUTING = "legacy_v7"


def get_routing(name: str | None = None) -> dict[str, Any]:
    key = name or DEFAULT_ROUTING
    if key not in ROUTINGS:
        raise KeyError(f"unknown routing {key}; choose from {list(ROUTINGS)}")
    return ROUTINGS[key].resolve()


def dump_protein_domains() -> list[dict[str, Any]]:
    return [asdict(domain_slice(n)) for n in PROTEIN_SCALE_DOMAINS]


def main() -> int:
    print("FSOT protein-scale domain interface (pin D1D38A)")
    print(f"  P_NEW={P_NEW:.6f}  C_EFF={C_EFF:.6f}  ETA_EFF={ETA_EFF:.6f}")
    print()
    print(f"{'domain':22s} {'D':>3s} {'δψ':>6s} {'|S|':>10s}  gate@D")
    for row in dump_protein_domains():
        g = long_range_gate(row["D_eff"])
        print(
            f"{row['name']:22s} {row['D_eff']:3d} {row['delta_psi']:6.3f} "
            f"{row['abs_S']:10.6f}  {g}"
        )
    print()
    print("Lawful routings:")
    out = {"domains": dump_protein_domains(), "routings": {}}
    for name, r in ROUTINGS.items():
        res = r.resolve()
        out["routings"][name] = res
        mark = " (DEFAULT)" if name == DEFAULT_ROUTING else ""
        print(
            f"  {name}{mark}: chem={res['chem']['name']}(D={res['chem']['D_eff']}) "
            f"ss={res['ss']['name']}(D={res['ss']['D_eff']}) "
            f"reg={res['region']['name']}(D={res['region']['D_eff']}) "
            f"gate={res['long_range_gate']} "
            f"chem_amp={res['chem_amp']:.5f} region_amp={res['region_amp']:.5f}"
        )
    path = ROOT / "data" / "domain_interface_table.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nWrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
