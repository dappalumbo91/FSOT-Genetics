#!/usr/bin/env python3
"""Python replay of SMT catalog bounds (no z3 required)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OBL = ROOT / "verification" / "obligations.json"


def main() -> int:
    doc = json.loads(OBL.read_text(encoding="utf-8"))
    failed = []
    for o in doc.get("obligations") or []:
        lhs, rhs, kind = int(o["lhs"]), int(o["rhs"]), o["kind"]
        ok = (
            (kind == "nat_eq" and lhs == rhs)
            or (kind == "nat_lt" and lhs < rhs)
            or (kind == "nat_le" and lhs <= rhs)
        )
        if not ok:
            failed.append(o["id"])
    print(f"smt_python  obligations={len(doc.get('obligations') or [])}  fail={len(failed)}")
    if failed:
        print("FAIL", failed)
        return 1
    print("FSOT_GENETICS_SMT_PYTHON_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
