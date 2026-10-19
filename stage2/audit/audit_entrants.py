#!/usr/bin/env python3
"""
Entrant audit CLI (ticket 02, re-run). This is the selection record, not a
benchmark runner.

    python audit_entrants.py

Writes `entrant_audit.csv` and `entrant_audit.md` next to this script, and fails
loudly if any verdict violates the capability gates, compute tier, or 16 h budget.
"""

import sys
from pathlib import Path

AUDIT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(AUDIT_DIR))
from entrant_audit import admitted, check_audit_gate, excluded, reference, write_audit  # noqa: E402


def main():
    check_audit_gate()
    write_audit(AUDIT_DIR / "entrant_audit.csv", AUDIT_DIR / "entrant_audit.md")
    print(f"[audit] admitted: {[r['model'] for r in admitted()]}")
    print(f"[audit] reference: {[r['model'] for r in reference()]}")
    print(f"[audit] excluded: {[r['model'] for r in excluded()]}")
    print("audit_entrants: OK")


if __name__ == "__main__":
    main()
