#!/usr/bin/env python3
"""
Pre-flight entrant audit CLI (ticket 02). Run before the benchmark.

    python audit_entrants.py

Writes `output/entrant_audit.csv` and `docs/entrant_audit.md`, and fails loudly
if any entrant violates the pretrained-weights-only constraint.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "_lib"))
from entrant_audit import admitted, check_audit_gate, excluded, write_audit  # noqa: E402
from utils import OUTPUT_DIR, WORKFLOW_ROOT  # noqa: E402


def main():
    check_audit_gate()
    write_audit(OUTPUT_DIR / "entrant_audit.csv", WORKFLOW_ROOT / "docs" / "entrant_audit.md")
    print(f"[audit] admitted: {[r['model'] for r in admitted()]}")
    print(f"[audit] excluded: {[r['model'] for r in excluded()]}")
    print("audit_entrants: OK")


if __name__ == "__main__":
    main()
