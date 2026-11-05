"""Offline contract for the per-model fidelity ledgers (ticket 09/17)."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
CODE = HERE.parent.parent / "script"  # execution/script
sys.path.insert(0, str(CODE))

from repo_root import BENCH, EXECUTION  # noqa: E402

import aggregate_fidelity  # noqa: E402

ALLOWED = {"config_only", "adapter_extension", "recipe_omission", "recipe_substitution", "reproducibility_break"}


class FidelityLedgerTest(unittest.TestCase):
    def test_every_model_has_the_three_baseline_rows(self) -> None:
        for model in aggregate_fidelity.MODELS:
            with self.subTest(model=model):
                rows, problems = aggregate_fidelity.parse(
                    EXECUTION / "models" / model / "FIDELITY.md", model
                )
                self.assertEqual(problems, [])
                self.assertGreaterEqual(len(rows), 10)
                components = " ".join(row["component"].lower() for row in rows)
                for required in ("gene space", "corpus scale", "initialisation"):
                    self.assertIn(required, components, f"{model}: missing {required}")
                for row in rows:
                    self.assertIn(row["class"], ALLOWED, row)
                    self.assertTrue(row["expected_effect"], row)
                    self.assertTrue(row["ticket"], row)


if __name__ == "__main__":
    unittest.main()
