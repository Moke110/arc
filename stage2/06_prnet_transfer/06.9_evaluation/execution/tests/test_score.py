"""Score-driver tests: the test-phase gate, the plan, and the resume predicate."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import _support  # noqa: E402
import score  # noqa: E402


class ScoreGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root = _support.clean_dir(_support.tmp_root() / "score")
        cls.pop = _support.write_synthetic_population(cls.root, drugs_per_fold=2, folds=(0,))

    def test_test_phase_requires_explicit_confirmation(self):
        code = score.main(
            [
                "--units", str(self.pop["units_path"]),
                "--phase", "test",
                "--out", str(self.root / "out"),
                "--plan",
            ]
        )
        self.assertEqual(code, 4)

    def test_plan_lists_units_and_scores_nothing(self):
        import contextlib
        import io

        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = score.main(
                [
                    "--units", str(self.pop["units_path"]),
                    "--phase", "test",
                    "--out", str(self.root / "out"),
                    "--plan",
                    "--confirm-frozen-selection",
                ]
            )
        self.assertEqual(code, 0)
        self.assertEqual(buffer.getvalue().count("e9_f"), 6)

    def test_subset_limit_is_honoured_in_plan(self):
        import contextlib
        import io

        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = score.main(
                [
                    "--units", str(self.pop["units_path"]),
                    "--phase", "test",
                    "--out", str(self.root / "out"),
                    "--plan",
                    "--confirm-frozen-selection",
                    "--limit", "4",
                ]
            )
        self.assertEqual(code, 0)
        self.assertEqual(buffer.getvalue().count("e9_f"), 4)

    def test_unknown_unit_id_is_rejected(self):
        code = score.main(
            [
                "--units", str(self.pop["units_path"]),
                "--phase", "test",
                "--out", str(self.root / "out"),
                "--plan",
                "--confirm-frozen-selection",
                "--only-unit", "does-not-exist",
            ]
        )
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
