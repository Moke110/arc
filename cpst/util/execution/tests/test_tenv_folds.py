"""Verify the released tenv partition identities at the shared reader boundary."""

from collections import Counter
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

from cpst_util.tenv_folds import COLD_DOSE_DRUG, COLD_DRUG, DEFAULT, open_tenv_folds


class TenvFoldsTest(unittest.TestCase):
    def setUp(self):
        self.release = open_tenv_folds()

    def test_parent_rotation_covers_every_condition(self):
        protocol = self.release.protocol(COLD_DRUG)
        self.assertEqual(protocol.seeds, (17, 42, 43, 2026, 2027))
        self.assertEqual(len(protocol.coordinates), 25)
        self.assertEqual(len(self.release.axis_inputs), 8)
        self.assertEqual(len(self.release.conditions), 1383)
        validation = Counter()
        test = Counter()
        sizes = []
        for fold in protocol.folds:
            positions = fold.positions(self.release.conditions)
            groups = {
                phase: {self.release.conditions[i].drug_parent for i in positions[phase]}
                for phase in ("train", "valid", "test")
            }
            self.assertEqual(sum(map(len, groups.values())), 459)
            self.assertEqual(len(set.union(*groups.values())), 459)
            self.assertEqual(sum(map(len, positions.values())), 1383)
            validation.update(groups["valid"])
            test.update(groups["test"])
            sizes.append(len(groups["test"]))
        self.assertEqual(test, validation)
        self.assertEqual(set(test.values()), {1})
        self.assertEqual(sorted(sizes), [91, 92, 92, 92, 92])

    def test_dose_crossing_covers_every_condition_and_excludes_both_strips(self):
        protocol = self.release.protocol(COLD_DOSE_DRUG)
        self.assertEqual(protocol.seeds, (17, 42, 2026))
        self.assertEqual(len(protocol.coordinates), 27)
        seen = Counter()
        held = {}
        for fold in protocol.folds:
            positions = fold.positions(self.release.conditions)
            self.assertEqual(set(positions), {
                "train", "valid", "test", "excluded-held-parent", "excluded-held-dose",
            })
            self.assertEqual(sum(map(len, positions.values())), 1383)
            self.assertEqual(len(fold.held_parents), 153)
            held.setdefault(fold.held_dose_um, set()).update(fold.held_parents)
            for phase, rows in positions.items():
                for i in rows:
                    condition = self.release.conditions[i]
                    expected = {
                        "test": (True, True),
                        "excluded-held-parent": (True, False),
                        "excluded-held-dose": (False, True),
                        "train": (False, False),
                        "valid": (False, False),
                    }[phase]
                    self.assertEqual((condition.drug_parent in fold.held_parents,
                                      condition.dose_um == fold.held_dose_um), expected)
            seen.update(positions["test"])
        self.assertEqual({dose: len(parents) for dose, parents in held.items()},
                         {0.2: 459, 1.0: 459, 10.0: 459})
        self.assertEqual(seen, dict.fromkeys(range(1383), 1))

    def test_align_caller_row_order_and_reject_population_drift(self):
        fold = self.release.protocol(COLD_DRUG).folds[0]
        conditions = self.release.conditions
        reordered = conditions[::-1]
        normal = fold.positions(conditions)
        observed = fold.positions(reordered)
        self.assertEqual(observed["test"], tuple(
            len(conditions) - index - 1 for index in reversed(normal["test"])
        ))
        altered = (replace(conditions[0], plate="unknown"), *conditions[1:])
        for invalid in (altered, conditions[:-1], (*conditions, conditions[0])):
            with self.assertRaisesRegex(ValueError, "condition population"):
                fold.positions(invalid)

    def test_reject_rewritten_ledger_with_pinned_manifest(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root)
            (path / "manifest.json").write_bytes((DEFAULT / "manifest.json").read_bytes())
            ledger = (DEFAULT / "conditions.csv").read_bytes()
            (path / "conditions.csv").write_bytes(ledger + b"\n")
            with self.assertRaisesRegex(ValueError, "condition ledger digest differs"):
                open_tenv_folds(path)
