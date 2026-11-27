"""Acceptance of publication's downstream report reader."""

import csv
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/06_tsg_shape/06.4_round4/execution"))

from publish_measured import verify_ledger
from measured_roster import coordinates, fit_coordinates


class FormalPublicationTest(unittest.TestCase):
    def test_score_ledger_rejects_a_mismatched_fit_manifest(self):
        roster = coordinates(False)
        new = set(fit_coordinates(False))
        def identity(index, key):
            return dict(shape=key[0], fold=key[1], seed=key[2],
                        registration_sha256=f"registration-{index}", run_manifest_sha256=f"manifest-{index}",
                        best_epoch=1, stopping_epoch=2)

        runs = [identity(index, key) for index, key in enumerate(roster)]
        fresh = {key: run for key, run in zip(roster, runs, strict=True) if key in new}
        retained = [run for key, run in zip(roster, runs, strict=True) if key not in new]
        with TemporaryDirectory() as directory:
            table = Path(directory) / "runs.csv"
            score_rows = [dict(shape=row["shape"], fold=row["fold"], seed=row["seed"],
                               run_sha256=row["registration_sha256"], manifest_sha256=row["run_manifest_sha256"],
                               best_epoch=row["best_epoch"], stopping_epoch=row["stopping_epoch"]) for row in runs]

            def write_rows():
                with table.open("w", newline="", encoding="utf-8") as stream:
                    writer = csv.DictWriter(stream, fieldnames=score_rows[0])
                    writer.writeheader()
                    writer.writerows(score_rows)

            write_rows()
            verify_ledger(table, fresh, retained)
            score_rows[18]["manifest_sha256"] = "wrong-fit-manifest"
            write_rows()
            with self.assertRaisesRegex(ValueError, "scored fit identity"):
                verify_ledger(table, fresh, retained)
