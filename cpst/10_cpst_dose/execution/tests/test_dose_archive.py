"""Exercise Step 10 receipt checks with disposable evidence."""

from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from dose_prediction_archive import acceptance, receipt_checker
from dose_runs import axis_run
from dose_formal import FIGURES_RECEIPT, RECEIPT, SCHEMA
from dose_matrix import run_paths
from dose_protocol import SCENARIO, coordinates, protocol


class ArchiveReceiptTest(unittest.TestCase):
    def test_only_a_pinned_run_can_omit_the_weight_and_its_other_outputs_remain_checked(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            run = root / "runs/example"
            run.mkdir(parents=True)
            (run / "trace.csv").write_bytes(b"trace")
            record = {"registration_sha256": "registration", "process_token": "another process",
                      "outputs": {"checkpoint.pt": "weight digest", "trace.csv": sha256(b"trace").hexdigest()}}
            (run / "manifest.json").write_text(json.dumps(record))
            files = ("checkpoint.pt", "trace.csv")
            with self.assertRaises(ValueError):
                axis_run.verified_receipt(run, "manifest.json", "registration", files)
            checker = receipt_checker(root, {"runs/example/manifest.json": axis_run.digest(run / "manifest.json")})
            self.assertEqual(checker(run, "manifest.json", "registration", files), record)
            with self.assertRaises(ValueError):
                checker(run, "manifest.json", "registration", ("checkpoint.pt",))
            (run / "trace.csv").write_bytes(b"changed")
            with self.assertRaises(ValueError):
                checker(run, "manifest.json", "registration", files)

    def test_untrusted_acceptance_pin_fails_before_reading_matrix_or_data(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            receipt = root / RECEIPT
            receipt.parent.mkdir(parents=True)
            receipt.write_text("{}")
            with self.assertRaisesRegex(ValueError, "acceptance receipt digest"):
                acceptance(root / "missing-bundle", root, sha256(b"other receipt").hexdigest())

    def test_disposition_cannot_admit_a_path_outside_the_accepted_roster(self):
        with TemporaryDirectory() as temporary:
            bundle = Path(temporary)
            publication = bundle / "readback"
            matrix = bundle / "matrix"
            report = bundle / "report"
            for directory in (publication / "log", matrix, report):
                directory.mkdir(parents=True)
            (matrix / "preflight.json").write_text("{}")
            (report / "manifest.json").write_text("{}")
            roster = coordinates()
            first = roster[0]
            remaining = [{"coordinate": list(key), "path": run.relative_to(matrix).as_posix(),
                          "manifest_sha256": sha256("/".join(map(str, key)).encode()).hexdigest()}
                         for key, run in run_paths(matrix, False).items()]
            self.assertEqual(len(remaining), 270)
            for escaped in (False, True):
                if escaped:
                    remaining[0]["path"] = "../different-fit"
                (matrix / "manifest.json").write_text(json.dumps({"fits": remaining}))
                receipt = publication / RECEIPT
                receipt.write_text(json.dumps({
                    "schema": SCHEMA, "status": "PASS", "ticket": 405, "scenario": SCENARIO,
                    "seeds": list(protocol().seeds), "coordinates": 270, "independent_exact_replays": 270,
                    "matrix_sha256": axis_run.digest(matrix / "manifest.json"),
                    "report_sha256": axis_run.digest(report / "manifest.json"),
                    "preflight_sha256": axis_run.digest(matrix / "preflight.json"),
                }))
                pin = axis_run.digest(receipt)
                (publication / FIGURES_RECEIPT).write_text(json.dumps({"scientific_receipt_sha256": pin}))
                if escaped:
                    with self.assertRaisesRegex(ValueError, "fit path differs"):
                        acceptance(bundle, publication, pin)
                else:
                    verified = acceptance(bundle, publication, pin)
                    first_run = run_paths(matrix, False)[first]
                    name = (first_run / "manifest.json").relative_to(matrix).as_posix()
                    self.assertEqual(verified[name], remaining[0]["manifest_sha256"])
                    self.assertEqual(len(verified), 270)
