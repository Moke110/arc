from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from report_publication import payload_inventory


class ReportPublication(unittest.TestCase):
    def test_report_inventory_rejects_symlinked_payload(self):
        with TemporaryDirectory() as tmp:
            target = Path(tmp)
            (target / "real.csv").write_text("content")
            (target / "alias.csv").symlink_to("real.csv")
            with self.assertRaisesRegex(ValueError, "symlink"):
                payload_inventory(target)

    def test_report_inventory_is_relative_and_omits_receipt(self):
        with TemporaryDirectory() as tmp:
            target = Path(tmp)
            (target / "log").mkdir()
            (target / "log/acceptance.json").write_text("{}")
            (target / "README.md").write_text("Report")
            observed = payload_inventory(target)
            self.assertEqual(set(observed), {"README.md"})
            self.assertEqual(len(observed["README.md"]), 64)

    def test_publication_refuses_to_replace_existing_report(self):
        from report_publication import install_report

        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            step = root / "step"
            (step / "figure").mkdir(parents=True)
            figure = step / "figure/metrics.png"
            figure.write_bytes(b"previous accepted figure")
            with self.assertRaises(FileExistsError):
                install_report(root / "candidate", step, {}, "Reviewed all three figures.")
            self.assertEqual(figure.read_bytes(), b"previous accepted figure")
