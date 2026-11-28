"""Publication evidence checks for the canonical architecture figures."""

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
STAGE = REPO / "cpst"
sys.path.insert(0, str(STAGE / "06_tsg_module/execution"))

from build_comparison import accepted_fit
from cpst_util.derivation.release import sha256, write_json


class ComparisonEvidenceTests(unittest.TestCase):
    def test_rejects_a_changed_replay_manifest_for_an_accepted_fit(self):
        with TemporaryDirectory() as directory:
            fit = Path(directory)
            registration = fit / "registration.json"
            manifest = fit / "manifest.json"
            write_json(registration, {"seed": 17})
            write_json(manifest, {"test_replay": "exact"})
            record = SimpleNamespace(shape="XΔ-XΔ", fold=0, seed=17,
                                     run_sha256=sha256(registration), manifest_sha256=sha256(manifest))
            paths = {(record.run_sha256, record.manifest_sha256): fit}

            write_json(manifest, {"test_replay": "different"})
            with self.assertRaisesRegex(ValueError, "XΔ-XΔ.*0.*17"):
                accepted_fit(paths, record)
