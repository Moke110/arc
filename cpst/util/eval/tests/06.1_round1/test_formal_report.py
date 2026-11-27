from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

import torch

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst/06_tsg_module/06.1_round1/execution"))
from report_shape import verify_checkpoint


class FormalReferenceTest(unittest.TestCase):
    def test_exact_selected_checkpoint_requires_the_same_state_and_objective(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.pt"
            fresh = root / "fresh.pt"
            checkpoint = {"epoch": 15, "validation_loss": 2.5, "validation_mse": 1.5,
                          "state_dict": {"weight": torch.tensor([1.0, 2.0])}}
            torch.save({**checkpoint, "registration_sha256": "older"}, reference)
            torch.save({**checkpoint, "registration_sha256": "fresh"}, fresh)
            verify_checkpoint(fresh, reference, "fresh", "older")
            torch.save({**checkpoint, "registration_sha256": "fresh",
                        "state_dict": {"weight": torch.tensor([1.0, 2.1])}}, fresh)
            with self.assertRaisesRegex(ValueError, "checkpoint tensor differs"):
                verify_checkpoint(fresh, reference, "fresh", "older")
