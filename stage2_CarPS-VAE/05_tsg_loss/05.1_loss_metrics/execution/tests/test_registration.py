"""The loss coordinate is part of the public registration command."""

from pathlib import Path
import subprocess
import sys
from unittest import TestCase

RUN = Path(__file__).resolve().parents[1] / "run_loss.py"


class RegistrationCommandTests(TestCase):
    def test_registration_advertises_the_registered_coordinate(self):
        result = subprocess.run([sys.executable, str(RUN), "register", "--help"],
                                capture_output=True, text=True, check=True)
        self.assertIn("--loss {mse,cos,pcc,l1}", result.stdout)
        self.assertIn("--alpha {0.75,1,1.25}", result.stdout)
        self.assertIn("--fold {0,1,2}", result.stdout)
