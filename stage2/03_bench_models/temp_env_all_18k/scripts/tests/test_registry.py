"""Offline registry contract for the real-model adapters (ticket 13).

No model import and no GPU: this only checks that every non-baseline registry
entry points at an existing adapter, config, run card, and fidelity ledger, and
that the config declares the 18,164-gene output.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
CODE = HERE.parents[1]
BENCH = CODE.parent
sys.path.insert(0, str(CODE))

import runner  # noqa: E402


class RealAdapterRegistryTest(unittest.TestCase):
    def test_every_real_adapter_is_wired(self) -> None:
        registry = runner.load_registry()
        real = {name: entry for name, entry in registry.items() if not entry.get("baseline")}
        self.assertTrue({"prnet", "cpa", "chemcpa"} <= set(real))

        for name, entry in real.items():
            with self.subTest(model=name):
                self.assertTrue((CODE / entry["adapter"]).exists(), entry["adapter"])
                self.assertIn(entry["env"]["kind"], {"python", "venv", "conda"})

                config_path = BENCH / entry["author_defaults"]["config"]
                self.assertTrue(config_path.exists(), config_path)
                config = json.loads(config_path.read_text())
                self.assertEqual(config["model"], name)

                gene_setting = json.dumps(config)
                self.assertIn("18164", gene_setting)

                model_dir = BENCH / "models" / name
                for artifact in ("README.md", "FIDELITY.md", "config.json"):
                    self.assertTrue((model_dir / artifact).exists(), model_dir / artifact)

                ledger = (model_dir / "FIDELITY.md").read_text()
                self.assertIn("18,164", ledger)
                self.assertIn("random init", ledger)


if __name__ == "__main__":
    unittest.main()
