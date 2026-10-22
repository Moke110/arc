"""Offline contract test for the Squidiff adapter (ticket 15).

Deliberately does **not** import the adapter or the vendored model (no torch
model construction, no GPU): it checks the registry entry, the author-default
config, and the fidelity ledger.  The functional smoke is a GPU run recorded in
`models/squidiff/README.md`.

Run:
    cd stage2/03_bench_models && python -m unittest tests.test_squidiff_adapter -v
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

CODE = Path(__file__).resolve().parents[1]
BENCH = CODE.parent
MODEL_DIR = BENCH / "models" / "squidiff"


class SquidiffAdapterContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads((MODEL_DIR / "config.json").read_text())
        cls.registry = json.loads((BENCH / "configs" / "models.json").read_text())
        cls.adapter_source = (CODE / "adapters" / "squidiff" / "run.py").read_text()

    def test_registry_entry(self) -> None:
        entry = self.registry["squidiff"]
        self.assertEqual(entry["adapter"], "adapters/squidiff/run.py")
        self.assertEqual(entry["env"], {"kind": "conda", "name": "arc-stage2-bench-squidiff"})
        self.assertEqual(entry["ideal_fit"], "full")
        self.assertFalse(entry["baseline"])
        self.assertEqual(entry["author_defaults"]["config"], "models/squidiff/config.json")
        self.assertTrue((CODE / entry["adapter"]).exists())

    def test_gene_space(self) -> None:
        self.assertEqual(self.config["model"], "squidiff")
        self.assertEqual(self.config["gene_size"], 18164)
        self.assertEqual(self.config["output_dim"], 18164)
        # The sampler uses gene_size, so the two must agree.
        self.assertEqual(self.config["gene_size"], self.config["output_dim"])
        self.assertTrue(self.config["use_drug_structure"])
        self.assertEqual(self.config["drug_dimension"], 1024)

    def test_author_recipe_defaults(self) -> None:
        self.assertEqual(self.config["batch_size"], 64)
        self.assertEqual(self.config["total_steps"], 100000)
        self.assertEqual(self.config["schedule"]["total_steps"], 100000)
        self.assertEqual(self.config["schedule"]["name"], "linear_anneal_to_zero")
        self.assertEqual(self.config["optimizer"], {"name": "adamw", "lr": 1e-4, "weight_decay": 0.0})
        self.assertEqual(self.config["ema_rate"], 0.9999)
        self.assertEqual(self.config["diffusion_steps"], 1000)
        self.assertEqual(self.config["num_layers"], 3)
        self.assertEqual(self.config["checkpoint_grid"], 10)
        self.assertEqual(self.config["fingerprint"]["dose_transform"], "log10(dose_uM + 1)")
        self.assertEqual(self.config["fingerprint"]["radius"], 2)

    def test_smoke_budget(self) -> None:
        smoke = self.config["smoke"]
        self.assertLess(smoke["total_steps"], self.config["total_steps"])
        self.assertLess(smoke["sample_steps"], self.config["diffusion_steps"])
        self.assertGreater(smoke["total_steps"], 0)

    def test_adapter_flags(self) -> None:
        for needle in ("BENCH_SMOKE", "dont_write_bytecode", "ddim", "update_ema", "100,000"):
            self.assertIn(needle, self.adapter_source)

    def test_fidelity_ledger(self) -> None:
        text = (MODEL_DIR / "FIDELITY.md").read_text()
        for needle in (
            "18,164",
            "random init",
            "90x",
            "paired control",
            "checkpoint grid",
            "EMA",
            "drug condition",
            "dose transform",
            "corpus scale",
            "recipe_substitution",
            "adapter_extension",
        ):
            self.assertIn(needle, text)
        self.assertTrue((MODEL_DIR / "README.md").exists())


if __name__ == "__main__":
    unittest.main()
