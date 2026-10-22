"""Offline contract test for the biolord adapter (ticket 14).

Deliberately does **not** import the adapter or the vendored model (no torch
model construction, no scvi/anndata import, no GPU): it checks the registry
entry, the author-default config, the adapter source for the required wiring,
and the fidelity ledger.  The functional smoke is a GPU run recorded in
`models/biolord/README.md`.

Run:
    cd stage2/03_bench_models && python -m unittest tests.test_biolord_adapter -v
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

CODE = Path(__file__).resolve().parents[1]
BENCH = CODE.parent
MODEL_DIR = BENCH / "models" / "biolord"


class BiolordAdapterContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads((MODEL_DIR / "config.json").read_text())
        cls.registry = json.loads((BENCH / "configs" / "models.json").read_text())
        cls.adapter_source = (CODE / "adapters" / "biolord" / "run.py").read_text()

    def test_registry_entry(self) -> None:
        entry = self.registry["biolord"]
        self.assertEqual(entry["adapter"], "adapters/biolord/run.py")
        self.assertEqual(
            entry["env"], {"kind": "venv", "path": "stage2/_envs/biolord"}
        )
        self.assertEqual(entry["ideal_fit"], "partial")
        self.assertFalse(entry["baseline"])
        self.assertEqual(entry["author_defaults"]["config"], "models/biolord/config.json")
        self.assertTrue((CODE / entry["adapter"]).exists())

    def test_gene_space(self) -> None:
        self.assertEqual(self.config["model"], "biolord")
        self.assertEqual(self.config["n_genes"], 18164)

    def test_author_recipe_defaults(self) -> None:
        author = self.config["author_config"]
        self.assertEqual(author["ordered_attributes_keys"], ["rdkit2d_dose"])
        self.assertEqual(author["categorical_attributes_keys"], ["cell_type"])
        self.assertEqual(author["batch_size"], 512)
        self.assertEqual(author["max_epochs"], 500)
        self.assertEqual(author["n_latent"], 256)
        self.assertEqual(author["n_latent_attribute_ordered"], 256)
        self.assertEqual(author["n_latent_attribute_categorical"], 3)
        self.assertEqual(author["decoder_width"], 4096)
        self.assertEqual(author["decoder_depth"], 4)
        self.assertEqual(author["attribute_nn_width"], 2048)
        self.assertEqual(author["attribute_nn_depth"], 2)
        self.assertEqual(author["latent_lr"], 1e-4)
        self.assertEqual(author["reconstruction_penalty"], 1e4)
        self.assertEqual(author["unknown_attribute_noise_param"], 20.0)
        self.assertEqual(author["unknown_attribute_penalty"], 0.5)
        self.assertTrue(author["cosine_scheduler"])
        # module_params is what the adapter actually passes to BiolordModule.
        params = self.config["module_params"]
        self.assertEqual(params["reconstruction_penalty"], 1e4)
        self.assertEqual(params["unknown_attribute_noise_param"], 20.0)
        self.assertEqual(params["decoder_width"], 4096)
        self.assertEqual(params["gene_likelihood"], "normal")
        self.assertEqual(params["n_latent_attribute_ordered"], 256)
        plan = self.config["plan_kwargs"]
        self.assertEqual(plan["latent_lr"], 1e-4)
        self.assertEqual(plan["attribute_nn_lr"], 1e-2)
        self.assertEqual(plan["scheduler_max_epochs"], 500)

    def test_attribute_encoding(self) -> None:
        self.assertEqual(self.config["descriptor"]["kind"], "rdkit2d")
        self.assertEqual(self.config["descriptor"]["n_raw"], 200)
        self.assertEqual(self.config["dose"]["transform"], "raw_uM")
        self.assertEqual(self.config["dose"]["position"], "appended_last")
        self.assertEqual(self.config["categorical"]["cell_type"], "ipsc_cm")

    def test_checkpoint_grid_and_smoke(self) -> None:
        self.assertEqual(self.config["checkpoint_grid"], 50)
        self.assertEqual(self.config["max_epochs"], 500)
        self.assertEqual(self.config["checkpoint_grid"] * 10, self.config["max_epochs"])
        smoke = self.config["smoke"]
        self.assertLess(smoke["max_epochs"], self.config["max_epochs"])
        self.assertLessEqual(smoke["checkpoint_grid"], smoke["max_epochs"])
        self.assertGreater(smoke["limit_train_batches"], 0)

    def test_adapter_wiring(self) -> None:
        for needle in (
            "BENCH_SMOKE",
            "dont_write_bytecode",
            "rdkit2d_dose",
            "get_dataset",
            "get_expression",
            "control_rows",
            "limit_val_batches",
            "best.pt",
            "train_adata.h5ad",
            "ensure_float32",
        ):
            self.assertIn(needle, self.adapter_source)
        # The categorical-only helper cannot target the ordered attribute.
        self.assertNotIn("compute_prediction_adata", self.adapter_source)

    def test_fidelity_ledger(self) -> None:
        text = (MODEL_DIR / "FIDELITY.md").read_text()
        for needle in (
            "18,164",
            "random init",
            "7,561",
            "4,999",
            "corpus scale",
            "transductive prediction",
            "ideal_fit: partial",
            "control handling",
            "attribute encoding",
            "checkpoint grid",
            "recipe_substitution",
            "recipe_omission",
            "adapter_extension",
            "config_only",
        ):
            self.assertIn(needle, text)
        self.assertTrue((MODEL_DIR / "README.md").exists())


if __name__ == "__main__":
    unittest.main()
