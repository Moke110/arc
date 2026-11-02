"""Config resolution/validation and the structure/recipe extension points (#10/#11)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import _support  # noqa: F401  (path setup)
from ptlib import config


def canonical() -> dict:
    return {
        "structure": {
            "n_genes": 18164,
            "hidden": [128],
            "z_dim": 64,
            "c_dim": 64,
            "n_dim": 10,
            "drug_dim": 1024,
            "adaptor": [128],
            "dropout": 0.05,
        },
        "recipe": {
            "batch_size": 512,
            "max_epochs": 100,
            "min_epochs": 20,
            "early_stop_patience": 20,
            "sched_factor": 0.5,
            "sched_patience": 10,
            "min_lr": 1e-8,
            "weight_decay": 1e-8,
        },
        "metric": {"name": "direction_auprc", "k": 200, "plateau_mode": "max"},
        "loss": "gaussian_nll",
        "optimizer": "adam",
        "prediction_noise": "zero",
        "n_folds": 5,
        "seeds": [20261021, 20261022, 20261023],
        "lrs": {"pretrain": 1e-3, "cm": [1e-4, 3e-4, 1e-3]},
    }


class NormalizeTests(unittest.TestCase):
    def test_canonical_passthrough(self):
        cfg, notes = config.normalize_locked_config(canonical())
        self.assertEqual(cfg["structure"]["n_genes"], 18164)
        self.assertEqual(cfg["metric"]["k"], 200)
        self.assertEqual(notes["mapping"], [])
        config.validate_config(cfg)

    def test_field_aliases_are_recorded(self):
        raw = {
            "structure": {
                "in_dim": 18164,
                "hidden_layer_sizes": 128,
                "state_latent": 64,
                "drug_latent": 64,
                "noise": 10,
                "n_bits": 1024,
                "adaptor_layer_sizes": [128],
                "dr_rate": 0.05,
            },
            "optim": {
                "batch": 512, "epochs": 100, "min_epochs": 20, "patience": 20,
                "plateau_factor": 0.5, "plateau_patience": 10, "wd": 1e-8,
            },
            "metric": {"metric": "direction AUPRC@200"},
            "learning_rates": {"pretrain": 1e-3, "cm": [1e-4, 3e-4, 1e-3]},
            "seeds": [1, 2, 3],
        }
        cfg, notes = config.normalize_locked_config(raw)
        self.assertEqual(cfg["structure"]["hidden"], [128])
        self.assertEqual(cfg["structure"]["z_dim"], 64)
        self.assertEqual(cfg["structure"]["drug_dim"], 1024)
        self.assertEqual(cfg["recipe"]["max_epochs"], 100)
        self.assertEqual(cfg["metric"]["k"], 200)
        self.assertTrue(any("in_dim" in item for item in notes["mapping"]))
        self.assertTrue(any("state_latent" in item for item in notes["mapping"]))

    def test_missing_recipe_values_default_and_are_recorded(self):
        raw = canonical()
        del raw["recipe"]["batch_size"]
        raw["structure"].pop("drug_dim")
        cfg, notes = config.normalize_locked_config(raw)
        self.assertEqual(cfg["recipe"]["batch_size"], 512)
        self.assertEqual(cfg["structure"]["drug_dim"], 1024)
        self.assertTrue(any("batch_size" in item for item in notes["defaulted"]))


class ValidationTests(unittest.TestCase):
    def _cfg(self):
        return canonical()

    def test_metric_k_above_axis_rejected(self):
        cfg = self._cfg()
        cfg["structure"]["n_genes"] = 100
        with self.assertRaises(config.ConfigError):
            config.validate_config(cfg)

    def test_plateau_mode_must_be_max(self):
        cfg = self._cfg()
        cfg["metric"]["plateau_mode"] = "min"
        with self.assertRaises(config.ConfigError):
            config.validate_config(cfg)

    def test_metric_must_be_direction_auprc(self):
        cfg = self._cfg()
        cfg["metric"]["name"] = "valid_loss"
        with self.assertRaises(config.ConfigError):
            config.validate_config(cfg)

    def test_loss_must_be_gaussian_nll(self):
        cfg = self._cfg()
        cfg["loss"] = "mse"
        with self.assertRaises(config.ConfigError):
            config.validate_config(cfg)

    def test_batch_size_one_rejected_for_batchnorm(self):
        cfg = self._cfg()
        cfg["recipe"]["batch_size"] = 1
        with self.assertRaises(config.ConfigError):
            config.validate_config(cfg)

    def test_min_epochs_above_max_rejected(self):
        cfg = self._cfg()
        cfg["recipe"]["min_epochs"] = 200
        with self.assertRaises(config.ConfigError):
            config.validate_config(cfg)

    def test_nonzero_prediction_noise_rejected(self):
        cfg = self._cfg()
        cfg["prediction_noise"] = "train"
        with self.assertRaises(config.ConfigError):
            config.validate_config(cfg)


class StructureTests(unittest.TestCase):
    def test_presets(self):
        base = canonical()["structure"]
        struct, notes = config.structure_from_spec("hidden128", base)
        self.assertEqual(struct["hidden"], [128])
        self.assertEqual(notes["source"], "preset")
        struct, _ = config.structure_from_spec("latent128", base)
        self.assertEqual(struct["z_dim"], 128)
        self.assertEqual(struct["hidden"], [128])
        struct, _ = config.structure_from_spec("hidden256", base)
        self.assertEqual(struct["hidden"], [256])

    def test_json_override(self):
        base = canonical()["structure"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "structure.json"
            path.write_text(json.dumps({"state_latent": 128, "hidden": [256]}))
            struct, notes = config.structure_from_spec(str(path), base)
        self.assertEqual(struct["z_dim"], 128)
        self.assertEqual(struct["hidden"], [256])
        self.assertEqual(notes["source"], "json")

    def test_unknown_spec_fails(self):
        with self.assertRaises(config.ConfigError):
            config.structure_from_spec("gigantic", canonical()["structure"])

    def test_compatibility(self):
        a = canonical()["structure"]
        ok, mismatches = config.structure_compatible(a, dict(a))
        self.assertTrue(ok)
        self.assertEqual(mismatches, [])
        b = dict(a, z_dim=128)
        ok, mismatches = config.structure_compatible(a, b)
        self.assertFalse(ok)
        self.assertTrue(any("z_dim" in item for item in mismatches))

    def test_capacity_parameter_expectations(self):
        base = canonical()["structure"]
        struct128, _ = config.structure_from_spec("latent128", base)
        struct256, _ = config.structure_from_spec("hidden256", base)
        from ptlib import model as model_mod

        # the audited capacity arithmetic (A §1.5 / design review)
        self.assertEqual(model_mod.parameter_count(model_mod.build_model(base)), 7_185_000)
        self.assertEqual(model_mod.parameter_count(model_mod.build_model(struct128)), 7_201_448)
        self.assertEqual(model_mod.parameter_count(model_mod.build_model(struct256)), 14_194_280)


class RecipeOverrideTests(unittest.TestCase):
    def test_overrides_recorded_and_validated(self):
        recipe, applied = config.apply_recipe_overrides(
            canonical()["recipe"], {"max_epochs": 2, "min_epochs": 1}
        )
        self.assertEqual(recipe["max_epochs"], 2)
        self.assertEqual(recipe["min_epochs"], 1)
        self.assertEqual(applied, {"max_epochs": 2, "min_epochs": 1})
        with self.assertRaises(config.ConfigError):
            config.apply_recipe_overrides(recipe, {"max_epochs": 0})

    def test_fingerprint_changes_with_structure_and_seed(self):
        import copy

        cfg = canonical()
        h1 = config.config_fingerprint(cfg)
        cfg2 = copy.deepcopy(cfg)
        cfg2["structure"]["z_dim"] = 128
        self.assertNotEqual(h1, config.config_fingerprint(cfg2))
        cfg3 = copy.deepcopy(cfg)
        cfg3["seed"] = 1
        self.assertNotEqual(h1, config.config_fingerprint(cfg3))


if __name__ == "__main__":
    unittest.main()
