"""Exercise the frozen training plan and formal autoencoder commands."""
import hashlib
import json
import subprocess
import sys
import unittest

import test_inputs


class FormalFitTest(unittest.TestCase):
    def setUp(self):
        fixture = test_inputs.InputReleaseTest()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        fixture.prepare("--seed", "1")
        self.release = fixture.release
        self.run_dir = fixture.root / "training"

    def cli(self, action, *args, error=None):
        command = [
            sys.executable, str(test_inputs.CLI), action, "--release", str(self.release),
            "--fixture", "--device", "cpu", *args,
        ]
        if action in {"freeze", "develop", "refit"}:
            command.extend(("--run-dir", str(self.run_dir)))
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if error:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(error, result.stderr)
            return None
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout)

    def freeze(self):
        return self.cli("freeze", "--batch-size", "2", "--seed", "17",
                        "--max-epochs", "2", "--max-seconds", "60")

    def test_freeze_binds_recipe_and_split_before_training(self):
        frozen = self.freeze()
        plan = self.run_dir / "plan.json"
        self.assertEqual(frozen["plan_sha256"], hashlib.sha256(plan.read_bytes()).hexdigest())
        record = json.loads(plan.read_text())
        self.assertEqual(record["recipe"]["optimizer"]["lr"], 0.001)
        self.assertEqual(record["recipe"]["seed"], 17)
        self.assertEqual(record["recipe"]["max_epochs"], 2)
        self.assertEqual(record["recipe"]["max_seconds"], 60)
        manifest = json.loads((self.release / "manifest.json").read_text())
        self.assertEqual(record["source_membership"], manifest["sources"])
        self.assertEqual(record["split_seed"], 1)
        self.cli("freeze", "--max-epochs", "2", "--max-seconds", "60", error="already exists")
        with plan.open("a") as file:
            file.write(" ")
        self.cli("develop", error="record digest mismatch: plan.json")
        self.assertFalse((self.run_dir / "development").exists())

    def test_development_selects_from_source_specific_validation(self):
        frozen = self.freeze()
        development = self.cli("develop")
        self.assertEqual(development["plan_sha256"], frozen["plan_sha256"])
        epochs = development["epochs"]
        self.assertEqual([epoch["epoch"] for epoch in epochs], [1, 2])
        for epoch in epochs:
            self.assertEqual(epoch["training_profiles"], {"tahoe": 4, "tenv": 3})
            self.assertEqual(epoch["validation_profiles"], {"tahoe": 1, "tenv": 1})
            self.assertEqual(epoch["steps"], 4)
            self.assertEqual(set(epoch["validation_mse"]), {"tahoe", "tenv"})
        best = min(epochs, key=lambda epoch: (epoch["validation_mse"]["tenv"], epoch["epoch"]))
        self.assertEqual(development["selected_epochs"], best["epoch"])
        self.assertEqual(development["stage"], "development")

    def test_refit_restarts_and_independent_reload_checks_provenance(self):
        self.freeze()
        development = self.cli("develop")
        output = self.release.parent / "autoencoder"
        refit = self.cli("refit", "--output", str(output))
        self.assertEqual(refit["fit_stage"], "full_data_refit")
        self.assertEqual(refit["initial_model_sha256"], development["initial_model_sha256"])
        self.assertEqual(refit["selected_epochs"], development["selected_epochs"])
        self.assertEqual(len(refit["epochs"]), development["selected_epochs"])
        for epoch in refit["epochs"]:
            self.assertEqual(epoch["training_profiles"], {"tahoe": 5, "tenv": 4})
            self.assertEqual(epoch["steps"], 5)
            self.assertNotIn("validation_mse", epoch)
        checkpoint_args = ("--checkpoint", str(output / "checkpoint.pt"))
        replay = self.cli("reload", *checkpoint_args)
        self.assertEqual(replay["checkpoint_sha256"], refit["checkpoint_sha256"])
        self.assertEqual(replay["reconstruction_shape"], [2, 3])
        self.assertEqual(replay["plan_sha256"], refit["plan_sha256"])
        self.assertLessEqual(replay["max_abs_error"], replay["atol"])
        with (output / "development.json").open("a") as file:
            file.write(" ")
        self.cli("reload", *checkpoint_args, error="record digest mismatch: development.json")
