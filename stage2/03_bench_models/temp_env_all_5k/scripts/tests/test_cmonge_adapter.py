"""Offline contract tests for the CMonge adapter (ticket 15).

Two layers, neither trains a model or touches CUDA:

- `ConfigContractTest` reads `models/cmonge/config.json` and the adapter source
  with the *harness* interpreter (the `arc` env has no JAX).
- `DescriptorHelperTest` runs the adapter's pure RDKit/condition helpers in the
  bench venv through a subprocess, with `JAX_PLATFORMS=cpu`; it is skipped when
  `envs/cmonge` has not been built.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

CODE = Path(__file__).resolve().parents[1]
BENCH = CODE.parent
CONFIG = BENCH / "models/cmonge/config.json"
ADAPTER = CODE / "adapters/cmonge/run.py"
REGISTRY = BENCH / "configs/models.json"
VENV_PY = BENCH.parents[1] / "_envs/cmonge/bin/python"


class ConfigContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.cfg = json.loads(CONFIG.read_text())

    def test_author_recipe_is_frozen(self) -> None:
        self.assertEqual(self.cfg["model"], "cmonge")
        self.assertEqual(self.cfg["n_genes"], 5000)
        self.assertEqual(self.cfg["ae"]["data_dim"], 5000)
        self.assertEqual(self.cfg["ae"]["latent_dim"], 50)
        self.assertEqual(self.cfg["ae"]["hidden_dims"], [512, 512])
        self.assertEqual(self.cfg["ae"]["n_epochs"], 50)
        self.assertEqual(self.cfg["ae"]["batch_size"], 512)
        self.assertEqual(self.cfg["ae"]["optimizer"]["lr"], 1e-4)
        self.assertEqual(self.cfg["ae"]["optimizer"]["weight_decay"], 1e-5)
        self.assertEqual(self.cfg["monge"]["dim_data"], 50)
        self.assertEqual(self.cfg["monge"]["hidden_dims"], [64, 64, 64, 64])
        self.assertEqual(self.cfg["monge"]["num_train_iters"], 10000)
        self.assertEqual(self.cfg["monge"]["batch_size"], 512)
        self.assertEqual(self.cfg["monge"]["fitting_loss"]["epsilon_fitting"], 1.0)
        self.assertEqual(self.cfg["monge"]["regularizer"]["epsilon_regularizer"], 1e-2)
        self.assertEqual(self.cfg["seed"], 20261021)

    def test_smoke_is_shorter_than_the_author_run(self) -> None:
        smoke = self.cfg["smoke"]
        self.assertLess(smoke["ae_epochs"], self.cfg["ae"]["n_epochs"])
        self.assertLess(smoke["cm_iters"], self.cfg["monge"]["num_train_iters"])
        self.assertLessEqual(smoke["eval_every"], smoke["cm_iters"])

    def test_registry_entry_points_at_this_adapter(self) -> None:
        registry = json.loads(REGISTRY.read_text())
        entry = registry["cmonge"]
        self.assertEqual(entry["adapter"], "adapters/cmonge/run.py")
        self.assertEqual(entry["env"], {"kind": "venv", "path": "stage2/_envs/cmonge"})
        self.assertEqual(entry["ideal_fit"], "full")
        self.assertEqual(entry["author_defaults"]["config"], "models/cmonge/config.json")

    def test_adapter_has_no_truth_access(self) -> None:
        source = ADAPTER.read_text()
        for forbidden in (
            "group_delta",
            "delta_log1p_cp10k.npy",
            "true_sets",
            "splits.parquet",
            "counts.npy",
        ):
            self.assertNotIn(forbidden, source)
        self.assertIn("sys.dont_write_bytecode = True", source)


@unittest.skipUnless(VENV_PY.exists(), "cmonge bench venv not built")
class DescriptorHelperTest(unittest.TestCase):
    def test_pure_helpers(self) -> None:
        script = r"""
import importlib.util
import sys
import numpy as np

spec = importlib.util.spec_from_file_location("cmonge_adapter", r"%s")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

smiles = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
rows, width = module._rdkit_rows(smiles)
assert rows.shape == (3, width), rows.shape
assert width >= 100, width
assert np.isfinite(rows).all()
stats = module._fit_descriptor_standardizer(rows, 0.01)
keep = np.asarray(stats["keep"], dtype=bool)
assert keep.shape == (width,)
normalized = module._apply_descriptor_standardizer(rows, stats)
assert normalized.shape == (3, int(keep.sum()))
assert np.isfinite(normalized).all()
order = ["d0", "d1"]
descriptor_rows = rows[:2]
conditions = module._condition_vectors(
    descriptor_rows, order, ["d0", "d1"], np.asarray([1.0, 10.0]), stats
)
assert conditions.shape == (2, int(keep.sum()) + 1), conditions.shape
assert abs(float(conditions[0, -1]) - 6.907755) < 1e-4, conditions[:, -1]
assert abs(float(conditions[1, -1]) - 9.210340) < 1e-4, conditions[:, -1]
print("ok", width, int(keep.sum()))
""" % str(ADAPTER)
        env = os.environ.copy()
        env.update(
            {
                "JAX_PLATFORMS": "cpu",
                "PYTHONDONTWRITEBYTECODE": "1",
                "CUDA_VISIBLE_DEVICES": "",
                "PYTHONPATH": str(ADAPTER.parents[2]),
            }
        )
        result = subprocess.run(
            [str(VENV_PY), "-c", script],
            capture_output=True,
            text=True,
            timeout=600,
            env=env,
            cwd=str(ADAPTER.parents[2]),
        )
        self.assertEqual(result.returncode, 0, result.stderr[-2000:])
        self.assertIn("ok", result.stdout)


if __name__ == "__main__":
    unittest.main()
