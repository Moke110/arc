import sys
from pathlib import Path
import unittest

REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / ".git").exists())
sys.path.insert(0, str(REPO / "cpst"))
from util.training import experiment_configuration


class ExperimentConfigurationTests(unittest.TestCase):
    def test_explicit_execution_identity_keeps_ordered_axis_and_objective(self):
        transform = {"mean": 0.48, "std": 0.23}
        model = {"method_identity": "native", "molecule_feature": "KPGT",
                 "transform": transform, "beta": 0.1}
        objective = {"pcc": 0.25, "delta_l1": 1.75}
        bound = {"model": model, "objective": objective, "shape": "XΔ-XΔ",
                 "axis_sha256": "ordered-genes", "scenario": "cold-drug", "seed": 17}
        arguments = {"axis": "denoise-5k", "fold": 1}

        self.assertEqual(experiment_configuration(arguments, bound), {
            "architecture": "XΔ-XΔ",
            "objective": objective,
            "representation": "KPGT",
            "dose": {"conditioning": "latent-numeric", "transform": transform},
            "axis": {"name": "denoise-5k", "sha256": "ordered-genes"},
            "split": "cold-drug",
            "fold": 1,
            "seed": 17,
        })


if __name__ == "__main__":
    unittest.main()
