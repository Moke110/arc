"""Exercise the shared expression comparison publication workflow."""

from hashlib import file_digest
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "script" / "compare.py"
MODELS = ("mean", "PCA_1024", "MLP_1024", "AAEc_1024", "AAEcs_1024", "AAEcs_1024+64")


class SharedCompareTest(unittest.TestCase):
    def test_fixture_fit_seal_score_and_replay(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            release = root / "release"
            self.make_release(release)
            split = root / "split.json"
            models = root / "models"
            validation = root / "validation"
            interface = models / "shared_expression"
            seal = root / "seal.json"
            scores = root / "scores.json"
            budget = root / "budget.json"
            budget.write_text(json.dumps({
                "train_batch_size": 1, "validation_batch_size": 1, "max_optimizer_steps": 2,
            }))
            self.cli("reconstruct.py", "freeze-split", release, split)
            for model in ("mean", "shared_expression", *MODELS[1:]):
                output = models / model
                arguments = ["--interface", interface] if model in MODELS[1:] else []
                fit_arguments = ["--budget", budget] if model in MODELS[2:] else []
                self.cli("reconstruct.py", "fit", release, split, "--model", model,
                         "--output", output, *arguments, *fit_arguments)
                self.cli("reconstruct.py", "score", release, split, "--artifact", output,
                         "--mean-artifact", models / "mean", *arguments,
                         "--output", validation / f"{model}-validation.json")
            arguments = ["--models", models, "--seal", seal]
            validation_path = validation / "PCA_1024-validation.json"
            valid_report = validation_path.read_text()
            mismatched = json.loads(valid_report)
            mismatched["interface_manifest_sha256"] = "0" * 64
            validation_path.write_text(json.dumps(mismatched))
            self.cli("compare.py", "seal", release, split, *arguments, "--shared-expression",
                     "--validation-dir", validation,
                     error="validation uses a different shared expression interface")
            self.assertFalse(seal.exists())
            validation_path.write_text(valid_report)
            self.cli("compare.py", "seal", release, split, *arguments, "--shared-expression",
                     "--validation-dir", validation)
            test_scores = self.cli("compare.py", "score", release, split, *arguments,
                                   "--allow-test", "--output", scores)
            replay = self.cli("compare.py", "reload", release, split, *arguments,
                              "--allow-test", "--expected", scores)
            self.assertEqual(test_scores, replay)
            self.assertEqual(set(test_scores["models"]), set(MODELS))
            for model in MODELS:
                self.assertEqual(test_scores["models"][model]["role"], "test")
            self.assertEqual(test_scores["shared_expression"]["validation"]["role"], "validation")
            self.assertEqual(test_scores["models"]["AAEcs_1024+64"]["latent"]["encoder_tokens"], 5)
            figure = root / "models_metric.png"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "figure", "--expected", str(scores), "--output", str(figure)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(figure.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

    def cli(self, script, action, release, split, *arguments, error=None):
        result = subprocess.run(
            [sys.executable, str(SCRIPT.with_name(script)), action, "--fixture", "--release",
             str(release), "--split", str(split), *map(str, arguments)],
            capture_output=True, text=True,
        )
        if error:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(error, result.stderr)
            return
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout)

    @staticmethod
    def make_release(release):
        release.mkdir()
        (release / "genes.csv").write_text("symbol,hgnc_id\nX,HGNC:1\nY,HGNC:2\nZ,HGNC:3\n")
        for source, offset in (("tahoe", 0), ("tenv", 1)):
            expression = np.array([[1, 2, 3], [2, 4, 8], [3, 1, 2], [4, 6, 8]], dtype=np.float32)
            np.save(release / f"{source}_expression.npy", expression + offset)
            profiles = ["source_row,profile_id,plate,is_control"]
            profiles.extend(f"{index},{source}-{index},{plate},False"
                            for index, plate in enumerate("ABCA"))
            (release / f"{source}_profiles.csv").write_text("\n".join(profiles) + "\n")
            (release / f"{source}_mapping.csv").write_text("fixture mapping\n")
        files = {}
        for path in release.iterdir():
            with path.open("rb") as stream:
                files[path.name] = file_digest(stream, "sha256").hexdigest()
        (release / "manifest.json").write_text(json.dumps({
            "schema_version": 1, "fixture": True,
            "normalization": "log1p(CP10K of raw counts projected to the ordered CM axis)",
            "inputs": {"genes.csv": files["genes.csv"]},
            "sources": {source: {"shape": [4, 3]} for source in ("tahoe", "tenv")},
            "files": files,
        }))

    def test_figure_accepts_exact_shared_expression_roster(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scores = root / "scores.json"
            figure = root / "models_metric.png"
            metrics = {"profiles": 1, "pcc": 0.5, "rmse": 0.4, "relative_squared_error": 0.75}
            entries = {}
            for name in MODELS:
                if name == "mean":
                    latent = {"tokens": None, "scalars": None}
                else:
                    encoder_tokens = 1088 if name.endswith("+64") else 1024
                    latent = {"tokens": 1024, "scalars": 262144, "channels": 256,
                              "encoder_tokens": encoder_tokens}
                entries[name] = {
                    "display_name": name, "latent": latent, "sources": {"tahoe": metrics, "tenv": metrics},
                }
            scores.write_text(json.dumps({
                "role": "test", "schema_version": 4, "gene_count": 18164, "models": entries,
            }))
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "figure", "--expected", str(scores), "--output", str(figure)],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(figure.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

    def test_real_score_requires_committed_seal_before_opening_corpus(self):
        with tempfile.TemporaryDirectory(dir=SCRIPT.parents[4]) as directory:
            root = Path(directory)
            seal = root / "uncommitted.json"
            seal.write_text(json.dumps({"models": dict.fromkeys(MODELS), "schema_version": 4}))
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "score", "--allow-test", "--seal", str(seal),
                 "--release", str(root / "absent"), "--output", str(root / "scores.json")],
                capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("seal must match a file committed in HEAD", result.stderr)
