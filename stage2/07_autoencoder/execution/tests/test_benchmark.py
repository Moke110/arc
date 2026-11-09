"""Bounded reconstruction fixtures exercise the fit and replay workflow."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "script"
sys.path.insert(0, str(SCRIPT))


class BenchmarkTest(unittest.TestCase):
    def setUp(self):
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.root = Path(work.name)
        self.release = self.root / "release"
        self.release.mkdir()
        self.split_path = self.root / "split.json"
        self.artifact = self.root / "mean"
        self.profile_plates = {"tahoe": ["A", "B", "C", "A"], "tenv": ["A", "B", "D"]}
        self.expressions = {
            "tahoe": np.array([[1, 2, 3], [2, 4, 8], [3, 1, 2], [4, 6, 8]], dtype=np.float32),
            "tenv": np.array([[8, 3, 0], [5, 5, 5], [0, 6, 4]], dtype=np.float32),
        }
        (self.release / "genes.csv").write_text("symbol,hgnc_id\nX,HGNC:1\nY,HGNC:2\nZ,HGNC:3\n")
        for source, expression in self.expressions.items():
            np.save(self.release / f"{source}_expression.npy", expression)
            with (self.release / f"{source}_profiles.csv").open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["source_row", "profile_id", "plate", "is_control"])
                for index, plate in enumerate(self.profile_plates[source]):
                    writer.writerow([index * 2, f"{source}-{index}", plate, index == 0])
            (self.release / f"{source}_mapping.csv").write_text("fixture mapping\n")
        self.write_manifest()

    def write_manifest(self):
        files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in self.release.iterdir() if p.name != "manifest.json"}
        manifest = {
            "schema_version": 1,
            "fixture": True,
            "normalization": "log1p(CP10K of raw counts projected to the ordered CM axis)",
            "inputs": {"genes.csv": files["genes.csv"]},
            "sources": {source: {"shape": list(expr.shape)}
                        for source, expr in self.expressions.items()},
            "files": files,
        }
        (self.release / "manifest.json").write_text(json.dumps(manifest))

    def cli(self, action, *args, success=True):
        result = subprocess.run(
            [sys.executable, str(SCRIPT / "reconstruct.py"), action, "--fixture",
             "--release", str(self.release), "--split", str(self.split_path), *map(str, args)],
            capture_output=True, text=True,
        )
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
            return json.loads(result.stdout)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        return result.stderr

    def test_mean_fit_and_independent_reload_use_frozen_roles(self):
        from benchmark import Corpus, Split

        split_record = self.cli("freeze-split")
        split = Split(Corpus(self.release, fixture=True), self.split_path)
        training = []
        for source, plates in self.profile_plates.items():
            rank = sorted(set(plates), key=lambda plate: hashlib.sha256(
                json.dumps([20261109, source, plate]).encode()).digest())
            expected = {"validation": rank[:1], "test": rank[1:2], "train": rank[2:]}
            for role, role_plates in expected.items():
                rows = [i for i, plate in enumerate(plates) if plate in role_plates]
                np.testing.assert_array_equal(split.rows[source][role], rows)
                record = split_record["sources"][source]["roles"][role]
                self.assertEqual(record["profiles"], len(rows))
                self.assertEqual(record["plates"], role_plates)
                self.assertEqual(len(record["membership_sha256"]), 64)
                if role == "train":
                    training.extend(self.expressions[source][rows])
        self.cli("fit", "--output", self.artifact, "--batch-size", 1)
        mean = np.load(self.artifact / "mean.npy")
        np.testing.assert_array_equal(mean, np.mean(np.array(training), axis=0))

        for role in ("train", "validation", "test"):
            path = self.root / f"{role}.json"
            args = ["--artifact", self.artifact, "--role", role, "--batch-size", 1]
            scores = self.cli("score", *args, "--output", path)
            replay = self.cli("reload", *args, "--expected", path)
            self.assertEqual(scores, replay)
            for source, result in scores["sources"].items():
                self.assertEqual(result["profiles"], len(split.rows[source][role]))
                self.assertEqual(result["relative_squared_error"], 1.0)

        with (self.artifact / "mean.npy").open("ab") as handle:
            handle.write(b"changed")
        self.assertIn("digest", self.cli("reload", "--artifact", self.artifact, "--expected", path, success=False))

    def test_identity_and_undefined_metrics(self):
        from benchmark import Corpus, Split, score

        self.cli("freeze-split")
        corpus = Corpus(self.release, fixture=True)
        split = Split(corpus, self.split_path)
        for role in ("train", "validation", "test"):
            scores = score(corpus, split, role, lambda values: values, np.zeros(3), 1)
            for source, result in scores["sources"].items():
                rows = split.rows[source][role]
                defined = int(np.count_nonzero(np.ptp(self.expressions[source][rows], axis=1)))
                self.assertEqual(result["profiles"], len(rows))
                self.assertEqual(result["defined_pcc_profiles"], defined)
                if defined:
                    self.assertAlmostEqual(result["pcc"], 1.0)
                else:
                    self.assertIsNone(result["pcc"])
                self.assertEqual(result["rmse"], 0.0)
                self.assertEqual(result["relative_squared_error"], 0.0)
        role = next(role for role, rows in split.rows["tenv"].items() if 1 in rows)
        scores = score(corpus, split, role, lambda values: values, np.full(3, 5), 2)
        self.assertIsNone(scores["sources"]["tenv"]["relative_squared_error"])

    def test_refuses_changed_split_and_release(self):
        self.cli("freeze-split")
        frozen = json.loads(self.split_path.read_text())
        frozen["sources"]["tahoe"]["roles"]["test"]["plates"] = ["other"]
        self.split_path.write_text(json.dumps(frozen))
        self.assertIn("split", self.cli("fit", "--output", self.artifact, success=False))
        with (self.release / "tenv_expression.npy").open("ab") as handle:
            handle.write(b"changed")
        self.assertIn("digest", self.cli("fit", "--output", self.artifact, success=False))


if __name__ == "__main__":
    unittest.main()
