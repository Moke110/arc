"""Exercise the Step 07 release through its command-line interface."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


CLI = next(path for path in Path(__file__).resolve().parents
           if path.name == "execution") / "script/lxpert.py"


class InputReleaseTest(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.release = self.root / "input"
        (self.root / "tables").mkdir()
        (self.root / "tables/genes.csv").write_text(
            "symbol,hgnc_id,ensembl_id\n"
            "A,HGNC:1,ENSG1\nFAM153B,HGNC:2,ENSG00000182230\nC,HGNC:3,ENSG3\n"
        )
        for source in ("tahoe", "tenv"):
            self.write_corpus(source)

    def write_corpus(self, source, *, missing_gene=False, zero_projected=False):
        obs = pd.DataFrame(
            {"plate": ["p3", "p1", "p1", "p2", "p3"],
             "is_control": [False, True, False, False, True]},
            index=pd.Index(["first", "control", "blank", "third", "last"],
                           name="profile_id" if source == "tahoe" else "well"),
        )
        if source == "tahoe":
            var = pd.DataFrame({"ensembl_id": ["ENSG3", "ENSG0", "ENSG1"]},
                               index=["renamed-C", "off-axis", "renamed-A"])
            counts = np.tile([3, 900, 1], (len(obs), 1)).astype(np.int32)
        else:
            obs["is_blank"] = [False, False, True, False, False]
            var = pd.DataFrame({"hgnc_id": ["HGNC:3", "HGNC:0", "HGNC:2", "HGNC:1"]},
                               index=["C", "off-axis", "FAM153B", "A"])
            counts = np.tile([3, 900, 4, 1], (len(obs), 1)).astype(np.int32)
        if missing_gene:
            var.iloc[0, 0] = "unknown-id"
        if zero_projected:
            counts[0, :] = 0
            counts[0, 1] = 900
        expression = ad.AnnData(X=csr_matrix(counts) if source == "tahoe" else counts,
                                obs=obs, var=var)
        expression.layers["log1p_cp10k"] = np.full(counts.shape, 42, dtype=np.float32)
        expression.write_h5ad(self.root / f"{source}.h5ad")

    def run_cli(self, action, *args, error=None):
        command = [sys.executable, str(CLI), action, "--release", str(self.release), *args]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if error:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(error, result.stderr)
        else:
            self.assertEqual(result.returncode, 0, result.stderr)

    def prepare(self, *args, **kwargs):
        self.run_cli("prepare", "--corpus-dir", str(self.root), "--fixture", *args, **kwargs)

    def test_projects_raw_counts_then_normalizes_and_keeps_source_identity(self):
        self.prepare()
        self.run_cli("verify")
        for source, expected, observed, rows in (
            ("tahoe", [2500, 0, 7500], [True, False, True], [0, 1, 2, 3, 4]),
            ("tenv", [1250, 5000, 3750], [True, True, True], [0, 1, 3, 4]),
        ):
            profiles = self.profiles(source)
            self.assertEqual([int(row["source_row"]) for row in profiles], rows)
            self.assertEqual([row["profile_id"] for row in profiles],
                             [["first", "control", "blank", "third", "last"][i] for i in rows])
            expression = np.load(self.release / f"{source}_expression.npy")
            self.assertEqual(expression.dtype, np.float32)
            np.testing.assert_allclose(expression, np.tile(np.log1p(expected), (len(rows), 1)))
            np.testing.assert_array_equal(np.load(self.release / f"{source}_target_mask.npy"),
                                          observed)


    def profiles(self, source):
        with (self.release / f"{source}_profiles.csv").open() as file:
            return list(csv.DictReader(file))

    def test_development_partition_is_reproducible_and_covers_refit(self):
        self.prepare("--seed", "17")
        first_release = self.release
        self.release = self.root / "repeated"
        self.prepare("--seed", "17")
        self.assertEqual((self.release / "manifest.json").read_bytes(),
                         (first_release / "manifest.json").read_bytes())
        manifest = json.loads((self.release / "manifest.json").read_text())
        self.assertEqual(manifest["split_seed"], 17)
        for source in ("tahoe", "tenv"):
            profiles = self.profiles(source)
            groups = {
                role: {row["plate"] for row in profiles if row["development"] == role}
                for role in ("train", "validation")
            }
            self.assertTrue(groups["train"])
            self.assertTrue(groups["validation"])
            self.assertTrue(groups["train"].isdisjoint(groups["validation"]))
            self.assertEqual(groups["train"] | groups["validation"], {"p1", "p2", "p3"})
            membership = manifest["sources"][source]
            self.assertEqual(sum(role["profiles"] for role in membership["development"].values()),
                             membership["refit"]["profiles"])
            self.assertEqual(membership["refit"]["profiles"], len(profiles))

    def test_verify_rejects_corrupted_expression_and_changed_source(self):
        self.prepare()
        self.run_cli("verify", "--corpus-dir", str(self.root))
        with (self.root / "tahoe.h5ad").open("ab") as file:
            file.write(b"changed")
        self.run_cli("verify", "--corpus-dir", str(self.root), error="source digest mismatch: tahoe")
        with (self.release / "tenv_expression.npy").open("ab") as file:
            file.write(b"changed")
        self.run_cli("verify", error="release file digest mismatch: tenv_expression.npy")

    def test_verify_rejects_mixed_plate_even_if_file_digest_was_updated(self):
        self.prepare()
        profiles = self.profiles("tahoe")
        profiles[0]["development"] = (
            "train" if profiles[0]["development"] == "validation" else "validation"
        )
        path = self.release / "tahoe_profiles.csv"
        with path.open("w") as file:
            writer = csv.DictWriter(file, fieldnames=profiles[0])
            writer.writeheader()
            writer.writerows(profiles)
        manifest_path = self.release / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["files"][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest))
        self.run_cli("verify", error="invalid development plate partition")

    def test_prepare_rejects_missing_axis_gene(self):
        self.write_corpus("tenv", missing_gene=True)
        self.prepare(error="tenv: unexpected gene mapping gaps")
        self.assertFalse(self.release.exists())

    def test_prepare_rejects_zero_projected_library(self):
        self.write_corpus("tahoe", zero_projected=True)
        self.prepare(error="tahoe: zero projected library size at source row 0")
        self.assertFalse(self.release.exists())

    def test_fixture_opt_in_is_required_for_noncanonical_axis(self):
        self.run_cli("prepare", "--corpus-dir", str(self.root),
                     error="gene axis differs from the frozen ordered CM axis")
        self.assertFalse(self.release.exists())
