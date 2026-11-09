"""The shared CM release is exercised through the Step 02 CLI."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


CLI = Path(__file__).resolve().parents[1] / "build_corpus.py"


class SharedCMTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        self.native = root / "native"
        self.release = root / "release"
        tables = self.native / "tables"
        tables.mkdir(parents=True)
        (tables / "genes.csv").write_text(
            "symbol,hgnc_id,ensembl_id\nA,HGNC:1,ENSG1\nFAM153B,HGNC:27323,ENSG00000182230\n"
        )
        self.crosswalk = tables / "gene_crosswalk.csv"
        self.crosswalk.write_text(
            "dataset_id,source_gene_id,hgnc_id\n"
            "tahoe,ENSG1,HGNC:1\ntahoe,1,HGNC:1\n"
            "tahoe,ENSG00000289731,HGNC:27323\ntahoe,61482,HGNC:27323\n"
            "tahoe,ENSG3,\ntahoe,3,\n"
        )
        self.write_native("tahoe", [[90, 6, 2], [100, 0, 8]])
        self.write_native("tenv", [[8, 900, 2], [0, 0, 0], [5, 100, 0]])

    def write_native(self, source, counts):
        if source == "tahoe":
            obs = pd.DataFrame({"plate": ["A", "B"], "is_control": [False, True]},
                               index=pd.Index(["h1", "h2"], name="profile_id"))
            var = pd.DataFrame(
                {"ensembl_id": ["ENSG3", "ENSG00000289731", "ENSG1"],
                 "token_id": [3, 61482, 1]}, index=["off-axis", "FAM153B", "A"])
            matrix = csr_matrix(np.array(counts, dtype=np.int32))
        else:
            obs = pd.DataFrame(
                {"plate": ["D", "E", "F"], "is_control": [True, False, False],
                 "is_blank": [False, True, False]},
                index=pd.Index(["w1", "w2", "w3"], name="well"))
            var = pd.DataFrame(
                {"ensembl_id": ["ENSG1", "ENSG3", "ENSG00000182230"],
                 "hgnc_id": ["HGNC:1", "HGNC:3", "HGNC:27323"]},
                index=["A", "off-axis", "FAM153B"])
            matrix = np.array(counts)
        ad.AnnData(matrix, obs=obs, var=var).write_h5ad(self.native / f"{source}.h5ad")

    def cli(self, mode):
        return subprocess.run(
            [sys.executable, str(CLI), "--shared-cm", mode, "--shared-cm-fixture",
             "--shared-cm-native", str(self.native), "--shared-cm-output", str(self.release)],
            capture_output=True, text=True,
        )

    def test_build_and_verify_uses_observed_columns_in_shared_order(self):
        build = self.cli("build")
        self.assertEqual(build.returncode, 0, build.stderr)
        self.assertEqual((self.release / "genes.csv").read_bytes(),
                         (self.native / "tables/genes.csv").read_bytes())
        expected = {
            "tahoe": [[math.log1p(2500), math.log1p(7500)], [math.log1p(10000), 0]],
            "tenv": [[math.log1p(8000), math.log1p(2000)], [math.log1p(10000), 0]],
        }
        for source, values in expected.items():
            expression = np.load(self.release / f"{source}_expression.npy", mmap_mode="r")
            self.assertEqual(expression.dtype, np.float32)
            self.assertTrue(expression.flags.c_contiguous)
            np.testing.assert_allclose(expression, values, rtol=1e-7)
        profiles = pd.read_csv(self.release / "tenv_profiles.csv").to_dict("records")
        self.assertEqual(profiles, [
            {"source_row": 0, "profile_id": "w1", "plate": "D", "is_control": True},
            {"source_row": 2, "profile_id": "w3", "plate": "F", "is_control": False},
        ])
        manifest = json.loads((self.release / "manifest.json").read_text())
        self.assertEqual(manifest["sources"]["tenv"]["population"],
                         {"raw": 3, "blank": 1, "treated": 1, "control": 1})
        verify = self.cli("verify")
        self.assertEqual(verify.returncode, 0, verify.stderr)
        spots = json.loads(verify.stdout)["sources"]["tahoe"]["fam153b"]
        self.assertEqual(spots["nonzero"]["raw_count"], 6)
        self.assertEqual(spots["nonzero"]["projected_library_total"], 8)
        self.assertEqual(spots["zero"]["raw_count"], 0)

    def test_exact_ensembl_takes_priority_over_crosswalk_aliases(self):
        original = self.crosswalk.read_text()
        self.crosswalk.write_text(original.replace("tahoe,ENSG3,\ntahoe,3,\n",
                                                   "tahoe,ENSG3,HGNC:1\ntahoe,3,HGNC:1\n"))
        build = self.cli("build")
        self.assertEqual(build.returncode, 0, build.stderr)
        mapping = pd.read_csv(self.release / "tahoe_mapping.csv")
        self.assertEqual(list(mapping.native_column), [2, 1])
        self.assertEqual(list(mapping.resolution), ["exact_ensembl", "crosswalk_hgnc"])

    def test_exact_ensembl_can_use_identity_corroborated_by_crosswalk(self):
        original = self.crosswalk.read_text()
        self.crosswalk.write_text(original.replace("tahoe,ENSG1,HGNC:1\ntahoe,1,HGNC:1\n",
                                                   "tahoe,ENSG1,\ntahoe,1,\nscreen_ctox,ENSG1,HGNC:1\n"))
        build = self.cli("build")
        self.assertEqual(build.returncode, 0, build.stderr)
        mapping = pd.read_csv(self.release / "tahoe_mapping.csv")
        self.assertEqual(list(mapping.resolution), ["exact_ensembl_corroborated", "crosswalk_hgnc"])

    def test_build_rejects_unresolved_or_conflicting_identity(self):
        original = self.crosswalk.read_text()
        fam = "tahoe,ENSG00000289731,HGNC:27323\ntahoe,61482,HGNC:27323\n"
        off_axis = "tahoe,ENSG3,\ntahoe,3,\n"
        for text, expected in [
            (original.replace(fam, fam.replace("HGNC:27323", "")), "missing HGNC"),
            (original.replace(fam, fam.replace("61482,HGNC:27323", "61482,HGNC:5")), "conflicting"),
            (original.replace(off_axis, "tahoe,ENSG3,HGNC:27323\ntahoe,3,HGNC:27323\n"), "ambiguous"),
            (original.replace("tahoe,ENSG1,HGNC:1\ntahoe,1,HGNC:1\n",
                              "tahoe,ENSG1,HGNC:2\ntahoe,1,HGNC:2\n"), "identity conflicts"),
        ]:
            with self.subTest(expected=expected):
                self.crosswalk.write_text(text)
                build = self.cli("build")
                self.assertNotEqual(build.returncode, 0)
                self.assertIn(expected, build.stderr)
                self.assertFalse(self.release.exists())

    def test_build_rejects_invalid_counts_and_zero_projected_libraries(self):
        for counts, expected in [
            ([[8, -1, 2], [0, 0, 0], [5, 100, 0]], "negative raw counts"),
            ([[0, 100, 0], [0, 0, 0], [5, 100, 0]], "zero projected library"),
        ]:
            with self.subTest(expected=expected):
                self.write_native("tenv", counts)
                build = self.cli("build")
                self.assertNotEqual(build.returncode, 0)
                self.assertIn(expected, build.stderr)
                self.assertFalse(self.release.exists())

    def test_verify_replays_raw_counts_even_if_changed_expression_is_rehashed(self):
        build = self.cli("build")
        self.assertEqual(build.returncode, 0, build.stderr)
        expression_path = self.release / "tahoe_expression.npy"
        expression = np.load(expression_path, mmap_mode="r+")
        # Swap the two nonzero values, keeping this row's library and all zeros.
        expression[0] = expression[0, ::-1]
        expression.flush()
        del expression
        manifest_path = self.release / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["files"][expression_path.name] = hashlib.sha256(expression_path.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest))
        verify = self.cli("verify")
        self.assertNotEqual(verify.returncode, 0)
        self.assertIn("normalization differs from raw counts", verify.stderr)


if __name__ == "__main__":
    unittest.main()
