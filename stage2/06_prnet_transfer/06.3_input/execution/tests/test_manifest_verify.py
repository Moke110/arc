"""End-to-end build/verify tests on a stub build (no frozen artifacts written).

Runs the real pipeline with `--stub-handoffs --skip-tahoe` into a temporary
output root, then verifies it, tampers with it and checks the failure modes.
Requires the frozen CM view (skipped on a tracked-files-only worktree).
"""

from __future__ import annotations

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))

import anndata as ad
import numpy as np
import pandas as pd

import build_input
import contract

HAVE_DATA = contract.CM_VIEW.exists()


def _writer_env_ok() -> bool:
    major, minor = (int(x) for x in contract.runtime_env()["anndata"].split(".")[:2])
    return (major, minor) >= (0, 11)


@unittest.skipUnless(HAVE_DATA, f"frozen data not under {contract.DATA_ROOT}")
@unittest.skipUnless(_writer_env_ok(), "the builder needs anndata >= 0.11")
class StubBuildVerifyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="06_3_stub_")
        cls.root = Path(cls.tmp.name)
        cls.rc = build_input.main([
            "--stub-handoffs", "--skip-tahoe", "--out-root", str(cls.root)])
        cls.out = contract.output_dir(cls.root)
        cls.manifest = json.loads((cls.out / "manifest.json").read_text())

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_build_succeeds_and_is_labelled_stub(self):
        self.assertEqual(self.rc, 0, "stub build failed")
        self.assertEqual(self.manifest["status"], "stub")
        self.assertTrue(self.manifest["stub_handoffs"])
        for name in ("cm_transfer.h5ad", "folds.parquet", "consensus_controls.npz",
                     "consensus_membership.parquet", "locked_config.json",
                     "counts.csv", "fold_sizes.csv", "consensus_summary.csv",
                     "var_checks.csv", "manifest.json"):
            self.assertTrue((self.out / name).exists(), name)
        # the 16 GB Tahoe artifact is not written in this mode
        self.assertFalse((self.out / "tahoe_transfer.h5ad").exists())

    def test_verify_accepts_a_stub_only_with_allow_stub(self):
        self.assertEqual(build_input.main([
            "--verify", "--out-root", str(self.root)]), 1)
        self.assertEqual(build_input.main([
            "--verify", "--allow-stub", "--out-root", str(self.root)]), 0)

    def test_verify_rejects_a_tampered_artifact(self):
        scratch = Path(tempfile.mkdtemp(prefix="06_3_tamper_"))
        try:
            target = contract.output_dir(scratch)
            shutil.copytree(self.out, target)
            with open(target / "folds.parquet", "ab") as fh:
                fh.write(b"tamper")
            self.assertEqual(build_input.main([
                "--verify", "--allow-stub", "--out-root", str(scratch)]), 1)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

    def test_verify_rejects_a_changed_locked_config(self):
        scratch = Path(tempfile.mkdtemp(prefix="06_3_config_"))
        try:
            target = contract.output_dir(scratch)
            shutil.copytree(self.out, target)
            path = target / "locked_config.json"
            config = json.loads(path.read_text())
            config["optimisation"]["batch_size"] = 256
            path.write_text(json.dumps(config))
            self.assertEqual(build_input.main([
                "--verify", "--allow-stub", "--out-root", str(scratch)]), 1)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

    def test_rebuild_without_force_is_refused(self):
        self.assertEqual(build_input.main([
            "--stub-handoffs", "--skip-tahoe", "--out-root", str(self.root)]), 1)

    def test_cm_artifact_contract(self):
        cm = ad.read_h5ad(self.out / "cm_transfer.h5ad")
        self.assertEqual(list(cm.obs.columns), contract.CM_OBS)
        self.assertEqual(list(cm.var.columns), contract.VAR_COLUMNS)
        self.assertEqual(int(cm.n_vars), contract.N_GENES)
        self.assertTrue((np.asarray(cm.X) >= 0).all())
        self.assertTrue(np.isfinite(np.asarray(cm.layers["log1p_cp10k"])).all())
        self.assertNotIn("control_mean", cm.uns)
        # main_set is the include-and-QC flag and only applies to treated rows
        main = cm.obs["main_set"].astype(bool)
        self.assertFalse(main[cm.obs["is_control"].astype(bool)].any())
        self.assertTrue((cm.obs.loc[main, "eligibility_decision"].astype(str)
                         == "include").all())

    def test_folds_cover_exactly_the_main_set(self):
        cm = ad.read_h5ad(self.out / "cm_transfer.h5ad", backed="r")
        main = set(cm.obs.loc[cm.obs["main_set"].astype(bool), "drug_key"].astype(str))
        cm.file.close()
        table = pd.read_parquet(self.out / "folds.parquet")
        self.assertEqual(set(table["drug_key"].astype(str)), main)
        self.assertEqual(set(table["role"]), {"train", "valid", "test"})
        per_drug = table.drop_duplicates(["drug_key", "fold"])["role"].value_counts()
        self.assertEqual(per_drug["test"], len(main))
        self.assertEqual(per_drug["valid"], len(main))
        self.assertEqual(per_drug["train"], 3 * len(main))

    def test_consensus_recomputes_and_weights_sum(self):
        cm = ad.read_h5ad(self.out / "cm_transfer.h5ad")
        membership = pd.read_parquet(self.out / "consensus_membership.parquet")
        npz = np.load(self.out / "consensus_controls.npz")
        layer = np.asarray(cm.layers["log1p_cp10k"])
        import consensus as consensus_mod
        recomputed = consensus_mod.recompute_c_star(layer, membership, n_folds=5)
        np.testing.assert_allclose(recomputed, npz["c_star"], atol=1e-6)
        self.assertEqual(npz["c_star"].shape, (5, contract.N_GENES))

    def test_counts_csv_matches_manifest(self):
        table = pd.read_csv(self.out / "counts.csv", keep_default_na=False)
        flat = dict(zip(table["key"].astype(str), table["value"].astype(str)))
        expected = build_input._flatten_dict(self.manifest["counts"])
        for key, value in expected.items():
            self.assertEqual(flat.get(key), value, key)


class MissingHandoffTest(unittest.TestCase):
    @unittest.skipUnless(HAVE_DATA, f"frozen data not under {contract.DATA_ROOT}")
    def test_build_fails_closed_without_handoffs(self):
        with tempfile.TemporaryDirectory(prefix="06_3_nohandoff_") as tmp:
            rc = build_input.main(["--out-root", tmp])
            self.assertEqual(rc, 1)


if __name__ == "__main__":
    unittest.main()
