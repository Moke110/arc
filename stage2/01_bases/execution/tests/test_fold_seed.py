"""Fold/seed index products: row-prefix selection and the re-homed assignment.

The seam is the published product contract (ADR 0063): ``cohorts.csv`` keeps the
accepted assignment byte-identical, ``folds.csv``/``seeds.csv`` order frozen
identities by row, and ``manifest.json`` binds every file digest.  An experiment
takes the first ``n`` rows and a request beyond the table is an error, not a
partial answer.

Fixtures are tiny worked examples; the real accepted assignment is verified by
``python stage2/01_bases/execution/build_fold_n_seed.py --verify``.
"""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_fold_n_seed as cli
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.fold_seed import (FoldSeedError, FoldSeedIndex, build_fold_seed,
                             sha256, verify_fold_seed)

COHORTS_HEADER = "sample_id,drug_parent,dose_um,group_id,role,fold\n"
FOLDS_HEADER = "row,fold,n_wells,n_train,n_validation,n_test,n_test_drugs\n"
SEEDS_HEADER = "row,seed\n"


def fixture_cohorts():
    """Two parents per test fold, one validation parent and three training.

    Every parent has two doses and two replicate wells, so each fold partitions the
    same 24 wells: 8 test + 4 validation + 12 training.
    """
    lines = [COHORTS_HEADER]
    for fold in range(3):
        test = {fold * 2, fold * 2 + 1}
        remaining = sorted(set(range(6)) - test)
        validation = {remaining[0]}
        for parent in range(6):
            role = ("test" if parent in test else
                    "validation" if parent in validation else "train")
            for dose in (1, 10):
                for replicate in (0, 1):
                    lines.append(
                        f"w{parent}-{dose}-{replicate},{parent},{dose},"
                        f"g{parent}-{dose},{role},{fold}\n")
    return "".join(lines).encode()


def fixture_product(root):
    """A hand-written product directory: files plus an independently hashed manifest."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    cohorts = fixture_cohorts()
    folds = (FOLDS_HEADER
             + "0,0,24,12,4,8,2\n"
             + "1,1,24,12,4,8,2\n"
             + "2,2,24,12,4,8,2\n").encode()
    seeds = (SEEDS_HEADER + "0,17\n1,42\n2,2026\n").encode()
    (root / "cohorts.csv").write_bytes(cohorts)
    (root / "folds.csv").write_bytes(folds)
    (root / "seeds.csv").write_bytes(seeds)
    manifest = {
        "schema": "stage2-fold-seed-index-v1",
        "files": {"folds.csv": {"sha256": hashlib.sha256(folds).hexdigest(),
                                "bytes": len(folds)},
                  "seeds.csv": {"sha256": hashlib.sha256(seeds).hexdigest(),
                                "bytes": len(seeds)}},
        "cohorts": {"file": "cohorts.csv",
                    "sha256": hashlib.sha256(cohorts).hexdigest(),
                    "n_wells": 24, "n_drugs": 6},
        "folds": {"file": "folds.csv", "count": 3, "order": [0, 1, 2]},
        "seeds": {"file": "seeds.csv", "count": 3, "order": [17, 42, 2026]},
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return root, manifest


class FoldSeedIndexTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root, self.manifest = fixture_product(Path(temporary.name) / "product")

    def test_selection_returns_an_in_order_row_prefix(self):
        index = FoldSeedIndex.load(self.root)
        self.assertEqual(index.select_folds(2)["fold"].tolist(), [0, 1])
        self.assertEqual(index.select_seeds(2)["seed"].tolist(), [17, 42])
        self.assertEqual(index.select_folds(3)["fold"].tolist(), [0, 1, 2])
        self.assertEqual(index.select_seeds(3)["seed"].tolist(), [17, 42, 2026])
        # A prefix is the head of the table, so appending rows never reshuffles.
        pd.testing.assert_frame_equal(index.select_folds(2),
                                      index.folds.iloc[:2].reset_index(drop=True))

    def test_selection_beyond_the_table_is_rejected(self):
        index = FoldSeedIndex.load(self.root)
        for call, kind in ((index.select_folds, "fold"), (index.select_seeds, "seed")):
            with self.assertRaisesRegex(FoldSeedError, f"requested 4 {kind}"):
                call(4)

    def test_appending_a_seed_never_reshuffles_the_earlier_rows(self):
        (self.root / "seeds.csv").write_bytes(
            (SEEDS_HEADER + "0,17\n1,42\n2,2026\n3,7\n").encode())
        manifest = json.loads((self.root / "manifest.json").read_text())
        manifest["files"]["seeds.csv"] = {
            "sha256": sha256(self.root / "seeds.csv"),
            "bytes": (self.root / "seeds.csv").stat().st_size,
        }
        manifest["seeds"]["count"] = 4
        manifest["seeds"]["order"] = [17, 42, 2026, 7]
        (self.root / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        index = FoldSeedIndex.load(self.root)
        self.assertEqual(index.select_seeds(2)["seed"].tolist(), [17, 42])
        self.assertEqual(index.select_seeds(4)["seed"].tolist(), [17, 42, 2026, 7])
        with self.assertRaisesRegex(FoldSeedError, "requested 5 seed"):
            index.select_seeds(5)

    def test_manifest_digest_binds_every_file(self):
        index = FoldSeedIndex.load(self.root)
        self.assertEqual(index.select_folds(1)["n_wells"].tolist(), [24])
        tampered = self.root / "folds.csv"
        tampered.write_bytes(tampered.read_bytes().replace(b"24,12,4,8,2", b"24,11,5,8,2"))
        with self.assertRaisesRegex(FoldSeedError, "folds.csv"):
            FoldSeedIndex.load(self.root)

    def test_table_order_must_match_the_declared_order(self):
        (self.root / "seeds.csv").write_bytes(
            (SEEDS_HEADER + "0,42\n1,17\n2,2026\n").encode())
        manifest = json.loads((self.root / "manifest.json").read_text())
        manifest["files"]["seeds.csv"] = {
            "sha256": sha256(self.root / "seeds.csv"),
            "bytes": (self.root / "seeds.csv").stat().st_size,
        }
        (self.root / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        with self.assertRaisesRegex(FoldSeedError, "order"):
            FoldSeedIndex.load(self.root)


class BuildFoldSeedTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "accepted-cohorts.csv"
        self.source.write_bytes(fixture_cohorts())
        self.digest = sha256(self.source)

    def build(self, output=None):
        return build_fold_seed(self.source, output or self.root / "product",
                               accepted_sha256=self.digest)

    def test_assignment_is_rehomed_byte_identically_with_ordered_index_tables(self):
        manifest = self.build()
        product = self.root / "product"
        self.assertEqual((product / "cohorts.csv").read_bytes(),
                         self.source.read_bytes())
        folds = pd.read_csv(product / "folds.csv")
        self.assertEqual(folds.columns.tolist(),
                         ["row", "fold", "n_wells", "n_train", "n_validation",
                          "n_test", "n_test_drugs"])
        self.assertEqual(folds["row"].tolist(), [0, 1, 2])
        self.assertEqual(folds["fold"].tolist(), [0, 1, 2])
        self.assertEqual(folds["n_wells"].tolist(), [24, 24, 24])
        self.assertEqual(folds["n_train"].tolist(), [12, 12, 12])
        self.assertEqual(folds["n_validation"].tolist(), [4, 4, 4])
        self.assertEqual(folds["n_test"].tolist(), [8, 8, 8])
        self.assertEqual(folds["n_test_drugs"].tolist(), [2, 2, 2])
        seeds = pd.read_csv(product / "seeds.csv")
        self.assertEqual(seeds.columns.tolist(), ["row", "seed"])
        self.assertEqual(seeds["row"].tolist(), [0, 1, 2])
        self.assertEqual(seeds["seed"].tolist(), [17, 42, 2026])
        self.assertEqual(manifest["schema"], "stage2-fold-seed-index-v1")
        self.assertEqual(manifest["cohorts"]["sha256"], self.digest)
        for name in ("folds.csv", "seeds.csv"):
            self.assertEqual(manifest["files"][name]["sha256"],
                             sha256(product / name))
        index = FoldSeedIndex.load(product)
        self.assertEqual(index.select_folds(2)["fold"].tolist(), [0, 1])
        self.assertEqual(index.select_seeds(2)["seed"].tolist(), [17, 42])

    def test_rebuild_is_byte_deterministic(self):
        first = self.build(self.root / "one")
        second = self.build(self.root / "two")
        self.assertEqual(first, second)
        for name in ("cohorts.csv", "folds.csv", "seeds.csv", "manifest.json"):
            self.assertEqual((self.root / "one" / name).read_bytes(),
                             (self.root / "two" / name).read_bytes())

    def test_assignment_that_breaks_the_frozen_fold_contract_is_rejected(self):
        frame = pd.read_csv(self.source, keep_default_na=False)
        broken = frame.copy()
        broken.loc[frame.index[frame.role.eq("validation") & frame.fold.eq(0)][0], "role"] = "test"  # crosses a parent-dose group
        path = self.root / "broken.csv"
        broken.to_csv(path, index=False)
        with self.assertRaisesRegex(FoldSeedError, "crosses roles"):
            build_fold_seed(path, self.root / "out", accepted_sha256=sha256(path))
        shifted = frame.copy()
        shifted.loc[shifted.fold.eq(2), "fold"] = 3  # fold identities must be 0..n-1
        path = self.root / "shifted.csv"
        shifted.to_csv(path, index=False)
        with self.assertRaisesRegex(FoldSeedError, "0..n-1"):
            build_fold_seed(path, self.root / "out2", accepted_sha256=sha256(path))

    def test_a_source_that_is_not_the_accepted_assignment_is_rejected(self):
        with self.assertRaisesRegex(FoldSeedError, "accepted"):
            build_fold_seed(self.source, self.root / "product",
                            accepted_sha256="0" * 64)
        # The default gate is the accepted digest, not whatever the caller supplies.
        with self.assertRaisesRegex(FoldSeedError, "accepted"):
            build_fold_seed(self.source, self.root / "product2")

    def test_the_pinned_accepted_digest_is_the_published_assignment(self):
        self.assertEqual(
            "c1dcce7954b3fd1757dd68fdef4a9e5bf80a0015be6354d29c7be37628a843cc",
            __import__("bases.fold_seed", fromlist=["ACCEPTED_COHORTS_SHA256"])
            .ACCEPTED_COHORTS_SHA256)



class VerifyFoldSeedTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "accepted-cohorts.csv"
        self.source.write_bytes(fixture_cohorts())
        self.digest = sha256(self.source)
        self.product = self.root / "product"
        build_fold_seed(self.source, self.product, accepted_sha256=self.digest)

    def verify(self, expected=None):
        return verify_fold_seed(self.product,
                                expected_cohorts_sha256=expected or self.digest)

    def restamp(self, *names):
        manifest = json.loads((self.product / "manifest.json").read_text())
        for name in names:
            manifest["files"][name] = {
                "sha256": sha256(self.product / name),
                "bytes": (self.product / name).stat().st_size,
                "rows": len(pd.read_csv(self.product / name)),
            }
        (self.product / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    def test_accepted_product_passes_every_check(self):
        report = self.verify()
        self.assertEqual(report["result"], "PASS", report["failures"])
        self.assertEqual([check["check"] for check in report["checks"] if check["status"] != "PASS"], [])
        self.assertEqual(report["cohorts"]["sha256"], self.digest)
        self.assertEqual(report["folds"][0]["n_wells"], 24)
        self.assertEqual([row["seed"] for row in report["seeds"]], [17, 42, 2026])
        self.assertEqual(report["files"]["folds.csv"]["sha256"],
                         sha256(self.product / "folds.csv"))

    def test_tampered_bytes_are_reported(self):
        folds = self.product / "folds.csv"
        folds.write_bytes(folds.read_bytes().replace(b"24,12,4,8,2", b"24,12,4,8,3"))
        report = self.verify()
        self.assertEqual(report["result"], "FAIL")
        self.assertTrue(any("digest" in name for name in report["failures"]))

    def test_wrong_fold_counts_are_reported(self):
        folds = self.product / "folds.csv"
        folds.write_bytes(folds.read_bytes().replace(b"24,12,4,8,2", b"24,11,5,8,2"))
        self.restamp("folds.csv")
        report = self.verify()
        self.assertEqual(report["result"], "FAIL")
        self.assertTrue(any("count" in name for name in report["failures"]))

    def test_assignment_that_breaks_drug_disjointness_is_reported(self):
        frame = pd.read_csv(self.product / "cohorts.csv", keep_default_na=False)
        mask = frame.fold.eq(1) & frame.drug_parent.eq(0)
        self.assertTrue(mask.any())
        frame.loc[mask, "role"] = "test"
        frame.to_csv(self.product / "cohorts.csv", index=False)
        folds = pd.read_csv(self.product / "folds.csv")
        folds.loc[folds.fold.eq(1), ["n_train", "n_test"]] = [8, 12]
        folds.to_csv(self.product / "folds.csv", index=False)
        self.restamp("folds.csv")
        manifest = json.loads((self.product / "manifest.json").read_text())
        manifest["cohorts"]["sha256"] = sha256(self.product / "cohorts.csv")
        (self.product / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        report = self.verify(expected=sha256(self.product / "cohorts.csv"))
        self.assertEqual(report["result"], "FAIL")
        self.assertTrue(any("disjoint" in name or "cover" in name
                            for name in report["failures"]), report["failures"])

    def test_wrong_seed_order_is_reported(self):
        seeds = self.product / "seeds.csv"
        seeds.write_bytes(SEEDS_HEADER.encode() + b"0,42\n1,17\n2,2026\n")
        self.restamp("seeds.csv")
        report = self.verify()
        self.assertEqual(report["result"], "FAIL")
        self.assertTrue(any("seed" in name for name in report["failures"]))



class CommandLineTests(unittest.TestCase):
    def run_cli(self, argv):
        with contextlib.redirect_stdout(io.StringIO()):
            return cli.main(argv)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "accepted-cohorts.csv"
        self.source.write_bytes(fixture_cohorts())
        self.digest = sha256(self.source)
        self.product = self.root / "product"
        self.record = self.root / "record.json"

    def test_build_then_verify_is_offline_and_fail_closed(self):
        code = self.run_cli(["--cohorts", str(self.source), "--output", str(self.product),
                             "--expected-cohorts-sha256", self.digest,
                             "--record", str(self.record)])
        self.assertEqual(code, 0)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["verification"]["result"], "PASS")
        self.assertEqual(record["manifest"]["cohorts"]["sha256"], self.digest)
        self.assertEqual(self.run_cli(["--verify", "--output", str(self.product),
                                       "--expected-cohorts-sha256", self.digest]), 0)
        folds = self.product / "folds.csv"
        folds.write_bytes(folds.read_bytes().replace(b"0,0,24", b"0,0,25"))
        self.assertEqual(self.run_cli(["--verify", "--output", str(self.product),
                                       "--expected-cohorts-sha256", self.digest]), 1)

    def test_a_source_that_is_not_the_accepted_assignment_fails(self):
        code = self.run_cli(["--cohorts", str(self.source), "--output", str(self.product),
                             "--expected-cohorts-sha256", "0" * 64])
        self.assertEqual(code, 1)
        self.assertFalse((self.product / "manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
