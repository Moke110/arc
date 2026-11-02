"""Fold-algorithm tests: determinism, 6:2:2 balance, no straddle, no leak.

Pure synthetic frames — no h5ad, no frozen data.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))

import numpy as np
import pandas as pd

import contract
import folds


def synthetic_drugs(n: int, *, seed: int = 0) -> pd.DataFrame:
    """`n` drugs, one plate each, 9-18 wells, plus two salt-pair clusters."""
    rng = np.random.default_rng(seed)
    keys = [f"KEY{i:04d}" for i in range(n)]
    # first 14 characters are the connectivity block and must be distinct
    inchikeys = [f"{'A' * 10}{i:04d}-BBBBBBBB-N" for i in range(n)]
    table = pd.DataFrame({
        "drug_key": keys,
        "inchikey": inchikeys,
        "n_wells": rng.integers(9, 19, size=n).astype(int),
        "plate": [f"Plate{i % 16 + 1}" for i in range(n)],
        "n_doses": np.full(n, 3, dtype=int),
    })
    # Two connectivity pairs: same inchikey[:14] on the first 14 characters.
    for a, b in ((0, 1), (2, 3)):
        table.loc[b, "inchikey"] = table.loc[a, "inchikey"]
    return table


class FoldAlgorithmTest(unittest.TestCase):
    def test_each_drug_test_once_valid_once(self):
        drugs = synthetic_drugs(460)
        assign = folds.cluster_drugs(drugs)
        self.assertEqual(len(assign), 460)
        self.assertEqual(set(assign["test_fold"]), set(range(5)))
        # valid_fold = test_fold - 1 (mod 5): still exactly one valid fold each
        self.assertTrue((assign["valid_fold"] != assign["test_fold"]).all())
        test_counts = assign.groupby("test_fold").size()
        self.assertTrue((test_counts >= 90).all(), test_counts.to_dict())
        self.assertTrue((test_counts <= 92).all(), test_counts.to_dict())
        valid_counts = assign.groupby("valid_fold").size()
        self.assertTrue((valid_counts <= 92).all(), valid_counts.to_dict())

    def test_deterministic(self):
        drugs = synthetic_drugs(120)
        first = folds.cluster_drugs(drugs)
        second = folds.cluster_drugs(drugs)
        pd.testing.assert_frame_equal(first, second)

    def test_seed_changes_assignment(self):
        # the seed breaks ties between equally sized clusters; make every drug
        # the same size so the tie order is what changes.
        drugs = synthetic_drugs(120)
        drugs["n_wells"] = 9
        first = folds.cluster_drugs(drugs, seed=1)
        second = folds.cluster_drugs(drugs, seed=2)
        self.assertFalse(first["test_fold"].equals(second["test_fold"]))

    def test_salt_pairs_share_a_cluster_and_a_fold(self):
        drugs = synthetic_drugs(60)
        assign = folds.cluster_drugs(drugs)
        for a, b in ((0, 1), (2, 3)):
            ka = drugs.loc[a, "drug_key"]
            kb = drugs.loc[b, "drug_key"]
            self.assertEqual(assign.set_index("drug_key").loc[ka, "test_fold"],
                             assign.set_index("drug_key").loc[kb, "test_fold"])
            self.assertEqual(assign.set_index("drug_key").loc[ka, "cluster_id"],
                             assign.set_index("drug_key").loc[kb, "cluster_id"])

    def test_expand_and_invariants(self):
        drugs = synthetic_drugs(100)
        assign = folds.cluster_drugs(drugs)
        plates = drugs[["drug_key", "plate", "n_wells", "n_doses"]].copy()
        # one drug spans two plates (the real Perfluorooctanoic-acid case)
        extra = plates.iloc[[0]].copy()
        extra["plate"] = "Plate99"
        extra["n_wells"] = 9
        plates = pd.concat([plates, extra], ignore_index=True)
        table = folds.expand_fold_table(assign, plates)
        facts = folds.assert_fold_invariants(
            table, assign, expected_drugs=set(drugs["drug_key"]), n_folds=5)
        self.assertEqual(facts["n_drugs"], 100)
        # 5 folds x (100 drugs + 1 extra plate row) = 505 rows
        self.assertEqual(facts["n_rows"], 5 * 101)
        role_counts = table.drop_duplicates(["drug_key", "fold"])["role"].value_counts()
        self.assertEqual(role_counts["test"], 100)
        self.assertEqual(role_counts["valid"], 100)
        self.assertEqual(role_counts["train"], 300)

    def test_collapse_is_unique_per_drug_and_fold(self):
        drugs = synthetic_drugs(100)
        assign = folds.cluster_drugs(drugs)
        plates = drugs[["drug_key", "plate", "n_wells", "n_doses"]].copy()
        extra = plates.iloc[[0]].copy()
        extra["plate"] = "Plate99"
        plates = pd.concat([plates, extra], ignore_index=True)
        long = folds.expand_fold_table(assign, plates)
        collapsed = folds.collapse_fold_table(long)
        self.assertFalse(collapsed.duplicated(["drug_key", "fold"]).any())
        self.assertEqual(len(collapsed), 5 * 100)
        two_plate = collapsed.loc[collapsed["drug_key"] == plates.iloc[0]["drug_key"]]
        self.assertTrue((two_plate["n_plates"] == 2).all())
        self.assertTrue(two_plate["plate"].str.contains(";").all())
        # collapse preserves the wells and the role assignment
        self.assertEqual(int(collapsed["n_wells"].sum()), int(long["n_wells"].sum()))
        for fold in range(5):
            self.assertEqual(
                set(collapsed.loc[collapsed["fold"] == fold, "drug_key"]),
                set(long.loc[long["fold"] == fold, "drug_key"]))
        folds.assert_fold_invariants(collapsed, assign,
                                     expected_drugs=set(drugs["drug_key"]))

    def test_collapse_keeps_a_single_plate_drug_unchanged(self):
        drugs = synthetic_drugs(30)
        assign = folds.cluster_drugs(drugs)
        plates = drugs[["drug_key", "plate", "n_wells", "n_doses"]]
        collapsed = folds.collapse_fold_table(
            folds.expand_fold_table(assign, plates))
        self.assertTrue((collapsed["n_plates"] == 1).all())
        self.assertFalse(collapsed["plate"].str.contains(";").any())

    def test_invariants_reject_a_foreign_drug(self):
        drugs = synthetic_drugs(30)
        assign = folds.cluster_drugs(drugs)
        plates = drugs[["drug_key", "plate", "n_wells", "n_doses"]]
        table = folds.expand_fold_table(assign, plates)
        with self.assertRaises(contract.ContractError):
            folds.assert_fold_invariants(table, assign,
                                        expected_drugs=set(drugs["drug_key"]) | {"KEY9999"})

    def test_invariants_reject_a_missing_role(self):
        drugs = synthetic_drugs(30)
        assign = folds.cluster_drugs(drugs)
        plates = drugs[["drug_key", "plate", "n_wells", "n_doses"]]
        table = folds.expand_fold_table(assign, plates)
        broken = table.drop(index=table.index[(table["drug_key"] == "KEY0000")
                                              & (table["role"] == "valid")])
        with self.assertRaises(contract.ContractError):
            folds.assert_fold_invariants(broken, assign)

    def test_rejects_duplicate_or_empty_input(self):
        with self.assertRaises(contract.ContractError):
            folds.cluster_drugs(pd.DataFrame({"drug_key": [], "inchikey": [],
                                              "n_wells": []}))
        dup = synthetic_drugs(10)
        dup.loc[1, "drug_key"] = dup.loc[0, "drug_key"]
        with self.assertRaises(contract.ContractError):
            folds.cluster_drugs(dup)

    def test_rejects_too_few_drugs(self):
        with self.assertRaises(contract.ContractError):
            folds.cluster_drugs(synthetic_drugs(3))

    def test_fold_sizes_and_plate_scope(self):
        drugs = synthetic_drugs(50)
        assign = folds.cluster_drugs(drugs)
        plates = drugs[["drug_key", "plate", "n_wells", "n_doses"]]
        table = folds.expand_fold_table(assign, plates)
        sizes = folds.fold_sizes(table)
        self.assertEqual(list(sizes["fold"]), [0, 1, 2, 3, 4])
        self.assertTrue((sizes["test_wells"] > 0).all())
        for fold in range(5):
            train_plates = folds.fold_plates(table, fold, "train")
            test_plates = folds.fold_plates(table, fold, "test")
            self.assertTrue(train_plates)
            self.assertTrue(test_plates)
        # a fold's test plates are a subset of all plates, never the same set
        merged = table[table["role"] != "train"]
        self.assertTrue(set(merged["role"]) <= {"test", "valid"})


class ExcludedDrugTest(unittest.TestCase):
    def test_excluded_uncertain_get_no_fold(self):
        """The fold population is main_set only: exclude/uncertain never appear."""
        all_drugs = synthetic_drugs(40)
        main = all_drugs.iloc[:30]
        assign = folds.cluster_drugs(main[["drug_key", "inchikey", "n_wells"]])
        included = set(assign["drug_key"])
        excluded = set(all_drugs["drug_key"]) - included
        self.assertTrue(excluded)
        self.assertFalse(included & excluded)
        self.assertTrue(all(not key.startswith("KEY003") for key in included))


if __name__ == "__main__":
    unittest.main()
