"""Consensus-control tests: equal plate weighting, permissible pool, leak gates.

Synthetic obs/log frames only — no h5ad, no frozen data.  The values are
constructed so the expected `c*` is exact arithmetic.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))

import numpy as np
import pandas as pd

import consensus
import contract

N_FOLDS_TEST = 2
N_GENES_TEST = 8

# Per-plate control level: c* is then a clean arithmetic mean of plate levels.
PLATE_LEVEL = {"P1": 1.0, "P2": 3.0, "P3": 7.0, "P4": 11.0}


def synthetic_obs(plates=("P1", "P2", "P3"), wells_per_plate=3):
    rows = []
    for plate in plates:
        for index in range(wells_per_plate):
            rows.append({"observation_id": f"{plate}-ctrl{index}", "plate": plate,
                         "is_control": True, "vehicle": "DMSO", "qc_pass": True})
        rows.append({"observation_id": f"{plate}-media", "plate": plate,
                     "is_control": True, "vehicle": "Media", "qc_pass": True})
        rows.append({"observation_id": f"{plate}-treated", "plate": plate,
                     "is_control": False, "vehicle": "DMSO", "qc_pass": True})
    obs = pd.DataFrame(rows)
    obs.index = pd.Index(obs["observation_id"].astype(str), name="well")
    return obs


def synthetic_log(obs):
    log1p = np.zeros((len(obs), N_GENES_TEST), dtype=np.float32)
    for row, (_, record) in enumerate(obs.iterrows()):
        if record["is_control"]:
            log1p[row, :] = PLATE_LEVEL[record["plate"]]
        else:
            log1p[row, :] = 100.0 + PLATE_LEVEL[record["plate"]]
    return log1p


def synthetic_folds():
    """fold 0 trains on P1+P2 (test drug on P2); fold 1 trains on P2+P3."""
    rows = []
    for fold, train_plates, test_plate in ((0, ("P1", "P2"), "P2"),
                                           (1, ("P2", "P3"), "P3")):
        for plate in train_plates:
            rows.append({"drug_key": f"key-{plate}", "fold": fold, "role": "train",
                         "plate": plate})
        rows.append({"drug_key": f"key-test-{fold}", "fold": fold, "role": "test",
                     "plate": test_plate})
        rows.append({"drug_key": "key-valid", "fold": fold, "role": "valid",
                     "plate": train_plates[0]})
    return pd.DataFrame(rows)


class ConsensusTest(unittest.TestCase):
    def setUp(self):
        self.obs = synthetic_obs()
        self.log1p = synthetic_log(self.obs)
        self.folds = synthetic_folds()
        self.c_star, self.membership, self.summary, _, self.plate_weights = (
            consensus.build_consensus(self.log1p, self.obs, self.folds,
                                      n_folds=N_FOLDS_TEST, n_genes=N_GENES_TEST))

    def test_equal_plate_weighting(self):
        self.assertEqual(self.summary["n_plates"].tolist(), [2, 2])
        self.assertAlmostEqual(float(self.plate_weights[0, 0]), 0.5, places=6)
        per_plate = self.membership.groupby(["fold", "plate"], observed=True)["weight"].sum()
        np.testing.assert_allclose(per_plate.to_numpy(), 0.5, atol=1e-6)
        totals = self.membership.groupby("fold")["weight"].sum()
        np.testing.assert_allclose(totals.to_numpy(), 1.0, atol=1e-6)

    def test_cstar_is_the_mean_of_plate_log_means(self):
        # fold 0 pools P1 (level 1.0) and P2 (level 3.0) -> 2.0
        np.testing.assert_allclose(self.c_star[0], 2.0, atol=1e-6)
        # fold 1 pools P2 (3.0) and P3 (7.0) -> 5.0
        np.testing.assert_allclose(self.c_star[1], 5.0, atol=1e-6)
        self.assertEqual(self.c_star.dtype, np.float32)
        self.assertEqual(self.c_star.shape, (N_FOLDS_TEST, N_GENES_TEST))

    def test_pool_is_dmso_only_and_never_media(self):
        vehicles = set(self.membership["vehicle"].astype(str))
        self.assertEqual(vehicles, {"DMSO"})
        idx = self.membership["row"].to_numpy()
        self.assertTrue(bool(self.obs["is_control"].to_numpy()[idx].all()))
        self.assertFalse((self.obs["vehicle"].astype(str).to_numpy()[idx] == "Media").any())

    def test_pool_excludes_held_out_plates(self):
        # fold 0 trains on P1+P2, so held-out P3 never contributes to it
        fold0 = self.membership.loc[self.membership["fold"] == 0]
        self.assertEqual(set(fold0["plate"]), {"P1", "P2"})
        fold1 = self.membership.loc[self.membership["fold"] == 1]
        self.assertEqual(set(fold1["plate"]), {"P2", "P3"})

    def test_treated_rows_are_never_in_the_pool(self):
        idx = self.membership["row"].to_numpy()
        treated = ~self.obs["is_control"].to_numpy()[idx]
        self.assertFalse(treated.any())
        # and the treated expression (100 + level) is absent from c*
        self.assertTrue((self.c_star < 50).all())

    def test_recompute_from_membership_matches(self):
        recomputed = consensus.recompute_c_star(
            self.log1p, self.membership, n_folds=N_FOLDS_TEST, n_genes=N_GENES_TEST)
        np.testing.assert_allclose(recomputed, self.c_star, atol=1e-6)

    def test_qc_failed_control_is_excluded(self):
        obs = synthetic_obs()
        bad_id = obs.index[0]
        obs.loc[bad_id, "qc_pass"] = False
        c_star, membership, summary, _, _ = consensus.build_consensus(
            synthetic_log(obs), obs, self.folds, n_folds=N_FOLDS_TEST,
            n_genes=N_GENES_TEST)
        self.assertNotIn(bad_id, set(membership["observation_id"]))
        self.assertEqual(summary["n_wells"].tolist(), [5, 6])
        np.testing.assert_allclose(c_star[0], 2.0, atol=1e-6)

    def test_no_qc_pass_control_fails_closed(self):
        obs = synthetic_obs()
        # both plates that carry a training drug of fold 0 lose every control
        obs.loc[obs["plate"].isin(["P1", "P2"]), "qc_pass"] = False
        with self.assertRaises(contract.ContractError):
            consensus.build_consensus(synthetic_log(obs), obs, self.folds,
                                      n_folds=N_FOLDS_TEST, n_genes=N_GENES_TEST)

    def test_a_plate_losing_all_controls_is_skipped_not_fabricated(self):
        obs = synthetic_obs()
        obs.loc[obs["plate"] == "P1", "qc_pass"] = False
        c_star, membership, summary, _, _ = consensus.build_consensus(
            synthetic_log(obs), obs, self.folds, n_folds=N_FOLDS_TEST,
            n_genes=N_GENES_TEST)
        self.assertEqual(set(membership.loc[membership["fold"] == 0, "plate"]), {"P2"})
        self.assertEqual(summary.loc[summary["fold"] == 0, "n_plates"].iloc[0], 1)
        np.testing.assert_allclose(c_star[0], 3.0, atol=1e-6)

    def test_provenance_gate_rejects_a_held_out_plate(self):
        poisoned = self.membership.copy()
        row = poisoned.index[poisoned["fold"] == 0][0]
        poisoned.loc[row, "plate"] = "P3"     # no training drug on P3 in fold 0
        obs = self.obs.copy()
        obs.iloc[poisoned.loc[row, "row"], obs.columns.get_loc("plate")] = "P3"
        with self.assertRaises(contract.ContractError):
            consensus.assert_consensus_provenance(poisoned, obs, self.folds,
                                                  n_folds=N_FOLDS_TEST)

    def test_provenance_gate_rejects_a_treated_row(self):
        poisoned = self.membership.copy()
        treated = int(np.flatnonzero(~self.obs["is_control"].to_numpy())[0])
        poisoned.loc[poisoned.index[0], "row"] = treated
        with self.assertRaises(contract.ContractError):
            consensus.assert_consensus_provenance(poisoned, self.obs, self.folds,
                                                  n_folds=N_FOLDS_TEST)

    def test_provenance_gate_rejects_unequal_plate_weights(self):
        poisoned = self.membership.copy()
        fold0 = poisoned["fold"] == 0
        poisoned.loc[fold0, "weight"] = 0.0
        poisoned.loc[poisoned.index[fold0][0], "weight"] = 1.0
        with self.assertRaises(contract.ContractError):
            consensus.assert_consensus_provenance(poisoned, self.obs, self.folds,
                                                  n_folds=N_FOLDS_TEST)

    def test_provenance_gate_rejects_qc_failed_member(self):
        obs = self.obs.copy()
        row_ix = int(self.membership.loc[self.membership["fold"] == 0, "row"].iloc[0])
        obs.iloc[row_ix, obs.columns.get_loc("qc_pass")] = False
        with self.assertRaises(contract.ContractError):
            consensus.assert_consensus_provenance(self.membership, obs, self.folds,
                                                  n_folds=N_FOLDS_TEST)

    def test_allowed_plates(self):
        self.assertEqual(consensus.allowed_plates(self.folds, 0), ["P1", "P2"])
        self.assertEqual(consensus.allowed_plates(self.folds, 1), ["P2", "P3"])

    def test_membership_columns_are_frozen(self):
        self.assertEqual(list(self.membership.columns), consensus.MEMBERSHIP_COLUMNS)


if __name__ == "__main__":
    unittest.main()
