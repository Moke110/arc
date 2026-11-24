"""Public publication, replay and reader examples for plate-matched training."""
import json
from hashlib import file_digest
from pathlib import Path
import sys
import tempfile
import unittest

from anndata import AnnData
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[2]
STAGE = STEP.parent
sys.path[:0] = [str(path) for path in (
    STEP / "execution", STEP / "output/api",
    STAGE / "01_data_raw/output/api", STAGE / "02_gene_axis/02.1_hvg-deg/output/api",
)]

from native_corpus import open_corpus, seal_corpus
from prediction_axis import COMMON_SHA256, PredictionAxis

from build_training import publish_training, replay_training
from training_data import open_training


def digest(path):
    with path.open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


class TrainingTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
        raw = self.root / "raw"
        raw.mkdir()
        ids = ["t1", "t2", "v1", "v2", "t3", "v3", "t4", "v4", "blank",
               "bad_t", "media", "zero_t", "zero_v", "missing_control", "blank_dose",
               "zero_control_t", "held_parent", "held_dose"]
        role = ["treated", "treated", "vehicle", "vehicle", "treated", "vehicle",
                "treated", "vehicle", "blank", "treated", "vehicle", "treated",
                "vehicle", "treated", "treated", "treated", "treated", "treated"]
        plate = ["A", "A", "A", "A", "B", "B", "A", "A", "A", "A",
                 "A", "A", "Z", "D", "A", "Z", "A", "A"]
        cell = ["CM"] * len(ids)
        cell[6:8] = ["other", "other"]
        cell[-1] = "other"
        group = ["A-CM", "A-CM", "A-CM", "A-CM", "B-CM", "B-CM", "A-other", "A-other",
                 "", "A-CM", "A-Media", "A-CM", "Z-CM", "", "A-CM", "Z-CM", "A-CM", "A-other"]
        parent = ["P", "P", "", "", "P", "", "Q", "", "", "", "", "Z", "", "R", "S", "S", "P", "Q"]
        dose = [1, 1, 0, 0, 1, 0, 2, 0, 0, 1, 0, 3, 0, 2, np.nan, 2, 2, 1]
        obs = pd.DataFrame({
            "role": role, "matched_control_group": group, "plate": plate, "cell_line": cell,
            "drug_parent_inchikey": parent, "dose_um": dose,
            "vehicle": ["Media" if name == "media" else "DMSO" for name in ids],
            "view_exclusion_reason": ["technical_blank" if name == "blank" else "" for name in ids],
        }, index=pd.Index(ids, name="profile_id"))
        counts = np.zeros((len(ids), 18164), dtype=np.int32)
        counts[:, 1] = 100000
        counts[:, [5, 2]] = [[1, 0], [0, 9], [3, 0], [0, 1], [2, 2], [0, 8],
                             [4, 0], [1, 1], [0, 0], [1, 1], [99, 0], [0, 0],
                             [0, 0], [2, 2], [2, 2], [2, 2], [1, 3], [99, 1]]
        AnnData(X=counts, obs=obs, var=pd.DataFrame(index=map(str, range(18164)))).write_h5ad(
            raw / "tenv.h5ad"
        )
        seal_corpus(raw, {"tenv": dict(enumerate(range(18164)))}, 18164, COMMON_SHA256, {}, {})
        self.corpus = self.enterContext(open_corpus(raw))
        self.axis = PredictionAxis((5, 2), ("five", "two"), "a" * 64)
        self.output = self.root / "output"
        self.staging = self.root / "staging"
        self.fold_root = self.root / "folds"
        self.fold_root.joinpath("tenv").mkdir(parents=True)
        self.receipt = self.root / "receipt.json"
        self.fold_rows = []
        for scenario in ("cold-drug", "cold-drug-dose"):
            for fold in range(3):
                held_dose = (1, 2, 3)[fold]
                held_parent = ("P", "Q", "Z")[fold]
                for row, (profile, drug, dose_um) in enumerate(
                    [("t1", "P", 1), ("t2", "P", 1), ("t3", "P", 1),
                     ("t4", "Q", 2), ("zero_t", "Z", 3),
                     ("zero_control_t", "S", 2), ("held_parent", "P", 2),
                     ("held_dose", "Q", 1)]
                ):
                    if drug == held_parent and (scenario == "cold-drug" or dose_um == held_dose):
                        label = "test"
                    elif scenario == "cold-drug-dose" and drug == held_parent:
                        label = "excluded-held-parent"
                    elif scenario == "cold-drug-dose" and dose_um == held_dose:
                        label = "excluded-held-dose"
                    else:
                        label = "train" if drug == "P" else "validation"
                    self.fold_rows.append([row, profile, drug, dose_um, scenario, fold, label])
        self.freeze()

    def freeze(self):
        source = self.fold_root / "tenv"
        pd.DataFrame(self.fold_rows, columns=[
            "row", "profile_id", "drug_parent", "dose_um", "scenario", "fold", "partition",
        ]).to_csv(source / "splits.csv", index=False)
        (source / "manifest.json").write_text(json.dumps({
            "outputs": {"splits.csv": digest(source / "splits.csv")},
        }))
        pair = self.fold_root / "manifest.json"
        pair.write_text(json.dumps({
            "schema": "tsg-dose-round1-axis-local-pair-v1",
            "view_generation": {"axis": {"sha256": self.axis.sha256}},
            "sources": {"tenv": digest(source / "manifest.json")},
        }))
        self.receipt.write_text(json.dumps({
            "pair_sha256": digest(pair),
            "input_source_evidence": {"tenv": {"split_seed": 226}},
        }))

    def publish(self):
        return publish_training(self.output, self.corpus, self.axis, self.fold_root, self.receipt,
                                sources=("tenv",), staging=self.staging)

    def test_plate_matched_integer_pooling_and_bounded_replay(self):
        self.publish()
        with open_training(self.output, source="tenv") as data:
            self.assertEqual(list(data.genes.gene_index), [5, 2])
            self.assertEqual(list(zip(data.controls.plate, data.controls.cell_line)),
                             [("A", "CM"), ("B", "CM"), ("A", "other")])
            self.assertEqual(list(data.treatments.plate), ["A", "B", "A", "A", "A"])
            self.assertEqual(list(data.treatments.control_row), [0, 1, 2, 0, 2])
            expected_controls = np.log1p([[7500, 2500], [0, 10000], [5000, 5000]])
            expected_treated = np.log1p([[1000, 9000], [5000, 5000], [10000, 0]])
            np.testing.assert_allclose(data.read_controls(0, 3), expected_controls)
            np.testing.assert_allclose(data.read_treatments(0, 3), expected_treated)
            np.testing.assert_allclose(data.read_treatments(0, 3, "delta"),
                                       expected_treated - expected_controls, rtol=1e-6)
            np.testing.assert_allclose(data.control(), expected_controls.mean(axis=0))
            self.assertEqual(data.treatments.cold_drug_0.iloc[:3].tolist(), ["test", "test", "valid"])
            self.assertEqual(data.treatments.cold_drug_dose_1.iloc[:3].tolist(), ["train", "train", "test"])
            self.assertEqual(data.manifest["inputs"]["folds"]["held_dose_um"], [1, 2, 3])
            memberships = data.membership().set_index("profile_id")
            self.assertEqual(memberships.loc["v1", "row"], 0)
            self.assertEqual(memberships.loc["t2", "row"], 0)
            self.assertEqual(memberships.loc["t4", "row"], 2)
            self.assertEqual(data.exclusions().set_index("profile_id").reason.to_dict(), {
                "blank": "technical_blank", "bad_t": "unresolved_drug",
                "media": "incompatible_vehicle", "zero_t": "zero_treatment_counts",
                "zero_v": "zero_control_counts", "missing_control": "no_matched_vehicle",
                "blank_dose": "unresolved_dose", "zero_control_t": "zero_control_counts",
            })
            with self.assertRaisesRegex(ValueError, "window"):
                data.read_treatments(0, 257)
        result = replay_training(self.output, self.corpus, self.axis, self.fold_root, self.receipt,
                                 sources=("tenv",), staging=self.staging)
        self.assertEqual(result["result"], "PASS")


    def test_both_excluded_strips_remain_readable_in_cold_drug_dose(self):
        self.publish()
        with open_training(self.output, source="tenv") as data:
            rows = data.treatments.iloc[-2:]
            self.assertEqual(rows.cold_drug_dose_0.tolist(), ["excluded", "excluded"])
            self.assertEqual(rows.cold_drug_dose_0_reason.tolist(), ["held_parent", "held_dose"])
            np.testing.assert_allclose(data.read_treatments(3, 5),
                                       np.log1p([[2500, 7500], [9900, 100]]))

    def test_failed_publication_keeps_previous_accepted_bytes(self):
        self.publish()
        accepted = digest(self.output / "manifest.json")
        self.fold_rows = [
            row for row in self.fold_rows if (row[1], row[4], row[5]) != ("t1", "cold-drug", 1)
        ]
        self.freeze()
        with self.assertRaisesRegex(ValueError, "missing frozen fold assignment"):
            self.publish()
        self.assertEqual(digest(self.output / "manifest.json"), accepted)
        with open_training(self.output, source="tenv") as data:
            self.assertEqual(data.treatments.control_row.iloc[0], 0)

    def test_disagreeing_member_assignments_fail_publication(self):
        row = next(row for row in self.fold_rows
                   if (row[1], row[4], row[5]) == ("t2", "cold-drug", 1))
        row[-1] = "validation"
        self.freeze()
        with self.assertRaisesRegex(ValueError, "assignments disagree"):
            self.publish()
        self.assertFalse(self.output.exists())
