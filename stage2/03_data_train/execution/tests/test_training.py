"""Public publication, replay and reader examples for plate-matched training."""
import gzip
import json
from hashlib import file_digest
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

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

from build_training import promote_training, publish_training, replay_training
from training_data import open_training


def digest(path):
    with path.open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


class TrainingTest(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory()))
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
        tahoe_obs = obs.copy()
        tahoe_obs["vehicle"] = np.where(
            obs.role.eq("vehicle"), obs.vehicle.map({"DMSO": "DMSO_TF", "Media": "DMSO"}), ""
        )
        self.counts = counts
        self.source_observations = {"tenv": obs, "tahoe_bulk": tahoe_obs}
        self.corpus = self.make_corpus("raw", self.source_observations)
        self.axis = PredictionAxis((5, 2), ("five", "two"), "a" * 64)
        self.output = self.root / "output"
        self.staging = self.root / "staging"
        self.fold_root = self.root / "folds"
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

    def make_corpus(self, name, source_observations):
        raw = self.root / name
        raw.mkdir()
        genes = pd.DataFrame(index=map(str, range(18164)))
        for source, observations in source_observations.items():
            AnnData(X=self.counts, obs=observations, var=genes).write_h5ad(raw / f"{source}.h5ad")
        native_columns = dict(enumerate(range(18164)))
        seal_corpus(raw, dict.fromkeys(source_observations, native_columns),
                    18164, COMMON_SHA256, {}, {})
        return self.enterContext(open_corpus(raw))

    def test_tenv_axis_release_has_three_keyed_products_and_replays_from_raw_counts(self):
        import h5py

        from build_tenv_training import publish_tenv, replay_tenv
        from training_data import open_tenv_training

        tenv = self.source_observations["tenv"].copy()
        tenv.loc[tenv.cell_line.eq("other"), "plate"] = "C"
        corpus = self.make_corpus("tenv-raw", {"tenv": tenv})
        args = (self.output, corpus, self.axis, self.axis, self.fold_root, self.receipt)
        publish_tenv(*args, staging=self.staging)
        with open_tenv_training(self.output) as data:
            self.assertEqual(data.treatments[["plate", "drug_parent", "dose_um"]].values.tolist(),
                             [["A", "P", 1], ["B", "P", 1], ["C", "Q", 2],
                              ["A", "P", 2], ["C", "Q", 1]])
            expected_controls = np.log1p([[7500, 2500], [0, 10000], [5000, 5000]]).astype(np.float32)
            expected_treated = np.log1p([[1000, 9000], [5000, 5000], [10000, 0]]).astype(np.float32)
            np.testing.assert_array_equal(data.read_controls(0, 3), expected_controls)
            np.testing.assert_array_equal(data.read_treatments(0, 3), expected_treated)
            np.testing.assert_array_equal(data.read_treatments(0, 3, "delta"),
                                          expected_treated - expected_controls)
            np.testing.assert_allclose(data.control(), expected_controls.mean(axis=0))
            self.assertEqual(data.treatments.cold_dose_drug_0.tolist(), ["test", "test", "valid", "", ""])
            self.assertEqual(data.treatments.cold_dose_drug_0_reason.tolist(),
                             ["", "", "", "held_parent", "held_dose"])
            with self.assertRaisesRegex(ValueError, "window"):
                data.read_treatments(0, 257)
        with h5py.File(self.output / "plate-ref-mtx.h5", "r") as handle:
            self.assertEqual(handle["plate"].asstr()[:].tolist(), ["A", "B", "C"])
        with h5py.File(self.output / "pert-mtx.h5", "r") as handle:
            self.assertEqual(handle["control_row"][:].tolist(), [0, 1, 2, 0, 2])
            self.assertEqual(handle["cold_dose_drug_0"].asstr()[:].tolist(),
                             ["test", "test", "valid", "", ""])
        self.assertTrue((self.output / "tenv-ref-vector.h5").is_file())
        self.assertEqual(replay_tenv(*args, staging=self.staging)["result"], "PASS")

    def test_tenv_axis_publication_requires_the_reference_population(self):
        from build_tenv_training import publish_tenv, replay_tenv
        from training_data import open_tenv_training

        tenv = self.source_observations["tenv"].copy()
        tenv.loc[tenv.cell_line.eq("other"), "plate"] = "C"
        self.counts[:2, 3] = [1, 2]
        corpus = self.make_corpus("tenv-raw", {"tenv": tenv})
        folds = (self.axis, self.fold_root, self.receipt)
        reference = self.root / "reference"
        publish_tenv(reference, corpus, self.axis, *folds, staging=self.staging)

        axes = (PredictionAxis((5, 2, 3), ("five", "two", "three"), "b" * 64),
                PredictionAxis((2,), ("two",), "c" * 64))
        args = self.output, corpus, axes[0], *folds
        publish_tenv(*args, staging=self.staging, reference=reference)
        self.assertEqual(replay_tenv(*args, staging=self.staging, reference=reference)["result"], "PASS")
        with open_tenv_training(self.output) as data:
            digest = data.manifest_sha256
        with self.assertRaisesRegex(ValueError, "treatment population or fold assignment"):
            publish_tenv(self.output, corpus, axes[1], *folds, staging=self.staging, reference=reference)
        with open_tenv_training(self.output) as data:
            self.assertEqual(data.manifest_sha256, digest)

    def freeze(self):
        rows = pd.DataFrame(self.fold_rows, columns=[
            "row", "profile_id", "drug_parent", "dose_um", "scenario", "fold", "partition",
        ])
        source_digests = {}
        for key in ("tenv", "tahoe_bulk"):
            source = self.fold_root / key
            source.mkdir(parents=True, exist_ok=True)
            rows.to_csv(source / "splits.csv", index=False)
            (source / "manifest.json").write_text(json.dumps({
                "outputs": {"splits.csv": digest(source / "splits.csv")},
            }))
            source_digests[key] = digest(source / "manifest.json")
        pair = self.fold_root / "manifest.json"
        pair.write_text(json.dumps({
            "schema": "tsg-dose-round1-axis-local-pair-v1",
            "view_generation": {"axis": {"sha256": self.axis.sha256}},
            "sources": source_digests,
        }))
        self.receipt.write_text(json.dumps({
            "pair_sha256": digest(pair),
            "input_source_evidence": {"tenv": {"split_seed": 226}, "tahoe_bulk": {"split_seed": 2026}},
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

    def test_pair_adds_tahoe_using_its_compatible_vehicle_and_frozen_strips(self):
        self.publish()
        tenv_digest = digest(self.output / "tenv/manifest.json")
        args = (self.output, self.corpus, self.axis, self.fold_root, self.receipt)
        publish_training(*args, staging=self.staging)
        self.assertEqual(digest(self.output / "tenv/manifest.json"), tenv_digest)
        with open_training(self.output, source="tahoe-bulk") as data:
            self.assertEqual(list(data.genes.gene_index), [5, 2])
            self.assertEqual(data.manifest["inputs"]["raw_source"], "tahoe_bulk")
            self.assertEqual(data.manifest["inputs"]["folds"]["split_seed"], 2026)
            self.assertEqual(list(zip(data.controls.plate, data.controls.cell_line)),
                             [("A", "CM"), ("B", "CM"), ("A", "other")])
            self.assertEqual(data.treatments.control_row.tolist(), [0, 1, 2, 0, 2])
            expected_controls = np.log1p([[7500, 2500], [0, 10000], [5000, 5000]]).astype(np.float32)
            expected_treatment = np.log1p([[1000, 9000], [5000, 5000], [10000, 0]]).astype(np.float32)
            np.testing.assert_array_equal(data.read_controls(0, 3), expected_controls)
            np.testing.assert_array_equal(data.read_treatments(0, 3), expected_treatment)
            np.testing.assert_array_equal(data.read_treatments(0, 3, "delta"),
                                          expected_treatment - expected_controls)
            np.testing.assert_allclose(data.control(), expected_controls.mean(axis=0))
            excluded = data.exclusions().set_index("profile_id")
            self.assertEqual(excluded.loc["media", "reason"], "incompatible_vehicle")
            strips = data.treatments.iloc[-2:]
            self.assertEqual(strips.cold_drug_dose_0.tolist(), ["excluded", "excluded"])
            self.assertEqual(strips.cold_drug_dose_0_reason.tolist(), ["held_parent", "held_dose"])
        result = replay_training(*args, staging=self.staging)
        self.assertEqual(result["sources"], ["tenv", "tahoe-bulk"])

    def test_verified_pair_can_be_promoted_with_the_public_api_intact(self):
        self.publish()
        api = self.output / "api"
        api.mkdir()
        reader = api / "reader.py"
        reader.write_text("PUBLISHED = True\n")
        candidate = self.root / "candidate"
        publish_training(candidate, self.corpus, self.axis, self.fold_root, self.receipt,
                         staging=self.staging)
        expected = digest(candidate / "manifest.json")

        promote_training(candidate, self.output, self.corpus, self.axis,
                             self.fold_root, self.receipt,
                         expected_sha256=expected, staging=self.staging)

        self.assertEqual(digest(self.output / "manifest.json"), expected)
        self.assertEqual(reader.read_text(), "PUBLISHED = True\n")
        candidate_counts = candidate / "tahoe-bulk/training.h5"
        accepted_counts = self.output / "tahoe-bulk/training.h5"
        self.assertFalse(candidate_counts.samefile(accepted_counts))
        self.assertEqual(digest(candidate_counts), digest(accepted_counts))
        with open_training(self.output, source="tahoe-bulk") as data:
            self.assertEqual(data.treatments.control_row.iloc[0], 0)

    def test_promotion_rejects_corrupt_payload_under_the_accepted_manifest(self):
        self.publish()
        before = digest(self.output / "manifest.json")
        candidate = self.root / "candidate"
        publish_training(candidate, self.corpus, self.axis, self.fold_root, self.receipt,
                         staging=self.staging)
        expected = digest(candidate / "manifest.json")
        with (candidate / "tahoe-bulk/training.h5").open("ab") as handle:
            handle.write(b"changed after verification")

        with self.assertRaisesRegex(ValueError, "training payload digest differs"):
            promote_training(candidate, self.output, self.corpus, self.axis,
                             self.fold_root, self.receipt,
                             expected_sha256=expected, staging=self.staging)

        self.assertEqual(digest(self.output / "manifest.json"), before)
        self.assertFalse((self.output / "tahoe-bulk").exists())

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

    def test_canonical_destination_requires_promotion(self):
        with patch("build_training.DEFAULT", self.output):
            with self.assertRaisesRegex(ValueError, "scratch destination"):
                self.publish()
        self.assertFalse(self.output.exists())

    def test_promotion_rejects_self_consistent_candidate_that_does_not_replay(self):
        self.publish()
        accepted = digest(self.output / "manifest.json")
        candidate = self.root / "candidate"
        publish_training(candidate, self.corpus, self.axis, self.fold_root, self.receipt,
                         staging=self.staging)
        source = candidate / "tenv"
        treatments = source / "treatments.csv"
        rows = pd.read_csv(treatments, keep_default_na=False)
        rows.loc[0, "cold_drug_0"] = "train"
        rows.to_csv(treatments, index=False)
        dataset_manifest = source / "manifest.json"
        record = json.loads(dataset_manifest.read_text())
        record["outputs"]["treatments.csv"] = digest(treatments)
        dataset_manifest.write_text(json.dumps(record))
        pair_manifest = candidate / "manifest.json"
        record = json.loads(pair_manifest.read_text())
        record["sources"]["tenv"] = digest(dataset_manifest)
        pair_manifest.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, "replay"):
            promote_training(candidate, self.output, self.corpus, self.axis,
                             self.fold_root, self.receipt, expected_sha256=digest(pair_manifest),
                             staging=self.staging)
        self.assertEqual(digest(self.output / "manifest.json"), accepted)

    def test_compressed_fold_archive_uses_original_decompressed_bytes(self):
        self.publish()
        accepted = digest(self.output / "manifest.json")
        path = self.fold_root / "tenv/splits.csv"
        original = path.read_bytes()
        archived = path.with_suffix(".csv.gz")
        archived.write_bytes(gzip.compress(original, mtime=0))
        path.unlink()
        self.publish()
        self.assertEqual(digest(self.output / "manifest.json"), accepted)
        archived.write_bytes(gzip.compress(original + b"modified", mtime=0))
        with self.assertRaisesRegex(ValueError, "fold digest differs"):
            self.publish()
        self.assertEqual(digest(self.output / "manifest.json"), accepted)
