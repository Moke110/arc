"""Public publication, replay and reader examples for plate-matched training."""
import gzip
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

from build_tenv_training import publish_tenv, replay_tenv
from tenv_axes import tenv_axis_output
from training_data import open_tenv_training


class TrainingAxisTest(unittest.TestCase):
    def test_named_tenv_release_uses_the_canonical_axis_id(self):
        step = Path("custom-training-root")
        self.assertEqual(tenv_axis_output("denoise-5k", step), step / "03.7_tenv_denoise-5k/output")
        with self.assertRaisesRegex(ValueError, "unknown tenv training axis"):
            tenv_axis_output("denoise-deg-5k", step)

    def test_both_2k_axes_have_public_training_release_paths(self):
        step = Path("custom-training-root")
        self.assertEqual(tenv_axis_output("hvg-2k", step), step / "03.8_tenv_hvg-2k/output")
        self.assertEqual(tenv_axis_output("denoise-2k", step), step / "03.9_tenv_denoise-2k/output")


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
        obs.loc[obs.cell_line.eq("other"), "plate"] = "C"
        self.counts = counts
        self.obs = obs
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
        self.corpus = self.make_corpus("raw", {"tenv": self.obs})

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

        corpus = self.make_corpus("tenv-raw", {"tenv": self.obs})
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
            members = data.membership().set_index("profile_id")
            self.assertEqual(members.loc["t2", "row"], 0)
            self.assertEqual(members.loc["t4", "row"], 2)
            self.assertEqual(data.exclusions().set_index("profile_id").reason.to_dict(), {
                "blank": "technical_blank", "bad_t": "unresolved_drug",
                "media": "incompatible_vehicle", "zero_t": "zero_treatment_counts",
                "zero_v": "zero_control_counts", "missing_control": "no_matched_vehicle",
                "blank_dose": "unresolved_dose", "zero_control_t": "zero_control_counts",
            })
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
        self.counts[:2, 3] = [1, 2]
        corpus = self.make_corpus("tenv-raw", {"tenv": self.obs})
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
        for key in ("tenv",):
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
            "input_source_evidence": {"tenv": {"split_seed": 226}},
        }))

    def publish(self):
        return publish_tenv(self.output, self.corpus, self.axis, self.axis, self.fold_root, self.receipt,
                            staging=self.staging)

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
        with open_tenv_training(self.output) as data:
            self.assertEqual(data.treatments.control_row.iloc[0], 0)

    def test_disagreeing_member_assignments_fail_publication(self):
        row = next(row for row in self.fold_rows
                   if (row[1], row[4], row[5]) == ("t2", "cold-drug", 1))
        row[-1] = "validation"
        self.freeze()
        with self.assertRaisesRegex(ValueError, "assignments disagree"):
            self.publish()
        self.assertFalse(self.output.exists())

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
