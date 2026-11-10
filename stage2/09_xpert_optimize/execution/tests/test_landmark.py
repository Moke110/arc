"""CLI seam for the Step 09 prepare, fit, predict and replay contract."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


EXECUTION = Path(__file__).resolve().parents[1]
REPO = next(path for path in EXECUTION.parents if (path / ".git").exists())
sys.path.insert(0, str(REPO / "stage2/03_modes_tenv/execution/script"))
from inputs import split_drugs


def run(relative, *arguments, check=True):
    return subprocess.run(
        [sys.executable, str(EXECUTION / relative), *map(str, arguments)],
        check=check, capture_output=True, text=True,
        env={**os.environ, "OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "2",
             "OPENBLAS_NUM_THREADS": "2", "NUMEXPR_NUM_THREADS": "2"},
    )


def make_fixture(root):
    native, tenv_input = root / "native", root / "frozen_tenv"
    (native / "tables").mkdir(parents=True)
    tenv_input.mkdir()
    axis = pd.DataFrame({
        "symbol": ["A", "B", "C", "HDGFL3"], "hgnc_id": ["H1", "H2", "H3", "H4"],
        "entrez_id": [1, 2, 3, 50810], "ensembl_id": ["E1", "E2", "E3", "E4"],
    })
    axis.iloc[:3].to_csv(native / "tables/genes.csv", index=False)
    landmarks = pd.DataFrame({"gene_idx": [0, 1, 2], "gene_id": [2, 50810, 1], "ensembl_id": ["E2", "E4", "E1"]})
    landmarks.to_csv(root / "landmarks.csv", index=False)

    drugs = pd.DataFrame({"drug_key": list("abcdefgh"), "parent_connectivity": ["LFQSCWFLJHTTHZ", *"RSTUVWX"],
                          "morgan_index": range(8)})
    drugs.to_csv(tenv_input / "drugs.csv", index=False)
    folds = split_drugs(drugs)
    folds.to_csv(tenv_input / "folds.csv", index=False)
    np.save(tenv_input / "morgan.npy", np.eye(8, 2048, dtype=np.uint8))
    tenv_obs = pd.DataFrame([
        {"well": f"tenv:{i}", "drug_key": drug, "dose_um": 1.0, "plate": "P1",
         "time_h": 24.0, "is_control": i >= 8, "is_blank": False,
         "compound_name": "DMSO" if i >= 8 else drug,
         "vehicle": "DMSO"}
        for i, drug in enumerate([*drugs.drug_key, "", ""])
    ]).set_index("well", drop=False)
    wells = pd.DataFrame({"case_id": tenv_obs.index[:8], "drug_key": drugs.drug_key,
                          "dose_um": 1.0, "control_index": 0, "group_id": list("abcdefgh")})
    wells.to_csv(tenv_input / "wells.csv", index=False)
    pd.DataFrame({"control_index": [0], "plate": ["P1"], "n_wells": [2]}).to_csv(
        tenv_input / "controls.csv", index=False)
    for fold in range(3):
        target = tenv_input / f"fold_{fold}"
        target.mkdir()
        for role, name in [("train", "train"), ("validation", "validation"), ("test", "test_cases")]:
            selected = folds.query("fold == @fold and role == @role").drug_key
            wells[wells.drug_key.isin(selected)].to_csv(target / f"{name}.csv", index=False)
    ad.AnnData(X=np.array([[2, 3, 1000, 5]] * 8 + [[0, 0, 99, 0], [2, 3, 3, 5]], dtype=np.int32),
               obs=tenv_obs, var=axis.set_index("symbol")).write_h5ad(native / "tenv.h5ad")

    tahoe_obs = pd.DataFrame([
        {"profile_id": f"{plate}:{well}:L1", "plate": plate, "well": well, "cell_line_id": "L1",
         "drug": drug, "dose_um": 1.0, "is_control": drug == "DMSO_TF"}
        for plate in ["P1", "P2", "P3"] for well, drug in [("A1", "a"), ("A2", "DMSO_TF")]
    ]).set_index("profile_id")
    tahoe_var = axis[["ensembl_id"]].assign(token_id=list("1234")).set_index("ensembl_id", drop=False)
    ad.AnnData(X=csr_matrix(np.array([[1, 2, 100, 7], [2, 3, 9, 5]] * 3, dtype=np.int32)),
               obs=tahoe_obs, var=tahoe_var).write_h5ad(native / "tahoe.h5ad")
    pd.DataFrame({"drug": ["a"], "canonical_smiles": ["CCO"]}).to_parquet(root / "tahoe_drugs.parquet")
    crosswalk = [
        {"dataset_id": source, "source_gene_id": native_id, "symbol": row.symbol, "hgnc_id": row.hgnc_id}
        for i, row in enumerate(axis.itertuples(index=False), 1)
        for source, native_id in [("tahoe", row.ensembl_id), ("tahoe", str(i)), ("temp_env", row.ensembl_id)]
    ]
    pd.DataFrame(crosswalk).to_csv(native / "tables/gene_crosswalk.csv", index=False)
    return native, tenv_input


class LandmarkTest(unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        tmp_path = self.tmp_path = Path(temp.name)
        self.native, self.tenv = make_fixture(tmp_path)
        self.release = tmp_path / "tahoe_tenv_bulk_l1000"
        self.prepare(self.release)

    def prepare(self, target):
        tmp_path = self.tmp_path
        run("01_corpus/build.py", "--native", self.native, "--tenv-release", self.tenv,
            "--landmarks", tmp_path / "landmarks.csv", "--tahoe-drugs", tmp_path / "tahoe_drugs.parquet",
            "--output", target, "--fixture")

    def test_prepare_preserves_folds_and_normalizes_after_projection(self):
        release = self.release
        predictor = release / "predictor"
        assert pd.read_csv(predictor / "genes.csv").symbol.tolist() == ["B", "HDGFL3", "A"]
        assert not list(predictor.rglob("test_expression.npy"))
        test_parents = set()
        for fold in range(3):
            root = predictor / "tenv" / f"fold_{fold}"
            roles = [pd.read_csv(root / f"{role}.csv") for role in ("train", "validation", "test")]
            parents = [set(role.parent_connectivity) for role in roles]
            assert sum(map(len, parents)) == len(set.union(*parents)) == 8
            assert "source_row" not in roles[-1]
            test_parents |= parents[-1]
            cases = pd.read_csv(root / "train.csv")
            assert len(cases)
            np.testing.assert_allclose(np.load(root / "train_expression.npy")[0],
                                       [8.00670085, 8.51739317, 7.60140233], rtol=1e-6)
        control = np.load(predictor / "tenv/controls.npy")[0]
        np.testing.assert_allclose(control, np.array([8.00670085, 8.51739317, 7.60140233]) / 2, rtol=1e-6)
        train = pd.read_csv(predictor / "tahoe/train.csv")
        valid = pd.read_csv(predictor / "tahoe/validation.csv")
        assert set(train.plate).isdisjoint(valid.plate)
        assert all(set(part.parent_connectivity) & test_parents for part in (train, valid))
        manifest = json.loads((predictor / "manifest.json").read_text())
        assert manifest["zero_libraries"]["tenv"]["controls"] == 1
        assert manifest["identity_seed"] not in manifest["training_seeds"]
        vectors = np.load(predictor / "identity.npy")
        assert vectors.shape == (3, 256) and len(np.unique(vectors, axis=0)) == 3
        assert np.load(predictor / "identity_hdgfl3.npy").shape == (256,)
        regenerated = self.tmp_path / "regenerated"
        self.prepare(regenerated)
        for side in ("predictor", "evaluator"):
            assert (release / side / "manifest.json").read_bytes() == (regenerated / side / "manifest.json").read_bytes()

    def fixture_policy(self):
        policy = self.tmp_path / "policy.json"
        policy.write_text(json.dumps({
            "kind": "qualification", "device": "cpu", "kernel": "math", "dtype": "float32",
            "microbatch": 1, "effective_batch": 1,
            "phases": {"tahoe": {"updates": 1, "validation_interval": 1, "learning_rate": 0.0001},
                       "tenv": {"updates": 1, "validation_interval": 1, "learning_rate": 0.0001}},
        }))
        return policy

    def test_two_stage_checkpoint_replay_uses_only_predictor(self):
        release, tmp_path = self.release, self.tmp_path
        predictor = release / "predictor"
        (release / "evaluator").rename(tmp_path / "unmounted_evaluator")
        policy = self.fixture_policy()
        target = tmp_path / "fit"
        run("script/xpert.py", "fit", "--predictor", predictor, "--policy", policy,
            "--fold", 0, "--seed", 17, "--output", target)
        report = json.loads((target / "run.json").read_text())
        for stage in ("tahoe", "tenv"):
            assert report["phases"][stage]["best_update"] == 1
            assert np.isfinite(report["phases"][stage]["history"][0]["validation_terms"]).all()
            assert len(report["phases"][stage]["history"][0]["validation_terms"]) == 4
        first = tmp_path / "prediction"
        run("script/xpert.py", "predict", "--predictor", predictor,
            "--checkpoint", target / "tenv.pt", "--role", "test", "--output", first)
        run("script/xpert.py", "replay", "--predictor", predictor,
            "--checkpoint", target / "tenv.pt", "--prediction", first)
        assert np.isfinite(np.load(first / "delta.npy")).all()
        assert np.load(first / "delta.npy").shape[1] == 3

    def test_fit_rejects_redirect_to_evaluator(self):
        predictor = self.release / "predictor"
        private = self.release / "evaluator/tenv/fold_0/test_expression.npy"
        relative = "tenv/fold_0/train_expression.npy"
        (predictor / relative).unlink()
        (predictor / relative).symlink_to(private)
        manifest_path = predictor / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["files"][relative] = hashlib.sha256(private.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest))
        rejected = run("script/xpert.py", "fit", "--predictor", predictor,
                       "--policy", self.fixture_policy(), "--fold", 0, "--seed", 17,
                       "--output", self.tmp_path / "rejected", check=False)
        assert rejected.returncode != 0
        assert f"predictor path escaped or is a symlink: {relative}" in rejected.stderr
