"""Gene-axis release and fit/replay contract on a bounded native-count fixture."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import anndata as ad
import numpy as np
import pandas as pd

from test_landmark import make_fixture, run


AXES = ("G2000", "G5000", "G10000", "G18164")


def build(root, native, tenv):
    landmark, axes = root / "landmark", root / "gene_axes"
    run("01_corpus/build.py", "--native", native, "--tenv-release", tenv,
        "--landmarks", root / "landmarks.csv", "--tahoe-drugs", root / "tahoe_drugs.parquet",
        "--output", landmark, "--fixture")
    run("01_corpus/gene_axes.py", "--native", native, "--landmark", landmark,
        "--output", axes, "--fixture-sizes", 1, 2, 3)
    return axes, landmark


def outside_train(native, tenv):
    tenv_fold = pd.read_csv(tenv / "folds.csv").query("fold == 0").set_index("drug_key")
    for source in ("tahoe", "tenv"):
        path = native / f"{source}.h5ad"
        data = ad.read_h5ad(path)
        counts = data.X if source == "tenv" else data.X.toarray()
        if source == "tahoe":
            validation = min(data.obs.plate.unique(), key=lambda p: hashlib.sha256(f"20261110|{p}".encode()).digest())
            outside = data.obs.plate.eq(validation) & ~data.obs.is_control
        else:
            outside = data.obs.drug_key.map(tenv_fold.role).isin(("validation", "test"))
        changed = np.flatnonzero(outside)
        assert len(changed)
        counts[changed, 0] = 4444
        counts[changed, 1] = 2 + 10 * np.arange(len(changed))
        data.X = counts
        data.write_h5ad(path)


def conflicting_dispersion(native, tenv):
    genes_path = native / "tables/genes.csv"
    cm = pd.read_csv(genes_path)
    gene = {"symbol": "Z", "hgnc_id": "H5", "entrez_id": 5, "ensembl_id": "E5"}
    pd.concat([cm, pd.DataFrame([gene])], ignore_index=True).to_csv(genes_path, index=False)
    crosswalk_path = native / "tables/gene_crosswalk.csv"
    crosswalk = pd.read_csv(crosswalk_path)
    extra = pd.DataFrame([{"dataset_id": source, "source_gene_id": name, "symbol": "Z", "hgnc_id": "H5"}
                          for source, name in (("tahoe", "E5"), ("tahoe", "5"), ("temp_env", "E5"))])
    pd.concat([crosswalk, extra], ignore_index=True).to_csv(crosswalk_path, index=False)
    tenv_fold = pd.read_csv(tenv / "folds.csv").query("fold == 0").set_index("drug_key")
    for source in ("tahoe", "tenv"):
        path = native / f"{source}.h5ad"
        data = ad.read_h5ad(path)
        counts = data.X if source == "tenv" else data.X.toarray()
        counts = np.column_stack((counts, np.full(len(counts), 40, dtype=np.int32)))
        if source == "tahoe":
            validation = min(data.obs.plate.unique(), key=lambda p: hashlib.sha256(f"20261110|{p}".encode()).digest())
            train = data.obs.plate.ne(validation) & ~data.obs.is_control
            patterns = ([8, 22, 30], [32, 18, 10])
            added = pd.DataFrame([{"ensembl_id": "E5", "token_id": "5"}], index=["E5"])
        else:
            train = data.obs.drug_key.map(tenv_fold.role).eq("train")
            patterns = ([22, 10, 28], [18, 30, 12])
            added = pd.DataFrame([gene]).set_index("symbol")
        rows = np.flatnonzero(train)
        counts[rows, :3] = [patterns[i % 2] if i < len(rows) - len(rows) % 2 else [20, 20, 20]
                            for i in range(len(rows))]
        var = pd.concat([data.var, added], ignore_index=False)
        var.index.name = data.var.index.name
        ad.AnnData(X=counts, obs=data.obs, var=var).write_h5ad(path)


class GeneAxisTest(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_nested_panels_keep_cohort_and_ignore_held_out_expression(self):
        original = self.root / "original"
        native, tenv = make_fixture(original)
        axes, landmark = build(original, native, tenv)
        variant = self.root / "changed"
        native2, tenv2 = make_fixture(variant)
        outside_train(native2, tenv2)
        changed, _ = build(variant, native2, tenv2)
        for name in ("ranks", *AXES, "shared"):
            assert (axes / f"fold_0/{name}.csv").read_bytes() == (changed / f"fold_0/{name}.csv").read_bytes()
        for fold in range(3):
            panels = [set(pd.read_csv(axes / f"fold_{fold}/{name}.csv").hgnc_id) for name in AXES]
            assert [len(panel) for panel in panels] == [1, 2, 3, 3]
            assert all(a <= b for a, b in zip(panels, panels[1:]))
        shared = pd.read_csv(axes / "fold_0/shared.csv")
        assert shared.symbol.tolist() == ["A"]
        release = json.loads((axes / "manifest.json").read_text())
        assert release["folds"]["0"]["shared_continuous_provisional"]

        base = landmark / "predictor"
        for name in AXES:
            model = axes / f"views/{name}/predictor"
            for path in ("tahoe/train.csv", "tahoe/validation.csv", "tenv/fold_0/train.csv",
                         "tenv/fold_0/validation.csv", "tenv/fold_0/test.csv"):
                assert (model / path).read_bytes() == (base / path).read_bytes()
            assert not list(model.rglob("test_expression.npy"))
        small = np.load(axes / "views/G2000/predictor/tenv/fold_0/train_expression.npy")[0]
        full = np.load(axes / "views/G18164/predictor/tenv/fold_0/train_expression.npy")[0]
        np.testing.assert_allclose(small, [9.21044037], rtol=1e-6)
        np.testing.assert_allclose(full, [3.03977296, 3.42916094, 9.20545333], rtol=1e-6)

    def test_full_axis_fit_predict_and_checkpoint_replay(self):
        native, tenv = make_fixture(self.root)
        axes, _ = build(self.root, native, tenv)
        predictor = axes / "views/G18164/predictor"
        policy = self.root / "policy.json"
        policy.write_text(json.dumps({
            "kind": "qualification", "device": "cpu", "kernel": "math", "dtype": "float32",
            "microbatch": 1, "effective_batch": 2,
            "phases": {"tahoe": {"updates": 1, "validation_interval": 1, "learning_rate": .0001},
                       "tenv": {"updates": 1, "validation_interval": 1, "learning_rate": .0001}},
        }))
        fit, prediction = self.root / "fit", self.root / "prediction"
        run("script/xpert.py", "fit", "--predictor", predictor, "--policy", policy,
            "--fold", 0, "--seed", 17, "--output", fit)
        run("script/xpert.py", "predict", "--predictor", predictor,
            "--checkpoint", fit / "tenv.pt", "--role", "validation", "--output", prediction)
        run("script/xpert.py", "replay", "--predictor", predictor,
            "--checkpoint", fit / "tenv.pt", "--prediction", prediction)
        report = json.loads((fit / "run.json").read_text())
        assert report["model"] == "XPert_G18164"
        assert np.load(prediction / "delta.npy").shape[1] == 3

    def test_equal_rank_weighting_breaks_combined_ties_by_cm_order(self):
        native, tenv = make_fixture(self.root)
        conflicting_dispersion(native, tenv)
        axes, _ = build(self.root, native, tenv)
        ranks = pd.read_csv(axes / "fold_0/ranks.csv").set_index("symbol")
        assert ranks.loc[["A", "B", "C", "Z"], "tahoe_rank"].tolist() == [1, 3, 2, 4]
        assert ranks.loc[["A", "B", "C", "Z"], "tenv_rank"].tolist() == [3, 1, 2, 4]
        assert ranks.loc[["A", "B", "C"], "rank_sum"].tolist() == [4, 4, 4]
        assert pd.read_csv(axes / "fold_0/G2000.csv").symbol.tolist() == ["A"]
