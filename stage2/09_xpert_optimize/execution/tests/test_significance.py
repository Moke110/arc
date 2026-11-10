"""Bounded real-DESeq2 check through the evaluator release CLI."""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


SCRIPT = Path(__file__).resolve().parents[1] / "script" / "build_significance.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_manifest(directory, **metadata):
    files = {str(path.relative_to(directory)): sha(path)
             for path in directory.rglob("*.csv")}
    path = directory / "manifest.json"
    path.write_text(json.dumps({**metadata, "files": files}))
    return sha(path)


def write_fixture(root):
    genes = pd.DataFrame({"hgnc_id": [f"H{i}" for i in range(64)],
                          "symbol": [f"G{i}" for i in range(63)] + ["HDGFL3"],
                          "ensembl_id": [f"E{i}" for i in range(64)]})
    obs, counts, cases, controls = [], [], [], []
    rng = np.random.default_rng(17)
    mean = rng.uniform(50, 300, len(genes))
    for plate in range(3):
        for treatment, count in (("", 4), ("a", 4), ("b", 2)) if plate < 2 else (("", 4),):
            for _ in range(count):
                idx = len(obs)
                obs.append({"case_id": f"well-{idx}", "plate": f"P{plate}",
                            "drug_key": treatment, "dose_um": 1.0, "is_blank": False,
                            "is_control": not treatment, "compound_name": treatment or "DMSO",
                            "vehicle": "DMSO"})
                value = rng.poisson(rng.gamma(35, mean / 35) * (plate + 1))
                value[2] = 0
                if treatment == "a":
                    value[0] *= 10
                    value[-1] //= 10
                elif treatment == "b":
                    value[:] = 0
                counts.append(value)
                if treatment:
                    cases.append({**{key: obs[-1][key] for key in ("case_id", "drug_key", "dose_um", "plate")},
                                  "group_id": treatment, "control_index": plate, "source_row": idx})
                else:
                    controls.append({"control_index": plate, "source_row": idx})

    native = root / "tenv.h5ad"
    ad.AnnData(csr_matrix(np.asarray(counts, dtype=np.int32)), obs=pd.DataFrame(obs).set_index("case_id"),
               var=genes.set_index("symbol", drop=False)).write_h5ad(native)
    axes = root / "axes"
    predictor = axes / "views/G18164/predictor"
    evaluator = axes / "views/G18164/evaluator"
    landmark = root / "landmark/predictor"
    (predictor / "tenv/fold_0").mkdir(parents=True)
    evaluator.mkdir()
    landmark.mkdir(parents=True)
    genes.iloc[:-1].to_csv(predictor / "genes.csv", index=False)
    genes.iloc[[-1]].to_csv(landmark / "genes.csv", index=False)
    pd.DataFrame(cases).to_csv(predictor / "tenv/fold_0/validation.csv", index=False)
    pd.DataFrame({"control_index": range(3), "plate": ["P0", "P1", "P2"],
                  "n_wells": [4, 4, 4]}).to_csv(predictor / "tenv/controls.csv", index=False)
    pd.DataFrame(controls).to_csv(evaluator / "tenv_control_rows.csv", index=False)
    pdigest = write_manifest(predictor, schema="step09-predictor-v1", axis="G18164",
                             ready_folds=[0], gene_count=63, fixture=True)
    edigest = write_manifest(evaluator, schema="step09-evaluator-v1", predictor_manifest_sha256=pdigest)
    ldigest = write_manifest(landmark, schema="step09-predictor-v1", axis="landmark", fixture=True)
    write_manifest(axes, schema="step09-gene-axes-v1", ready_folds=[0], fixture=True,
                   source_hashes={"tenv.h5ad": sha(native)}, landmark_predictor_sha256=ldigest,
                   views={"G18164": {"predictor": pdigest, "evaluator": edigest}})
    return native, axes, landmark.parent


class SignificanceTest(unittest.TestCase):
    def test_all_contrasts_failed_systemically_does_not_publish_release(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            native, axes, landmark = write_fixture(root)
            worker = root / "failed-rscript"
            worker.write_text(f"#!{sys.executable}\n"
                              "import pathlib, sys\n"
                              "path = pathlib.Path(sys.argv[3]) / 'results.csv'\n"
                              "path.write_text('status,reason\\nuntestable,DESeq2: configuration error\\n'"
                              " 'untestable,zero_library\\n')\n")
            worker.chmod(0o755)
            output = root / "validation"
            process = subprocess.run(
                [sys.executable, SCRIPT, "--role", "validation", "--native", native,
                 "--gene-axes", axes, "--landmark", landmark, "--rscript", worker, "--output", output],
                text=True, capture_output=True,
            )
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("DESeq2 produced no eligible contrast", process.stderr)
            self.assertFalse(output.exists())

    def test_multiplate_wald_contrast_and_untestable_group(self):
        rscript = os.environ.get("STEP09_RSCRIPT")
        if not rscript:
            rscript = shutil.which("Rscript")
            if not rscript or subprocess.run(
                [rscript, "--vanilla", "-e", 'quit(status=!requireNamespace("DESeq2", quietly=TRUE))'],
                capture_output=True,
            ).returncode:
                self.skipTest("set STEP09_RSCRIPT to an Rscript with DESeq2")
        with TemporaryDirectory() as directory:
            root = Path(directory)
            native, axes, landmark = write_fixture(root)
            output = root / "validation"
            process = subprocess.run(
                [sys.executable, SCRIPT, "--role", "validation", "--native", native,
                 "--gene-axes", axes, "--landmark", landmark, "--rscript", rscript, "--output", output],
                text=True, capture_output=True,
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual((manifest["schema"], manifest["role"], manifest["fold"], manifest["fdr"]),
                             ("step09-significance-v1", "validation", 0, 0.05))
            for name, digest in manifest["files"].items():
                self.assertEqual(sha(output / name), digest)
            genes = pd.read_csv(output / "genes.csv")
            self.assertEqual(genes.symbol.iloc[-1], "HDGFL3")
            groups = pd.read_csv(output / "groups.csv").set_index("group_id")
            self.assertEqual(groups.index.tolist(), ["a", "b"])
            self.assertEqual(groups.loc["a", "n_controls"], 8)
            self.assertEqual(groups.loc["a", "n_plates"], 2)
            self.assertEqual(groups.loc["a", "design"], "~ plate + condition")
            self.assertEqual(groups.loc["a", "status"], "eligible")
            self.assertEqual(groups.loc["b", "status"], "untestable")
            samples = pd.read_csv(output / "contrast_samples.csv")
            self.assertEqual(set(samples.plate), {"P0", "P1"})
            pvalue = np.load(output / "pvalue.npy")
            padj = np.load(output / "padj.npy")
            effect = np.load(output / "log2_fold_change.npy")
            self.assertEqual(padj.shape, (2, 64))
            self.assertTrue(np.isnan(padj[1]).all() and np.isnan(padj[0, 2]))
            self.assertTrue(padj[0, 0] < .05 and padj[0, -1] < .05)
            self.assertTrue(effect[0, 0] > 2 and effect[0, -1] < -2)
            # The most significant gene has rank 1 and adjustment uses the
            # complete 64-gene family, including the unexpressed gene.
            idx = np.nanargmin(pvalue[0])
            self.assertAlmostEqual(padj[0, idx] / pvalue[0, idx], 64.0)

    def test_test_role_requires_committed_promotion_before_reading_native_counts(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "test"
            process = subprocess.run(
                [sys.executable, SCRIPT, "--role", "test", "--output", output,
                 "--native", Path(directory) / "missing.h5ad"],
                text=True, capture_output=True,
            )
            self.assertNotEqual(process.returncode, 0)
            self.assertIn("test labels require a committed validation promotion", process.stderr)
            self.assertFalse(output.exists())
