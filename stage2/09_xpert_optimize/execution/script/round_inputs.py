"""Resolve the sealed inputs and fit evidence used by Round 1 evaluation."""

import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from round_contract import checked_payload, read_json, role_result, sha256, verify_committed


STEP = Path(__file__).resolve().parents[2]
CASE_COLUMNS = ["case_id", "group_id", "drug_key", "dose_um", "control_index"]
GENE_COLUMNS = ["hgnc_id", "symbol"]
REFERENCES = ("no_change", "delta_mean")
MODEL_VIEWS = {
    "XPert_base": "tahoe_tenv_bulk_l1000/predictor",
    **{f"XPert_{name}": f"gene_axes/views/{name}/predictor"
       for name in ("G2000", "G5000", "G10000", "G18164")},
    **{f"XPert_MLP{width}": "gene_axes/views/G18164/predictor" for width in (512, 1024, 2048, 4096)},
}
SHARED_VIEW = "gene_axes/views/shared/predictor"
FULL_VIEW = MODEL_VIEWS["XPert_G18164"]


def table(path):
    return pd.read_csv(path, keep_default_na=False, dtype={"case_id": str, "group_id": str, "hgnc_id": str})


def matrix(path):
    return np.load(path, mmap_mode="r", allow_pickle=False)


def read_manifest(root):
    path = root / "manifest.json"
    if path.resolve() != path.absolute():
        raise ValueError(f"manifest is a symlink: {path}")
    return read_json(path), sha256(path)


def fit_attempts(work, model):
    """Measure each launch once, separately from checkpointed phase minutes."""
    attempts = []
    for path in sorted((work / "resources" / model).glob("fit-*")):
        if not path.is_dir() or not re.fullmatch(r"fit-[0-9]+", path.name):
            continue
        clock, timing = path / "wall-clock.tsv", path / "time.log"
        seconds, source = None, None
        if clock.exists():
            start, end = map(float, clock.read_text().split())
            if not (np.isfinite(start) and np.isfinite(end) and end >= start):
                raise ValueError(f"invalid fit attempt clock: {clock}")
            seconds, source = end - start, clock
        elif timing.exists():
            match = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:.]+)",
                              timing.read_text(errors="replace"))
            if match:
                fields = map(float, reversed(match[1].split(":")))
                seconds, source = sum(value * 60**i for i, value in enumerate(fields)), timing
        attempts.append({"attempt": path.name, "seconds": seconds,
                         "timing_sha256": sha256(source) if source else None})
    return attempts


class Release:
    def __init__(self, root):
        self.root = root
        self.manifest, self.digest = read_manifest(root)
        self.verified = {}

    def path(self, relative):
        if relative not in self.verified:
            self.verified[relative] = checked_payload(self.root, self.manifest, relative)
        return self.verified[relative]

    def table(self, relative):
        return table(self.path(relative))

    def matrix(self, relative):
        return matrix(self.path(relative))


class Cohort:
    def __init__(self, view, role, evaluator_digest):
        self.genes = view.table("genes.csv")[GENE_COLUMNS]
        self.cases = view.table(f"tenv/fold_0/{role}.csv")[CASE_COLUMNS]
        self.controls = view.matrix("tenv/controls.npy")
        self.control_index = self.cases.control_index.to_numpy()
        target = f"tenv/fold_0/{role}_expression.npy"
        self.evaluator = None
        if role == "test":
            self.evaluator = Release(view.root.parent / "evaluator")
            if (self.evaluator.manifest["predictor_manifest_sha256"] != view.digest
                    or (evaluator_digest is not None and self.evaluator.digest != evaluator_digest)):
                raise ValueError("test expression evaluator differs from its sealed predictor lineage")
        self.treated = (self.evaluator or view).matrix(target)
        if (self.cases.case_id.duplicated().any() or self.genes.hgnc_id.duplicated().any()
                or self.controls.ndim != 2 or self.controls.shape[1] != len(self.genes)
                or np.any(self.control_index < 0) or np.any(self.control_index >= len(self.controls))
                or self.treated.shape != (len(self.cases), len(self.genes))):
            raise ValueError(f"unaligned expression, genes or cases: {view.root}, {role}")

    def binding(self):
        if self.evaluator:
            return self.evaluator.digest
        return None


class Significance:
    def __init__(self, root, role, expected_genes, expected_groups):
        self.release = Release(root)
        manifest = self.release.manifest
        expected = ("step09-significance-v1", role, 0, .05)
        if tuple(manifest[key] for key in ("schema", "role", "fold", "fdr")) != expected:
            raise ValueError("invalid significance truth contract")
        self.genes = self.release.table("genes.csv")[GENE_COLUMNS]
        pd.testing.assert_frame_equal(self.genes, expected_genes)
        groups = self.release.table("groups.csv")
        if (groups.group_id.duplicated().any() or set(groups.group_id) != set(expected_groups)
                or not groups.status.isin(("eligible", "untestable")).all()):
            raise ValueError("significance groups differ from the role cohort")
        self.groups = groups.set_index("group_id")
        self.index = {name: i for i, name in enumerate(groups.group_id)}
        self.padj = self.release.matrix("padj.npy")
        if self.padj.shape != (len(groups), len(self.genes)):
            raise ValueError("significance table has wrong shape")
        for name in ("pvalue.npy", "log2_fold_change.npy", "contrast_samples.csv"):
            self.release.path(name)

    def row(self, group):
        metadata = self.groups.loc[group]
        padj = self.padj[self.index[group]]
        finite = np.isfinite(padj)
        if (np.any(np.isinf(padj)) or np.any(padj[finite] < 0) or np.any(padj[finite] > 1)
                or (metadata.status == "eligible") != bool(finite.any())
                or ("n_genes_testable" in metadata and metadata.n_genes_testable != finite.sum())):
            raise ValueError(f"invalid adjusted p-values for {group}")
        return metadata, padj < .05, finite


class RoundInputs:
    def __init__(self, root, work, corpus, plan_path):
        self.root = root
        self.work = work
        plan_bytes = verify_committed(plan_path)
        self.plan = plan = json.loads(plan_bytes)
        self.plan_sha = hashlib.sha256(plan_bytes).hexdigest()
        if tuple(plan[key] for key in ("schema", "fold", "seed")) != ("step09-round1-plan-v1", 0, 17):
            raise ValueError("invalid Round 1 plan identity")
        if tuple(plan["references"]) != REFERENCES or [(m["name"], m["view"]) for m in plan["models"]] != list(MODEL_VIEWS.items()):
            raise ValueError("Round 1 plan does not enumerate the complete frozen cohort")
        policy_path = plan_path.with_name("policy.json")
        policy_bytes = verify_committed(policy_path)
        self.policy = json.loads(policy_bytes)
        if hashlib.sha256(policy_bytes).hexdigest() != plan["policy_sha256"]:
            raise ValueError("policy differs from the Round 1 plan")
        self.axis_manifest, digest = read_manifest(corpus / "gene_axes")
        if digest != plan["gene_axes_manifest_sha256"]:
            raise ValueError("gene axes differ from the Round 1 plan")
        self.axis_sha = digest
        self.views = {relative: Release(corpus / relative)
                      for relative in dict.fromkeys((*MODEL_VIEWS.values(), SHARED_VIEW))}
        self.evaluator_digests = {}
        for entry in plan["models"]:
            if self.views[entry["view"]].digest != entry["predictor_sha256"]:
                raise ValueError(f"{entry['name']}: predictor differs from the Round 1 plan")
        for name, digests in self.axis_manifest["views"].items():
            relative = f"gene_axes/views/{name}/predictor"
            if self.views[relative].digest != digests["predictor"]:
                raise ValueError(f"{name}: predictor differs from the frozen axis release")
            self.evaluator_digests[relative] = digests["evaluator"]
        landmark = self.views[MODEL_VIEWS["XPert_base"]]
        self.fixture = landmark.manifest["fixture"]
        if not self.fixture:
            record = read_json(STEP / "execution/landmark_record/corpus.json")
            self.evaluator_digests[MODEL_VIEWS["XPert_base"]] = record["manifest_sha256"]["evaluator"]
        self.full_genes = self.views[FULL_VIEW].table("genes.csv")[GENE_COLUMNS]
        base_genes = landmark.table("genes.csv")[GENE_COLUMNS]
        supplement = base_genes.loc[~base_genes.hgnc_id.isin(self.full_genes.hgnc_id)]
        if supplement.symbol.tolist() != ["HDGFL3"]:
            raise ValueError("full significance family requires exactly one HDGFL3 supplement")
        self.significance_genes = pd.concat([self.full_genes, supplement], ignore_index=True)
        shared = self.views[SHARED_VIEW].table("genes.csv")[GENE_COLUMNS]
        intersection = set(base_genes.hgnc_id)
        for view in MODEL_VIEWS.values():
            intersection.intersection_update(self.views[view].table("genes.csv").hgnc_id)
        if (len(shared) != plan["shared_gene_count"]
                or shared.hgnc_id.tolist() != self.full_genes.loc[self.full_genes.hgnc_id.isin(intersection), "hgnc_id"].tolist()
                or plan["continuous_metrics_provisional"] != (len(shared) < 200)):
            raise ValueError("shared continuous panel differs from the frozen model intersection")

    def cohort(self, relative, role):
        return Cohort(self.views[relative], role, self.evaluator_digests.get(relative))

    def statuses(self, role):
        statuses = {}
        pending = []
        for model in MODEL_VIEWS:
            path = self.work / "status" / f"{model}.json"
            status = read_json(path) if path.exists() else {}
            if status and status.get("model") != model:
                raise ValueError(f"status record belongs to a different model: {path}")
            statuses[model] = status
            validation, _ = role_result(status, "validation")
            current, _ = role_result(status, role)
            if validation is None or (role == "test" and validation == "complete" and current is None):
                pending.append(model)
        if pending and not self.fixture:
            raise ValueError(f"Round 1 {role} requires terminal model statuses: {', '.join(pending)}")
        return statuses

    def prediction(self, model, role, cohort, status):
        root = self.root / "predictions" / role / model
        run_path = self.work / "fits" / model / "run.json"
        validation, validation_failure = role_result(status, "validation")
        outcome, failure = role_result(status, role)
        if role == "test" and validation == "failed":
            outcome, failure = "failed", validation_failure
        expected = {"model": model, "fold": 0, "seed": 17, "policy_sha256": self.plan["policy_sha256"],
                    "predictor_sha256": self.views[MODEL_VIEWS[model]].digest}
        binding = {"fit_attempts": fit_attempts(self.work, model)}
        run = None
        if run_path.exists():
            try:
                candidate = read_json(run_path)
                if (candidate["status"] != "complete" or any(candidate[key] != value for key, value in expected.items())
                        or candidate["policy"] != self.policy or not np.isfinite(candidate["training_minutes"])
                        or candidate["training_minutes"] < 0
                        or (not self.fixture and candidate["architecture_sha256"] != self.plan["fit_sources"]["model.py"])):
                    raise ValueError(f"{model}: fit result differs from the frozen Round 1 plan")
                run = candidate
                binding["run_sha256"] = sha256(run_path)
            except (OSError, ValueError, KeyError, TypeError):
                if outcome != "failed":
                    raise
        if outcome == "failed":
            binding["failure"] = {key: failure.get(key) for key in ("stage", "exit_code")}
            reason = f"{failure['stage']} failed (exit {failure.get('exit_code')})"
            return None, run, binding, reason
        replay = status.get(f"{role}_replay_sha256")
        if not run or not (root / "manifest.json").exists() or (not replay and not self.fixture):
            if outcome == "complete":
                raise ValueError(f"{model}: terminal {role} replay is missing its fit or prediction")
            reason = "fit unavailable" if not run else f"{role} prediction or checkpoint replay unavailable"
            return None, run, binding, reason
        release = Release(root)
        manifest = release.manifest
        if (manifest["schema"] != "step09-xpert-v1" or manifest["role"] != role
                or "qualification_limit" in manifest
                or any(manifest["run_key"][key] != value for key, value in expected.items())
                or manifest["run_key"]["policy"] != self.policy
                or manifest["checkpoint_sha256"] != run["phases"]["tenv"]["checkpoint_sha256"]
                or (replay is not None and replay != release.digest)):
            raise ValueError(f"{model}: prediction does not match the verified selected checkpoint")
        if not self.fixture and (run["policy"]["kind"] != "formal" or manifest["run_key"]["policy"]["kind"] != "formal"):
            raise ValueError(f"{model}: qualification predictions cannot enter Round 1")
        cases = release.table("cases.csv")[CASE_COLUMNS]
        genes = release.table("genes.csv")[GENE_COLUMNS]
        pd.testing.assert_frame_equal(cases, cohort.cases)
        pd.testing.assert_frame_equal(genes, cohort.genes)
        values = release.matrix("delta.npy")
        if values.shape != cohort.treated.shape:
            raise ValueError(f"{model}: prediction shape differs from cohort")
        binding["prediction_manifest_sha256"] = release.digest
        binding["checkpoint_replay_sha256"] = replay
        return values, run, binding, ""
