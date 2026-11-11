"""Build fold 0 evaluator-only tenv significance labels from native counts."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd

from round_contract import promotion_digest, verify_committed_promotion

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "01_corpus"))
from build import DEFAULT_NATIVE, DEFAULT_OUTPUT as DEFAULT_LANDMARK, STEP
from build import file_hashes, gene_mapping, raw_batches, sha256, write_json
from gene_axes import DEFAULT_OUTPUT as DEFAULT_AXES, hashed_payload


FDR = 0.05
FOLD = 0
FULL_GENE_COUNT = 18164
WORKER = Path(__file__).with_name("deseq2_contrasts.R")
STATISTICS = ("log2_fold_change", "pvalue", "padj")


def manifest(path):
    return json.loads(path.read_text())


def pinned_manifest(directory, digest):
    path = directory / "manifest.json"
    if path.resolve() != path.absolute() or sha256(path) != digest:
        raise ValueError(f"release manifest changed: {path}")
    return manifest(path)


def load_release(axes, landmark, role):
    axes_path = axes / "manifest.json"
    axis_manifest = manifest(axes_path)
    if axis_manifest["schema"] != "step09-gene-axes-v1" or FOLD not in axis_manifest["ready_folds"]:
        raise ValueError("expected a fold 0 gene-axis release")
    if not axis_manifest["fixture"]:
        accepted = STEP / "execution/gene_axes_record/release-manifest.json"
        if sha256(axes_path) != sha256(accepted):
            raise ValueError("gene-axis manifest differs from the accepted release")
    view = axes / "views/G18164"
    predictor, evaluator = view / "predictor", view / "evaluator"
    landmark_predictor = landmark / "predictor"
    digests = axis_manifest["views"]["G18164"]
    pred = pinned_manifest(predictor, digests["predictor"])
    eval_ = pinned_manifest(evaluator, digests["evaluator"])
    land = pinned_manifest(landmark_predictor, axis_manifest["landmark_predictor_sha256"])
    if (pred["axis"] != "G18164" or eval_["predictor_manifest_sha256"] != digests["predictor"]
            or land["axis"] != "landmark" or pred["fixture"] != axis_manifest["fixture"]):
        raise ValueError("significance release views do not share the frozen cohort")
    fields = ["hgnc_id", "symbol", "ensembl_id"]
    cm = pd.read_csv(hashed_payload(predictor, pred, "genes.csv"), dtype=str)[fields]
    supplemental = pd.read_csv(hashed_payload(landmark_predictor, land, "genes.csv"), dtype=str)[fields]
    supplemental = supplemental.loc[supplemental.symbol.eq("HDGFL3")]
    if (len(supplemental) != 1 or len(cm) != pred["gene_count"]
            or (not axis_manifest["fixture"] and len(cm) != FULL_GENE_COUNT)):
        raise ValueError("invalid CM family or HDGFL3 supplement")
    family = pd.concat([cm, supplemental], ignore_index=True)
    if family.hgnc_id.isna().any() or family.hgnc_id.duplicated().any():
        raise ValueError("the significance family must contain unique measured genes")

    cases = pd.read_csv(hashed_payload(predictor, pred, f"tenv/fold_{FOLD}/{role}.csv"))
    if role == "test":
        if "source_row" in cases:
            raise ValueError("test native row joins must remain evaluator-only")
        rows = pd.read_csv(hashed_payload(evaluator, eval_, f"tenv/fold_{FOLD}/test_rows.csv"))
        cases = cases.merge(rows, on="case_id", how="left", validate="one_to_one")
        if len(cases) != len(rows):
            raise ValueError("evaluator test rows do not match the predictor cases")
    controls = pd.read_csv(hashed_payload(predictor, pred, "tenv/controls.csv"))
    rows = pd.read_csv(hashed_payload(evaluator, eval_, "tenv_control_rows.csv"))
    control_rows = rows.merge(controls, on="control_index", how="left", validate="many_to_one")
    if (controls.plate.duplicated().any() or control_rows.plate.isna().any()
            or not np.array_equal(rows.groupby("control_index").size().reindex(controls.control_index),
                                  controls.n_wells)):
        raise ValueError("frozen DMSO controls do not reproduce matched plate counts")

    provenance = {
        "gene_axes_manifest_sha256": sha256(axes_path),
        "predictor_manifest_sha256": digests["predictor"],
        "evaluator_manifest_sha256": digests["evaluator"],
        "landmark_predictor_manifest_sha256": axis_manifest["landmark_predictor_sha256"],
        "native_sha256": axis_manifest["source_hashes"]["tenv.h5ad"],
    }
    return family, cases, control_rows, axis_manifest["fixture"], provenance


def check_native_rows(obs, cases, controls):
    if (cases.empty or cases.case_id.duplicated().any() or obs.index.has_duplicates
            or any(not pd.api.types.is_integer_dtype(frame.source_row)
                   or frame.source_row.duplicated().any()
                   or not frame.source_row.between(0, len(obs) - 1).all()
                   for frame in (cases, controls))):
        raise ValueError("invalid native tenv row joins")
    treated = obs.iloc[cases.source_row]
    matched = obs.iloc[controls.source_row]
    if (treated.index.tolist() != cases.case_id.tolist() or treated.is_control.any()
            or treated.is_blank.any() or not np.array_equal(treated.plate, cases.plate)
            or not np.array_equal(treated.drug_key, cases.drug_key)
            or not np.allclose(treated.dose_um, cases.dose_um, rtol=1e-14, atol=0)):
        raise ValueError("native treated wells do not match the frozen role")
    if (not matched.is_control.all() or matched.is_blank.any()
            or not matched.compound_name.str.upper().eq("DMSO").all()
            or not matched.vehicle.str.upper().eq("DMSO").all()
            or not np.array_equal(matched.plate, controls.plate)):
        raise ValueError("native control wells do not match same-plate DMSO")


def write_shard(matrix, columns, contrasts, controls, workspace):
    groups, samples = [], []
    matched_plate = controls.drop_duplicates("control_index").set_index("control_index").plate
    with (workspace / "counts.bin").open("wb") as handle:
        for group_id, treated in contrasts:
            plates = sorted(treated.plate.unique())
            matched = controls.loc[controls.control_index.isin(treated.control_index.unique())]
            if (set(plates) != set(matched.plate) or treated.drug_key.nunique() != 1
                    or treated.dose_um.nunique() != 1
                    or not np.array_equal(treated.control_index.map(matched_plate), treated.plate)):
                raise ValueError(f"invalid plate-matched drug–dose group: {group_id}")
            batch = pd.concat([treated.assign(condition="treated"),
                               matched.assign(condition="control")]).sort_values("source_row")
            groups.append({
                "group_id": group_id, "drug_key": treated.drug_key.iloc[0],
                "dose_um": float(treated.dose_um.iloc[0]), "n_treated": len(treated),
                "n_controls": len(matched), "n_plates": len(plates), "plates": "|".join(plates),
                "design": "~ plate + condition" if len(plates) > 1 else "~ condition",
            })
            samples.append(batch[["case_id", "source_row", "plate", "condition"]].assign(group_id=group_id))
            for _, counts in raw_batches(matrix, batch.source_row.to_numpy(), columns):
                if np.max(counts) > np.iinfo(np.int32).max:
                    raise ValueError("native count is too large for DESeq2 integers")
                counts.astype("<i4", copy=False).tofile(handle)
    groups = pd.DataFrame(groups)
    groups.to_csv(workspace / "groups.csv", index=False)
    pd.concat(samples, ignore_index=True).to_csv(workspace / "samples.csv", index=False)
    return groups


def write_contrasts(native, family, cases, controls, fixture, workspace, workers):
    with h5py.File(native, "r") as data:
        obs = read_elem(data["obs"])
        check_native_rows(obs, cases, controls)
        mapping = gene_mapping("tenv", family, read_elem(data["var"]), None, fixture)
        columns = mapping.native_column.to_numpy()
        controls = controls.assign(case_id=obs.index[controls.source_row])
        groups = list(cases.groupby("group_id", sort=True))
        if workers == 1 or len(groups) == 1:
            return [(workspace, write_shard(data["X"], columns, groups, controls, workspace))]
        split = (len(groups) + 1) // 2
        shards = []
        for index, part in enumerate((groups[:split], groups[split:])):
            path = workspace / f"shard-{index}"
            path.mkdir()
            shards.append((path, write_shard(data["X"], columns, part, controls, path)))
        return shards


def fit_shards(rscript, shards, n_genes):
    commands = [[rscript, "--vanilla", str(WORKER), str(path), str(n_genes)] for path, _ in shards]
    if len(commands) == 1:
        subprocess.run(commands[0], check=True)
        return
    sizes = " + ".join(str(len(groups)) for _, groups in shards)
    print(f"Starting two DESeq2 workers with {sizes} contrasts", flush=True)
    workers = []
    try:
        for command in commands:
            workers.append(subprocess.Popen(command))
        active = workers.copy()
        while active:
            for worker in active[:]:
                code = worker.poll()
                if code is not None:
                    if code:
                        raise subprocess.CalledProcessError(code, worker.args)
                    active.remove(worker)
            if active:
                time.sleep(0.2)
    finally:
        for worker in workers:
            if worker.poll() is None:
                worker.terminate()
        for worker in workers:
            try:
                worker.wait(timeout=5)
            except subprocess.TimeoutExpired:
                worker.kill()
                worker.wait()


def merge_samples(shards, destination):
    with destination.open("wb") as handle:
        header = None
        for path, _ in shards:
            with (path / "samples.csv").open("rb") as source:
                columns = source.readline()
                if header is None:
                    header = columns
                    handle.write(header)
                elif header != columns:
                    raise ValueError("DESeq2 sample columns differ across shards")
                shutil.copyfileobj(source, handle)


def write_statistic(shards, field, destination, n_genes):
    shape = (sum(len(groups) for _, groups in shards), n_genes)
    result = np.lib.format.open_memmap(destination, mode="w+", dtype="<f8", shape=shape)
    offset = 0
    for path, groups in shards:
        source = path / f"{field}.bin"
        if source.stat().st_size != len(groups) * n_genes * 8:
            raise ValueError(f"invalid DESeq2 result length: {source}")
        with source.open("rb") as handle:
            for row in range(len(groups)):
                result[offset + row] = np.fromfile(handle, dtype="<f8", count=n_genes)
        offset += len(groups)
    result.flush()


def build(native, axes, landmark, role, output, rscript, promotion=None, workers=1):
    if workers not in (1, 2):
        raise ValueError("DESeq2 requires one or two workers")
    if output.exists():
        raise FileExistsError(output)
    if role == "test":
        if promotion is None:
            raise ValueError("test labels require a committed validation promotion")
        verify_committed_promotion(promotion)
    sources = {Path(__file__).name: sha256(Path(__file__)), WORKER.name: sha256(WORKER)}
    family, cases, controls, fixture, provenance = load_release(axes, landmark, role)
    if role == "test":
        provenance["promotion_sha256"] = promotion_digest(promotion)
    if sha256(native) != provenance["native_sha256"]:
        raise ValueError("native tenv counts differ from the frozen gene-axis release")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=f".{output.name}-") as temporary:
        workspace = Path(temporary)
        result = workspace / "release"
        result.mkdir()
        shards = write_contrasts(native, family, cases, controls, fixture, workspace, workers)
        fit_shards(rscript, shards, len(family))
        groups = pd.concat([part for _, part in shards], ignore_index=True)
        result_tables = []
        for path, part in shards:
            table = pd.read_csv(path / "results.csv", keep_default_na=False)
            if len(part) != len(table):
                raise ValueError("DESeq2 did not return every contrast")
            result_tables.append(table)
        metadata = pd.concat(result_tables, ignore_index=True)
        eligible = int(metadata.status.eq("eligible").sum())
        errors = metadata.reason[metadata.reason.str.startswith("DESeq2: ") | metadata.reason.eq("no_finite_pvalues")]
        if not eligible and not errors.empty:
            raise ValueError(f"DESeq2 produced no eligible contrast. Reasons: {errors.value_counts().to_dict()}")
        pd.concat([groups, metadata], axis=1).to_csv(result / "groups.csv", index=False)
        family[["hgnc_id", "symbol"]].to_csv(result / "genes.csv", index=False)
        merge_samples(shards, result / "contrast_samples.csv")
        for field in STATISTICS:
            write_statistic(shards, field, result / f"{field}.npy", len(family))
        versions = {(path / "versions.txt").read_text() for path, _ in shards}
        if len(versions) != 1:
            raise ValueError("DESeq2 runtime versions differ across shards")
        r_version, deseq_version = versions.pop().splitlines()
        write_json(result / "manifest.json", {
            "schema": "step09-significance-v1", "role": role, "fold": FOLD, "fdr": FDR,
            "groups": len(groups), "eligible_groups": eligible,
            "gene_count": len(family), "r_version": r_version, "deseq2_version": deseq_version,
            "method": "two-sided DESeq2 Wald, per-contrast BH over CM genes plus HDGFL3",
            "independent_filtering": False, "size_factor_fallback": "poscounts",
            "sources": sources,
            **provenance, "files": file_hashes(result),
        })
        result.rename(output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--role", choices=("validation", "test"), required=True)
    parser.add_argument("--native", type=Path, default=DEFAULT_NATIVE / "tenv.h5ad")
    parser.add_argument("--gene-axes", type=Path, default=DEFAULT_AXES)
    parser.add_argument("--landmark", type=Path, default=DEFAULT_LANDMARK)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1,
                        help="split into up to two independent DESeq2 processes")
    parser.add_argument("--promotion", type=Path, help="committed validation seal required for test")
    args = parser.parse_args()
    print(build(args.native, args.gene_axes, args.landmark, args.role, args.output, args.rscript,
                args.promotion, args.workers),
          flush=True)


if __name__ == "__main__":
    main()
