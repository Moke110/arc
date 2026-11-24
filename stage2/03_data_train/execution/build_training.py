"""Publish and replay pooled Stage 2 training products on frozen G_5000."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parent
sys.path[:0] = [str(path) for path in (
    STEP / "output/api", STAGE / "util/output/api", STAGE / "01_data_raw/output/api",
    STAGE / "02_gene_axis/02.1_hvg-deg/output/api", STAGE / "07_tsg_dose/output/api",
)]

from native_corpus import open_corpus
from prediction_axis import open_axes
from stage2_util.derivation import release
from stage2_util.publication import publication_lock, publish_candidate
from training_data import (DEFAULT, DATASET_SCHEMA, FILES, FOLD_COLUMNS, SCHEMA, STRIP_COLUMNS,
                           open_training, sha256)

from tsg_round1.fold_ledger import FOLD_RECEIPT, FOLD_ROOT, load_folds

from frozen_folds import group_labels

SOURCES = {
    "tenv": ("tenv", "DMSO"),
    "tahoe-bulk": ("tahoe_bulk", "DMSO_TF"),
}
STAGING = STAGE.parent / "tmp/stage2-training-staging"
TREATMENT_KEY = ("plate", "cell_line", "drug_parent", "dose_um")
CONTROL_KEY = TREATMENT_KEY[:2]
RECIPE = {
    "normalization": "G_5000 integer sum, CP10K on selected-axis total, natural log1p in float32",
    "C": "equal-weight mean of normalized M_C rows",
    "delta": "normalized M_T treatment minus matched normalized M_C",
}


@dataclass(frozen=True)
class PooledGroup:
    key: tuple
    pooled_row: int
    members: list[int]
    library_total: int

    @property
    def control_key(self):
        return self.key[:2]


def _treatment_reason(row):
    if row.view_exclusion_reason:
        return row.view_exclusion_reason
    if not isinstance(row.drug_parent_inchikey, str) or not row.drug_parent_inchikey:
        return "unresolved_drug"
    if pd.isna(row.dose_um) or not np.isfinite(row.dose_um) or row.dose_um <= 0:
        return "unresolved_dose"
    if not row.plate or not row.cell_line:
        return "unresolved_control_context"
    return ""


def _population(source, solvent):
    text_columns = ["view_exclusion_reason", "plate", "cell_line", "drug_parent_inchikey", "vehicle"]
    obs = source.obs[[*text_columns, "role", "dose_um"]].copy()
    obs[text_columns] = obs[text_columns].astype("string").fillna("")
    controls = defaultdict(list)
    treatments = defaultdict(list)
    exclusions = {}
    keys = {}
    for i, row in enumerate(obs.itertuples(index=False)):
        if row.role == "vehicle":
            if row.vehicle != solvent:
                exclusions[i] = "incompatible_vehicle"
            elif not row.plate or not row.cell_line:
                exclusions[i] = "unresolved_control_context"
            else:
                key = (row.plate, row.cell_line)
                controls[key].append(i)
                keys[i] = key
        elif row.role != "treated":
            exclusions[i] = row.view_exclusion_reason or f"role_{row.role}"
    for i, row in enumerate(obs.itertuples(index=False)):
        if row.role != "treated":
            continue
        reason = _treatment_reason(row)
        key = (row.plate, row.cell_line)
        if not reason and not any(keys.get(member) == key for member in source.matched_control_rows(i)):
            reason = "no_matched_vehicle"
        if reason:
            exclusions[i] = reason
        else:
            treatments[(*key, row.drug_parent_inchikey, row.dose_um)].append(i)
    return obs, controls, treatments, exclusions


def _pool_counts(path, projection, n_rows, groups):
    group_index = np.full(n_rows, -1, dtype=np.int64)
    for index, members in enumerate(groups):
        group_index[members] = index
    pooled = np.lib.format.open_memmap(
        path, mode="w+", dtype=np.int64, shape=(len(groups), len(projection.columns))
    )
    pooled[:] = 0
    for start, block in projection.batches():
        indices = group_index[start:start + len(block)]
        valid = indices >= 0
        np.add.at(pooled, indices[valid], block[valid])
    return pooled


def _select_groups(pooled, controls, treatments, exclusions):
    active_controls = {}
    active_treatments = []
    for group_index, (key, members) in enumerate(controls.items()):
        total = pooled[group_index].sum(dtype=np.int64)
        if total > 0:
            active_controls[key] = PooledGroup(key, group_index, members, total)
        else:
            exclusions.update(dict.fromkeys(members, "zero_control_counts"))
    for group_index, (key, members) in enumerate(treatments.items(), len(controls)):
        total = pooled[group_index].sum(dtype=np.int64)
        if key[:2] not in active_controls:
            exclusions.update(dict.fromkeys(members, "zero_control_counts"))
        elif total <= 0:
            exclusions.update(dict.fromkeys(members, "zero_treatment_counts"))
        else:
            active_treatments.append(PooledGroup(key, group_index, members, total))
    if not active_controls or not active_treatments:
        raise ValueError("training release requires usable controls and treatments")
    return active_controls, active_treatments


def _normalize(counts, totals):
    return np.log1p(counts.astype(np.float64) * (10000.0 / np.asarray(totals)[:, None])).astype(np.float32)


def _write_tables(root, obs, controls, treatments, exclusions, folds, control_rows):
    control_table = [(*group.key, len(group.members), group.library_total)
                     for group in controls.values()]
    pd.DataFrame(control_table, columns=[*CONTROL_KEY, "member_count", "library_total"]).to_csv(
        root / "controls.csv", index=False
    )
    sizes = [len(group.members) for group in treatments]
    members = [member for group in treatments for member in group.members]
    groups = np.repeat(np.arange(len(treatments)), sizes)
    labels = group_labels(folds, obs.iloc[members], groups)
    treatment_table = [
        (*group.key, control_rows[group.control_key], size, group.library_total, *label)
        for group, size, label in zip(treatments, sizes, labels, strict=True)
    ]
    pd.DataFrame(treatment_table, columns=[
        *TREATMENT_KEY, "control_row", "member_count", "library_total", *FOLD_COLUMNS, *STRIP_COLUMNS,
    ]).to_csv(root / "treatments.csv", index=False)
    with (root / "members.csv").open("w", newline="") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(("profile_id", "role", "row"))
        for role, groups in (("control", controls.values()), ("treatment", treatments)):
            for row, group in enumerate(groups):
                writer.writerows((obs.index[member], role, row) for member in group.members)
    with (root / "exclusions.csv").open("w", newline="") as file:
        writer = csv.writer(file, lineterminator="\n")
        writer.writerow(("profile_id", "role", "reason"))
        writer.writerows((obs.index[index], obs.role.iloc[index], reason)
                         for index, reason in sorted(exclusions.items()))


def _write_matrices(root, pooled, controls, treatments, control_rows):
    width = pooled.shape[1]
    control_groups = list(controls.values())
    matrix = np.empty((len(controls), width), dtype=np.float32)
    for start in range(0, len(controls), 64):
        groups = control_groups[start:start + 64]
        rows = [group.pooled_row for group in groups]
        totals = [group.library_total for group in groups]
        matrix[start:start + len(groups)] = _normalize(pooled[rows], totals)
    with h5py.File(root / "training.h5", "w") as handle:
        handle.create_dataset("C", data=matrix.mean(axis=0, dtype=np.float64).astype(np.float32))
        handle.create_dataset("M_C", data=matrix)
        shape = (len(treatments), width)
        chunk = (min(64, shape[0]), width)
        treatment = handle.create_dataset("M_T/treatment", shape, np.float32,
                                          chunks=chunk, compression="gzip", compression_opts=1)
        delta = handle.create_dataset("M_T/delta", shape, np.float32,
                                      chunks=chunk, compression="gzip", compression_opts=1)
        for start in range(0, len(treatments), 64):
            block = treatments[start:start + 64]
            rows = [group.pooled_row for group in block]
            totals = [group.library_total for group in block]
            values = _normalize(pooled[rows], totals)
            matched = matrix[[control_rows[group.control_key] for group in block]]
            treatment[start:start + len(block)] = values
            delta[start:start + len(block)] = values - matched


def _build_dataset(root, key, corpus, axis, fold_root, fold_receipt):
    native_source, solvent = SOURCES[key]
    folds = load_folds(fold_root, fold_receipt, native_source, axis.sha256)
    projection = axis.project(corpus, native_source)
    obs, controls, treatments, exclusions = _population(projection.source, solvent)
    if not controls or not treatments:
        raise ValueError(f"{key} has no candidate controls or treatments")
    root.mkdir()
    cache = root / "pooled.npy"
    pooled = _pool_counts(cache, projection, len(obs), [*controls.values(), *treatments.values()])
    try:
        active_controls, active_treatments = _select_groups(pooled, controls, treatments, exclusions)
        control_rows = {key: row for row, key in enumerate(active_controls)}
        _write_tables(root, obs, active_controls, active_treatments, exclusions, folds, control_rows)
        _write_matrices(root, pooled, active_controls, active_treatments, control_rows)
    finally:
        del pooled
        cache.unlink()
    pd.DataFrame({"gene_index": axis.gene_indices, "symbol": axis.symbols}).to_csv(
        root / "genes.csv", index=False
    )
    source_record = corpus.manifest["sources"][native_source]
    release.write_json(root / "manifest.json", {
        "schema": DATASET_SCHEMA,
        "source": key,
        "axis_size": len(axis.gene_indices),
        "recipe": {**RECIPE, "vehicle": solvent},
        "inputs": {
            "raw_source": native_source,
            "raw_sha256": corpus.manifest_sha256,
            "raw_dataset_sha256": source_record["manifest_sha256"],
            "raw_counts_sha256": source_record["sha256"],
            "axis_sha256": axis.sha256,
            "folds": folds.evidence,
        },
        "counts": {
            "raw_rows": len(obs), "candidate_controls": len(controls),
            "candidate_treatments": len(treatments), "controls": len(active_controls),
            "treatments": len(active_treatments), "excluded_raw_rows": len(exclusions),
        },
        "outputs": {name: sha256(root / name) for name in FILES},
    })


def _build(root, corpus, axis, fold_root, fold_receipt, sources):
    for key in sources:
        _build_dataset(root / key, key, corpus, axis, Path(fold_root), Path(fold_receipt))
    release.write_json(root / "manifest.json", {
        "schema": SCHEMA,
        "raw_sha256": corpus.manifest_sha256,
        "axis_sha256": axis.sha256,
        "sources": {key: sha256(root / key / "manifest.json") for key in sources},
    })


def _verify_dataset(data, obs, axis):
    if (tuple(data.genes.gene_index) != axis.gene_indices
            or tuple(data.genes.symbol) != axis.symbols):
        raise ValueError("training gene order differs")
    members = data.membership()
    excluded = data.exclusions()
    identities = pd.concat([members.profile_id, excluded.profile_id], ignore_index=True)
    if (identities.duplicated().any() or set(identities) != set(obs.index)
            or not members.role.isin({"control", "treatment"}).all()
            or excluded.reason.eq("").any()
            or len(excluded) != data.manifest["counts"]["excluded_raw_rows"]):
        raise ValueError("training membership and exclusions do not partition the raw observations")
    for role, table in (("control", data.controls), ("treatment", data.treatments)):
        counts = members.loc[members.role.eq(role), "row"].value_counts().sort_index()
        if (list(counts.index) != list(table.index)
                or not np.array_equal(counts.to_numpy(), table.member_count.to_numpy())):
            raise ValueError(f"{role} member counts differ from the ledger")
    control_sum = np.zeros(len(data.genes), dtype=np.float64)
    controls = []
    for start in range(0, len(data.controls), 64):
        block = data.read_controls(start, min(start + 64, len(data.controls)))
        if not np.isfinite(block).all() or (block < 0).any():
            raise ValueError("normalized control values must be finite and nonnegative")
        control_sum += block.sum(axis=0, dtype=np.float64)
        controls.append(block)
    if not np.array_equal(data.control(), (control_sum / len(data.controls)).astype(np.float32)):
        raise ValueError("C is not the equal-weight mean of M_C")
    control_matrix = np.concatenate(controls)
    for start in range(0, len(data.treatments), 64):
        stop = min(start + 64, len(data.treatments))
        treatment = data.read_treatments(start, stop)
        control = control_matrix[data.treatments.control_row.iloc[start:stop].to_numpy()]
        delta = data.read_treatments(start, stop, "delta")
        if (not np.isfinite(treatment).all() or (treatment < 0).any()
                or not np.array_equal(delta, treatment - control)):
            raise ValueError("M_T contains invalid expression or unmatched Delta")


def verify_training(root, corpus, axis, *, sources=None):
    """Verify release pins, complete provenance ledgers, and bounded scientific reads."""
    root = Path(root)
    manifest_path = root / "manifest.json"
    record = release.read_json(manifest_path)
    expected = tuple(SOURCES if sources is None else sources)
    if (record["schema"] != SCHEMA or set(record["sources"]) != set(expected)
            or record["raw_sha256"] != corpus.manifest_sha256 or record["axis_sha256"] != axis.sha256):
        raise ValueError("training release inventory or generation differs")
    names = {path.name for path in root.iterdir()}
    if names - {"api"} != {"manifest.json", *expected}:
        raise ValueError("unexpected files in training release")
    for key in expected:
        native_source, solvent = SOURCES[key]
        with open_training(root, source=key) as dataset:
            source_record = corpus.manifest["sources"][native_source]
            inputs = dataset.manifest["inputs"]
            if (inputs["raw_source"] != native_source
                    or inputs["raw_dataset_sha256"] != source_record["manifest_sha256"]
                    or inputs["raw_counts_sha256"] != source_record["sha256"]
                    or dataset.manifest["recipe"] != {**RECIPE, "vehicle": solvent}):
                raise ValueError(f"training dataset inputs or recipe differ: {key}")
            _verify_dataset(dataset, corpus.source(native_source).obs, axis)
    return {"result": "PASS", "manifest_sha256": sha256(manifest_path), "sources": list(expected)}


def _sources(sources):
    selected = tuple(SOURCES if sources is None else sources)
    if not selected or len(set(selected)) != len(selected) or not set(selected).issubset(SOURCES):
        raise ValueError("unknown or duplicate training release sources")
    return selected


def _install_training(destination, corpus, axis, selected, build, staging):
    with publication_lock(destination) as output:
        stage = Path(staging or STAGING).resolve()
        if stage.is_relative_to(output):
            raise ValueError("training staging must be outside the accepted output")

        def verify_existing(root):
            sources = tuple(release.read_json(root / "manifest.json")["sources"])
            return verify_training(root, corpus, axis, sources=sources)

        if (output / "manifest.json").exists():
            existing = release.read_json(output / "manifest.json")
            if not set(existing["sources"]).issubset(selected):
                raise ValueError("publication cannot drop a previously accepted training source")
        elif output.exists() and {path.name for path in output.iterdir()} - {"api"}:
            raise ValueError("refusing to replace unrecognized files in training output")

        def prepare(candidate):
            build(candidate)
            if (output / "api").is_dir():
                shutil.copytree(output / "api", candidate / "api",
                                ignore=shutil.ignore_patterns("__pycache__"))

        return publish_candidate(
            release, output, prepare,
            lambda root: verify_training(root, corpus, axis, sources=selected),
            staging=stage, verify_existing=verify_existing,
        )


def publish_training(destination, corpus, axis, fold_root, fold_receipt, *, sources=None, staging=None):
    """Build and verify a release at a caller-specified scratch destination."""
    if Path(destination).resolve() == DEFAULT.resolve():
        raise ValueError("publish requires a scratch destination; use promote for the stable output")
    selected = _sources(sources)

    def build(candidate):
        _build(candidate, corpus, axis, fold_root, fold_receipt, selected)

    return _install_training(destination, corpus, axis, selected, build, staging)


def promote_training(candidate, destination, corpus, axis, fold_root=FOLD_ROOT,
                     fold_receipt=FOLD_RECEIPT, *, expected_sha256, staging=None):
    """Rebuild and install a complete pair only when it matches the pinned candidate."""
    source = Path(candidate).resolve()
    output = Path(destination).resolve()
    stage = Path(staging or STAGING).resolve()
    if source.is_relative_to(output) or output.is_relative_to(source) or stage.is_relative_to(source):
        raise ValueError("training candidate, output and staging must keep the candidate separate")
    if output == DEFAULT.resolve() and sha256(Path(fold_receipt)) != sha256(FOLD_RECEIPT):
        raise ValueError("stable promotion requires the original Step 07 fold receipt")
    selected = _sources(None)

    def build(prepared):
        for directory in (source, *(source / key for key in selected)):
            if directory.is_symlink() or not directory.is_dir():
                raise ValueError(f"training candidate directory must be a real directory: {directory.name}")
            for path in directory.iterdir():
                if path.name != "api" and path.is_symlink():
                    raise ValueError(f"training candidate asset must not be a symlink: {path.name}")
        if sha256(source / "manifest.json") != expected_sha256:
            raise ValueError("training candidate manifest digest differs from the accepted digest")
        accepted = verify_training(source, corpus, axis, sources=selected)
        if accepted["manifest_sha256"] != expected_sha256:
            raise ValueError("training candidate manifest changed during verification")
        _build(prepared, corpus, axis, fold_root, fold_receipt, selected)
        if sha256(prepared / "manifest.json") != expected_sha256:
            raise ValueError("replayed training manifest differs from the pinned candidate")

    return _install_training(destination, corpus, axis, selected, build, stage)


def replay_training(destination, corpus, axis, fold_root, fold_receipt, *, sources=None, staging=None):
    """Rebuild separately and compare both the release and dataset manifest digests."""
    selected = _sources(sources)
    with publication_lock(destination) as output:
        stage = Path(staging or STAGING).resolve()
        if stage.is_relative_to(output):
            raise ValueError("training staging must be outside the accepted output")
        accepted = verify_training(output, corpus, axis, sources=selected)
        candidate = release.candidate(output, stage)
        try:
            _build(candidate, corpus, axis, fold_root, fold_receipt, selected)
            replayed = verify_training(candidate, corpus, axis, sources=selected)
            if replayed != accepted:
                raise ValueError("replayed training manifest differs from the accepted release")
        finally:
            shutil.rmtree(candidate)
        return accepted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("publish", "verify", "promote"))
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--staging", type=Path, default=STAGING)
    parser.add_argument("--raw", type=Path, default=STAGE / "01_data_raw/output")
    parser.add_argument("--axis", type=Path, default=STAGE / "02_gene_axis/02.1_hvg-deg/output/release")
    parser.add_argument("--fold-root", type=Path, default=FOLD_ROOT)
    parser.add_argument("--fold-receipt", type=Path, default=FOLD_RECEIPT)
    parser.add_argument("--source", choices=tuple(SOURCES), action="append")
    parser.add_argument("--candidate", type=Path, help="accepted paired candidate for promotion")
    parser.add_argument("--expected-sha256", help="paired manifest SHA-256 from verification")
    args = parser.parse_args()
    if args.command == "publish" and args.output.resolve() == DEFAULT.resolve():
        parser.error("publish requires a scratch destination via --output; use promote for the stable output")
    if args.command == "promote":
        if args.candidate is None or args.expected_sha256 is None or args.source is not None:
            parser.error("promote requires --candidate and --expected-sha256 for the complete pair")
    elif args.candidate is not None or args.expected_sha256 is not None:
        parser.error("--candidate and --expected-sha256 apply to promote only")
    with open_corpus(args.raw) as corpus:
        axis = open_axes(args.axis).axis(5000)
        arguments = (args.output, corpus, axis, args.fold_root, args.fold_receipt)
        kwargs = {"sources": args.source, "staging": args.staging}
        if args.command == "promote":
            output = promote_training(args.candidate, *arguments,
                                      expected_sha256=args.expected_sha256, staging=args.staging)
            result = {"result": "PASS", "manifest_sha256": sha256(output / "manifest.json")}
        elif args.command == "publish":
            output = publish_training(*arguments, **kwargs)
            result = {"result": "PASS", "manifest_sha256": sha256(output / "manifest.json")}
        else:
            result = replay_training(*arguments, **kwargs)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
