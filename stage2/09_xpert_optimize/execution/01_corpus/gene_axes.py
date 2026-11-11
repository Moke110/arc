"""Freeze the fold-local HVG panels and build the Round 1 model-ready views."""

import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import shutil
import tempfile

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd

from build import (DEFAULT_NATIVE, DEFAULT_OUTPUT as LANDMARK, STEP, file_hashes, fingerprint,
                   gene_mapping, input_files, normalize_counts, raw_batches, read_csv, sha256,
                   write_json)


AXES = ("G2000", "G5000", "G10000", "G18164")
DEFAULT_SIZES = (2000, 5000, 10000)
DEFAULT_OUTPUT = LANDMARK.with_name("gene_axes")
N_MEAN_BINS = 20


def hashed_payload(root, manifest, relative):
    path = root / relative
    if path.resolve() != path.absolute() or sha256(path) != manifest["files"][relative]:
        raise ValueError(f"landmark payload differs from its manifest: {relative}")
    return path


def source_dispersion(matrix, rows, columns):
    """Accumulate per-gene moments of train-only linear CP10K in bounded slabs."""
    if len(rows) < 2:
        raise ValueError("HVG ranking requires at least two training wells per source")
    mean, m2 = np.zeros(len(columns)), np.zeros(len(columns))
    total, zeros = 0, 0
    for _, raw in raw_batches(matrix, np.sort(rows), columns):
        linear = raw.astype(np.float32)
        library = raw.sum(axis=1, dtype=np.float64)
        scale = np.divide(10_000, library, out=np.zeros_like(library), where=library != 0)
        linear *= scale.astype(np.float32)[:, None]
        batch_mean = linear.mean(axis=0, dtype=np.float64)
        centered = linear - batch_mean
        n = len(linear)
        delta = batch_mean - mean
        merged = total + n
        m2 += np.einsum("ij,ij->j", centered, centered) + delta**2 * (total * n / merged)
        mean += delta * (n / merged)
        total = merged
        zeros += int(np.count_nonzero(library == 0))

    safe_mean = np.where(mean == 0, 1e-12, mean)
    log_mean = np.log1p(safe_mean)
    dispersion = m2 / ((total - 1) * safe_mean)
    log_dispersion = np.full_like(dispersion, np.nan)
    np.log(dispersion, out=log_dispersion, where=dispersion > 0)
    groups = pd.Series(log_dispersion).groupby(pd.cut(log_mean, N_MEAN_BINS), observed=True)
    avg, sd = groups.transform("mean").to_numpy(), groups.transform("std").to_numpy()
    score = np.full_like(dispersion, np.nan)
    np.divide(log_dispersion - avg, sd, out=score, where=sd > 0)
    score[(groups.transform("count").to_numpy() == 1) & np.isfinite(log_dispersion)] = 1
    order = np.argsort(-np.nan_to_num(score, nan=-np.inf), kind="stable")
    rank = np.empty_like(order)
    rank[order] = np.arange(1, len(order) + 1)
    return {"log_mean": log_mean, "log_dispersion": log_dispersion, "normalized_dispersion": score,
            "rank": rank}, {"train_wells": total, "zero_libraries": zeros}


def source_metadata(landmark, predictor_manifest, evaluator_manifest, source, fold=0):
    predictor, evaluator = landmark / "predictor", landmark / "evaluator"

    def read_predictor(relative):
        return pd.read_csv(hashed_payload(predictor, predictor_manifest, relative), keep_default_na=False)

    prefix = "tahoe" if source == "tahoe" else f"tenv/fold_{fold}"
    case_files = {name: f"{prefix}/{name}.csv" for name in ("train", "validation")}
    cases = {name: read_predictor(path) for name, path in case_files.items()}
    if source == "tenv":
        case_files["test"] = f"{prefix}/test.csv"
        row_file = f"tenv/fold_{fold}/test_rows.csv"
        test_rows = pd.read_csv(hashed_payload(evaluator, evaluator_manifest, row_file), keep_default_na=False)
        cases["test"] = read_predictor(case_files["test"]).merge(
            test_rows, on="case_id", validate="one_to_one", how="left"
        )
        if cases["test"].source_row.isna().any():
            raise ValueError("test case missing native row")
    controls_relative = f"{source}/controls.csv"
    controls = read_predictor(controls_relative)
    control_rows_relative = f"{source}_control_rows.csv"
    control_rows = pd.read_csv(hashed_payload(evaluator, evaluator_manifest, control_rows_relative))
    if (not np.array_equal(controls.control_index, np.arange(len(controls)))
            or not np.array_equal(control_rows.control_index,
                                  np.repeat(controls.control_index, controls.n_wells))):
        raise ValueError(f"{source}: invalid control row identities")
    files = {**case_files, "controls": controls_relative, "morgan": f"{source}/morgan.npy",
             "control_rows": control_rows_relative}
    if source == "tenv":
        files["test_rows"] = row_file
    return cases, controls, control_rows, files


def panel_tables(directory, cm, landmark_genes, ranked, fold, sizes, axes=AXES):
    directory.mkdir()
    ranks = cm.copy()
    for source, data in ranked.items():
        for field, values in data.items():
            ranks[f"{source}_{field}"] = values
    ranks["rank_sum"] = ranks.tahoe_rank + ranks.tenv_rank
    ordered = ranks.rank_sum.sort_values(kind="stable").index.to_numpy()
    ranks.to_csv(directory / "ranks.csv", index=False)
    panels = {name: cm.iloc[np.sort(ordered[:size])] for name, size in zip(axes[:3], sizes, strict=True)}
    panels["G18164"] = cm
    panels["shared"] = panels[axes[0]].loc[panels[axes[0]].hgnc_id.isin(landmark_genes.hgnc_id)]
    for name, panel in panels.items():
        panel.to_csv(directory / f"{name}.csv", index=False)
    shared_count = len(panels["shared"])
    print(json.dumps({"fold": fold, "shared_gene_count": shared_count,
                      "shared_continuous_provisional": shared_count < 200}), flush=True)
    return panels


def open_expression(stack, path, rows, genes):
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = stack.enter_context(path.open("wb"))
    np.lib.format.write_array_header_2_0(handle, {
        "descr": np.dtype(np.float32).str, "fortran_order": False, "shape": (rows, genes),
    })
    return handle


def stream_views(matrix, source, columns, view_roots, panels, metadata, landmark, manifests):
    cases, controls, control_rows, files = metadata
    names = tuple(cases)
    for role, table in cases.items():
        if np.any(np.diff(table.source_row) <= 0):
            raise ValueError(f"{source} {role} cases are not in native order")
    requested = np.concatenate([*(cases[role].source_row for role in names), control_rows.source_row])
    roles = np.repeat(np.arange(len(names) + 1), [*(len(cases[role]) for role in names), len(control_rows)])
    indices = np.concatenate([np.zeros(sum(map(len, cases.values())), dtype=int), control_rows.control_index])
    order = np.argsort(requested)
    requested, roles, indices = requested[order], roles[order], indices[order]
    zeros = {name: {role: 0 for role in (*names, "controls")} for name in view_roots}
    sums = {name: np.zeros((len(controls), len(panels[name])), dtype=np.float64) for name in view_roots}

    predictor, evaluator = landmark / "predictor", landmark / "evaluator"
    with ExitStack() as stack:
        sinks = {}
        for name, root in view_roots.items():
            destination = root / "predictor" / source
            destination.mkdir(parents=True)
            panel = panels[name]
            cols = panel.identity_index.to_numpy(dtype=int)
            mapping = columns.iloc[cols]
            mapping.to_csv(destination / "mapping.csv", index=False)
            for key in (*names, "controls", "morgan"):
                relative = files[key]
                path = root / "predictor" / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(hashed_payload(predictor, manifests["predictor"], relative), path)
            evaluator_root = root / "evaluator"
            for key in ("control_rows", "test_rows"):
                if key in files:
                    relative = files[key]
                    path = evaluator_root / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(hashed_payload(evaluator, manifests["evaluator"], relative), path)
            outputs = []
            for role in names:
                side = "evaluator" if role == "test" else "predictor"
                path = root / side / Path(files[role]).with_name(f"{role}_expression.npy")
                outputs.append(open_expression(stack, path, len(cases[role]), len(panel)))
            sinks[name] = (cols, outputs)

        for part, full_counts in raw_batches(matrix, requested, columns.native_column.to_numpy()):
            role_indices, stratum = roles[part], indices[part]
            masks = [role_indices == i for i in range(len(names) + 1)]
            for name, (cols, outputs) in sinks.items():
                counts = full_counts[:, cols]
                for role, mask, out in zip(names, masks[:-1], outputs, strict=True):
                    if mask.any():
                        values, zero = normalize_counts(counts[mask])
                        values.tofile(out)
                        zeros[name][role] += zero
                mask = masks[-1]
                if mask.any():
                    values, zero = normalize_counts(counts[mask])
                    np.add.at(sums[name], stratum[mask], values)
                    zeros[name]["controls"] += zero
    for name, root in view_roots.items():
        values = (sums[name] / controls.n_wells.to_numpy()[:, None]).astype(np.float32)
        np.save(root / "predictor" / source / "controls.npy", values)
    return zeros


def manifest_record(root, panels, zeros, landmark_sha, frozen, fold=0):
    manifests = {}
    for name, panel in panels.items():
        path = root / "views" / name
        predictor, evaluator = path / "predictor", path / "evaluator"
        view = {"schema": "step09-predictor-v1", "axis": name, "gene_count": len(panel),
                "ready_folds": [fold], "fixture": frozen["fixture"], "fit_eligible": name != "shared",
                "training_seeds": frozen["training_seeds"], "landmark_predictor_sha256": landmark_sha,
                "panel_sha256": sha256(root / f"fold_{fold}/{name}.csv"),
                "zero_libraries": {source: {role: n for role, n in views[name].items() if role != "test"}
                                   for source, views in zeros.items() if name in views},
                "files": file_hashes(predictor)}
        write_json(predictor / "manifest.json", view)
        predictor_sha = sha256(predictor / "manifest.json")
        write_json(evaluator / "manifest.json", {
            "schema": "step09-evaluator-v1", "predictor_manifest_sha256": predictor_sha,
            "test_zero_libraries": {"tenv": {f"fold_{fold}/test": zeros["tenv"][name]["test"]}},
            "files": file_hashes(evaluator),
        })
        manifests[name] = {"predictor": predictor_sha, "evaluator": sha256(evaluator / "manifest.json")}
    return manifests


def build(native, landmark, output, fixture_sizes, rerun=False, materialized_fold=0):
    if output.exists():
        raise FileExistsError(output)
    manifests = {side: json.loads((landmark / side / "manifest.json").read_text())
                 for side in ("predictor", "evaluator")}
    frozen = manifests["predictor"]
    landmark_sha = sha256(landmark / "predictor/manifest.json")
    if (frozen["schema"] != "step09-predictor-v1" or frozen["axis"] != "landmark"
            or manifests["evaluator"]["predictor_manifest_sha256"] != landmark_sha):
        raise ValueError("landmark predictor and evaluator manifests differ")
    if frozen["fixture"] != (fixture_sizes is not None):
        raise ValueError("fixture panel sizes require a fixture landmark release")
    if not frozen["fixture"]:
        corpus_record = json.loads((STEP / "execution/landmark_record/corpus.json").read_text())
        for side, digest in corpus_record["manifest_sha256"].items():
            if sha256(landmark / side / "manifest.json") != digest:
                raise ValueError(f"landmark {side} manifest differs from the accepted release")

    inputs = input_files(native)
    starting_stats = fingerprint(inputs)
    for name, path in inputs.items():
        if sha256(path) != frozen["source_hashes"][name]:
            raise ValueError(f"source differs from the landmark release: {name}")
    cm = read_csv(inputs["genes.csv"])
    cm["identity_index"] = np.arange(len(cm))
    axes = ("G512", "G1024", "G2048", "G18164") if rerun else AXES
    sizes = tuple(fixture_sizes or ((512, 1024, 2048) if rerun else DEFAULT_SIZES))
    if not (0 < sizes[0] < sizes[1] < sizes[2] <= len(cm)):
        raise ValueError("invalid nested panel sizes")
    crosswalk = read_csv(inputs["gene_crosswalk.csv"])
    predictor = landmark / "predictor"
    landmarks = read_csv(hashed_payload(predictor, frozen, "genes.csv"))
    metadata = {source: source_metadata(landmark, frozen, manifests["evaluator"], source, materialized_fold)
                for source in ("tahoe", "tenv")}
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=f".{output.name}-") as temp, ExitStack() as stack:
        root = Path(temp) / "release"
        root.mkdir()
        datasets = {source: stack.enter_context(h5py.File(inputs[f"{source}.h5ad"], "r"))
                    for source in metadata}
        mappings = {source: gene_mapping(source, cm, read_elem(ds["var"]), crosswalk, frozen["fixture"])
                    for source, ds in datasets.items()}
        tahoe, tahoe_record = source_dispersion(
            datasets["tahoe"]["X"], metadata["tahoe"][0]["train"].source_row.to_numpy(),
            mappings["tahoe"].native_column.to_numpy(),
        )
        folds, fold_0 = {}, None
        for fold in range(3):
            train_path = hashed_payload(predictor, frozen, f"tenv/fold_{fold}/train.csv")
            train_rows = pd.read_csv(train_path).source_row.to_numpy()
            tenv, tenv_record = source_dispersion(datasets["tenv"]["X"], train_rows,
                                                   mappings["tenv"].native_column.to_numpy())
            panels = panel_tables(root / f"fold_{fold}", cm, landmarks, {"tahoe": tahoe, "tenv": tenv}, fold, sizes, axes)
            shared_count = len(panels["shared"])
            folds[str(fold)] = {"shared_gene_count": shared_count,
                                "shared_continuous_provisional": shared_count < 200,
                                "train_rank_cohorts": {"tahoe": tahoe_record, "tenv": tenv_record}}
            if fold == materialized_fold:
                fold_0 = panels

        view_roots = {name: root / "views" / name for name in fold_0}
        for name, view in view_roots.items():
            view_predictor = view / "predictor"
            view_predictor.mkdir(parents=True)
            shutil.copyfile(root / f"fold_{materialized_fold}/{name}.csv", view_predictor / "genes.csv")
            for relative in ("identity.npy", "identity_hdgfl3.npy", "cm_genes.csv", "tenv_folds.csv"):
                shutil.copyfile(hashed_payload(predictor, frozen, relative), view_predictor / relative)
        zeros = {}
        for source, ds in datasets.items():
            roots = {name: view for name, view in view_roots.items() if source == "tenv" or name != "shared"}
            zeros[source] = stream_views(ds["X"], source, mappings[source], roots, fold_0,
                                         metadata[source], landmark, manifests)
        views = manifest_record(root, fold_0, zeros, landmark_sha, frozen, materialized_fold)
        if starting_stats != fingerprint(inputs):
            raise ValueError("native input changed during panel build")
        write_json(root / "manifest.json", {
            "schema": "step09-gene-axes-v1", "fixture": frozen["fixture"], "ready_folds": [materialized_fold],
            "landmark_predictor_sha256": landmark_sha, "builder_sha256": sha256(Path(__file__)),
            "source_hashes": {name: frozen["source_hashes"][name] for name in inputs},
            "rank_recipe": {"source": "Seurat normalized dispersion on train-only linear CP10K",
                            "mean_bins": N_MEAN_BINS, "source_weights": {"tahoe": .5, "tenv": .5},
                            "tie_break": "frozen CM order"},
            "axis_sizes": {name: len(fold_0[name]) for name in axes},
            "folds": folds, "views": views,
            "panels_sha256": {str(path.relative_to(root)): sha256(path)
                              for path in sorted(root.glob("fold_*/*.csv"))},
        })
        root.rename(output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, default=DEFAULT_NATIVE)
    parser.add_argument("--landmark", type=Path, default=LANDMARK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--fixture-sizes", type=int, nargs=3)
    parser.add_argument("--rerun", action="store_true", help="Build the ADR0060 smaller gene axes")
    parser.add_argument("--fold", type=int, choices=(0, 1, 2), default=0)
    args = parser.parse_args()
    print(build(args.native, args.landmark, args.output, args.fixture_sizes, args.rerun, args.fold))


if __name__ == "__main__":
    main()
