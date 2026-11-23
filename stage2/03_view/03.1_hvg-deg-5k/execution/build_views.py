"""Build and deterministically replay the raw-first HVG-DEG full 5k views."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from functools import lru_cache
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem

STEP = Path(__file__).resolve().parents[1]
STAGE = STEP.parents[1]
sys.path[:0] = [str(path) for path in (
    STEP / "output/api",
    STAGE / "util/output/api",
    STAGE / "01_corpus/output/api",
    STAGE / "02_gene_axis/02.1_hvg-deg/output/api",
)]

import normalization
from stage2_util.derivation import release
from expression_view import (DEFAULT, LAYERS, PAYLOAD_SCHEMA, RECIPE, SCHEMA, SOURCES,
                             check_inventory, layer_path, open_view, sha256)
from native_corpus import open_corpus
from prediction_axis import open_axes
from stage2_util.publication import publication_lock, publish_candidate

CONDITION = ["cell_line", "drug_parent_inchikey", "dose_um"]


def _population(source):
    obs = source.obs.fillna({"view_exclusion_reason": ""})
    valid = []
    excluded = []
    for index, row in enumerate(obs.itertuples()):
        reason = row.view_exclusion_reason
        if not reason and row.role == "treated" and row.match_status != "matched":
            reason = "no_matched_vehicle"
        if not reason and row.role not in ("treated", "vehicle"):
            reason = f"role_{row.role}"
        if reason:
            excluded.append((str(row.Index), row.role, reason))
        else:
            valid.append(index)
    return obs, valid, excluded


def _check_conditions(name, obs, valid):
    treated = obs.iloc[valid].loc[lambda rows: rows.role == "treated", CONDITION]
    if treated.empty:
        raise ValueError(f"{name} has no conditions eligible for paired views")
    if treated.isna().any().any() or treated.drug_parent_inchikey.eq("").any():
        raise ValueError("a retained treated observation has no condition identity")


def _controls(source, obs, valid):
    groups = {}
    for index in valid:
        if obs.role.iloc[index] == "vehicle":
            group = obs.matched_control_group.iloc[index]
            groups.setdefault(group, []).append(index)
    accepted = set(valid)
    return {
        index: tuple(row for row in source.matched_control_rows(index) if row in accepted)
        if obs.role.iloc[index] == "treated"
        else tuple(groups[obs.matched_control_group.iloc[index]])
        for index in valid
    }


def _frame(obs, positions, totals):
    rows = obs.iloc[positions].copy()
    rows.index = pd.Index([str(identity) for identity in rows.index], name="profile_id")
    if rows.index.has_duplicates:
        raise ValueError("expression view identities are not unique")
    if "profile_id" in rows:
        rows["profile_id"] = rows.index
    rows["library_total"] = totals[positions]
    rows["member_count"] = 1
    return rows


def _normalized(counts, members):
    raw = counts[list(members)].sum(axis=0, dtype=np.int64, keepdims=True)
    return normalization.log1p_cp10k(raw, normalization.library_totals(raw))[0].astype(np.float32)


def _write_full(root, name, obs, positions, counts, totals, controls, tags, var):
    rows = _frame(obs, positions, totals)

    @lru_cache(maxsize=512)
    def norm_controls(control_members):
        if not control_members:
            raise ValueError(f"retained {name}/full profile has no matched vehicle")
        return _normalized(counts, control_members)

    with (root / f"{name}.members.csv").open("w", newline="") as ledger, h5py.File(
        root / f"{name}.h5ad", "w"
    ) as handle:
        writer = csv.writer(ledger, lineterminator="\n")
        writer.writerow(("profile_id", "member_profile_id", "role"))
        handle.attrs.update(tags, family="full", source=name)
        handle.attrs["encoding-type"], handle.attrs["encoding-version"] = "anndata", "0.1.0"
        write_elem(handle, "obs", rows)
        write_elem(handle, "var", var)
        layers = handle.create_group("layers")
        layers.attrs["encoding-type"], layers.attrs["encoding-version"] = "dict", "0.1.0"
        datasets = {}
        width = len(var)
        for layer in LAYERS:
            dataset = handle.create_dataset(layer_path(layer), (len(rows), width), dtype=np.float32,
                                            chunks=(min(64, len(rows)), width), compression="gzip",
                                            compression_opts=1)
            dataset.attrs["encoding-type"], dataset.attrs["encoding-version"] = "array", "0.2.0"
            datasets[layer] = dataset
        for start in range(0, len(positions), 64):
            block = positions[start:start + 64]
            arrays = {layer: np.empty((len(block), width), dtype=np.float32) for layer in LAYERS}
            for offset, position in enumerate(block):
                identity = obs.index[position]
                vehicle = tuple(sorted(set(controls[position])))
                arrays["X"][offset] = _normalized(counts, (position,))
                arrays["control"][offset] = norm_controls(vehicle)
                arrays["delta"][offset] = arrays["X"][offset] - arrays["control"][offset]
                writer.writerow((identity, identity, "expression"))
                writer.writerows((identity, obs.index[index], "matched_vehicle") for index in vehicle)
            for layer, dataset in datasets.items():
                dataset[start:start + len(block)] = arrays[layer]


def _build_source(root, name, corpus, axis, tags, var):
    projection = axis.project(corpus, name)
    obs, valid, excluded = _population(projection.source)
    controls = _controls(projection.source, obs, valid)
    _check_conditions(name, obs, valid)
    pd.DataFrame(excluded, columns=["profile_id", "role", "reason"]).to_csv(
        root / f"{name}.exclusions.csv", index=False
    )
    cache = root / "selected-counts.npy"
    counts = np.lib.format.open_memmap(cache, mode="w+", dtype=np.int64, shape=(len(obs), len(var)))
    try:
        totals = np.empty(len(obs), dtype=np.int64)
        for start, block in projection.batches():
            stop = start + len(block)
            counts[start:stop] = block
            totals[start:stop] = block.sum(axis=1, dtype=np.int64)
        _write_full(root, name, obs, valid, counts, totals, controls, tags, var)
    finally:
        del counts
        cache.unlink()


def _build(root, corpus, axis, sources, corpus_sha256):
    tags = {"schema": PAYLOAD_SCHEMA, "corpus_sha256": corpus_sha256, "axis_sha256": axis.sha256}
    var = pd.DataFrame({"gene_index": axis.gene_indices, "symbol": axis.symbols},
                       index=pd.Index([str(gene) for gene in axis.gene_indices], name="gene_index_label"))
    for name in sources:
        _build_source(root, name, corpus, axis, tags, var)
    manifest = {**tags, **RECIPE, "schema": SCHEMA, "axis_size": len(var), "sources": list(sources),
                "outputs": {path.relative_to(root).as_posix(): sha256(path)
                            for path in sorted(root.rglob("*")) if path.is_file()}}
    release.write_json(root / "manifest.json", manifest)


def _verify(root, axis, sources, corpus_sha256):
    manifest = release.read_json(root / "manifest.json")
    release.check_relative_paths(manifest)
    if tuple(check_inventory(manifest)) != tuple(sources):
        raise ValueError("expression view release sources differ")
    if {path.name for path in root.iterdir()} != {"manifest.json", *manifest["outputs"]}:
        raise ValueError("expression view release inventory differs")
    for name in sources:
        with open_view(root, source=name, corpus_sha256=corpus_sha256,
                       axis_sha256=axis.sha256) as view:
            if tuple(view.var.gene_index) != axis.gene_indices:
                raise ValueError("expression view axis order differs")
        membership = f"{name}.members.csv"
        if sha256(root / membership) != manifest["outputs"][membership]:
            raise ValueError("expression view membership digest differs")


def publish_views(destination, corpus, axis, *, sources=SOURCES, staging=None):
    if not sources or len(set(sources)) != len(sources) or not set(sources).issubset(SOURCES):
        raise ValueError("unknown or duplicate expression sources")
    with publication_lock(destination) as output:
        return publish_candidate(
            release, output, lambda root: _build(root, corpus, axis, sources, corpus.manifest_sha256),
            lambda root: _verify(root, axis, sources, corpus.manifest_sha256), staging=staging,
            verify_existing=lambda root: release.check_manifest(root, release.read_json(root / "manifest.json")),
        )


def replay_views(destination, corpus, axis, *, sources=SOURCES, staging=None,
                 historical_manifest=None):
    with publication_lock(destination) as output:
        selected = release.read_json(output / "manifest.json")["corpus_sha256"]
        verified = corpus.resolve_generation(selected, sources, historical_manifest=historical_manifest)
        return release.replay(output, lambda root: _build(root, corpus, axis, sources, verified),
                              lambda root: _verify(root, axis, sources, verified), root=staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--corpus", type=Path, default=STAGE / "01_corpus/output")
    parser.add_argument("--axis", type=Path)
    parser.add_argument("--historical-corpus", type=Path,
                        help="historical v1 corpus directory for replay of a sealed view")
    parser.add_argument("--destination", type=Path, default=DEFAULT)
    parser.add_argument("--staging", type=Path, default=STEP / "execution/checkpoints")
    parser.add_argument("--source", choices=SOURCES, action="append")
    args = parser.parse_args()
    sources = SOURCES if not args.source else tuple(args.source)
    if sources != SOURCES and args.destination.resolve() == DEFAULT.resolve():
        parser.error("single-source rehearsals require an explicit non-default destination")
    operation = publish_views if args.action == "publish" else replay_views
    options = {"sources": sources, "staging": args.staging}
    if args.historical_corpus is not None:
        if args.action != "verify":
            parser.error("--historical-corpus is only used to replay sealed views")
        options["historical_manifest"] = args.historical_corpus / "manifest.json"
    with open_corpus(args.corpus) as corpus:
        result = operation(args.destination, corpus, open_axes(args.axis).axis(5000), **options)
    print(json.dumps(result) if isinstance(result, dict) else result)


if __name__ == "__main__":
    main()
