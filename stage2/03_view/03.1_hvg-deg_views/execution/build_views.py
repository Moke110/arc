"""Build and deterministically replay the raw-first HVG-DEG view families."""
from __future__ import annotations

import argparse
import csv
import hashlib
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
from expression_view import (ARMS, DEFAULT, LAYERS, RECIPE, SCHEMA, SOLO_SEED,
                             layer_path, open_view, sha256)
from native_corpus import open_corpus
from prediction_axis import open_axes
from stage2_util.publication import publication_lock, publish_candidate

CONDITION = ["cell_line", "drug_parent_inchikey", "dose_um"]
FULL_FAMILIES = ("full", "full-pooled-label")


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


def _groups(obs, valid):
    positions = [index for index in valid if obs.role.iloc[index] == "treated"]
    conditions = obs.iloc[positions][CONDITION].copy()
    conditions["source_position"] = positions
    conditions["identity"] = conditions.index.astype(str)
    if conditions[CONDITION].isna().any().any() or conditions.drug_parent_inchikey.eq("").any():
        raise ValueError("a retained treated observation has no condition identity")
    conditions.sort_values(CONDITION + ["identity"], inplace=True)
    return {
        key: tuple(members.source_position)
        for key, members in conditions.groupby(CONDITION, sort=False, observed=True)
    }


def _identity(key):
    cell, parent, dose = key
    return f"{cell}|{parent}|{float(dose):g}"


def _solo(key, members, name):
    encoded = f"{SOLO_SEED}|{name}|{_identity(key)}".encode()
    number = int.from_bytes(hashlib.sha256(encoded).digest()[:8], "big")
    return members[number % len(members)]


def _plan(name, family, obs, valid, conditions, totals):
    if family in FULL_FAMILIES:
        return [(str(obs.index[index]), (index,)) for index in valid]
    if family == "merge":
        return [(f"merge::{name}::{_identity(key).replace('|', '::')}", rows)
                for key, rows in conditions.items()]
    if family == "solo":
        chosen = [_solo(key, rows, name) for key, rows in conditions.items()]
    else:
        chosen = [min(rows, key=lambda index: (-totals[index], str(obs.index[index])))
                  for rows in conditions.values()]
    return sorted((str(obs.index[index]), (index,)) for index in chosen)


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


def _control_union(members, controls):
    return tuple(sorted({position for member in members for position in controls[member]}))


def _frame(obs, plan, totals, merged):
    rows = obs.iloc[[members[0] for _, members in plan]].copy()
    rows.index = pd.Index([identity for identity, _ in plan], name="profile_id")
    if rows.index.has_duplicates:
        raise ValueError("expression view identities are not unique")
    if "profile_id" in rows:
        rows["profile_id"] = rows.index
    rows["library_total"] = [sum(int(totals[index]) for index in members) for _, members in plan]
    rows["member_count"] = [len(members) for _, members in plan]
    if merged:
        if "source_row" in rows:
            rows["source_row"] = -1
        for field in ("plate", "matched_control_group"):
            if field in rows:
                source = obs[field].astype(str).to_numpy()
                rows[field] = [";".join(sorted({source[index] for index in members}))
                               for _, members in plan]
    return rows


def _normalized(counts, members):
    raw = counts[list(members)].sum(axis=0, dtype=np.int64, keepdims=True)
    return normalization.log1p_cp10k(raw, normalization.library_totals(raw))[0].astype(np.float32)


def _write_arm(root, name, family, obs, plan, conditions, counts, totals, controls, tags, var):
    directory = root / family
    directory.mkdir(exist_ok=True)
    rows = _frame(obs, plan, totals, family == "merge")
    origins = {member: members for members in conditions.values() for member in members}
    pooled = family == "full-pooled-label"

    @lru_cache(maxsize=512)
    def norm_controls(control_members):
        if not control_members:
            raise ValueError(f"retained {name}/{family} profile has no matched vehicle")
        return _normalized(counts, control_members)

    @lru_cache(maxsize=512)
    def pooled_delta(members):
        return _normalized(counts, members) - norm_controls(_control_union(members, controls))

    with (directory / f"{name}.members.csv").open("w", newline="") as ledger, h5py.File(
        directory / f"{name}.h5ad", "w"
    ) as handle:
        writer = csv.writer(ledger, lineterminator="\n")
        writer.writerow(("profile_id", "member_profile_id", "role"))
        handle.attrs.update(tags, family=family, source=name)
        handle.attrs["encoding-type"], handle.attrs["encoding-version"] = "anndata", "0.1.0"
        write_elem(handle, "obs", rows)
        write_elem(handle, "var", var)
        layers = handle.create_group("layers")
        layers.attrs["encoding-type"], layers.attrs["encoding-version"] = "dict", "0.1.0"
        layer_names = ("delta",) if pooled else LAYERS
        datasets = {}
        width = len(var)
        for layer in layer_names:
            dataset = handle.create_dataset(layer_path(layer), (len(rows), width), dtype=np.float32,
                                            chunks=(min(64, len(rows)), width), compression="gzip",
                                            compression_opts=1)
            dataset.attrs["encoding-type"], dataset.attrs["encoding-version"] = "array", "0.2.0"
            datasets[layer] = dataset
        for start in range(0, len(plan), 64):
            block = plan[start:start + 64]
            arrays = {layer: np.empty((len(block), width), dtype=np.float32) for layer in layer_names}
            for offset, (identity, members) in enumerate(block):
                vehicle = _control_union(members, controls)
                label = origins.get(members[0], members) if pooled else members
                if pooled:
                    arrays["delta"][offset] = pooled_delta(label)
                else:
                    arrays["X"][offset] = _normalized(counts, members)
                    arrays["control"][offset] = norm_controls(vehicle)
                    arrays["delta"][offset] = arrays["X"][offset] - arrays["control"][offset]
                writer.writerows((identity, obs.index[index], "expression") for index in members)
                writer.writerows((identity, obs.index[index], "matched_vehicle") for index in vehicle)
                if pooled:
                    writer.writerows((identity, obs.index[index], "label_expression") for index in label)
                    writer.writerows((identity, obs.index[index], "label_vehicle")
                                     for index in _control_union(label, controls))
            for layer, dataset in datasets.items():
                dataset[start:start + len(block)] = arrays[layer]


def _build_source(root, name, corpus, axis, tags, var):
    projection = axis.project(corpus, name)
    obs, valid, excluded = _population(projection.source)
    controls = _controls(projection.source, obs, valid)
    conditions = _groups(obs, valid)
    if not conditions:
        raise ValueError(f"{name} has no conditions eligible for paired views")
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
        for family in ARMS[name]:
            plan = _plan(name, family, obs, valid, conditions, totals)
            _write_arm(root, name, family, obs, plan, conditions, counts, totals, controls, tags, var)
    finally:
        del counts
        cache.unlink()


def _build(root, corpus, axis, sources):
    tags = {"schema": SCHEMA, "corpus_sha256": corpus.manifest_sha256, "axis_sha256": axis.sha256}
    var = pd.DataFrame({"gene_index": axis.gene_indices, "symbol": axis.symbols},
                       index=pd.Index([str(gene) for gene in axis.gene_indices], name="gene_index_label"))
    for name in sources:
        _build_source(root, name, corpus, axis, tags, var)
    manifest = {**tags, **RECIPE, "axis_size": len(var), "sources": list(sources),
                "outputs": {path.relative_to(root).as_posix(): sha256(path)
                            for path in sorted(root.rglob("*")) if path.is_file()}}
    release.write_json(root / "manifest.json", manifest)


def _verify(root, corpus, axis, sources):
    manifest = release.read_json(root / "manifest.json")
    release.check_relative_paths(manifest)
    for name in sources:
        for family in ARMS[name]:
            with open_view(root, source=name, family=family, corpus_sha256=corpus.manifest_sha256,
                           axis_sha256=axis.sha256) as view:
                if tuple(view.var.gene_index) != axis.gene_indices:
                    raise ValueError("expression view axis order differs")
            name_in_manifest = f"{family}/{name}.members.csv"
            if sha256(root / name_in_manifest) != manifest["outputs"][name_in_manifest]:
                raise ValueError("expression view membership digest differs")


def publish_views(destination, corpus, axis, *, sources=tuple(ARMS), staging=None):
    if not sources or len(set(sources)) != len(sources) or not set(sources).issubset(ARMS):
        raise ValueError("unknown or duplicate expression sources")
    with publication_lock(destination) as output:
        return publish_candidate(
            release, output, lambda root: _build(root, corpus, axis, sources),
            lambda root: _verify(root, corpus, axis, sources), staging=staging,
            verify_existing=lambda root: release.check_manifest(root, release.read_json(root / "manifest.json")),
        )


def replay_views(destination, corpus, axis, *, sources=tuple(ARMS), staging=None):
    with publication_lock(destination) as output:
        return release.replay(output, lambda root: _build(root, corpus, axis, sources),
                              lambda root: _verify(root, corpus, axis, sources), root=staging)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--corpus", type=Path, default=STAGE / "01_corpus/output/native")
    parser.add_argument("--axis", type=Path)
    parser.add_argument("--destination", type=Path, default=DEFAULT)
    parser.add_argument("--staging", type=Path, default=STEP / "execution/checkpoints")
    parser.add_argument("--source", choices=tuple(ARMS), action="append")
    args = parser.parse_args()
    sources = tuple(ARMS) if not args.source else tuple(args.source)
    if sources != tuple(ARMS) and args.destination.resolve() == DEFAULT.resolve():
        parser.error("single-source rehearsals require an explicit non-default destination")
    operation = publish_views if args.action == "publish" else replay_views
    with open_corpus(args.corpus) as corpus:
        result = operation(args.destination, corpus, open_axes(args.axis).axis(5000),
                           sources=sources, staging=args.staging)
    print(json.dumps(result) if isinstance(result, dict) else result)


if __name__ == "__main__":
    main()
