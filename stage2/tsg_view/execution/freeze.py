#!/usr/bin/env python3
"""Freeze one view product into the Step 04-shaped `tsg-view-freeze-v1` format.

The freeze keeps the Step 04 file set (plus ``control_sources.csv``) but indexes
``expression.npy``/``control.npy`` by population row and carries no duplicated vehicle
tail (plan C7.6). Partitions come only from the sealed Step 04 ``(parent, dose)`` mapping,
so a view cannot move a condition between folds; for ``full`` the ``splits.csv`` is
byte-equal to the sealed freeze, and every view is cross-checked against the sealed
parent/dose population, parents and Morgan recipe.
"""
import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from views import (BATCH, CONDITION_KEYS, SCHEMA as PRODUCT_SCHEMA, HarnessError,  # noqa: E402
                   SOURCES, VIEWS, bind, relative, repo_root, sha256, write_json)
from inputs import (COMPANION, FREEZE_SCHEMA, OUTPUTS, POPULATION_COLUMNS, check_splits,
                    read_sealed_splits, sealed_binding, sealed_reference)
from inputs import Inputs  # noqa: E402

CONTROL_RULE = {
    "full": "mean-of-normalized-vehicles",
    "solo": "mean-of-normalized-vehicles",
    "solo-deepest": "mean-of-normalized-vehicles",
    "merge": "pooled-sum",
    "full-pooled-label": "full-view-mean-of-normalized-vehicles",
}


def _read_rows(dataset, positions):
    positions = np.asarray(positions, dtype=np.int64)
    order = np.argsort(positions, kind="stable")
    block = np.asarray(dataset[list(positions[order])], dtype=np.float32)
    restored = np.empty_like(block)
    restored[order] = block
    return restored


class _Product:
    """The product row space plus a bounded expression/control provider."""

    def __init__(self, manifest, manifest_path, views_root, root, accepted_dir=None):
        import h5py
        from anndata.io import read_elem

        self.manifest = manifest
        self.manifest_path = Path(manifest_path)
        self.mode = "delta" if manifest.get("payload", {}) and manifest["payload"].get(
            "kind") == "delta" else ("payload" if manifest.get("payload") else "alias")
        self.root = root
        accepted_manifest = (Path(accepted_dir) if accepted_dir else Path(root)) / (
            "manifest.json" if accepted_dir else
            manifest["source_manifests"]["accepted_view"])
        self.accepted = accepted_manifest.parent
        if not accepted_manifest.is_file():
            raise HarnessError(
                f"the accepted release is unavailable: {accepted_manifest}; freeze where "
                "the payload exists (the main checkout) or pass --repo-root")
        if self.mode == "payload":
            self.payload_path = Path(views_root) / manifest["view"] / manifest["payload"]["file"]
            if not self.payload_path.is_file():
                raise HarnessError(f"product payload is missing: {self.payload_path}")
            with h5py.File(self.payload_path, "r") as handle:
                obs = read_elem(handle["obs"])
        else:
            view = bind(root)["view"]
            with view.open_view(self.accepted) as release:
                with release.source(manifest["source"]) as source:
                    obs = source.obs.copy()
            self.payload_path = self.accepted / f"{manifest['source']}.h5ad"
        obs = pd.DataFrame(obs)
        obs["profile_id"] = obs.index.astype(str) if "profile_id" not in obs else obs[
            "profile_id"].astype(str)
        self.rows = obs.reset_index(drop=True)
        self.rows["profile_id"] = self.rows["profile_id"].astype(str)
        self.rows["view_row"] = np.arange(len(self.rows), dtype=np.int64)
        missing = [name for name in POPULATION_COLUMNS if name not in self.rows]
        if missing:
            raise HarnessError(f"product rows lack the freeze columns: {missing}")
        self.population_positions = np.flatnonzero(
            self.rows["role"].astype(str).to_numpy() == "treated")
        if not len(self.population_positions):
            raise HarnessError("a view product carries no treated row")
        self.delta = None
        if self.mode == "delta":
            self.delta = np.load(Path(views_root) / manifest["view"]
                                 / manifest["payload"]["file"], mmap_mode="r")
            if self.delta.shape != (len(self.rows), int(manifest["payload"]["cols"])):
                raise HarnessError("the delta payload shape differs from the product manifest")

    def __call__(self, start, stop):
        """Population rows ``[start, stop)`` as ``(expression, control)`` float32."""
        import h5py

        positions = self.population_positions[start:stop]
        with h5py.File(self.payload_path, "r") as handle:
            control = _read_rows(handle["layers/control"], positions)
            if self.mode == "delta":
                delta = np.asarray(self.delta[positions], dtype=np.float64)
                expression = (control.astype(np.float64) + delta).astype(np.float32)
            else:
                expression = _read_rows(handle["X"], positions)
        return expression, control


def _groups(manifest, population, views_root):
    """Per-population-row matched control group lists."""
    source_record = manifest.get("control_sources")
    if source_record:
        table = pd.read_csv(Path(views_root) / manifest["view"] / source_record["file"])
        lookup = {str(profile): [group for group in str(groups).split(";") if group]
                  for profile, groups in zip(table["profile_id"], table["control_groups"])}
        groups = [lookup[str(profile)] for profile in population["profile_id"]]
        totals = {str(profile): total for profile, total in
                  zip(table["profile_id"], table["pooled_retained_total"])}
    else:
        groups = [[group for group in str(value).split(";") if group]
                  for value in population["matched_control_group"]]
        totals = {}
    if any(not group for group in groups):
        raise HarnessError("a population row has no matched control group")
    return groups, totals


def freeze_product(*, view, source, views_root, destination, repo_root_path=None,
                   sealed=None, axis_dir=None, accepted_dir=None, record=None):
    """Write and self-verify one view freeze; return its receipt."""
    if view not in VIEWS or source not in SOURCES:
        raise HarnessError(f"unknown view freeze identity: {view}/{source}")
    root = repo_root(repo_root_path)
    modules = bind(root)
    views_root, destination = Path(views_root), Path(destination)
    manifest_path = views_root / view / f"{source}.json"
    if not manifest_path.is_file():
        raise HarnessError(f"view product manifest is missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if (manifest.get("schema") != PRODUCT_SCHEMA or manifest.get("view") != view
            or manifest.get("source") != source):
        raise HarnessError("view product manifest identity mismatch")
    if destination.exists() and any(destination.iterdir()):
        raise HarnessError("a view freeze requires a new empty candidate")
    sealed_dir = Path(sealed) if sealed else sealed_reference(root, source)
    binding = sealed_binding(sealed_dir, root)

    product = _Product(manifest, manifest_path, views_root, root,
                       accepted_dir=accepted_dir)
    population = product.rows.iloc[product.population_positions][list(
        POPULATION_COLUMNS)].copy()
    population = population.reset_index(drop=True)
    population.index = pd.Index(population["profile_id"].astype(str), name="profile_id")
    groups, pooled_totals = _groups(manifest, population, views_root)

    sealed_splits = read_sealed_splits(sealed_dir)
    protocol = modules["protocol"]
    splits = protocol.make_splits(population)
    check_splits(splits, population, sealed_splits, source)

    # -- parents and the Morgan recipe come from the sealed freeze ---------- #
    parents = population[["drug_parent_inchikey", "drug_parent_smiles"]].drop_duplicates(
        ).sort_values("drug_parent_inchikey").reset_index(drop=True)
    sealed_parents = pd.read_csv(sealed_dir / "parents.csv")
    if not parents.equals(sealed_parents):
        raise HarnessError(
            "the view parent population differs from the sealed Step 04 freeze")

    # -- control-source vehicle metadata ----------------------------------- #
    view_api = modules["view"]
    with view_api.open_view(product.accepted) as release:
        with release.source(source) as accepted_source:
            accepted_obs = accepted_source.obs.copy()
            control_groups = accepted_source.control_groups()
    accepted_obs["profile_id"] = accepted_obs.index.astype(str)
    accepted_obs["view_row"] = np.arange(len(accepted_obs), dtype=np.int64)
    referenced = sorted({group for row_groups in groups for group in row_groups})
    missing = [group for group in referenced if group not in control_groups]
    if missing:
        raise HarnessError(f"matched control groups absent from the accepted view: {missing[:3]}")
    vehicle_positions = np.unique(np.concatenate(
        [np.asarray(control_groups[group], dtype=np.int64) for group in referenced]))
    tail = accepted_obs.iloc[vehicle_positions]
    if not tail["role"].eq("vehicle").all():
        raise HarnessError("matched control resolution picked a non-vehicle row")
    tail_index = {int(position): len(population) + offset
                  for offset, position in enumerate(vehicle_positions)}
    matched = {}
    for row, row_groups in enumerate(groups):
        offsets = sorted({tail_index[int(position)] for group in row_groups
                          for position in control_groups[group]})
        matched[row] = offsets

    # -- write the freeze --------------------------------------------------- #
    destination.mkdir(parents=True, exist_ok=True)
    population.to_csv(destination / "population.csv", index=False, lineterminator="\n")
    pd.concat([population[list(POPULATION_COLUMNS)], tail[list(POPULATION_COLUMNS)]],
              ignore_index=True).to_csv(destination / "all_obs.csv", index=False,
                                        lineterminator="\n")
    splits.to_csv(destination / "splits.csv", index=False, lineterminator="\n")
    axis_dir = Path(axis_dir) if axis_dir else (
        Path(root) / Path(manifest["axis"]["file"]).parent)
    frame = modules["gene_axis"].GeneAxis(axis_dir).axis_frame(5000)
    frame[["gene_index", "symbol"]].to_csv(destination / "genes.csv", index=False,
                                           lineterminator="\n")
    shutil.copyfile(sealed_dir / "parents.csv", destination / "parents.csv")
    shutil.copyfile(sealed_dir / "morgan.npy", destination / "morgan.npy")
    write_json(destination / "matched_controls.json", matched)
    controls = pd.DataFrame({
        "profile_id": population["profile_id"].astype(str),
        "control_rule": CONTROL_RULE[view],
        "control_groups": [";".join(sorted(set(row_groups))) for row_groups in groups],
        "pooled_retained_total": [pooled_totals.get(str(profile), "")
                                  for profile in population["profile_id"]],
    })
    controls.to_csv(destination / COMPANION, index=False, lineterminator="\n")

    rows = len(population)
    expression = np.lib.format.open_memmap(destination / "expression.npy", mode="w+",
                                           dtype=np.float32, shape=(rows, 5000))
    control = np.lib.format.open_memmap(destination / "control.npy", mode="w+",
                                        dtype=np.float32, shape=(rows, 5000))
    for start in range(0, rows, BATCH):
        stop = min(start + BATCH, rows)
        x, c = product(start, stop)
        if x.shape != (stop - start, 5000) or c.shape != (stop - start, 5000):
            raise HarnessError("the product provider returned the wrong shape")
        expression[start:stop] = x
        control[start:stop] = c
    expression.flush()
    control.flush()
    del expression, control

    generation = {
        "schema": manifest["schema"],
        "product": manifest["product"],
        "view": manifest["view"],
        "source": manifest["source"],
        "recipe": manifest["recipe"],
        "normalization": manifest["normalization"],
        "control": {"rule": manifest["control_rule"], "delta": "expression - control"},
        "axis": manifest["axis"],
    }
    config = {
        "schema": FREEZE_SCHEMA,
        "source": source,
        "view": view,
        "view_generation": generation,
        "view_product": {"file": relative(manifest_path, root),
                         "sha256": sha256(manifest_path),
                         "mode": product.mode},
        "sealed": binding,
        "delta_kind": "pooled" if view == "full-pooled-label"
        else "expression_minus_control",
        "control_rule": CONTROL_RULE[view],
        "split_seed": int(protocol.SPLIT_SEED),
        "training_seed": int(protocol.TRAINING_SEED),
        "n_treated": int(rows),
        "n_conditions": int(population[list(CONDITION_KEYS)].drop_duplicates().shape[0]),
        "n_parents": int(len(parents)),
        "n_control_vehicles": int(len(vehicle_positions)),
        "morgan": {"source": "sealed Step 04 freeze", "radius": 2, "bits": 2048,
                   "include_chirality": True},
        "members": manifest.get("members"),
        "preparation_sha256": sha256(Path(__file__)),
    }
    write_json(destination / "config.json", config)
    outputs = {name: sha256(destination / name) for name in sorted(OUTPUTS)}
    write_json(destination / "manifest.json", {
        "schema": FREEZE_SCHEMA,
        "outputs": outputs,
        "view_product": config["view_product"],
        "sealed": binding,
        "delta_kind": config["delta_kind"],
    })
    verified = Inputs(destination, sealed=sealed_dir, repo_root=root)
    receipt = {
        "schema": "tsg-view-freeze-receipt-v1",
        "view": view,
        "source": source,
        "destination": relative(destination, root),
        "manifest_sha256": sha256(destination / "manifest.json"),
        "outputs": outputs,
        "rows": int(rows),
        "conditions": config["n_conditions"],
        "splits_sha256": outputs["splits.csv"],
        "sealed": binding,
        "splits_byte_equal_to_sealed": outputs["splits.csv"] == binding["splits_sha256"],
        "delta_kind": config["delta_kind"],
        "verified": True,
        "verified_rows": int(len(verified.population)),
    }
    if record is not None:
        write_json(record, receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--view", required=True, choices=VIEWS)
    parser.add_argument("--source", required=True, choices=SOURCES)
    parser.add_argument("--views", type=Path, default=HERE.parent / "output/views",
                        help="views root holding <view>/<source>.json")
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=None)
    parser.add_argument("--sealed", type=Path, default=None,
                        help="sealed Step 04 freeze directory for this source")
    parser.add_argument("--axis", type=Path, default=None,
                        help="gene-axis release directory (default: from the product manifest)")
    parser.add_argument("--accepted", type=Path, default=None,
                        help="accepted view-5k release directory (default: from the manifest)")
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    receipt = freeze_product(view=args.view, source=args.source, views_root=args.views,
                            destination=args.destination,
                            repo_root_path=args.repo_root, sealed=args.sealed,
                            axis_dir=args.axis, accepted_dir=args.accepted,
                            record=args.record)
    sys.stdout.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
