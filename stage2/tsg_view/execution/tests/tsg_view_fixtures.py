"""Hermetic, multi-condition fixtures for the `tsg_view` tests.

The Base step's own fixture writer can only express one condition per source, so this
module writes a contract-shaped integer corpus with explicit per-row identity columns
(parent, dose, plate, cell line, compound) while keeping the same h5ad layout and manifest
schema. Every count is a literal in this file, so a test can hand-compute an expected
normalization from the fixture instead of re-running the harness.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem

AXIS_SIZE = 18164
DEFAULT_PARENT = {"tenv": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "tahoe_bulk": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"}
DEFAULT_SMILES = {"tenv": "CCO", "tahoe_bulk": "CCO"}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _row(name, index, spec):
    treated = spec["role"] == "treated"
    return {
        "profile_id": f"{name}:p{index:02d}",
        "source_item": "fixture",
        "source_item_manifest_sha256": "0" * 64,
        "source_row": index,
        "role": spec["role"],
        "plate": spec["plate"],
        "cell_line": spec["cell_line"],
        "cell_line_id": spec["cell_line_id"],
        "compound_name": spec.get("compound", "compound" if treated else "vehicle"),
        "vehicle": spec.get("vehicle", "" if treated else "DMSO"),
        "dose_um": float(spec.get("dose", 0.0)),
        "dose_original_value": float(spec.get("dose", 0.0)),
        "dose_original_unit": "uM",
        "drug_parent_inchikey": spec.get("parent", "") if treated else "",
        "drug_parent_smiles": spec.get("smiles", "") if treated else "",
        "source_inchikey": spec.get("parent", "") if treated else "",
        "n_cells": 0,
        "matched_control_group": spec["group"],
    }


def write_corpus(root, sources):
    """``root/corpus``: one h5ad per source plus a contract-shaped manifest."""
    from bases import corpus as corpus_module

    release = Path(root) / "corpus"
    release.mkdir(parents=True, exist_ok=True)
    outputs, summaries, frames = {}, {}, {}
    axis_size = None
    for name in sorted(sources):
        entry = sources[name]
        counts = np.asarray(entry["counts"], dtype=np.int32)
        if axis_size is None:
            axis_size = counts.shape[1]
        elif counts.shape[1] != axis_size:
            raise ValueError("fixture sources disagree on the frozen axis size")
        rows = [_row(name, index, spec) for index, spec in enumerate(entry["identity"])]
        obs = pd.DataFrame(rows).set_index("profile_id", drop=False)
        obs.index.name = "profile_id"
        if len(obs) != len(counts):
            raise ValueError("fixture identity and counts disagree on row count")
        var = pd.DataFrame({"gene_index": np.arange(counts.shape[1], dtype=np.int64),
                            "symbol": [f"S{gene}" for gene in range(counts.shape[1])]})
        var.index = pd.Index([str(gene) for gene in var.gene_index],
                             name="gene_index_label")
        path = release / f"{name}.h5ad"
        with h5py.File(path, "w") as handle:
            handle.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
            write_elem(handle, "obs", obs)
            write_elem(handle, "var", var)
            write_elem(handle, "uns", {"schema": corpus_module.SCHEMA, "source": name})
            write_elem(handle, "layers", {})
            matrix = handle.create_dataset("X", data=counts)
            matrix.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})
        frames[name] = (obs, counts)
        roles = obs.role.to_numpy()
        groups = obs.matched_control_group.to_numpy()
        treated = roles == "treated"
        outputs[f"{name}.h5ad"] = {"sha256": sha256(path), "bytes": path.stat().st_size,
                                   "rows": int(len(obs)), "columns": int(axis_size),
                                   "dtype": "int32"}
        summaries[name] = {
            "counts": {"profiles": int(len(obs)), "treated": int(treated.sum()),
                       "vehicle": int((~treated).sum()), "raw_profiles": int(len(obs))},
            "excluded_by_reason": {},
            "controls": {"rule": "fixture",
                         "vehicle_rows": int((~treated).sum()),
                         "vehicle_groups": int(len({str(group) for group in groups[~treated]})),
                         "matched_control_rows": int((~treated).sum()),
                         "matched_groups": int(len({str(group) for group in groups[treated]})),
                         "treated_without_control": 0},
            "output": f"{name}.h5ad",
        }
    pd.DataFrame(columns=["source", "source_row", "profile_id", "reason"]).to_csv(
        release / "exclusions.csv", index=False, lineterminator="\n")
    manifest = {
        "schema": corpus_module.SCHEMA,
        "fixture": True,
        "axis": {"size": int(axis_size), "path": "genes.csv", "sha256": "0" * 64},
        "outputs": outputs,
        "sources": summaries,
        "controls": {name: summaries[name]["controls"] for name in summaries},
        "counts": {"axis_size": int(axis_size),
                   "tenv_profiles": summaries["tenv"]["counts"]["profiles"],
                   "tahoe_bulk_profiles": summaries["tahoe_bulk"]["counts"]["profiles"]},
    }
    (release / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return release


def write_axis(root, axes):
    """``root/gene-axis``: one ``gene_axis_<N>.csv`` per size plus a manifest."""
    from bases import gene_axis as axis_module

    release = Path(root) / "gene-axis"
    release.mkdir(parents=True, exist_ok=True)
    outputs, sizes = {}, {}
    for size, genes in sorted(axes.items()):
        genes = np.asarray(genes, dtype=np.int64)
        frame = pd.DataFrame({"gene_index": genes,
                              "symbol": [f"S{gene}" for gene in genes]})
        path = release / f"gene_axis_{int(size)}.csv"
        frame.to_csv(path, index=False, lineterminator="\n")
        outputs[f"stage2/01_bases/output/gene-axis/gene_axis_{int(size)}.csv"] = {
            "sha256": sha256(path), "rows": int(len(frame))}
        sizes[str(int(size))] = {"n": int(len(frame)),
                                 "file": f"gene_axis_{int(size)}.csv"}
    manifest = {"schema": axis_module.SCHEMA, "fixture": True, "axes": sizes,
                "outputs": outputs,
                "release": {"algorithm": axis_module.ALGORITHM,
                            "axis_sizes": [int(size) for size in sorted(axes)]}}
    (release / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return release


def write_release(root, sources, axes):
    return {"corpus": write_corpus(root, sources), "axis": write_axis(root, axes)}


def _counts(index, scale):
    """Six-gene fixture counts with a per-row ratio, plus large off-axis counts.

    Members of one condition deliberately carry different selected-gene proportions, so a
    normalize-of-sum pool can never silently equal a mean of normalized member values.
    """
    signal_a = 3 + (index * 3 + scale) % 7
    signal_b = 4 + (index * 5 + 2 * scale) % 11
    return [signal_a, 100, 100, signal_b, 100, 100]


def grid_sources():
    """A 6-gene, 9-condition (tenv) / 10-condition (tahoe_bulk) fixture.

    Tenv carries two same-plate 2-member conditions, a cross-plate 2-member condition and
    two singleton conditions; tahoe_bulk mixes singletons with one cross-plate pair and one
    second-cell-line condition. Vehicle wells in one group have different totals, so a mean
    of normalized controls can never equal a normalize-of-sum pool.
    """
    tenv_identity = []
    tenv_counts = []
    grid = [("P1", 1.0, "Plate1"), ("P1", 1.0, "Plate1"),
            ("P1", 10.0, "Plate1"), ("P1", 10.0, "Plate1"),
            ("P1", 100.0, "Plate2"), ("P1", 100.0, "Plate1"),
            ("P2", 1.0, "Plate1"), ("P2", 1.0, "Plate1"),
            ("P2", 10.0, "Plate2"), ("P2", 10.0, "Plate2"),
            ("P2", 100.0, "Plate1"), ("P2", 100.0, "Plate1"),
            ("P3", 1.0, "Plate2"),
            ("P3", 10.0, "Plate1"), ("P3", 10.0, "Plate1"),
            ("P3", 100.0, "Plate2")]
    for index, (parent, dose, plate) in enumerate(grid):
        tenv_identity.append({"role": "treated", "plate": plate, "cell_line": "hiPSC-CM",
                              "cell_line_id": "", "dose": dose, "parent": parent,
                              "smiles": DEFAULT_SMILES["tenv"],
                              "group": f"tenv::{plate}::DMSO"})
        tenv_counts.append(_counts(index, 1))
    for index, (plate, counts) in enumerate(
            [("Plate1", [6, 100, 100, 9, 100, 100]),
             ("Plate1", [18, 100, 100, 3, 100, 100]),
             ("Plate2", [3, 100, 100, 12, 100, 100]),
             ("Plate2", [9, 100, 100, 7, 100, 100])]):
        tenv_identity.append({"role": "vehicle", "plate": plate, "cell_line": "hiPSC-CM",
                              "cell_line_id": "", "group": f"tenv::{plate}::DMSO"})
        tenv_counts.append(counts)

    tahoe_grid = [("T1", 0.05, "pl1", "L1"), ("T1", 0.05, "pl2", "L1"),
                  ("T1", 0.5, "pl1", "L1"), ("T1", 5.0, "pl1", "L1"),
                  ("T2", 0.05, "pl2", "L1"), ("T2", 0.5, "pl2", "L1"),
                  ("T2", 5.0, "pl1", "L1"), ("T3", 0.05, "pl1", "L1"),
                  ("T3", 0.5, "pl2", "L1"), ("T3", 5.0, "pl3", "L2")]
    tahoe_identity, tahoe_counts = [], []
    for index, (parent, dose, plate, line) in enumerate(tahoe_grid):
        tahoe_identity.append({"role": "treated", "plate": plate, "cell_line": f"line{line}",
                               "cell_line_id": line, "dose": dose, "parent": parent,
                               "smiles": DEFAULT_SMILES["tahoe_bulk"],
                               "group": f"tahoe::{plate}::{line}"})
        tahoe_counts.append(_counts(index, 2))
    for plate, line, counts in [("pl1", "L1", [1, 0, 0, 2, 0, 0]),
                                ("pl1", "L1", [5, 0, 0, 1, 0, 0]),
                                ("pl2", "L1", [2, 0, 0, 7, 0, 0]),
                                ("pl2", "L1", [7, 0, 0, 3, 0, 0]),
                                ("pl3", "L2", [3, 0, 0, 4, 0, 0]),
                                ("pl3", "L2", [9, 0, 0, 2, 0, 0])]:
        tahoe_identity.append({"role": "vehicle", "plate": plate,
                               "cell_line": f"line{line}", "cell_line_id": line,
                               "group": f"tahoe::{plate}::{line}", "vehicle": "DMSO_TF"})
        tahoe_counts.append(counts)
    return {
        "sources": {"tenv": {"counts": tenv_counts, "identity": tenv_identity},
                    "tahoe_bulk": {"counts": tahoe_counts, "identity": tahoe_identity}},
        "axes": {2: [0, 3]},
    }


def grid_fixture(root):
    """Write the grid corpus and axis release under ``root``; return paths plus arrays."""
    spec = grid_sources()
    written = write_release(root, spec["sources"], spec["axes"])
    return {**written, **spec}


def dense_fixture(root, axis_size=5000):
    """A 5,000-gene view fixture: the native model, the axis and the freeze are 5000-wide.

    The integer corpus still spans the frozen 18,164-gene axis, so the released-size view
    genuinely selects a 5,000-gene axis rather than running off the end of the corpus.
    """
    rng = np.random.default_rng(192)
    spec = grid_sources()
    sources = {}
    for name, entry in spec["sources"].items():
        counts = [(1 + rng.integers(0, 40, AXIS_SIZE)).tolist()
                  for _ in entry["identity"]]
        sources[name] = {"counts": counts, "identity": entry["identity"]}
    genes = np.sort(rng.choice(AXIS_SIZE, axis_size, replace=False)).astype(np.int64)
    written = write_release(root, sources, {axis_size: genes})
    return {**written, "sources": sources, "axes": {axis_size: genes}}


def build_view(root, size=2, name="view-5k"):
    """Publish a fixture expression view next to ``root`` and return its path + manifest."""
    from bases import view

    written = root if isinstance(root, dict) else grid_fixture(root)
    accepted = Path(written["corpus"]).parent / f"fixture-view-{size}"
    manifest = view.build(accepted, corpus_dir=written["corpus"],
                          axis_dir=written["axis"], size=size, name=name,
                          repo_root=None, fixture=True)
    return accepted, manifest, written


def treated_population(accepted, source):
    """The accepted fixture view's treated population frame, index = profile_id."""
    from bases import view

    with view.open_view(accepted) as release:
        with release.source(source) as reader:
            obs = reader.obs.copy()
    obs["profile_id"] = obs.index.astype(str)
    obs["view_row"] = np.arange(len(obs), dtype=np.int64)
    population = obs[obs.role == "treated"].copy()
    population.index = pd.Index(population["profile_id"].astype(str), name="profile_id")
    return population


def write_sealed(root, source, population):
    """A sealed Step 04-shaped fold reference for ``population`` (parents + Morgan)."""
    from views import modules

    protocol = modules()["protocol"]
    destination = Path(root) / "sealed" / source
    destination.mkdir(parents=True, exist_ok=True)
    splits = protocol.make_splits(population)
    splits.to_csv(destination / "splits.csv", index=False)
    population.to_csv(destination / "population.csv", index=False)
    parents = population[["drug_parent_inchikey", "drug_parent_smiles"]].drop_duplicates(
        ).sort_values("drug_parent_inchikey").reset_index(drop=True)
    parents.to_csv(destination / "parents.csv", index=False)
    rng = np.random.default_rng(7)
    np.save(destination / "morgan.npy",
            rng.integers(0, 2, size=(len(parents), 2048), dtype=np.uint8))
    return destination
