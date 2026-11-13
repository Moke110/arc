"""Hermetic fixtures for the Stage 2 model-facing view tests (issue #191, ADR 0066).

Writes a contract-shaped integer corpus (one h5ad per source plus the manifest)
and a released-axis directory (``gene_axis_<N>.csv`` plus the manifest) directly,
so the view build/replay boundary is exercised without the multi-GB accepted
products. Every count is hand-worked in the tests; the fixture only drives the
public build and reader.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem

from bases import corpus as corpus_module
from bases import gene_axis as axis_module

AXIS_SIZE = 18164


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _obs_frame(name, counts, roles, groups):
    """Corpus-shaped obs for one fixture source; identity is derived from the group."""
    rows = []
    for index in range(len(counts)):
        group = str(groups[index])
        parts = group.split("::")
        plate = parts[1] if len(parts) > 1 else "plate1"
        treated = roles[index] == "treated"
        rows.append({
            "profile_id": f"{name}:p{index}",
            "source_item": "fixture",
            "source_item_manifest_sha256": "0" * 64,
            "source_row": index,
            "role": roles[index],
            "plate": plate,
            "cell_line": "hiPSC-CM" if name == "tenv" else "line1",
            "cell_line_id": "" if name == "tenv" else "L1",
            "compound_name": "compound" if treated else "DMSO",
            "vehicle": "" if treated and name != "tenv" else ("DMSO" if name == "tenv" else "DMSO_TF"),
            "dose_um": 1.0 if treated else 0.0,
            "dose_original_value": 1.0 if treated else 0.0,
            "dose_original_unit": "uM",
            "drug_parent_inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N" if treated else "",
            "drug_parent_smiles": "CCO" if treated else "",
            "source_inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N" if treated else "",
            "n_cells": 0 if name == "tenv" else 50,
            "matched_control_group": group,
        })
    frame = pd.DataFrame(rows).set_index("profile_id", drop=False)
    frame.index.name = "profile_id"
    return frame


def _write_source(path, name, counts, roles, groups):
    counts = np.asarray(counts, dtype=np.int32)
    obs = _obs_frame(name, counts, roles, groups)
    var = pd.DataFrame({"gene_index": np.arange(counts.shape[1], dtype=np.int64),
                        "symbol": [f"S{gene}" for gene in range(counts.shape[1])]})
    var.index = pd.Index([str(gene) for gene in var.gene_index],
                         name="gene_index_label")
    with h5py.File(path, "w") as handle:
        handle.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        write_elem(handle, "obs", obs)
        write_elem(handle, "var", var)
        write_elem(handle, "uns", {"schema": corpus_module.SCHEMA, "source": name})
        write_elem(handle, "layers", {})
        matrix = handle.create_dataset("X", data=counts)
        matrix.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})
    return obs


def write_corpus(root, sources):
    """Write ``root/corpus``: one h5ad per source plus a contract-shaped manifest.

    ``sources`` maps a source name to ``{"counts": [[...]], "role": [...],
    "matched_control_group": [...]}``. Every vehicle row is assumed to belong to
    a group that holds it; the view builder resolves controls through the group.
    """
    release = Path(root) / "corpus"
    release.mkdir(parents=True, exist_ok=True)
    outputs, summaries = {}, {}
    axis_size = None
    for name in sorted(sources):
        entry = sources[name]
        counts = np.asarray(entry["counts"], dtype=np.int32)
        if axis_size is None:
            axis_size = counts.shape[1]
        elif counts.shape[1] != axis_size:
            raise ValueError("fixture sources disagree on the frozen axis size")
        path = release / f"{name}.h5ad"
        obs = _write_source(path, name, counts, entry["role"],
                            entry["matched_control_group"])
        roles = obs.role.to_numpy()
        groups = obs.matched_control_group.to_numpy()
        treated = roles == "treated"
        matched_groups = sorted({str(group) for group in groups[treated]})
        outputs[f"{name}.h5ad"] = {"sha256": sha256(path), "bytes": path.stat().st_size,
                                   "rows": int(len(obs)), "columns": int(axis_size),
                                   "dtype": "int32"}
        summaries[name] = {
            "counts": {"profiles": int(len(obs)),
                       "treated": int(treated.sum()),
                       "vehicle": int((~treated).sum()),
                       "raw_profiles": int(len(obs))},
            "excluded_by_reason": {},
            "controls": {
                "rule": "fixture",
                "vehicle_rows": int((~treated).sum()),
                "vehicle_groups": int(len({str(group) for group in groups[~treated]})),
                "matched_control_rows": int(sum(
                    int((groups == group).sum()) for group in matched_groups)),
                "matched_groups": len(matched_groups),
                "treated_without_control": 0,
            },
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


def write_axis_release(root, axes):
    """Write ``root/gene-axis``: one ``gene_axis_<N>.csv`` per size plus a manifest."""
    release = Path(root) / "gene-axis"
    release.mkdir(parents=True, exist_ok=True)
    outputs, sizes = {}, {}
    for size, genes in sorted(axes.items()):
        genes = np.asarray(genes, dtype=np.int64)
        frame = pd.DataFrame({"gene_index": genes,
                              "symbol": [f"S{gene}" for gene in genes]})
        path = release / f"gene_axis_{int(size)}.csv"
        frame.to_csv(path, index=False, lineterminator="\n")
        key = f"stage2/01_bases/output/gene-axis/gene_axis_{int(size)}.csv"
        outputs[key] = {"sha256": sha256(path), "rows": int(len(frame))}
        sizes[str(int(size))] = {"n": int(len(frame)),
                                 "file": f"gene_axis_{int(size)}.csv"}
    manifest = {
        "schema": axis_module.SCHEMA,
        "fixture": True,
        "axes": sizes,
        "outputs": outputs,
        "release": {"axis_sizes": [int(size) for size in sorted(axes)]},
    }
    (release / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return release


def write_release(root, sources, axes):
    """The two pinned inputs of a view build, under ``root``."""
    root = Path(root)
    return {"corpus": write_corpus(root, sources),
            "axis": write_axis_release(root, axes)}


def normalization_fixture(root):
    """Six frozen genes, one treated tenv row, two same-plate DMSO controls.

    Only genes 0 and 3 are carried by the view; genes 1, 2, 4 and 5 exist so the
    full-axis library total differs from the view-local one. A Media vehicle
    well forms a singleton group, so its control is itself.
    """
    sources = {
        "tenv": {
            "counts": [
                [12, 100, 100, 24, 100, 100],   # treated
                [6, 100, 100, 6, 100, 100],     # DMSO vehicle 1
                [18, 100, 100, 18, 100, 100],   # DMSO vehicle 2
                [3, 100, 100, 3, 100, 100],     # Media vehicle (singleton group)
            ],
            "role": ["treated", "vehicle", "vehicle", "vehicle"],
            "matched_control_group": ["tenv::P1::DMSO"] * 3 + ["tenv::P1::Media"],
        },
        "tahoe_bulk": {
            "counts": [
                [5, 0, 0, 5, 0, 0],             # treated
                [1, 0, 0, 1, 0, 0],             # DMSO_TF vehicle 1
                [3, 0, 0, 3, 0, 0],             # DMSO_TF vehicle 2
            ],
            "role": ["treated", "vehicle", "vehicle"],
            "matched_control_group": ["tahoe::pl1::L1"] * 3,
        },
    }
    return {**write_release(root, sources, {2: [0, 3]}), "sources": sources}


def full_size_fixture(root):
    """A 5,010-gene corpus with released-size axes of exactly 2,000 and 5,000 genes."""
    axis_size = 5010
    columns = np.arange(axis_size, dtype=np.int32)
    tenv_counts = [
        (1 + (columns + 3 * row) % 7).tolist() for row in range(3)
    ]
    tahoe_counts = [
        (2 + (columns + 5 * row) % 11).tolist() for row in range(2)
    ]
    sources = {
        "tenv": {
            "counts": tenv_counts,
            "role": ["treated", "vehicle", "vehicle"],
            "matched_control_group": ["tenv::P1::DMSO"] * 3,
        },
        "tahoe_bulk": {
            "counts": tahoe_counts,
            "role": ["treated", "vehicle"],
            "matched_control_group": ["tahoe::pl1::L1"] * 2,
        },
    }
    rng = np.random.default_rng(191)
    axes = {
        2000: np.sort(rng.choice(axis_size, 2000, replace=False)).astype(np.int64),
        5000: np.sort(rng.choice(axis_size, 5000, replace=False)).astype(np.int64),
    }
    fixture = write_release(root, sources, axes)
    return fixture, axes
