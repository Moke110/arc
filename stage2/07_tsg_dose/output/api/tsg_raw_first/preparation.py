"""Prepare versioned dose inputs from verified raw-first views and shared vectors."""
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from expression_view import open_view, sha256
from prediction_axis import open_axes
from stage2_util.representations import RepresentationTable
from stage2_util.seeds import read_seeds
from tsg_round1.protocol import SPLIT_SEED, TRAINING_SEED, make_splits

from .inputs import (GENE_COLUMNS, INPUT_SCHEMA, KPGT_SCHEMA, PAIR_SCHEMA,
                     REQUIRED, SOURCES, Inputs, KPGTInputs, check_generation,
                     match_controls, ordered_parents)


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _representation(source, specification):
    path, digest = specification if specification is not None else (None, None)
    return RepresentationTable(path, source=source, manifest_sha256=digest)


def _fresh_destination(destination, inputs):
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(f"raw-first preparation requires a fresh destination: {destination}")
    if any(Path(source).resolve() in destination.resolve().parents for source in inputs):
        raise ValueError("candidate cannot reside inside a frozen input")
    destination.parent.mkdir(parents=True, exist_ok=True)
    return destination


def _freeze(view, output, generation, table):
    obs = view.obs.copy()
    if not obs.index.astype(str).equals(pd.Index(obs.profile_id.astype(str))):
        raise ValueError("raw-first profile index differs from source observation identities")
    if not obs.role.isin(("treated", "vehicle")).all():
        raise ValueError("full expression view contains an unsupported profile role")
    obs["view_row"] = np.arange(len(obs))
    treated = np.flatnonzero(obs.role.eq("treated"))
    vehicles = np.flatnonzero(obs.role.eq("vehicle"))
    controls = match_controls(obs.iloc[treated], obs.iloc[vehicles])
    population = obs.iloc[treated]
    if population.empty:
        raise ValueError("full expression view contains no treated rows")
    target = np.empty(len(obs), dtype=np.intp)
    target[np.r_[treated, vehicles]] = np.arange(len(obs))
    output.mkdir()
    for name, layer in (("expression", "X"), ("control", "control")):
        matrix = np.lib.format.open_memmap(output / f"{name}.npy", mode="w+", dtype=np.float32,
                                           shape=view.shape)
        for start, block in view.batches(layer):
            if not np.isfinite(block).all():
                raise ValueError(f"nonfinite values in verified {layer} view")
            matrix[target[start:start + len(block)]] = block
        matrix.flush()
        del matrix

    parents = ordered_parents(population)
    fingerprints = table.vectors("ECFP4C", table.verified_parent_keys(parents)).astype(np.uint8)
    np.save(output / "morgan.npy", fingerprints)
    parents.to_csv(output / "parents.csv", index=False)
    view.var[GENE_COLUMNS].to_csv(output / "genes.csv", index=False)
    population.to_csv(output / "population.csv", index=False)
    obs.iloc[np.r_[treated, vehicles]].to_csv(output / "all_obs.csv", index=False)
    make_splits(population).to_csv(output / "splits.csv", index=False)
    _json(output / "matched_controls.json", {str(row): members for row, members in controls.items()})
    payload = view.manifest["outputs"][f"full/{view.source}.h5ad"]
    if TRAINING_SEED not in read_seeds():
        raise ValueError("shared seed table no longer contains the study training seed")
    _json(output / "config.json", {
        "source": view.source, "view_generation": generation, "view_payload_sha256": payload,
        "split_seed": SPLIT_SEED, "training_seed": TRAINING_SEED,
        "representation_sha256": table.manifest_sha256, "n_treated": len(treated),
        "n_matched_vehicles": len(vehicles), "n_parents": len(parents),
    })
    _json(output / "manifest.json", {
        "schema": INPUT_SCHEMA, "source": view.source, "view_generation": generation,
        "outputs": {name: sha256(output / name) for name in sorted(REQUIRED)},
    })


def prepare(destination, *, view=None, representations=None, sources=SOURCES):
    """Freeze the paired accepted full view; an explicit source allows rehearsal."""
    if not sources or set(sources) - set(SOURCES) or len(set(sources)) != len(sources):
        raise ValueError("choose distinct raw-first source names")
    representations = representations or {}
    input_paths = ([] if view is None else [view]) + [spec[0] for spec in representations.values()]
    destination = _fresh_destination(destination, input_paths)
    axis = open_axes().axis(5000)
    with TemporaryDirectory(prefix=".raw-first-", dir=destination.parent) as temp:
        staging = Path(temp) / "pair"
        staging.mkdir()
        generation = None
        source_digests = {}
        for source in sources:
            with open_view(view, source=source, axis_sha256=axis.sha256) as full:
                source_generation = {"schema": full.manifest["schema"], "family": full.family,
                                     "view_sha256": full.manifest_sha256}
                source_generation.update({key: full.manifest[key] for key in (
                    "corpus_sha256", "axis_sha256", "axis_size", "normalization", "control"
                )})
                check_generation(source_generation, axis)
                if (full.manifest.get("default_family") != "full"
                    or not set(sources) <= set(full.manifest.get("sources", ()))) or (
                    generation is not None and source_generation != generation
                ):
                    raise ValueError("raw-first source arms differ from one default release")
                generation = source_generation
                _freeze(full, staging / source, generation, _representation(source, representations.get(source)))
            source_digests[source] = sha256(staging / source / "manifest.json")
        _json(staging / "manifest.json", {"schema": PAIR_SCHEMA, "view_generation": generation,
                                          "sources": source_digests})
        for source in sources:
            Inputs(staging / source)
        staging.rename(destination)
    return destination


def prepare_kpgt(tenv, destination, *, representation=None):
    """Seal KPGT vectors from the same shared representation bound to the tenv input."""
    tenv = Path(tenv)
    input_paths = [tenv.parent] + ([representation[0]] if representation is not None else [])
    destination = _fresh_destination(destination, input_paths)
    data = Inputs(tenv)
    if data.source != "tenv":
        raise ValueError("KPGT wrapper requires tenv")
    table = _representation("tenv", representation)
    if table.manifest_sha256 != data.config["representation_sha256"]:
        raise ValueError("KPGT representation release differs from the input freeze")
    with TemporaryDirectory(prefix=".raw-first-kpgt-", dir=destination.parent) as temp:
        output = Path(temp) / "kpgt"
        pair_dir = output / "pair"
        pair_dir.mkdir(parents=True)
        shutil.copytree(tenv, pair_dir / "tenv")
        shutil.copyfile(tenv.parent / "manifest.json", pair_dir / "manifest.json")
        data.parents.to_csv(output / "representation_parents.csv", index=False)
        np.save(output / "KPGT.npy", table.vectors("KPGT", table.verified_parent_keys(data.parents)))
        files = ("pair/manifest.json", "pair/tenv/manifest.json", "KPGT.npy", "representation_parents.csv")
        _json(output / "manifest.json", {"schema": KPGT_SCHEMA, "representation_sha256": table.manifest_sha256,
                                         "files": {name: sha256(output / name) for name in files}})
        KPGTInputs(output)
        output.rename(destination)
    return destination
