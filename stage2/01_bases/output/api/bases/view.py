"""Model-facing ``X`` / ``control`` / ``delta`` views of the integer raw corpus (ADR 0066).

A view is the model-facing product of the Stage 2 foundation. It is built from
the integer corpus in three deterministic steps:

1. normalize every corpus profile on the **full frozen 18,164-gene library
   total**: ``CP10K``, then natural ``log1p`` (ADR 0061's frozen normalization);
2. subset to one released axis's genes, in ascending original ``gene_index``;
3. attach ``control`` = the mean of that row's matched corpus control rows
   normalized the same way, and ``delta = X - control``.

The denominator is the **full axis, never the view**. Re-normalizing on the
view's own gene subset was rejected on measured evidence (ADR 0066 *Corpus and
view semantics*): the released ``G_2000`` carries 12.93% and ``G_5000`` 37.21%
of the tenv full-axis CP10K, so a view-local denominator would change the input
baseline by roughly 2.9x between the two axes and turn an artifact of the
denominator into a difference between models. There is one normalization and
two gene subsets.

Layout (the step's established matrix-and-metadata pattern, one h5ad per corpus
source): ``X`` is the dense ``float32`` matrix, ``layers/control`` and
``layers/delta`` carry the matched-control mean and the exact difference, ``obs``
carries the corpus row identity (``role``, dose/structure identity,
``matched_control_group``) plus the row's full-axis ``library_total``, and
``var`` carries the axis genes as original integers (never renumbered).
``library_totals.csv`` records every profile's full-axis library total.

The public reader (:class:`View` / :func:`open_view`) is the only supported
access path; it exposes bounded dense windows and never materializes the dense
observation x full-axis matrix. ``{N}k_eval_truth/`` is deliberately not created:
the evaluation redesign owns each view's truth later.
"""
from __future__ import annotations

import platform
from importlib import metadata
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

from . import corpus as corpus_release
from .gene_axis import ALGORITHM as AXIS_ALGORITHM, GeneAxis
from .release import (ReleaseError, check_manifest, output_digests, publish,
                      read_json, replay, sha256, write_json)

SCHEMA = "stage2-view-v2"
BATCH = 256
AXIS_SIZE = corpus_release.AXIS_SIZE
LAYERS = ("control", "delta")
PRODUCT_FILES = ("tenv.h5ad", "tahoe_bulk.h5ad", "library_totals.csv")
TOTALS_COLUMNS = ("source", "profile_id", "role", "library_total")
VIEW_NAMES = {2000: "view-2k", 3000: "view-3k", 4000: "view-4k", 5000: "view-5k"}
REQUIRED_OBS = (
    "role", "source_item", "source_item_manifest_sha256", "source_row", "plate",
    "cell_line", "cell_line_id", "compound_name", "vehicle", "dose_um",
    "dose_original_value", "dose_original_unit", "drug_parent_inchikey",
    "drug_parent_smiles", "source_inchikey", "n_cells", "matched_control_group",
    "library_total",
)

NORMALIZATION = {
    "library_total": "sum of the profile's full 18,164-gene integer counts",
    "formula": "log1p(count * 10000 / full_axis_library_total)",
    "denominator": ("full frozen 18,164-gene axis; the view is never normalized on "
                    "its own gene subset"),
    "source": "ADR 0061 frozen normalization, amended by ADR 0066 view semantics",
}

CONTROL = {
    "rule": ("control = arithmetic mean over the row's matched_control_group vehicle "
             "rows, each normalized on its own full-axis library total"),
    "delta": "delta = X - control, elementwise",
    "correspondence": ("controls are corpus rows; matched_control_group is carried "
                       "verbatim and the reader resolves it to row positions"),
}

RECIPE = (
    "Stream each corpus source in bounded batches; compute each profile's full-axis "
    "library total; log1p(CP10K) on that total; subset the released axis's genes in "
    "ascending gene_index; control is the mean of the row's matched vehicle rows "
    "normalized the same way; delta = X - control. float32; no view-local re-normalization."
)

BASIS_DESIGN_RECORD = {
    "source": "ADR 0066 Corpus and view semantics",
    "shares": {"2000": 0.129, "5000": 0.363},
    "note": ("the ADR figures came from an interim pre-#190 axis; the released-axis "
             "measured shares recorded above are authoritative for this product"),
}


class ViewError(ReleaseError):
    """A view release violates its contract."""


# --------------------------------------------------------------------------- #
# pure normalization helpers (the one live implementation)
# --------------------------------------------------------------------------- #
def library_totals(counts):
    """Full-axis library total of every profile (int64); a zero total fails closed."""
    counts = np.asarray(counts)
    if counts.ndim != 2:
        raise ViewError("library totals need a two-dimensional count block")
    if counts.size and counts.min() < 0:
        raise ViewError("corpus counts must be non-negative")
    totals = counts.sum(axis=1, dtype=np.int64)
    if (totals <= 0).any():
        raise ViewError("a corpus profile has a non-positive full-axis library total")
    return totals


def log1p_cp10k(counts, totals):
    """``log1p(count * 10000 / library_total)`` on the supplied full-axis totals."""
    counts = np.asarray(counts)
    totals = np.asarray(totals, dtype=np.float64)
    if counts.ndim != 2:
        raise ViewError("normalization needs a two-dimensional count block")
    if totals.shape != (counts.shape[0],):
        raise ViewError("one full-axis library total per profile is required")
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise ViewError("library totals must be finite and positive")
    return np.log1p(counts.astype(np.float64) * (10000.0 / totals[:, None]))


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def _relative(path, repo_root, label):
    path = Path(path)
    if repo_root is None:
        return f"{label}/{path.name}"
    try:
        return path.absolute().relative_to(Path(repo_root).absolute()).as_posix()
    except ValueError as error:
        raise ViewError(f"view input outside the repository: {path}") from error


def _add_input(path, inputs, repo_root, label, role):
    path = Path(path)
    key = _relative(path, repo_root, label)
    if key in inputs:
        raise ViewError(f"duplicate view input path: {key}")
    inputs[key] = {"role": role, "sha256": sha256(path), "bytes": path.stat().st_size}
    return key


def _axis_genes(axis_reader, size):
    declared = axis_reader.axis_sizes()
    if int(size) not in declared:
        raise ViewError(f"axis size {size} is not declared by the gene-axis release {declared}")
    frame = axis_reader.axis_frame(int(size))
    genes = np.asarray(frame.gene_index, dtype=np.int64)
    if len(genes) != int(size):
        raise ViewError(f"axis {size} holds {len(genes)} genes")
    if len(np.unique(genes)) != len(genes):
        raise ViewError(f"axis {size} repeats a gene_index")
    if genes.size and genes.min() < 0:
        raise ViewError(f"axis {size} carries a negative gene_index")
    if not np.array_equal(genes, np.sort(genes)):
        raise ViewError(f"axis {size} is not ordered by ascending gene_index")
    return frame, genes


def _var_frame(axis_frame, genes):
    var = pd.DataFrame({"gene_index": np.asarray(genes, dtype=np.int64),
                        "symbol": axis_frame.symbol.to_numpy()})
    var.index = pd.Index([str(gene) for gene in genes], name="gene_index_label")
    return var


def _control_means(source, genes, totals):
    """First pass: per matched control group, the mean of its normalized vehicle rows."""
    obs = source.obs
    role = obs.role.to_numpy()
    groups = obs.matched_control_group.to_numpy()
    vehicle = role == "vehicle"
    sums, counts_per = {}, {}
    for start, block in source.batches(BATCH):
        local_vehicle = np.flatnonzero(vehicle[start:start + len(block)])
        if not len(local_vehicle):
            continue
        block_totals = totals[start:start + len(block)]
        normalized = log1p_cp10k(block[:, genes], block_totals)
        for local in local_vehicle:
            group = str(groups[start + local])
            target = sums.get(group)
            if target is None:
                target = sums[group] = np.zeros(len(genes), dtype=np.float64)
            target += normalized[local]
            counts_per[group] = counts_per.get(group, 0) + 1
    means = {group: (sums[group] / counts_per[group]).astype(np.float32)
             for group in sums}
    missing = sorted({str(group) for group in groups if str(group) not in means})
    if missing:
        raise ViewError(f"matched control groups without a vehicle row: {missing[:5]}")
    return means


def _write_view(path, source, obs, var, genes, size, axis_size, totals, means, name):
    temporary = path.with_suffix(path.suffix + ".tmp")
    rows, columns = len(obs), len(genes)
    groups = obs.matched_control_group.to_numpy()
    with h5py.File(temporary, "w") as handle:
        handle.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        write_elem(handle, "obs", obs)
        write_elem(handle, "var", var)
        write_elem(handle, "uns", {
            "schema": SCHEMA,
            "source": source.source,
            "view": {"name": name, "size": int(size)},
            "axis_size": int(axis_size),
            "normalization": NORMALIZATION,
            "control": CONTROL,
            "counting": "one normalization (full-axis denominator), then the axis subset",
        })
        layers = handle.create_group("layers")
        layers.attrs.update({"encoding-type": "dict", "encoding-version": "0.1.0"})
        arrays = {}
        for key, parent in (("X", handle), ("control", layers), ("delta", layers)):
            data = parent.create_dataset(
                key, shape=(rows, columns), dtype="float32",
                chunks=(min(BATCH, rows), columns), compression="gzip",
                compression_opts=1, shuffle=True)
            data.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})
            arrays[key] = data
        position = 0
        for start, block in source.batches(BATCH):
            stop = start + len(block)
            x = log1p_cp10k(block[:, genes], totals[start:stop]).astype(np.float32)
            control = np.empty((len(x), columns), dtype=np.float32)
            for index, group in enumerate(groups[start:stop]):
                control[index] = means[str(group)]
            arrays["X"][start:stop] = x
            arrays["control"][start:stop] = control
            arrays["delta"][start:stop] = x - control
            position = stop
        if position != rows:
            raise ViewError(f"{source.source}: view row accounting differs")
    temporary.replace(path)


def _build_source(source, genes):
    obs = source.obs.copy()
    role = obs.role.to_numpy()
    if not set(role) <= {"treated", "vehicle"}:
        raise ViewError(f"{source.source}: invalid corpus role values")
    groups = obs.matched_control_group.to_numpy()
    if any(str(group) == "" for group in groups):
        raise ViewError(f"{source.source}: corpus row without a matched_control_group")
    rows = len(obs)
    totals = np.empty(rows, dtype=np.int64)
    treated = role == "treated"
    full_mass = axis_mass = treated_full_mass = treated_axis_mass = 0
    for start, block in source.batches(BATCH):
        stop = start + len(block)
        totals[start:stop] = library_totals(block)
        axis_sum = block[:, genes].sum(axis=1, dtype=np.int64)
        full_mass += int(totals[start:stop].sum())
        axis_mass += int(axis_sum.sum())
        treated_local = treated[start:stop]
        treated_full_mass += int(totals[start:stop][treated_local].sum())
        treated_axis_mass += int(axis_sum[treated_local].sum())
    means = _control_means(source, genes, totals)
    obs["library_total"] = totals
    return obs, totals, means, {
        "rows": int(rows),
        "treated": int(treated.sum()),
        "vehicle": int((~treated).sum()),
        "raw_profiles": int(rows),
        "full_axis_library_total": {
            "min": int(totals.min()),
            "median": float(np.median(totals)),
            "max": int(totals.max()),
            "sum": int(totals.sum()),
        },
        "axis_share_of_full_axis_cp10k": float(axis_mass / full_mass),
        "treated_axis_share_of_full_axis_cp10k": float(treated_axis_mass / treated_full_mass)
        if treated_full_mass else 0.0,
        "matched_groups": int(len({str(g) for g in groups[treated]})),
        "matched_control_rows": int(sum(
            int((groups == group).sum()) for group in {str(g) for g in groups[treated]})),
    }


def _environment():
    def version(package):
        try:
            return metadata.version(package)
        except metadata.PackageNotFoundError:
            return ""

    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "h5py": version("h5py"),
        "anndata": version("anndata"),
        "data_access": ("offline; reads only output/corpus and output/gene-axis; "
                        "streams 256-row batches, never a dense observation x axis matrix"),
        "resource_contract": ("one heavy job at a time; 10 GiB RLIMIT_AS fallback; "
                              "OMP_NUM_THREADS=2; peak RSS recorded by the driver"),
    }


def build(candidate, *, corpus_dir, axis_dir, size, name=None, repo_root=None,
          fixture=False):
    """Write a complete view candidate outside ``output/``."""
    candidate = Path(candidate)
    candidate.mkdir(parents=True, exist_ok=True)
    size = int(size)
    if size <= 0:
        raise ViewError("a view size must be positive")
    if not fixture and size not in VIEW_NAMES:
        raise ViewError(f"a released view size must be one of {sorted(VIEW_NAMES)}")
    axis_reader = GeneAxis(axis_dir)
    axis_reader.manifest()
    axis_frame, genes = _axis_genes(axis_reader, size)
    corpus = corpus_release.open_corpus(corpus_dir)
    corpus_manifest = corpus.manifest
    axis_size = int(corpus.axis_size)
    if not fixture and axis_size != AXIS_SIZE:
        raise ViewError(f"corpus axis size {axis_size} is not the frozen {AXIS_SIZE}")
    if genes.size and int(genes.max()) >= axis_size:
        raise ViewError("a view gene is outside the corpus's frozen axis")
    inputs = {}
    _add_input(Path(corpus_dir) / "manifest.json", inputs, repo_root, "corpus",
               "corpus_release_manifest")
    for source_name in corpus.sources:
        payload = Path(corpus_dir) / f"{source_name}.h5ad"
        declared = corpus_manifest["outputs"].get(f"{source_name}.h5ad")
        if not isinstance(declared, dict) or len(str(declared.get("sha256", ""))) != 64:
            raise ViewError(f"{source_name}: corpus manifest does not pin its payload")
        _add_input(payload, inputs, repo_root, "corpus", f"corpus:{source_name}")
        if sha256(payload) != declared["sha256"]:
            raise ViewError(f"{source_name}: corpus payload digest differs from its manifest")
    _add_input(Path(axis_dir) / "manifest.json", inputs, repo_root, "gene-axis",
               "gene_axis_release_manifest")
    axis_key = _add_input(Path(axis_dir) / f"gene_axis_{size}.csv", inputs, repo_root,
                          "gene-axis", f"gene_axis:{size}")
    summaries = {}
    frames = []
    for source_name in corpus.sources:
        source = corpus.source(source_name)
        try:
            obs, totals, means, summary = _build_source(source, genes)
            var = _var_frame(axis_frame, genes)
            _write_view(candidate / f"{source_name}.h5ad", source, obs, var, genes,
                        size, axis_size, totals, means, name)
        finally:
            source.close()
        summaries[source_name] = summary
        frames.append(pd.DataFrame({
            "source": source_name,
            "profile_id": list(obs.index),
            "role": list(obs.role),
            "library_total": totals,
        }))
    totals_frame = pd.concat(frames, ignore_index=True)[list(TOTALS_COLUMNS)]
    totals_frame.to_csv(candidate / "library_totals.csv", index=False,
                        lineterminator="\n")
    outputs = {}
    for name_ in PRODUCT_FILES:
        path = candidate / name_
        record = {"sha256": sha256(path), "bytes": path.stat().st_size}
        if name_.endswith(".h5ad"):
            with h5py.File(path) as handle:
                record.update({"rows": int(handle["X"].shape[0]),
                               "columns": int(handle["X"].shape[1]),
                               "dtype": str(handle["X"].dtype)})
        else:
            record.update({"rows": int(len(totals_frame)),
                           "columns": list(TOTALS_COLUMNS)})
        outputs[name_] = record
    manifest = {
        "schema": SCHEMA,
        "generated_by": "stage2/01_bases/execution/build_view.py",
        "fixture": bool(fixture),
        "view": {"name": name, "size": size, "file_layout": "one h5ad per corpus source",
                 "matrix": "X", "layers": list(LAYERS), "dtype": "float32"},
        "recipe": RECIPE,
        "normalization": NORMALIZATION,
        "control": CONTROL,
        "inputs": dict(sorted(inputs.items())),
        "outputs": dict(sorted(outputs.items())),
        "axis": {"size": size, "file": axis_key, "algorithm": AXIS_ALGORITHM,
                 "sha256": inputs[axis_key]["sha256"],
                 "gene_index_min": int(genes.min()),
                 "gene_index_max": int(genes.max()),
                 "ascending": True, "renumbered": False},
        "sources": summaries,
        "basis": {
            "definition": ("share of the full-axis CP10K mass carried by the view genes: "
                           "sum(counts[axis]) / sum(counts[full 18,164-gene axis]) pooled "
                           "over the source profiles"),
            "treated_definition": "the same ratio restricted to treated profiles",
            "sources": {source_name: {
                "share": summary["axis_share_of_full_axis_cp10k"],
                "treated_share": summary["treated_axis_share_of_full_axis_cp10k"],
                "n_axis": size,
            } for source_name, summary in summaries.items()},
            "design_record": dict(BASIS_DESIGN_RECORD),
        },
        "counts": {"axis_size": axis_size, "view_size": size,
                   "profiles_total": int(sum(item["rows"] for item in summaries.values())),
                   "sources": {name_: item["rows"] for name_, item in summaries.items()}},
        "environment": _environment(),
    }
    write_json(candidate / "manifest.json", manifest)
    return manifest


# --------------------------------------------------------------------------- #
# verification
# --------------------------------------------------------------------------- #
def verify(candidate, fixture=False):
    """Re-check a view candidate's schema, index space, controls and exact delta."""
    candidate = Path(candidate)
    manifest_path = candidate / "manifest.json"
    if not manifest_path.is_file():
        raise ViewError("candidate has no manifest.json")
    manifest = read_json(manifest_path)
    if manifest.get("schema") != SCHEMA:
        raise ViewError("unexpected view schema")
    if manifest.get("axis", {}).get("algorithm") != AXIS_ALGORITHM:
        raise ViewError("unexpected view axis algorithm")
    if bool(manifest.get("fixture")) is not bool(fixture):
        raise ViewError("view fixture flag differs from the verification mode")
    size = int(manifest["view"]["size"])
    if not fixture and size not in VIEW_NAMES:
        raise ViewError(f"released view size must be one of {sorted(VIEW_NAMES)}")
    check_manifest(candidate, manifest)
    if set(manifest["outputs"]) != set(PRODUCT_FILES):
        raise ViewError("view release does not declare exactly its product files")
    for path in candidate.rglob("*"):
        if "eval_truth" in path.name:
            raise ViewError(f"a view release must not ship {path.name}")
    totals = pd.read_csv(candidate / "library_totals.csv")
    if list(totals.columns) != list(TOTALS_COLUMNS):
        raise ViewError("library_totals.csv columns differ from the contract")
    if totals.profile_id.duplicated().any():
        raise ViewError("library_totals.csv repeats a profile identity")
    seen = 0
    for source_name, summary in sorted(manifest["sources"].items()):
        path = candidate / f"{source_name}.h5ad"
        if not path.is_file():
            raise ViewError(f"{source_name}: missing view payload")
        with h5py.File(path, "r") as handle:
            if set(handle["layers"]) != set(LAYERS):
                raise ViewError(f"{source_name}: unexpected view layers")
            if handle["X"].dtype != np.dtype("float32"):
                raise ViewError(f"{source_name}: view matrices must be float32")
            obs = read_elem(handle["obs"])
            var = read_elem(handle["var"])
            if handle["X"].shape != (len(obs), len(var)):
                raise ViewError(f"{source_name}: matrix shape differs from metadata")
            if handle["X"].shape != (int(summary["rows"]), size):
                raise ViewError(f"{source_name}: matrix shape differs from the manifest")
            genes = np.asarray(var.gene_index, dtype=np.int64)
            if len(genes) != size:
                raise ViewError(f"{source_name}: view does not hold exactly {size} genes")
            if not np.array_equal(genes, np.sort(genes)) or len(np.unique(genes)) != size:
                raise ViewError(f"{source_name}: view genes are not a distinct ascending set")
            if genes.size and int(genes.min()) < 0:
                raise ViewError(f"{source_name}: a view gene is negative")
            if not fixture and genes.size and int(genes.max()) >= AXIS_SIZE:
                raise ViewError(f"{source_name}: a view gene is outside the frozen axis")
            missing = [column for column in REQUIRED_OBS if column not in obs.columns]
            if missing:
                raise ViewError(f"{source_name}: obs lacks {missing}")
            if not set(obs.role) <= {"treated", "vehicle"}:
                raise ViewError(f"{source_name}: invalid role values")
            if (obs.library_total <= 0).any():
                raise ViewError(f"{source_name}: a profile has a non-positive library total")
            if obs.matched_control_group.eq("").any():
                raise ViewError(f"{source_name}: a row lacks matched_control_group")
            for start in range(0, len(obs), BATCH):
                stop = min(start + BATCH, len(obs))
                x = np.asarray(handle["X"][start:stop], dtype=np.float32)
                control = np.asarray(handle["layers/control"][start:stop], dtype=np.float32)
                delta = np.asarray(handle["layers/delta"][start:stop], dtype=np.float32)
                if not np.isfinite(x).all() or (x < 0).any():
                    raise ViewError(f"{source_name}: invalid X values")
                if not np.isfinite(control).all() or (control < 0).any():
                    raise ViewError(f"{source_name}: invalid control values")
                if not np.array_equal(delta, x - control):
                    raise ViewError(f"{source_name}: delta is not exactly X - control")
        part = totals[totals.source == source_name]
        if len(part) != int(summary["rows"]):
            raise ViewError(f"{source_name}: library_totals.csv row count differs")
        if list(part.profile_id) != list(obs.index):
            raise ViewError(f"{source_name}: library_totals.csv profile order differs")
        if list(part.role) != list(obs.role):
            raise ViewError(f"{source_name}: library_totals.csv role differs")
        if list(part.library_total.astype(np.int64)) != list(
                obs.library_total.astype(np.int64)):
            raise ViewError(f"{source_name}: recorded library total differs from obs")
        seen += len(obs)
    if len(totals) != seen:
        raise ViewError("library_totals.csv covers rows the manifest does not declare")
    basis = manifest.get("basis", {})
    for source_name in manifest["sources"]:
        share = basis.get("sources", {}).get(source_name, {}).get("share")
        if not isinstance(share, (int, float)) or not 0.0 <= float(share) <= 1.0:
            raise ViewError(f"{source_name}: missing or invalid full-axis CP10K share")
    return {"result": "PASS", "manifest_sha256": sha256(manifest_path),
            "outputs": sorted(output_digests(manifest)),
            "view": manifest["view"], "sources": sorted(manifest["sources"])}


def publish_release(output, *, root=None, **kwargs):
    return publish(output, lambda stage: build(stage, **kwargs),
                   lambda stage: verify(stage, fixture=kwargs.get("fixture", False)),
                   root=root)


def replay_release(output, *, root=None, **kwargs):
    return replay(output, lambda stage: build(stage, **kwargs),
                  lambda stage: verify(stage, fixture=kwargs.get("fixture", False)),
                  root=root)


# --------------------------------------------------------------------------- #
# public reader
# --------------------------------------------------------------------------- #
class ViewSource:
    """One corpus source's view with bounded ``X``/``control``/``delta`` windows."""

    def __init__(self, root, source):
        self.root = Path(root)
        self.source = source
        self.path = self.root / f"{source}.h5ad"
        self._handle = h5py.File(self.path, "r")
        self.obs = read_elem(self._handle["obs"])
        self.var = read_elem(self._handle["var"])
        self._X = self._handle["X"]
        self._control = self._handle["layers/control"]
        self._delta = self._handle["layers/delta"]
        if self._X.dtype != np.dtype("float32"):
            raise ViewError(f"{source}: view matrices must be float32")
        if set(self._handle["layers"]) != set(LAYERS):
            raise ViewError(f"{source}: unexpected view layers")
        self.shape = (len(self.obs), len(self.var))
        genes = np.asarray(self.var.gene_index, dtype=np.int64)
        if len(genes) and not np.array_equal(genes, np.sort(genes)):
            raise ViewError(f"{source}: view genes are not ascending")
        self._groups = {}
        for position in np.flatnonzero(self.obs.role.eq("vehicle").to_numpy()):
            key = str(self.obs.matched_control_group.iloc[int(position)])
            self._groups.setdefault(key, []).append(int(position))

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def gene_index(self):
        """The axis genes as original frozen indices, ascending (never renumbered)."""
        return np.asarray(self.var.gene_index, dtype=np.int64)

    @property
    def axis_size(self):
        return self.shape[1]

    def _window(self, dataset, start, stop):
        stop = self.shape[0] if stop is None else int(stop)
        return np.asarray(dataset[int(start):stop], dtype=np.float32)

    def X(self, start=0, stop=None):
        return self._window(self._X, start, stop)

    def control(self, start=0, stop=None):
        return self._window(self._control, start, stop)

    def delta(self, start=0, stop=None):
        return self._window(self._delta, start, stop)

    def batches(self, batch_size=BATCH):
        if batch_size <= 0:
            raise ViewError("batch_size must be positive")
        for start in range(0, self.shape[0], batch_size):
            stop = min(start + batch_size, self.shape[0])
            yield start, self.X(start, stop), self.control(start, stop), self.delta(start, stop)

    def library_totals(self):
        """Every profile's full-axis library total (int64)."""
        return np.asarray(self.obs.library_total, dtype=np.int64)

    def control_groups(self):
        return {key: np.asarray(rows, dtype=np.int64)
                for key, rows in self._groups.items()}

    def control_rows(self, group):
        if group not in self._groups:
            raise ViewError(f"unknown matched control group: {group}")
        return np.asarray(self._groups[group], dtype=np.int64)

    def matched_control_rows(self, index):
        """Vehicle row positions whose mean is this row's ``control``."""
        if not 0 <= int(index) < self.shape[0]:
            raise ViewError(f"view row out of range: {index}")
        group = str(self.obs.matched_control_group.iloc[int(index)])
        return self.control_rows(group)


class View:
    """Read-only view release: one released axis over the two corpus row sets."""

    def __init__(self, root):
        self.root = Path(root)
        self._manifest = None
        self._sources = {}

    def manifest(self):
        if self._manifest is None:
            path = self.root / "manifest.json"
            if not path.is_file():
                raise FileNotFoundError(f"no view manifest at {path}")
            manifest = read_json(path)
            if manifest.get("schema") != SCHEMA:
                raise ViewError(f"{path}: unexpected view schema")
            if manifest.get("axis", {}).get("algorithm") != AXIS_ALGORITHM:
                raise ViewError(f"{path}: unexpected view axis algorithm")
            self._manifest = manifest
        return self._manifest

    def digests(self):
        """Declared sha256 of every view output."""
        return output_digests(self.manifest())

    @property
    def sources(self):
        return tuple(sorted(self.manifest()["sources"]))

    @property
    def size(self):
        return int(self.manifest()["view"]["size"])

    @property
    def axis_size(self):
        return int(self.manifest()["counts"]["axis_size"])

    def source(self, name):
        if name not in self.sources:
            raise ViewError(f"unknown view source: {name}")
        if name not in self._sources:
            self._sources[name] = ViewSource(self.root, name)
        return self._sources[name]

    def library_totals(self):
        return pd.read_csv(self.root / "library_totals.csv")

    def close(self):
        for source in self._sources.values():
            source.close()
        self._sources = {}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def open_view(root):
    """Open a published view directory (the only supported access path)."""
    return View(root)
