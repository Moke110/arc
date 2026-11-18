"""5k-local model-facing views of the unchanged integer corpus (ADR 0077).

Select G_5000 before calculating each profile's int64 library total, then use
natural log1p(CP10K). Normalize each vehicle independently before matched-group
averaging; delta is X minus that control. Payloads are float32 and accessed in
bounded windows. Four named gene axes remain separate scientific products.
"""
from __future__ import annotations

import hashlib
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

SCHEMA = "stage2-view-v3"
BATCH = 256
AXIS_SIZE = corpus_release.AXIS_SIZE
LAYERS = ("control", "delta")
PRODUCT_FILES = ("tenv.h5ad", "tahoe_bulk.h5ad", "library_totals.csv")
TOTALS_COLUMNS = ("source", "profile_id", "role", "library_total",
                  "normalization_schema", "denominator_axis_sha256")
VIEW_NAMES = {5000: "view-5k"}
REQUIRED_OBS = (
    "role", "source_item", "source_item_manifest_sha256", "source_row", "plate",
    "cell_line", "cell_line_id", "compound_name", "vehicle", "dose_um",
    "dose_original_value", "dose_original_unit", "drug_parent_inchikey",
    "drug_parent_smiles", "source_inchikey", "n_cells", "matched_control_group",
    "library_total",
)

NORMALIZATION = {
    "library_total": "sum of the profile's selected G_5000 integer counts (int64)",
    "formula": "log1p(count * 10000 / selected_axis_library_total)",
    "denominator": "selected gene axis G_5000, before CP10K and natural log1p",
    "source": "ADR 0077 five-k-only axis-local normalization",
}
CONTROL = {
    "rule": ("control = arithmetic mean over matched_control_group vehicle rows, "
             "each normalized on its own selected-axis library total"),
    "delta": "delta = X - control, elementwise",
    "correspondence": "matched_control_group is carried verbatim from the integer corpus",
}
RECIPE = (
    "Stream bounded integer corpus batches; select G_5000 in ascending original "
    "gene_index; sum selected counts in int64; normalize every row independently "
    "with natural log1p(CP10K); average normalized matched vehicles for control; "
    "delta = X - control. Publish float32."
)


class ViewError(ReleaseError):
    """A view release violates its contract."""


# --------------------------------------------------------------------------- #
# pure normalization helpers (the one live implementation)
# --------------------------------------------------------------------------- #
def library_totals(counts):
    """Selected-axis library total of every profile (int64); a zero total fails closed."""
    counts = np.asarray(counts)
    if counts.ndim != 2:
        raise ViewError("library totals need a two-dimensional count block")
    if counts.size and counts.min() < 0:
        raise ViewError("corpus counts must be non-negative")
    totals = counts.sum(axis=1, dtype=np.int64)
    if (totals <= 0).any():
        raise ViewError("a corpus profile has a non-positive selected-axis library total")
    return totals


def log1p_cp10k(counts, totals):
    """``log1p(count * 10000 / library_total)`` on the supplied selected-axis totals."""
    counts = np.asarray(counts)
    totals = np.asarray(totals, dtype=np.float64)
    if counts.ndim != 2:
        raise ViewError("normalization needs a two-dimensional count block")
    if totals.shape != (counts.shape[0],):
        raise ViewError("one selected-axis library total per profile is required")
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


def _write_view(path, source, obs, var, genes, size, axis_size, totals, means, name, axis):
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
            "axis": axis,
            "normalization": NORMALIZATION,
            "control": CONTROL,
            "counting": "select axis, then normalize each profile on its selected-axis total",
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
    for start, block in source.batches(BATCH):
        stop = start + len(block)
        selected = block[:, genes]
        bad = np.flatnonzero(selected.sum(axis=1, dtype=np.int64) <= 0)
        if len(bad):
            profile = str(obs.index[start + int(bad[0])])
            raise ViewError(f"{source.source}: profile {profile} has non-positive selected-axis total")
        totals[start:stop] = library_totals(selected)
    means = _control_means(source, genes, totals)
    obs["library_total"] = totals
    return obs, totals, means, {
        "rows": int(rows),
        "treated": int(treated.sum()),
        "vehicle": int((~treated).sum()),
        "raw_profiles": int(rows),
        "selected_axis_library_total": {
            "min": int(totals.min()), "median": float(np.median(totals)),
            "max": int(totals.max()), "sum": int(totals.sum()),
        },
        "matched_groups": int(len({str(g) for g in groups[treated]})),
        "matched_control_rows": int(sum(
            int((groups == group).sum()) for group in {str(g) for g in groups[treated]})),
    }


def _columns_digest(var):
    return hashlib.sha256(var[["gene_index", "symbol"]].to_csv(
        index=False, lineterminator="\n").encode()).hexdigest()


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
    if not fixture:
        name = name or VIEW_NAMES[size]
        if name != VIEW_NAMES[size]:
            raise ViewError("only the canonical view-5k expression release name is supported")
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
    axis_binding = {"size": size, "file": axis_key, "algorithm": AXIS_ALGORITHM,
                    "sha256": inputs[axis_key]["sha256"],
                    "columns_sha256": _columns_digest(_var_frame(axis_frame, genes)),
                    "gene_index_min": int(genes.min()), "gene_index_max": int(genes.max()),
                    "ascending": True, "renumbered": False}
    summaries = {}
    frames = []
    for source_name in corpus.sources:
        source = corpus.source(source_name)
        try:
            obs, totals, means, summary = _build_source(source, genes)
            var = _var_frame(axis_frame, genes)
            _write_view(candidate / f"{source_name}.h5ad", source, obs, var, genes,
                        size, axis_size, totals, means, name, axis_binding)
        finally:
            source.close()
        summaries[source_name] = summary
        frames.append(pd.DataFrame({
            "source": source_name,
            "profile_id": list(obs.index),
            "role": list(obs.role),
            "library_total": totals,
            "normalization_schema": SCHEMA,
            "denominator_axis_sha256": axis_binding["sha256"],
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
        "axis": axis_binding,
        "sources": summaries,
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
    _check_contract(manifest)
    check_manifest(candidate, manifest)
    if set(manifest["outputs"]) != set(PRODUCT_FILES):
        raise ViewError("view release does not declare exactly its product files")
    for path in candidate.rglob("*"):
        if "eval_truth" in path.name:
            raise ViewError(f"a view release must not ship {path.name}")
    totals = pd.read_csv(candidate / "library_totals.csv")
    if list(totals.columns) != list(TOTALS_COLUMNS):
        raise ViewError("library_totals.csv columns differ from the contract")
    if (set(totals.normalization_schema) != {SCHEMA}
            or set(totals.denominator_axis_sha256) != {manifest["axis"]["sha256"]}):
        raise ViewError("library totals normalization/axis binding differs")
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
            _check_payload(handle, manifest, source_name)
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
            if (int(obs.role.eq("treated").sum()) != int(summary["treated"])
                    or int(obs.role.eq("vehicle").sum()) != int(summary["vehicle"])):
                raise ViewError(f"{source_name}: role populations differ from manifest")
            vehicle_sums, vehicle_counts = {}, {}
            groups = obs.matched_control_group.astype(str).to_numpy()
            for start in range(0, len(obs), BATCH):
                stop = min(start + BATCH, len(obs))
                local = np.flatnonzero(obs.role.iloc[start:stop].eq("vehicle").to_numpy())
                if not len(local):
                    continue
                normalized = np.asarray(handle["X"][start:stop], dtype=np.float64)
                for row in local:
                    group = groups[start + row]
                    if group not in vehicle_sums:
                        vehicle_sums[group] = np.zeros(size, dtype=np.float64)
                        vehicle_counts[group] = 0
                    vehicle_sums[group] += normalized[row]
                    vehicle_counts[group] += 1
            vehicle_means = {key: value / vehicle_counts[key]
                             for key, value in vehicle_sums.items()}
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
                if not np.allclose(np.expm1(x.astype(np.float64)).sum(axis=1),
                                   10000, rtol=1e-6, atol=0.01):
                    raise ViewError(f"{source_name}: X is not selected-axis CP10K")
                for row, group in enumerate(groups[start:stop]):
                    if group not in vehicle_means or not np.allclose(
                            control[row], vehicle_means[group], rtol=0, atol=1e-6):
                        raise ViewError(f"{source_name}: control differs from matched normalized vehicles")
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
    return {"result": "PASS", "manifest_sha256": sha256(manifest_path),
            "outputs": sorted(output_digests(manifest)),
            "view": manifest["view"], "sources": sorted(manifest["sources"])}


def _check_contract(manifest):
    if manifest.get("schema") != SCHEMA:
        raise ViewError("unexpected view schema; full-axis generations are historical")
    if manifest.get("normalization") != NORMALIZATION or manifest.get("control") != CONTROL:
        raise ViewError("view normalization/control generation differs")
    if manifest.get("recipe") != RECIPE:
        raise ViewError("view recipe generation differs")
    if not manifest.get("fixture") and manifest.get("view", {}).get("size") != 5000:
        raise ViewError("only the 5k expression release is supported")
    if not manifest.get("fixture") and manifest.get("view", {}).get("name") != "view-5k":
        raise ViewError("only the canonical view-5k expression release name is supported")
    axis = manifest.get("axis", {})
    if axis.get("algorithm") != AXIS_ALGORITHM or axis.get("size") != manifest["view"]["size"]:
        raise ViewError("view axis binding differs")
    pinned = manifest.get("inputs", {}).get(axis.get("file"), {})
    if pinned.get("sha256") != axis.get("sha256") or len(str(axis.get("columns_sha256", ""))) != 64:
        raise ViewError("view axis digest binding differs")
    if set(manifest.get("sources", {})) != {"tenv", "tahoe_bulk"}:
        raise ViewError("a view must bind both corpus sources")


def _check_payload(handle, manifest, source):
    uns = read_elem(handle["uns"])
    if (uns.get("schema") != SCHEMA or uns.get("source") != source
            or uns.get("normalization") != NORMALIZATION or uns.get("control") != CONTROL):
        raise ViewError(f"{source}: mixed payload normalization generation")
    axis = uns.get("axis", {})
    if axis != manifest["axis"] or uns.get("view") != {"name": manifest["view"]["name"], "size": manifest["view"]["size"]}:
        raise ViewError(f"{source}: payload axis binding differs")
    var = read_elem(handle["var"])
    if _columns_digest(var) != manifest["axis"]["columns_sha256"]:
        raise ViewError(f"{source}: payload gene indices/symbols differ from the pinned axis")


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

    def __init__(self, root, source, manifest):
        self.root = Path(root)
        self.source = source
        self.path = self.root / f"{source}.h5ad"
        self._handle = h5py.File(self.path, "r")
        try:
            _check_payload(self._handle, manifest, source)
        except BaseException:
            self._handle.close()
            raise
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
        """Every profile's selected-axis denominator (int64)."""
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
            _check_contract(manifest)
            check_manifest(self.root, manifest)
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
            self._sources[name] = ViewSource(self.root, name, self.manifest())
        return self._sources[name]

    def library_totals(self):
        self.manifest()
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
