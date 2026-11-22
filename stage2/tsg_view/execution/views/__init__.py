"""Shared harness for the ADR 0082 expression-view products.

Every view family is a thin module under ``views/`` that calls this harness:

* :class:`Context` resolves the pinned inputs (integer corpus, G_5000 axis, accepted
  ``view-5k`` release), binds the owning Step 01 implementation and carries the staged
  product directory and installer.
* :func:`write_product` is the single manifest-driven product writer. It owns the
  ``tsg-view-product-v1`` schema, the h5ad layout, companion tables (``profile_ids``,
  ``members``, ``control_sources``), the digest bookkeeping and the alias/delta-only
  modes, so a variant module never re-implements them.
* The pooling / drawing primitives (:func:`condition_table`,
  :func:`pooled_condition_counts`, :func:`vehicle_group_sums`, :func:`pooled_control`,
  :func:`solo_choice`, :func:`condition_blocks`, :func:`selected_blocks`) implement the
  locked ADR 0082 decisions C7.1-C7.5: one row per condition with a materialized
  control, normalize-of-sum pooling, a frozen solo seed and a member-choice row model.

Derived views never pass through ``bases.view.verify`` or ``tsg_round1.inputs.Inputs``:
both hard-assert the released full population. Normalization itself is never re-derived;
this module calls ``bases.view.library_totals`` and ``bases.view.log1p_cp10k``.
"""
from __future__ import annotations

import hashlib
import os
import platform
import sys
from contextlib import contextmanager
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

SCHEMA = "tsg-view-product-v1"
VIEWS = ("full", "merge", "solo", "solo-deepest", "full-pooled-label")
SOURCES = ("tenv", "tahoe_bulk")
SIZE = 5000
NAME = "view-5k"
SOLO_SEED = 20261122
BATCH = 256
CONDITION_KEYS = ("cell_line", "drug_parent_inchikey", "dose_um")
REQUIRED_OBS = ("profile_id", "role", "cell_line", "cell_line_id", "dose_um",
                "drug_parent_inchikey", "drug_parent_smiles", "library_total")
MODES = {"full": "alias", "merge": "payload", "solo": "payload",
         "solo-deepest": "payload", "full-pooled-label": "delta"}
RECIPE = {
    "full": ("ADR 0077 unchanged: select G_5000 integer counts, normalize every row on its "
             "own selected-axis total with natural log1p(CP10K), average normalized matched "
             "vehicles, delta = X - control; published as a digest-equal alias."),
    "merge": ("Sum the members' retained G_5000 integer counts per condition (int64), "
              "normalize the sum on its own selected-axis total; pool the union of the "
              "members' matched control groups' vehicle counts the same way; delta = X - control."),
    "solo": ("Draw one treated profile per condition once per source from a frozen seed "
             "and keep that member's full-view X and control."),
    "solo-deepest": ("Keep the member with the largest retained selected-axis total per "
                     "condition, tie-breaking to the smallest profile_id; keep its full-view "
                     "X and control."),
    "full-pooled-label": ("Keep the full view's rows, X and control; replace each treated "
                          "profile's delta by its condition's merge-recipe pooled delta; "
                          "publish only the delta payload."),
}
CONTROL_RULE = {
    "full": ("mean of the matched_control_group vehicle rows, each normalized on its own "
             "selected-axis total"),
    "merge": ("normalize-of-sum over the union of the members' matched_control_group vehicle "
              "counts (never a mean of normalized values)"),
    "solo": ("mean of the matched_control_group vehicle rows, each normalized on its own "
             "selected-axis total (the full-view control)"),
    "solo-deepest": ("mean of the matched_control_group vehicle rows, each normalized on its "
                     "own selected-axis total (the full-view control)"),
    "full-pooled-label": ("full-view control; only the label is pooled, delta_kind=pooled"),
}
PRODUCT_FILES = ("profile_ids.csv", "members.csv", "control_sources.csv")


class HarnessError(ValueError):
    """A view-product build violates the harness contract."""


# --------------------------------------------------------------------------- #
# pinned Step 01 implementation
# --------------------------------------------------------------------------- #
_MODULES = {}


def repo_root(explicit=None):
    """The checkout that owns the pinned inputs and payloads."""
    if explicit is not None:
        root = Path(explicit).resolve()
        if not (root / "stage2/01_bases/output/api/bases/view.py").is_file():
            raise HarnessError(f"repo root does not carry the Bases API: {root}")
        return root
    for path in Path(__file__).resolve().parents:
        if (path / ".git").exists():
            return path
    raise HarnessError("cannot locate the repository root; pass --repo-root")


def bind(root):
    """Insert the owning steps' API roots and import the shared implementations."""
    if _MODULES and _MODULES["root"] == Path(root):
        return _MODULES
    root = Path(root)
    for relative in ("stage2/01_bases/output/api", "stage2/02_models_tenv-5k/output/api",
                    "stage2/04_tsg_dose/output/api"):
        path = str(root / relative)
        if path not in sys.path:
            sys.path.insert(0, path)
    from bases import corpus, gene_axis, release, view
    from tsg_round1 import protocol
    _MODULES.clear()
    _MODULES.update(root=root, corpus=corpus, gene_axis=gene_axis, release=release,
                    view=view, protocol=protocol)
    return _MODULES


def modules():
    if not _MODULES:
        raise HarnessError("the Bases implementation is not bound; build a Context first")
    return _MODULES


def relative(path, root):
    """Repo-relative posix path, else the bare file name (never an absolute path)."""
    try:
        return Path(path).absolute().relative_to(Path(root).absolute()).as_posix()
    except ValueError:
        return Path(path).name


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    """The repository's canonical deterministic JSON form."""
    modules()["release"].write_json(path, value)
    return path


def read_json(path):
    return modules()["release"].read_json(path)


def environment():
    def version(package):
        try:
            return metadata.version(package)
        except metadata.PackageNotFoundError:
            return ""

    return {"python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "h5py": version("h5py"),
            "anndata": version("anndata"),
            "data_access": ("offline; reads the pinned integer corpus, the G_5000 axis and "
                            "the accepted view-5k release in bounded windows"),
            "resource_contract": ("one heavy job at a time under the repository caps; "
                                  "OMP_NUM_THREADS=2; peak RSS recorded by the driver")}


# --------------------------------------------------------------------------- #
# build context
# --------------------------------------------------------------------------- #
class Context:
    """Resolved inputs, staged destination and installer for one product build."""

    def __init__(self, *, view, source, output_dir, staging_dir, repo_root_path=None,
                 size=SIZE, fixture=False, corpus_dir=None, axis_dir=None,
                 accepted_dir=None, name=None, sealed_dir=None, provenance_root="repo"):
        if view not in VIEWS:
            raise HarnessError(f"unknown view: {view}")
        if source not in SOURCES:
            raise HarnessError(f"unknown source: {source}")
        if int(size) <= 0:
            raise HarnessError("a view size must be positive")
        self.view = view
        self.source = source
        self.size = int(size)
        self.fixture = bool(fixture)
        self.name = name or NAME
        self.root = repo_root(explicit=repo_root_path)
        # `provenance_root` is what `bases.view.build` records its input keys against;
        # hermetic fixtures live outside the checkout and pass None.
        self.provenance_root = self.root if provenance_root == "repo" else provenance_root
        self.modules = bind(self.root)
        bases = self.root / "stage2/01_bases/output"
        self.corpus_dir = Path(corpus_dir or bases / "corpus")
        self.axis_dir = Path(axis_dir or bases / "gene-axis")
        self.accepted_dir = Path(accepted_dir or bases / f"view-{self.size // 1000}k")
        self.sealed_dir = Path(sealed_dir) if sealed_dir else None
        self.output_dir = Path(output_dir)
        self.staging_dir = Path(staging_dir)
        self.evidence = {}
        if not self.fixture:
            if self.size != SIZE:
                raise HarnessError(f"a released view size must be {SIZE}")
            if self.name != NAME:
                raise HarnessError(f"only the canonical {NAME} expression name is supported")

    # -- pinned inputs ---------------------------------------------------- #
    @property
    def axis_frame(self):
        frame = self.modules["gene_axis"].GeneAxis(self.axis_dir).axis_frame(self.size)
        if len(frame) != self.size:
            raise HarnessError(f"axis {self.size} holds {len(frame)} genes")
        return frame

    @property
    def genes(self):
        return np.asarray(self.axis_frame.gene_index, dtype=np.int64)

    @property
    def symbols(self):
        return self.axis_frame.symbol.to_numpy()

    def axis_binding(self):
        frame = self.axis_frame
        genes = np.asarray(frame.gene_index, dtype=np.int64)
        if len(np.unique(genes)) != len(genes) or not np.array_equal(genes, np.sort(genes)):
            raise HarnessError("the G_5000 axis is not distinct ascending gene indices")
        columns = hashlib.sha256(frame[["gene_index", "symbol"]].to_csv(
            index=False, lineterminator="\n").encode()).hexdigest()
        path = self.axis_dir / f"gene_axis_{self.size}.csv"
        return {"size": self.size, "file": relative(path, self.root),
                "algorithm": self.modules["gene_axis"].ALGORITHM,
                "sha256": sha256(path), "columns_sha256": columns,
                "gene_index_min": int(genes.min()), "gene_index_max": int(genes.max()),
                "ascending": True, "renumbered": False}

    def accepted_manifest(self):
        path = self.accepted_dir / "manifest.json"
        manifest = self.modules["release"].read_json(path)
        if manifest.get("schema") != "stage2-view-v3":
            raise HarnessError(f"accepted view schema differs: {path}")
        if int(manifest.get("view", {}).get("size", -1)) != self.size:
            raise HarnessError(f"accepted view size differs: {path}")
        return manifest

    def source_manifests(self):
        accepted = self.accepted_dir / "manifest.json"
        axis = self.axis_dir / f"gene_axis_{self.size}.csv"
        return {"accepted_view": relative(accepted, self.root),
                "accepted_view_sha256": sha256(accepted),
                "gene_axis": relative(axis, self.root),
                "gene_axis_sha256": sha256(axis)}

    @contextmanager
    def corpus(self, source=None):
        corpus = self.modules["corpus"].open_corpus(self.corpus_dir)
        reader = corpus.source(source or self.source)
        try:
            yield reader
        finally:
            reader.close()

    @contextmanager
    def accepted_view(self):
        release = self.modules["view"].open_view(self.accepted_dir)
        try:
            yield release
        finally:
            release.close()

    def accepted_source(self):
        """A live ``ViewSource`` for this source; close it with the returned context."""
        return _AcceptedSource(self)

    # -- output ----------------------------------------------------------- #
    def write_product(self, **kwargs):
        """Write the staged product and return its manifest."""
        kwargs.setdefault("view", self.view)
        kwargs.setdefault("source", self.source)
        kwargs.setdefault("size", self.size)
        axis = self.axis_binding()
        normalization = dict(self.modules["view"].NORMALIZATION)
        manifest = write_product(self.staging_dir, axis=axis,
                                 normalization=normalization,
                                 source_manifests=self.source_manifests(),
                                 genes=self.genes, symbols=self.symbols, **kwargs)
        return manifest

    def install(self):
        """Install the staged source files under ``<output>/<view>/``."""
        return install(self.staging_dir, self.output_dir / self.view, self.source)


class _AcceptedSource:
    """Context manager binding one accepted view source reader."""

    def __init__(self, context):
        self.context = context
        self.release = None

    def __enter__(self):
        self.release = self.context.modules["view"].open_view(self.context.accepted_dir)
        return self.release.source(self.context.source)

    def __exit__(self, *exc):
        if self.release is not None:
            self.release.close()
        return False


def install(staging_dir, output_dir, source):
    """Move one source's staged product files into the accepted product directory."""
    staging_dir, output_dir = Path(staging_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    installed = []
    for path in sorted(staging_dir.iterdir()):
        if not path.is_file() or not path.name.startswith(f"{source}."):
            continue
        target = output_dir / path.name
        if target.exists():
            target.unlink()
        os.replace(path, target)
        installed.append(target)
    if not installed:
        raise HarnessError(f"no staged product files for {source} in {staging_dir}")
    return installed


# --------------------------------------------------------------------------- #
# product writer (schema tsg-view-product-v1)
# --------------------------------------------------------------------------- #
def _record(path):
    return {"file": Path(path).name, "sha256": sha256(path),
            "bytes": Path(path).stat().st_size}


def _profile_ids(rows):
    frame = pd.DataFrame({"row": np.arange(len(rows), dtype=np.int64),
                          "profile_id": rows["profile_id"].astype(str).to_numpy()})
    return frame


def _retained(rows):
    totals = np.asarray(rows["library_total"], dtype=np.int64)
    if len(totals) and totals.min() < 0:
        raise HarnessError("a product row has a negative retained total")
    if not len(totals):
        return {"min": None, "median": None, "max": None, "sum": 0}
    return {"min": int(totals.min()), "median": float(np.median(totals)),
            "max": int(totals.max()), "sum": int(totals.sum(dtype=np.int64))}


def _condition_count(rows):
    treated = rows[rows["role"].astype(str) == "treated"]
    if not len(treated):
        return 0
    keys = treated[list(CONDITION_KEYS)].astype(str)
    keys["dose_um"] = treated["dose_um"].to_numpy(dtype=float)
    return int(len(keys.drop_duplicates()))


def _obs_frame(rows):
    frame = rows.copy()
    frame["profile_id"] = frame["profile_id"].astype(str)
    frame = frame.set_index("profile_id", drop=False)
    frame.index.name = "profile_id"
    return frame


def _var_frame(genes, symbols):
    var = pd.DataFrame({"gene_index": np.asarray(genes, dtype=np.int64),
                        "symbol": np.asarray(symbols)})
    var.index = pd.Index([str(gene) for gene in genes], name="gene_index_label")
    return var


def _h5ad(path, obs, var, uns, blocks, rows):
    import h5py
    from anndata.io import write_elem

    temporary = Path(path).with_suffix(Path(path).suffix + ".tmp")
    columns = len(var)
    with h5py.File(temporary, "w") as handle:
        handle.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        write_elem(handle, "obs", obs)
        write_elem(handle, "var", var)
        write_elem(handle, "uns", uns)
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
        expected = 0
        for start, values in blocks:
            if int(start) != expected:
                raise HarnessError("product blocks must cover the rows in order")
            if set(values) - {"X", "control", "delta"}:
                raise HarnessError("a product block may only carry X/control/delta")
            length = len(values.get("X", values.get("control", values.get("delta"))))
            if length <= 0:
                raise HarnessError("a product block is empty")
            x = values.get("X")
            control = values.get("control")
            delta = values.get("delta")
            if x is None or control is None:
                raise HarnessError("a materialized product block needs X and control")
            x = np.asarray(x, dtype=np.float32)
            control = np.asarray(control, dtype=np.float32)
            if x.shape != (length, columns) or control.shape != (length, columns):
                raise HarnessError(f"a product block has the wrong shape at row {start}")
            exact = x - control
            if delta is None:
                delta = exact
            else:
                delta = np.asarray(delta, dtype=np.float32)
                if delta.shape != (length, columns) or not np.array_equal(delta, exact):
                    raise HarnessError("a product delta must equal X - control exactly")
            arrays["X"][start:start + length] = x
            arrays["control"][start:start + length] = control
            arrays["delta"][start:start + length] = delta
            expected = start + length
        if expected != rows:
            raise HarnessError(f"product blocks cover {expected} of {rows} rows")
    temporary.replace(path)


def _delta_payload(path, blocks, rows, columns):
    temporary = Path(path).with_suffix(Path(path).suffix + ".tmp")
    matrix = np.lib.format.open_memmap(temporary, mode="w+", dtype=np.float32,
                                       shape=(rows, columns))
    expected = 0
    for start, values in blocks:
        if int(start) != expected or "delta" not in values:
            raise HarnessError("a delta payload needs ordered delta blocks")
        delta = np.asarray(values["delta"], dtype=np.float32)
        if delta.shape[1] != columns:
            raise HarnessError("a delta block has the wrong shape")
        matrix[start:start + len(delta)] = delta
        expected = start + len(delta)
    if expected != rows:
        raise HarnessError(f"delta blocks cover {expected} of {rows} rows")
    matrix.flush()
    del matrix
    temporary.replace(path)


def write_product(destination, *, view, source, size, rows, blocks, recipe, control_rule,
                  mode, axis, normalization, source_manifests, alias=None, members=None,
                  retained_totals=None, seed=None, selection=None, control_sources=None,
                  delta_kind=None, genes=None, symbols=None):
    """Write one ``tsg-view-product-v1`` manifest plus its declared files.

    ``mode`` is ``payload`` (h5ad with X/control/delta), ``delta`` (only a float32
    ``.npy`` delta payload, X/control aliased to the accepted release) or ``alias``
    (manifests only, no payload).
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    if view not in VIEWS:
        raise HarnessError(f"unknown view: {view}")
    if source not in SOURCES:
        raise HarnessError(f"unknown source: {source}")
    if mode not in ("payload", "delta", "alias"):
        raise HarnessError(f"unknown product mode: {mode}")
    if mode != MODES[view]:
        raise HarnessError(f"view {view} must publish in {MODES[view]} mode, not {mode}")
    if int(size) != int(axis["size"]):
        raise HarnessError("the product size differs from its axis binding")
    rows = pd.DataFrame(rows)
    missing = [name for name in REQUIRED_OBS if name not in rows]
    if missing:
        raise HarnessError(f"product rows are missing required columns: {missing}")
    if rows["profile_id"].astype(str).duplicated().any():
        raise HarnessError("product profile_id values must be unique")
    role = set(rows["role"].astype(str))
    if role - {"treated", "vehicle"}:
        raise HarnessError(f"unexpected product roles: {sorted(role)}")

    files = {}
    profile = _profile_ids(rows)
    profile_path = destination / f"{source}.profile_ids.csv"
    profile.to_csv(profile_path, index=False, lineterminator="\n")
    files["profile_ids"] = _record(profile_path)

    members_record = None
    if members is not None:
        members = pd.DataFrame(members)
        if list(members.columns) != ["profile_id", "member_profile_id"]:
            raise HarnessError("members must carry profile_id,member_profile_id")
        members_path = destination / f"{source}.members.csv"
        members.to_csv(members_path, index=False, lineterminator="\n")
        members_record = {**_record(members_path), "rows": int(len(members)),
                          "columns": list(members.columns)}
    elif view in ("merge", "solo", "solo-deepest"):
        raise HarnessError(f"view {view} must declare its member profile ids")

    controls_record = None
    if control_sources is not None:
        control_sources = pd.DataFrame(control_sources)
        required = ["profile_id", "control_rule", "control_groups", "pooled_retained_total"]
        if list(control_sources.columns) != required:
            raise HarnessError(f"control_sources must carry {required}")
        if len(control_sources) != len(rows):
            raise HarnessError("control_sources must declare one row per product row")
        controls_path = destination / f"{source}.control_sources.csv"
        control_sources.to_csv(controls_path, index=False, lineterminator="\n")
        controls_record = {**_record(controls_path), "rows": int(len(control_sources)),
                           "columns": required}

    payload = None
    if mode == "payload":
        if genes is None or symbols is None:
            raise HarnessError("a payload product needs the axis gene/symbol vectors")
        genes = np.asarray(genes, dtype=np.int64)
        symbols = np.asarray(symbols)
        if len(genes) != int(axis["size"]) or len(symbols) != len(genes):
            raise HarnessError("the axis gene/symbol vectors disagree with the binding")
        obs = _obs_frame(rows)
        var = _var_frame(genes, symbols)
        uns = _uns(view, source, size, axis, normalization, recipe, control_rule)
        path = destination / f"{source}.h5ad"
        _h5ad(path, obs, var, uns, blocks, len(rows))
        payload = {**_record(path), "rows": int(len(rows)), "cols": int(len(var)),
                   "dtype": "float32", "layers": ["control", "delta"]}
    elif mode == "delta":
        if delta_kind != "pooled":
            raise HarnessError("a delta product must declare delta_kind=pooled")
        path = destination / f"{source}.delta.npy"
        _delta_payload(path, blocks, len(rows), int(axis["size"]))
        payload = {**_record(path), "rows": int(len(rows)), "cols": int(axis["size"]),
                   "dtype": "float32", "layers": ["delta"], "kind": "delta"}
    if mode == "alias" and not alias:
        raise HarnessError("an alias product must declare the accepted release digests")
    if alias is not None:
        alias = dict(alias)
        release_path = alias.get("release", "")
        if not release_path or release_path.startswith("/"):
            raise HarnessError("an alias must record the accepted release repo-relative")
        for name, digest in alias.get("payload_sha256", {}).items():
            if not isinstance(digest, str) or len(digest) != 64:
                raise HarnessError(f"alias payload digest is invalid for {name}")
        if len(str(alias.get("manifest_sha256", ""))) != 64:
            raise HarnessError("an alias must record the accepted manifest sha256")
    if mode == "alias" and blocks is not None:
        raise HarnessError("an alias product does not consume expression blocks")

    manifest = {
        "schema": SCHEMA,
        "generated_by": "stage2/tsg_view/execution/build.py",
        "product": f"{source}-{view}-5k-view",
        "view": view,
        "source": source,
        "rows": int(len(rows)),
        "treated": int((rows["role"].astype(str) == "treated").sum()),
        "conditions": _condition_count(rows),
        "axis": dict(axis),
        "normalization": dict(normalization),
        "recipe": recipe,
        "control_rule": control_rule,
        "seed": seed,
        "selection": selection,
        "retained_totals": retained_totals or _retained(rows),
        "profile_ids": files["profile_ids"],
        "members": members_record,
        "control_sources": controls_record,
        "source_manifests": dict(source_manifests),
        "payload": payload,
        "alias": alias,
        "delta_kind": delta_kind,
        "environment": environment(),
    }
    write_json(destination / f"{source}.json", manifest)
    return manifest


def _uns(view, source, size, axis, normalization, recipe, control_rule):
    return {
        "schema": SCHEMA,
        "view_variant": view,
        "source": source,
        "view": {"name": NAME, "size": int(size)},
        "axis": axis,
        "normalization": normalization,
        "control": {"rule": control_rule, "delta": "delta = X - control, elementwise"},
        "recipe": recipe,
    }


# --------------------------------------------------------------------------- #
# derived-view primitives (C7.1-C7.5)
# --------------------------------------------------------------------------- #
def condition_key(cell_line, parent, dose):
    return f"{cell_line}|{parent}|{float(dose):g}"


def merge_profile_id(source, cell_line, parent, dose):
    """The deterministic derived-row identity locked by plan C7.4."""
    return f"merge::{source}::{cell_line}::{parent}::{float(dose):g}"


def condition_table(obs, source):
    """Canonical condition table plus a view-row -> condition index vector.

    Treated rows only. Conditions are ordered by ``(cell_line, parent, dose)`` and their
    members by ``profile_id``, so the table is independent of corpus row order, of folds
    and of training seeds.
    """
    obs = pd.DataFrame(obs)
    if "library_total" not in obs:
        raise HarnessError(
            "condition_table needs view obs rows (with library_total), not corpus obs")
    obs["profile_id"] = obs["profile_id"].astype(str) if "profile_id" in obs else obs.index.astype(str)
    treated = np.flatnonzero(obs["role"].astype(str).to_numpy() == "treated")
    frame = obs.iloc[treated].copy().reset_index(drop=True)
    frame["row"] = treated
    ordered = frame.sort_values(list(CONDITION_KEYS) + ["profile_id"], kind="mergesort")
    conditions = []
    row_to_condition = np.full(len(obs), -1, dtype=np.int64)
    position = 0
    for key, members in ordered.groupby(list(CONDITION_KEYS), sort=False, dropna=False):
        rows = members["row"].to_numpy(dtype=np.int64)
        groups = sorted({str(group) for group in members["matched_control_group"]})
        conditions.append({
            "condition_index": len(conditions),
            "condition_key": condition_key(*key),
            "cell_line": key[0],
            "drug_parent_inchikey": key[1],
            "dose_um": float(key[2]),
            "members": rows,
            "profile_ids": members["profile_id"].to_numpy(),
            "member_count": int(len(rows)),
            "control_groups": groups,
            "library_totals": members["library_total"].to_numpy(dtype=np.int64),
            "cell_line_id": str(members["cell_line_id"].iloc[0]),
            "drug_parent_smiles": str(members["drug_parent_smiles"].iloc[0]),
            "plate": str(members["plate"].iloc[0]),
            "plates": sorted({str(plate) for plate in members["plate"]}),
            "compound_names": sorted({str(name) for name in members["compound_name"]}),
        })
        row_to_condition[rows] = position
        position += 1
    if not conditions:
        raise HarnessError(f"{source}: no treated condition")
    return pd.DataFrame(conditions), row_to_condition


def merge_rows(conditions, source, retained_totals, source_manifest_sha256):
    """Canonical freeze-shaped obs rows for the replicate-pooled view (C7.1/C7.4)."""
    rows = []
    for condition in conditions.itertuples():
        rows.append({
            "profile_id": merge_profile_id(source, condition.cell_line,
                                           condition.drug_parent_inchikey,
                                           condition.dose_um),
            "source_item": "tsg_view:merge",
            "source_item_manifest_sha256": source_manifest_sha256,
            "source_row": -1,
            "role": "treated",
            "plate": ";".join(condition.plates),
            "cell_line": condition.cell_line,
            "cell_line_id": condition.cell_line_id,
            "compound_name": ";".join(condition.compound_names),
            "vehicle": "",
            "dose_um": float(condition.dose_um),
            "dose_original_value": float(condition.dose_um),
            "dose_original_unit": "uM",
            "drug_parent_inchikey": condition.drug_parent_inchikey,
            "drug_parent_smiles": condition.drug_parent_smiles,
            "source_inchikey": "",
            "n_cells": 0,
            "matched_control_group": ";".join(sorted(set(condition.control_groups))),
            "library_total": int(retained_totals[condition.condition_index]),
            "member_count": int(condition.member_count),
        })
    return pd.DataFrame(rows)


def member_table(conditions, profile_id=None):
    """Long-form ``profile_id,member_profile_id`` table for a derived product."""
    rows = []
    for condition in conditions.itertuples():
        identity = profile_id(condition) if profile_id else condition.condition_key
        for member in condition.profile_ids:
            rows.append({"profile_id": identity, "member_profile_id": str(member)})
    return pd.DataFrame(rows, columns=["profile_id", "member_profile_id"])


def pooled_condition_counts(reader, genes, row_to_condition, conditions):
    """Sum member retained integer counts per condition in one bounded corpus pass.

    ``float64`` accumulation is exact for at most 16 addends of small integer counts and
    keeps the accumulator at 1.06 GB for 52,777 x 5,000, instead of a 2.1 GB int64 array.
    """
    totals = np.zeros((len(conditions), len(genes)), dtype=np.float64)
    for start, block in reader.batches(BATCH):
        stop = start + len(block)
        index = row_to_condition[start:stop]
        keep = index >= 0
        if keep.any():
            np.add.at(totals, index[keep], block[keep][:, genes])
    return totals


def pooled_control(conditions, group_sums):
    """Normalize-of-sum pooled control per condition (C7.2)."""
    view = modules()["view"]
    control = np.empty((len(conditions), len(next(iter(group_sums.values()))[0])),
                       dtype=np.float32)
    retained = np.empty(len(conditions), dtype=np.int64)
    for position, condition in enumerate(conditions.itertuples()):
        groups = sorted(set(condition.control_groups))
        missing = [group for group in groups if group not in group_sums]
        if missing:
            raise HarnessError(f"matched control groups without vehicles: {missing[:3]}")
        counts = np.zeros_like(next(iter(group_sums.values()))[0], dtype=np.int64)
        for group in groups:
            counts += group_sums[group][0]
        total = int(counts.sum(dtype=np.int64))
        control[position] = view.log1p_cp10k(counts.reshape(1, -1),
                                             np.asarray([total], dtype=np.int64))[0].astype(
                                                 np.float32)
        retained[position] = total
    return control, retained


def vehicle_group_sums(reader, genes):
    """Sum retained integer counts and totals per matched control group (one pass)."""
    obs = reader.obs
    role = obs["role"].astype(str).to_numpy()
    groups = obs["matched_control_group"].astype(str).to_numpy()
    keys = sorted({group for group, is_vehicle in zip(groups, role == "vehicle") if is_vehicle})
    position = {key: index for index, key in enumerate(keys)}
    counts = np.zeros((len(keys), len(genes)), dtype=np.int64)
    for start, block in reader.batches(BATCH):
        stop = start + len(block)
        local = np.flatnonzero(role[start:stop] == "vehicle")
        for offset in local:
            counts[position[groups[start + offset]]] += block[offset][genes].astype(np.int64)
    totals = counts.sum(axis=1, dtype=np.int64)
    if len(totals) and (totals <= 0).any():
        raise HarnessError("a matched control group has a non-positive retained total")
    return {key: (counts[index], int(totals[index])) for index, key in enumerate(keys)}


def pooled_condition_blocks(pooled, control):
    """Stream ``X/control/delta`` blocks for materialized pooled rows."""
    view = modules()["view"]
    rows = len(pooled)
    for start in range(0, rows, BATCH):
        stop = min(start + BATCH, rows)
        counts = np.asarray(pooled[start:stop]).astype(np.int64)
        totals = view.library_totals(counts)
        x = view.log1p_cp10k(counts, totals).astype(np.float32)
        block_control = np.asarray(control[start:stop], dtype=np.float32)
        yield start, {"X": x, "control": block_control, "delta": x - block_control}


def solo_choice(conditions, rule="draw", seed=SOLO_SEED, source="tenv"):
    """Choose one member per condition, independent of folds and training seeds."""
    if rule not in ("draw", "deepest"):
        raise HarnessError(f"unknown solo selection rule: {rule}")
    chosen = np.empty(len(conditions), dtype=np.int64)
    for position, condition in enumerate(conditions.itertuples()):
        members = np.asarray(condition.members, dtype=np.int64)
        if rule == "deepest":
            order = sorted(range(len(members)),
                           key=lambda index: (-int(condition.library_totals[index]),
                                              str(condition.profile_ids[index])))
            chosen[position] = members[order[0]]
            continue
        if len(members) == 1:
            chosen[position] = members[0]
            continue
        digest = hashlib.sha256(
            f"{seed}|{source}|{condition.condition_key}".encode()).digest()
        chosen[position] = members[int.from_bytes(digest[:8], "big") % len(members)]
    return chosen


def selected_blocks(reader, positions, keys=("X", "control", "delta")):
    """Stream the chosen full-view rows in target order from bounded h5py reads."""
    import h5py

    positions = np.asarray(positions, dtype=np.int64)
    datasets = {"X": "X", "control": "layers/control", "delta": "layers/delta"}
    for key in keys:
        if key not in datasets:
            raise HarnessError(f"unknown full-view layer: {key}")
    with h5py.File(reader.path, "r") as handle:
        for start in range(0, len(positions), BATCH):
            chunk = positions[start:start + BATCH]
            order = np.argsort(chunk, kind="stable")
            sorted_positions = chunk[order]
            values = {}
            for key in keys:
                block = np.asarray(handle[datasets[key]][list(sorted_positions)],
                                   dtype=np.float32)
                restored = np.empty_like(block)
                restored[order] = block
                values[key] = restored
            yield start, values


def selected_rows(reader, positions):
    """The chosen full-view obs rows verbatim (so ``solo`` keeps its member identity)."""
    obs = pd.DataFrame(reader.obs)
    obs["profile_id"] = obs.index.astype(str)
    return obs.iloc[np.asarray(positions, dtype=np.int64)].copy()


def selected_groups(reader, positions):
    """Per-row matched control group strings for the chosen full-view rows."""
    obs = pd.DataFrame(reader.obs)
    return obs.iloc[np.asarray(positions, dtype=np.int64)]["matched_control_group"].astype(
        str).tolist()


def control_source_table(conditions, profile_ids, rule, retained=None):
    """Per-row ``control_sources.csv`` provenance rows (C7.6)."""
    rows = []
    for position, condition in enumerate(conditions.itertuples()):
        rows.append({
            "profile_id": str(profile_ids[position]),
            "control_rule": rule,
            "control_groups": ";".join(sorted(set(condition.control_groups))),
            "pooled_retained_total": (""
                                      if retained is None else int(retained[position])),
        })
    return pd.DataFrame(rows, columns=["profile_id", "control_rule", "control_groups",
                                       "pooled_retained_total"])
