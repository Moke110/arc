"""Prepare and verify the ordered CM-axis transcriptome input release."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


SOURCES = ("tahoe", "tenv")
COUNTS = {"tahoe": {"raw": 67018, "blank": 0, "treated": 65618, "control": 1400},
          "tenv": {"raw": 4510, "blank": 32, "treated": 4206, "control": 272}}
AXIS_SHA256 = "19ce44acfd615a601d5bc6617074941e13d51cd6e527c8d519174e5085385f42"
N_GENES = 18164
FAM153B = "ENSG00000182230"
PROFILE_FIELDS = ["source_row", "profile_id", "plate", "is_control", "development"]
ROW_BATCH = 64
SCHEMA_VERSION = 1
NORMALIZATION = "log1p(CP10K of raw counts projected to the ordered CM axis)"


def sha256(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def source_files(corpus):
    return {"axis": corpus / "tables/genes.csv",
            **{source: corpus / f"{source}.h5ad" for source in SOURCES}}


def read_axis(path, fixture):
    if not fixture and sha256(path) != AXIS_SHA256:
        raise ValueError("gene axis differs from the frozen ordered CM axis")
    axis = pd.read_csv(path, dtype=str, keep_default_na=False)
    for field in ("symbol", "hgnc_id", "ensembl_id"):
        if field not in axis or axis[field].eq("").any() or axis[field].duplicated().any():
            raise ValueError(f"invalid ordered CM gene axis: {field}")
    if (not fixture and len(axis) != N_GENES) or (
        axis.loc[axis.symbol.eq("FAM153B"), "ensembl_id"].tolist() != [FAM153B]
    ):
        raise ValueError("invalid ordered CM gene axis or FAM153B identity")
    return axis


def gene_projection(axis, var, source):
    field = "ensembl_id" if source == "tahoe" else "hgnc_id"
    targets = set(axis[field])
    native = {}
    for col, identity in enumerate(var[field]):
        if identity in targets:
            if identity in native:
                raise ValueError(f"{source}: ambiguous gene mapping for {identity}")
            native[identity] = col
    expected_missing = {FAM153B} if source == "tahoe" else set()
    if targets - native.keys() != expected_missing:
        raise ValueError(f"{source}: unexpected gene mapping gaps: {sorted(targets - native.keys())[:10]}")
    observed = axis[field].isin(native).to_numpy()
    columns = [native[gene] for gene in axis.loc[observed, field]]
    return observed, columns


def plate_partition(plates, source, seed):
    # A source-specific hash rank is stable across Python and NumPy versions.
    plates = sorted(set(plates))
    if len(plates) < 2:
        raise ValueError(f"{source}: need at least two plates for development")
    ranked = sorted(plates, key=lambda plate: hashlib.sha256(
        json.dumps([seed, source, plate]).encode()).digest())
    validation = set(ranked[:max(1, len(plates) // 5)])
    return {"train": sorted(set(plates) - validation), "validation": sorted(validation)}


def source_profiles(obs, source, seed, fixture):
    index_name = "profile_id" if source == "tahoe" else "well"
    if obs.index.name != index_name or obs.index.has_duplicates or obs.index.hasnans:
        raise ValueError(f"{source}: invalid source profile identities")
    if obs["plate"].isna().any() or obs["plate"].eq("").any():
        raise ValueError(f"{source}: missing plate identity")
    for field in ("is_control", "is_blank"):
        if field not in obs and field == "is_blank" and source == "tahoe":
            continue  # Tahoe's Step 02 schema has no technical blanks.
        if field not in obs or obs[field].dtype != np.dtype(bool):
            raise ValueError(f"{source}: missing or invalid {field}")
    blank = obs.get("is_blank", pd.Series(False, index=obs.index)).to_numpy()
    rows = np.flatnonzero(~blank)
    profiles = pd.DataFrame({
        "source_row": rows,
        "profile_id": obs.index[rows].astype(str),
        "plate": obs["plate"].iloc[rows].astype(str).to_numpy(),
        "is_control": obs["is_control"].to_numpy()[rows],
    })
    counts = population_counts(profiles, len(obs))
    if not fixture and counts != COUNTS[source]:
        raise ValueError(f"{source}: unexpected source population {counts}")
    split = plate_partition(profiles.plate, source, seed)
    profiles["development"] = np.where(profiles.plate.isin(split["validation"]), "validation", "train")
    return profiles


def population_counts(profiles, raw):
    controls = int(profiles.is_control.sum())
    return {"raw": raw, "blank": raw - len(profiles), "treated": len(profiles) - controls,
            "control": controls}


def normalized_batches(matrix, shape, rows, observed, columns, source):
    sparse = isinstance(matrix, h5py.Group)
    if sparse:
        if matrix.attrs.get("encoding-type") != "csr_matrix" or tuple(matrix.attrs["shape"]) != shape:
            raise ValueError(f"{source}: expected native-axis CSR raw counts")
        ptr = matrix["indptr"][:]
        data, indices = matrix["data"], matrix["indices"]
        if (ptr.shape != (shape[0] + 1,) or ptr[0] != 0 or np.any(np.diff(ptr) < 0)
                or ptr[-1] != len(data) or len(indices) != len(data)):
            raise ValueError(f"{source}: invalid CSR pointers")
        dtype = data.dtype
    else:
        if matrix.shape != shape:
            raise ValueError(f"{source}: raw counts shape differs from metadata")
        dtype = matrix.dtype
    if dtype.kind not in "iu":
        raise ValueError(f"{source}: expected integer raw counts")

    # Window by source row, so even long stretches of blanks cannot expand a read.
    boundaries = np.flatnonzero(np.diff(rows // ROW_BATCH)) + 1
    for requested in np.split(rows, boundaries):
        start, stop = int(requested[0]), int(requested[-1]) + 1
        if sparse:
            counts = data[ptr[start]:ptr[stop]]
            col_ids = indices[ptr[start]:ptr[stop]]
            if np.any(col_ids < 0) or np.any(col_ids >= shape[1]):
                raise ValueError(f"{source}: invalid CSR column indices")
            native = csr_matrix((counts, col_ids, ptr[start:stop + 1] - ptr[start]),
                                shape=(stop - start, shape[1]))
            projected = native[requested - start][:, columns].toarray()
        else:
            counts = matrix[start:stop]
            projected = counts[requested - start][:, columns]
        if np.any(counts < 0):
            raise ValueError(f"{source}: negative raw counts in rows {start}:{stop}")
        total = projected.sum(axis=1, dtype=np.float64)
        if np.any(total <= 0):
            row = requested[np.flatnonzero(total <= 0)[0]]
            raise ValueError(f"{source}: zero projected library size at source row {row}")
        values = np.zeros((len(requested), len(observed)), dtype=np.float32)
        values[:, observed] = projected * (10000 / total[:, None])
        np.log1p(values, out=values)
        if not np.isfinite(values).all():
            raise ValueError(f"{source}: non-finite normalized expression in rows {start}:{stop}")
        yield values


def release_files(directory, source):
    return (directory / f"{source}_expression.npy", directory / f"{source}_target_mask.npy",
            directory / f"{source}_profiles.csv")


def membership_summary(profiles, source, seed):
    plates = plate_partition(profiles.plate, source, seed)
    development = {}
    for role, role_plates in plates.items():
        members = profiles.loc[profiles.plate.isin(role_plates)]
        if not members.development.eq(role).all():
            raise ValueError(f"{source}: invalid development plate partition")
        development[role] = {
            "plates": role_plates, "profiles": len(members),
            "membership_sha256": membership_digest(members),
        }
    return {"development": development, "refit": {"profiles": len(profiles),
                                                   "membership_sha256": membership_digest(profiles)}}


def membership_digest(profiles):
    rows = profiles[["source_row", "profile_id", "plate"]]
    return hashlib.sha256(rows.to_csv(index=False, lineterminator="\n").encode()).hexdigest()


def write_source(path, source, directory, axis, seed, fixture):
    expression_path, mask_path, profiles_path = release_files(directory, source)
    with h5py.File(path, "r") as file:
        obs, var = read_elem(file["obs"]), read_elem(file["var"])
        profiles = source_profiles(obs, source, seed, fixture)
        observed, columns = gene_projection(axis, var, source)
        with expression_path.open("wb") as expression:
            # A standard .npy header followed by contiguous row batches gives
            # consumers a single memory-mappable array without caching it in RAM.
            np.lib.format.write_array_header_2_0(expression, {
                "descr": np.dtype(np.float32).str, "fortran_order": False,
                "shape": (len(profiles), len(axis)),
            })
            for values in normalized_batches(file["X"], (len(obs), len(var)),
                                             profiles.source_row.to_numpy(), observed, columns, source):
                values.tofile(expression)
    np.save(mask_path, observed, allow_pickle=False)
    profiles.to_csv(profiles_path, index=False, lineterminator="\n")
    return {"population": population_counts(profiles, len(obs)),
            **membership_summary(profiles, source, seed)}


def prepare(corpus, release, seed, fixture, source_pins):
    if release.exists():
        raise ValueError(f"release already exists: {release}")
    files = source_files(corpus)
    before = {name: path.stat() for name, path in files.items()}
    axis = read_axis(files["axis"], fixture)
    hashes = {name: sha256(path) for name, path in files.items()}
    if not fixture and any(hashes[source] != source_pins[source] for source in SOURCES):
        raise ValueError("source corpus differs from the frozen Step 02 file digests")
    release.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=release.parent, prefix=".input-") as staging:
        directory = Path(staging)
        shutil.copyfile(files["axis"], directory / "genes.csv")
        sources = {source: write_source(files[source], source, directory, axis, seed, fixture)
                   for source in SOURCES}
        for name, path in files.items():
            after = path.stat()
            if (after.st_ino, after.st_size, after.st_mtime_ns) != (
                before[name].st_ino, before[name].st_size, before[name].st_mtime_ns
            ):
                raise ValueError(f"source changed while preparing: {path}")
        manifest = {
            "schema_version": SCHEMA_VERSION, "fixture": fixture, "split_seed": seed,
            "normalization": NORMALIZATION, "gene_count": len(axis),
            "axis_sha256": hashes["axis"],
            "source_sha256": {source: hashes[source] for source in SOURCES},
            "sources": sources,
            "files": {path.name: sha256(path) for path in sorted(directory.iterdir())},
        }
        (directory / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        os.rename(directory, release)
    return manifest


def load_release(release):
    """Read the release manifest. Call verify() to check the payload."""
    manifest = json.loads((release / "manifest.json").read_text())
    if manifest["schema_version"] != SCHEMA_VERSION or manifest["normalization"] != NORMALIZATION:
        raise ValueError("unsupported input release schema or normalization")
    return manifest


def verify(release, corpus, source_pins):
    manifest = load_release(release)
    fixture = manifest["fixture"]
    if type(fixture) is not bool or type(manifest["split_seed"]) is not int:
        raise ValueError("invalid release mode or split seed")
    axis = read_axis(release / "genes.csv", fixture)
    expected_files = {"genes.csv", *(path.name for source in SOURCES
                                    for path in release_files(release, source))}
    if set(manifest["files"]) != expected_files:
        raise ValueError("invalid release file list")
    for name, expected in manifest["files"].items():
        if sha256(release / name) != expected:
            raise ValueError(f"release file digest mismatch: {name}")
    if len(axis) != manifest["gene_count"] or manifest["files"]["genes.csv"] != manifest["axis_sha256"]:
        raise ValueError("release gene axis differs from its manifest")
    source_hashes = manifest["source_sha256"]
    if set(source_hashes) != set(SOURCES) or any(
        not re.fullmatch(r"[0-9a-f]{64}", digest) for digest in source_hashes.values()
    ):
        raise ValueError("invalid release source digests")
    if not fixture and source_hashes != source_pins:
        raise ValueError("release source digests differ from the frozen Step 02 file digests")
    if corpus is not None:
        expected = {"axis": manifest["axis_sha256"], **manifest["source_sha256"]}
        for name, path in source_files(corpus).items():
            if sha256(path) != expected[name]:
                raise ValueError(f"source digest mismatch: {name}")
    for source in SOURCES:
        verify_source(release, axis, source, manifest)
    return manifest


def verify_source(release, axis, source, manifest):
    expression_path, mask_path, profiles_path = release_files(release, source)
    profiles = pd.read_csv(profiles_path, dtype={"profile_id": str, "plate": str},
                           keep_default_na=False)
    summary = manifest["sources"][source]
    raw = summary["population"]["raw"]
    rows = profiles.source_row.to_numpy()
    if (type(raw) is not int or list(profiles) != PROFILE_FIELDS
            or profiles.empty or rows.dtype.kind not in "iu"
            or rows[0] < 0 or rows[-1] >= raw or np.any(np.diff(rows) <= 0)
            or profiles.profile_id.eq("").any() or profiles.profile_id.duplicated().any()
            or profiles.plate.eq("").any() or profiles.is_control.dtype != np.dtype(bool)):
        raise ValueError(f"{source}: invalid release profile identities")
    population = population_counts(profiles, raw)
    if population != summary["population"] or (not manifest["fixture"] and population != COUNTS[source]):
        raise ValueError(f"{source}: release population mismatch")
    if summary != {"population": population,
                   **membership_summary(profiles, source, manifest["split_seed"])}:
        raise ValueError(f"{source}: release membership mismatch")

    mask = np.load(mask_path, allow_pickle=False)
    expected_mask = (axis.ensembl_id.ne(FAM153B).to_numpy() if source == "tahoe"
                     else np.ones(len(axis), dtype=bool))
    if mask.dtype != np.dtype(bool) or not np.array_equal(mask, expected_mask):
        raise ValueError(f"{source}: invalid reconstruction target mask")
    with expression_path.open("rb") as expression:
        if np.lib.format.read_magic(expression) != (2, 0):
            raise ValueError(f"{source}: expected expression in NPY 2.0 format")
        shape, fortran_order, dtype = np.lib.format.read_array_header_2_0(expression)
        if shape != (len(profiles), len(axis)) or fortran_order or dtype != np.dtype(np.float32):
            raise ValueError(f"{source}: invalid expression shape, order or dtype")
        for start in range(0, len(profiles), ROW_BATCH):
            size = min(len(profiles) - start, ROW_BATCH) * len(axis)
            values = np.fromfile(expression, dtype=dtype, count=size)
            if len(values) != size:
                raise ValueError(f"{source}: truncated expression at release row {start}")
            batch = values.reshape(-1, len(axis))
            if (not np.isfinite(batch).all() or np.any(batch < 0) or np.any(batch[:, ~mask] != 0)
                    or not np.allclose(np.expm1(batch).sum(axis=1), 10000, rtol=1e-5)):
                raise ValueError(f"{source}: invalid normalized expression at release row {start}")
        if expression.read(1):
            raise ValueError(f"{source}: unexpected bytes after expression values")
