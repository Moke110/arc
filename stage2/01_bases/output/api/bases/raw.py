"""Build and independently replay the ADR 0052 shared-axis expression release."""

from collections import defaultdict
import json
import os
from pathlib import Path
import shutil
import tempfile

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from .native import digest as file_sha256
from .identity import AXIS_SHA256


SOURCES = ("tahoe", "tenv")
POPULATIONS = {
    "tahoe": {"raw": 67018, "blank": 0, "treated": 65618, "control": 1400},
    "tenv": {"raw": 4510, "blank": 32, "treated": 4206, "control": 272},
}
AXIS_SIZE = 18164
FAM153B = "HGNC:27323"
ROW_BATCH = 64
SCHEMA_VERSION = 2
NORMALIZATION = "log1p(CP10K of raw counts projected to the ordered CM axis)"
INPUT_NAMES = ("genes.csv", "gene_crosswalk.csv", "tahoe.h5ad", "tenv.h5ad")
PAYLOAD_NAMES = (
    "genes.csv",
    *(f"{source}_{suffix}" for source in SOURCES
      for suffix in ("expression.npy", "profiles.csv", "mapping.csv")),
)


def read_csv(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def frozen_digests():
    """The shared gene axis is frozen; rebuilt native files have fresh source seals."""
    return {"genes.csv": AXIS_SHA256}


def input_files(native, sources=SOURCES):
    """Resolve explicit per-source files, retaining the legacy directory interface."""
    names = ("genes.csv", "gene_crosswalk.csv", *[f"{source}.h5ad" for source in sources])
    if isinstance(native, dict):
        return {name: Path(native[name]) for name in names}
    native = Path(native)
    return {name: native / ("tables" if name.endswith(".csv") else "") / name
            for name in names}


def default_bindings(root, *, native=None, tables=None, tenv=None, tahoe=None):
    """Bind independently owned raw releases without copying their payloads."""
    output = Path(root) / "stage2/01_bases/output"
    legacy = Path(native) if native is not None else output / "tenv-raw-counts"
    if tables is not None:
        tables = Path(tables)
    elif native is not None:
        tables = legacy / "tables"
    else:
        tables = output / "gene-identity/tables"
    return {"genes.csv": tables / "genes.csv",
            "gene_crosswalk.csv": tables / "gene_crosswalk.csv",
            "tenv.h5ad": Path(tenv) if tenv is not None else legacy / "tenv.h5ad",
            "tahoe.h5ad": Path(tahoe) if tahoe is not None else
                (legacy / "tahoe.h5ad" if native is not None else
                 output / "tahoe-raw-counts/tahoe.h5ad")}


def add_binding_arguments(parser):
    parser.add_argument("--native", type=Path, help="legacy directory containing both raw h5ads and tables")
    parser.add_argument("--tables", type=Path, help="independently owned genes/crosswalk directory")
    parser.add_argument("--tenv", type=Path, help="raw tenv h5ad")
    parser.add_argument("--tahoe", type=Path, help="raw projected Tahoe h5ad")


def cli_bindings(args, root):
    return default_bindings(root, native=args.native, tables=args.tables,
                            tenv=args.tenv, tahoe=args.tahoe)


def source_provenance(files, sources, fixture=False):
    evidence = {}
    for source in sources:
        raw = files[f"{source}.h5ad"]
        manifest = raw.parent / ("tenv_manifest.json" if source == "tenv" else "manifest.json")
        if not fixture and not manifest.is_file():
            raise ValueError(f"{source}: raw release manifest is missing: {manifest}")
        evidence[source] = {"raw_file": str(raw.resolve()),
                            "release_manifest": str(manifest.resolve()),
                            "release_manifest_sha256": file_sha256(manifest) if manifest.is_file() else None}
    return evidence


def fingerprint(files):
    return {name: (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
            for name, path in files.items() for stat in [path.stat()]}


def read_inputs(files, fixture):
    hashes = {name: file_sha256(path) for name, path in files.items()}
    if not fixture:
        for name, digest in frozen_digests().items():
            if hashes[name] != digest:
                raise ValueError(f"input differs from frozen bases digest: {name}")
    axis = read_csv(files["genes.csv"])
    if axis.empty or (not fixture and len(axis) != AXIS_SIZE):
        raise ValueError("invalid ordered CM gene count")
    for field in ("symbol", "hgnc_id", "ensembl_id"):
        if axis[field].eq("").any() or axis[field].duplicated().any():
            raise ValueError(f"invalid ordered CM axis: {field}")
    crosswalk = read_csv(files["gene_crosswalk.csv"])
    return hashes, axis, crosswalk


def unique_native(values, label):
    ids = pd.Series(values, dtype="string")
    if ids.isna().any() or ids.eq("").any() or ids.duplicated().any():
        raise ValueError(f"{label}: missing or ambiguous native identities")
    return {identity: col for col, identity in enumerate(ids)}


def tahoe_columns(axis, var, crosswalk):
    native = unique_native(var["ensembl_id"], "Tahoe Ensembl")
    unique_native(var["token_id"], "Tahoe token")
    tahoe = crosswalk[crosswalk.dataset_id.eq("tahoe")].set_index("source_gene_id")
    if tahoe.index.has_duplicates:
        raise ValueError("Tahoe crosswalk has ambiguous source identities")

    # Use both native namespaces, with one native column per fallback candidate.
    candidates = defaultdict(list)
    identities = {}
    for col, (ensembl, token) in enumerate(zip(var.ensembl_id, var.token_id.astype(str), strict=True)):
        if ensembl not in tahoe.index or token not in tahoe.index:
            raise ValueError(f"Tahoe crosswalk missing native identity: {ensembl}, {token}")
        hgnc = tahoe.at[ensembl, "hgnc_id"]
        if hgnc != tahoe.at[token, "hgnc_id"]:
            raise ValueError(f"Tahoe crosswalk has conflicting Ensembl and token identities: {ensembl}, {token}")
        identities[ensembl] = hgnc
        if hgnc:
            candidates[hgnc].append(col)

    # An exact native ID may lack an HGNC annotation in Tahoe (e.g. ANKRD11-AS2).
    # Confirm the identity using the same Ensembl ID in the pinned crosswalk.
    evidence = defaultdict(set)
    known = crosswalk[crosswalk.hgnc_id.ne("")]
    for ensembl, hgnc in zip(known.source_gene_id, known.hgnc_id, strict=True):
        evidence[ensembl].add(hgnc)
    columns, methods = [], []
    for ensembl, hgnc in zip(axis.ensembl_id, axis.hgnc_id, strict=True):
        if ensembl in native:
            if evidence[ensembl] - {hgnc}:
                raise ValueError(f"Tahoe exact Ensembl identity conflicts with {hgnc}: {ensembl}")
            if not evidence[ensembl]:
                raise ValueError(f"Tahoe exact Ensembl match lacks HGNC identity evidence: {ensembl}")
            columns.append(native[ensembl])
            methods.append("exact_ensembl" if identities[ensembl] else "exact_ensembl_corroborated")
        else:
            matches = candidates[hgnc]
            if len(matches) != 1:
                reason = "missing" if not matches else "ambiguous"
                raise ValueError(f"Tahoe {reason} HGNC crosswalk mapping: {hgnc}")
            match = matches[0]
            mapped_id = var.ensembl_id.iloc[match]
            if evidence[mapped_id] != {hgnc}:
                raise ValueError(f"Tahoe crosswalk Ensembl identity conflicts with {hgnc}: {mapped_id}")
            columns.append(match)
            methods.append("crosswalk_hgnc")
    return columns, methods


def gene_mapping(source, axis, var, crosswalk, fixture):
    if source == "tahoe":
        columns, methods = tahoe_columns(axis, var, crosswalk)
    else:
        native = unique_native(var["hgnc_id"], "tenv HGNC")
        if missing := set(axis.hgnc_id) - native.keys():
            raise ValueError(f"tenv missing native HGNC identity: {sorted(missing)[:5]}")
        columns = [native[hgnc] for hgnc in axis.hgnc_id]
        methods = ["exact_hgnc"] * len(columns)
    if len(set(columns)) != len(axis):
        raise ValueError(f"{source}: multiple shared genes resolve to one native column")

    selected = var.iloc[columns]
    mapping = axis[["symbol", "hgnc_id", "ensembl_id"]].copy()
    mapping["native_column"] = columns
    mapping["native_symbol"] = selected.index.astype(str)
    mapping["native_ensembl_id"] = selected.ensembl_id.to_numpy()
    mapping["native_token_id"] = selected.token_id.astype(str).to_numpy() if source == "tahoe" else ""
    mapping["resolution"] = methods
    if not fixture:
        expected = ("ENSG00000289731", "61482") if source == "tahoe" else ("ENSG00000182230", "")
        actual = mapping.loc[mapping.hgnc_id.eq(FAM153B), ["native_ensembl_id", "native_token_id"]]
        if list(actual.itertuples(index=False, name=None)) != [expected]:
            raise ValueError(f"{source}: unexpected FAM153B native identity")
    return mapping


def source_contract(file, source, axis, crosswalk, fixture):
    obs = read_elem(file["obs"])
    var = read_elem(file["var"])
    index_name = "profile_id" if source == "tahoe" else "well"
    if (obs.index.name != index_name or obs.index.has_duplicates or obs.index.hasnans
            or (obs.index == "").any() or obs.plate.isna().any() or obs.plate.eq("").any()):
        raise ValueError(f"{source}: missing or duplicate profile identities or missing plate")
    if obs.is_control.dtype != np.dtype(bool):
        raise ValueError(f"{source}: invalid is_control flags")
    if "is_blank" in obs or source == "tenv":
        if obs.is_blank.dtype != np.dtype(bool):
            raise ValueError(f"{source}: invalid is_blank flags")
        blank = obs.is_blank.to_numpy()
    else:
        blank = np.zeros(len(obs), dtype=bool)
    if np.any(blank & obs.is_control.to_numpy()):
        raise ValueError(f"{source}: blank profiles marked as controls")
    rows = np.flatnonzero(~blank)
    profiles = pd.DataFrame({
        "source_row": rows,
        "profile_id": obs.index[rows].astype(str),
        "plate": obs.plate.iloc[rows].astype(str).to_numpy(),
        "is_control": obs.is_control.to_numpy()[rows],
    })
    controls = int(profiles.is_control.sum())
    population = {"raw": len(obs), "blank": int(blank.sum()),
                  "treated": len(profiles) - controls, "control": controls}
    if profiles.empty or (not fixture and population != POPULATIONS[source]):
        raise ValueError(f"{source}: unexpected source population {population}")

    mapping = gene_mapping(source, axis, var, crosswalk, fixture)
    summary = {
        "shape": [len(profiles), len(axis)],
        "population": population,
        "plates": profiles.plate.nunique(),
        "mapping_resolution_counts": {
            key: int(count) for key, count in mapping.resolution.value_counts().items()
        },
    }
    return profiles, mapping, (len(obs), len(var)), summary


def projected_batches(matrix, shape, rows, columns, source):
    """Read bounded native raw-count windows. The normalized native layer is unused."""
    sparse = isinstance(matrix, h5py.Group)
    if sparse:
        if matrix.attrs.get("encoding-type") != "csr_matrix" or tuple(matrix.attrs["shape"]) != shape:
            raise ValueError(f"{source}: expected native-axis CSR raw counts")
        ptr = matrix["indptr"][:]
        data, indices = matrix["data"], matrix["indices"]
        if (ptr.dtype.kind not in "iu" or ptr.shape != (shape[0] + 1,)
                or ptr[0] != 0 or np.any(ptr[1:] < ptr[:-1])
                or ptr[-1] != len(data) or len(indices) != len(data)
                or indices.dtype.kind not in "iu"):
            raise ValueError(f"{source}: invalid CSR pointers or indices")
        dtype = data.dtype
    else:
        if matrix.shape != shape:
            raise ValueError(f"{source}: raw count shape differs from metadata")
        dtype = matrix.dtype
    if dtype.kind not in "iu":
        raise ValueError(f"{source}: expected integer raw counts")

    for start in range(0, shape[0], ROW_BATCH):
        stop = min(start + ROW_BATCH, shape[0])
        low, high = rows.searchsorted([start, stop])
        selected = rows[low:high] - start
        if sparse:
            counts = data[ptr[start]:ptr[stop]]
            col_ids = indices[ptr[start]:ptr[stop]]
            if np.any(col_ids < 0) or np.any(col_ids >= shape[1]):
                raise ValueError(f"{source}: invalid CSR column indices")
            native = csr_matrix((counts, col_ids, ptr[start:stop + 1] - ptr[start]),
                                shape=(stop - start, shape[1]))
            if not native.has_canonical_format:
                raise ValueError(f"{source}: native CSR is unsorted or has duplicate columns")
        else:
            counts = native = matrix[start:stop]
        if np.any(counts < 0):
            raise ValueError(f"{source}: negative raw counts in source rows {start}:{stop}")
        if len(selected):
            projected = native[selected][:, columns]
            if sparse:
                projected = projected.toarray()
            totals = projected.sum(axis=1, dtype=np.float64)
            if np.any(totals <= 0):
                row = int(rows[low:high][np.flatnonzero(totals <= 0)[0]])
                raise ValueError(f"{source}: zero projected library at source row {row}")
            yield low, projected, totals


def csv_text(frame):
    return frame.to_csv(index=False, lineterminator="\n")


def write_source(path, directory, source, axis, crosswalk, fixture):
    with h5py.File(path, "r") as file:
        profiles, mapping, shape, summary = source_contract(file, source, axis, crosswalk, fixture)
        with (directory / f"{source}_expression.npy").open("wb") as expression:
            np.lib.format.write_array_header_2_0(expression, {
                "descr": np.dtype(np.float32).str, "fortran_order": False,
                "shape": tuple(summary["shape"]),
            })
            for _, counts, totals in projected_batches(
                file["X"], shape, profiles.source_row.to_numpy(), mapping.native_column.to_numpy(), source
            ):
                values = counts.astype(np.float64)
                values /= totals[:, None]
                values *= 10000
                np.log1p(values, out=values)
                values.astype(np.float32).tofile(expression)
        (directory / f"{source}_profiles.csv").write_text(csv_text(profiles))
        (directory / f"{source}_mapping.csv").write_text(csv_text(mapping))
    return summary


def build(native, release, fixture=False, sources=SOURCES):
    if release.exists():
        raise ValueError(f"release already exists: {release}")
    files = input_files(native, sources)
    before = fingerprint(files)
    digests, axis, crosswalk = read_inputs(files, fixture)
    release.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=release.parent, prefix=".shared-cm-") as stage:
        directory = Path(stage)
        shutil.copyfile(files["genes.csv"], directory / "genes.csv")
        sources = {source: write_source(files[f"{source}.h5ad"], directory, source,
                                        axis, crosswalk, fixture) for source in sources}
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "fixture": fixture,
            "normalization": NORMALIZATION,
            "inputs": digests,
            "raw_provenance": source_provenance(files, sources, fixture),
            "files": {name: file_sha256(directory / name) for name in ("genes.csv", *[f"{source}_{suffix}" for source in sources for suffix in ("expression.npy", "profiles.csv", "mapping.csv")])},
            "sources": sources,
        }
        if fingerprint(files) != before:
            raise ValueError("native inputs changed during the build")
        (directory / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        os.rename(directory, release)
    return {"manifest_sha256": file_sha256(release / "manifest.json"), **manifest}


def verify_source(path, directory, source, axis, crosswalk, fixture, summary):
    with h5py.File(path, "r") as file:
        profiles, mapping, shape, expected = source_contract(file, source, axis, crosswalk, fixture)
        if summary != expected:
            raise ValueError(f"{source}: release summary differs from native sources")
        for suffix, frame in (("profiles", profiles), ("mapping", mapping)):
            if (directory / f"{source}_{suffix}.csv").read_text() != csv_text(frame):
                raise ValueError(f"{source}: {suffix} differ from native sources and crosswalk")

        positions = np.flatnonzero(axis.hgnc_id.eq(FAM153B))
        spots = {}
        max_error = 0.0
        with (directory / f"{source}_expression.npy").open("rb") as expression:
            if np.lib.format.read_magic(expression) != (2, 0):
                raise ValueError(f"{source}: unsupported NPY format")
            actual_shape, fortran, dtype = np.lib.format.read_array_header_2_0(expression)
            if actual_shape != tuple(expected["shape"]) or fortran or dtype != np.dtype(np.float32):
                raise ValueError(f"{source}: unexpected expression shape, dtype or order")

            for offset, counts, totals in projected_batches(
                file["X"], shape, profiles.source_row.to_numpy(), mapping.native_column.to_numpy(), source
            ):
                observed = np.fromfile(expression, dtype=dtype, count=counts.size)
                if len(observed) != counts.size:
                    raise ValueError(f"{source}: truncated expression")
                observed = observed.reshape(counts.shape)
                # Replay from raw counts without using the writer's normalization.
                reference = np.log1p(counts * (10000 / totals[:, None])).astype(np.float32)
                if not np.isfinite(observed).all() or np.any(observed[counts == 0] != 0):
                    raise ValueError(f"{source}: non-finite expression or changed observed zero")
                max_error = max(max_error, float(np.max(np.abs(observed - reference))))
                if not np.allclose(observed, reference, atol=0, rtol=1e-7):
                    raise ValueError(f"{source}: normalization differs from raw counts at release row {offset}")
                if len(positions):
                    col = int(positions[0])
                    for label, mask in (("zero", counts[:, col] == 0), ("nonzero", counts[:, col] > 0)):
                        if label not in spots and mask.any():
                            row = int(np.flatnonzero(mask)[0])
                            profile = profiles.iloc[offset + row]
                            spots[label] = {
                                "native_column": int(mapping.native_column.iloc[col]),
                                "source_row": int(profile.source_row),
                                "profile_id": profile.profile_id,
                                "raw_count": int(counts[row, col]),
                                "projected_library_total": int(totals[row]),
                                "expression": float(observed[row, col]),
                            }
            if expression.read(1):
                raise ValueError(f"{source}: unexpected bytes after expression")
    return {"rows_replayed": len(profiles), "max_abs_error": max_error, "fam153b": spots}


def verify(native, release, fixture=False):
    manifest_path = release / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if (manifest["schema_version"] != SCHEMA_VERSION or manifest["fixture"] is not fixture
            or manifest["normalization"] != NORMALIZATION):
        raise ValueError("unexpected release schema, normalization or fixture mode")
    sources = tuple(manifest["sources"])
    if set(manifest["files"]) != {"genes.csv", *[f"{source}_{suffix}" for source in sources for suffix in ("expression.npy", "profiles.csv", "mapping.csv")]}:
        raise ValueError("unexpected release payload list")

    files = input_files(native, sources)
    before = fingerprint(files)
    digests, axis, crosswalk = read_inputs(files, fixture)
    if manifest["inputs"] != digests:
        raise ValueError("release input digests differ from native corpora and crosswalk")
    current_provenance = source_provenance(files, sources, fixture)
    for source, evidence in current_provenance.items():
        if (manifest["raw_provenance"][source]["release_manifest_sha256"] !=
                evidence["release_manifest_sha256"]):
            raise ValueError(f"{source}: raw release provenance changed")
    for name, digest in manifest["files"].items():
        if file_sha256(release / name) != digest:
            raise ValueError(f"release file digest mismatch: {name}")
    if manifest["files"]["genes.csv"] != digests["genes.csv"]:
        raise ValueError("release gene order differs from the pinned axis")
    checks = {
        source: verify_source(files[f"{source}.h5ad"], release, source, axis, crosswalk,
                              fixture, manifest["sources"][source])
        for source in sources
    }
    if fingerprint(files) != before:
        raise ValueError("native inputs changed during verification")
    return {"result": "PASS", "manifest_sha256": file_sha256(manifest_path), "sources": checks}


def tahoe_axis_mapping(axis, crosswalk, gene_metadata):
    """Return one evidenced source token per frozen gene, without expression access.

    ``native_column`` refers to immutable gene_metadata rows; gene_index is the
    frozen axis row. The returned original ENSG/token IDs remain source evidence.
    """
    if not isinstance(axis, pd.DataFrame):
        axis = read_csv(axis)
    if not isinstance(crosswalk, pd.DataFrame):
        crosswalk = read_csv(crosswalk)
    if not isinstance(gene_metadata, pd.DataFrame):
        gene_metadata = pd.read_parquet(gene_metadata)
    var = gene_metadata.rename(columns={'gene_symbol': 'native_symbol'}).copy()
    var.index = pd.Index(var.native_symbol.astype(str), name='symbol')
    mapping = gene_mapping('tahoe', axis, var, crosswalk, False)
    mapping.insert(0, 'gene_index', np.arange(len(mapping), dtype='int32'))
    return mapping
