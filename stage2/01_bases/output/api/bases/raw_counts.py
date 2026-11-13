"""Raw-count release on the frozen gene axis (ADR 0066).

This is the ``source/`` -> axis bridge layer of the Stage 2 rebuild, *not* a
model-facing product: integer counts only, no normalization, no delta and no
axis subsetting at this layer. Two payloads share one release area:

* ``tenv.h5ad`` — 4,510 native-gene wells x 19,613 native HGNC genes, with
  ``tenv_probe_mapping.csv`` (probe evidence) and ``tenv_axis_mapping.csv`` (the
  native-column bridge from each frozen ``gene_index`` to this payload);
* ``tahoe.h5ad`` — 67,018 projected wells x 18,164 frozen-axis genes, with
  ``token_mapping.csv`` (the original Ensembl/token evidence per frozen gene).

Provenance model (ADR 0066 *Removals*: the pre-0066 raw-count payloads' "content
moves into ``raw_counts/``"). The default publication mode **rehomes** the
already source-accepted payloads: it re-checks their bytes against the recorded
source-replay acceptance records, re-derives the Tahoe token mapping from the
immutable Tahoe ``source/`` metadata, and re-seals everything under one manifest
whose ``provenance.mode`` is ``source-accepted-replay`` with a pointer to the
acceptance record. The stronger ``source-rebuild`` mode rebuilds both payloads
from ``source/`` through the same builders and replays them from source; it stays
available and is documented in the step README/AGENTS.

The public readers in this module are the only supported access path to the raw
payloads. They expose bounded row windows (dense ``int32`` for tenv, CSR windows
for Tahoe) and never materialise an observation x 18,164 dense matrix.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy.sparse import csr_matrix

from . import native, raw, tahoe
from .identity import AXIS_SHA256
from .release import (ReleaseError, check_manifest, output_digests, publish,
                      read_json, replay, sha256, write_json)

SCHEMA = "stage2-bases-raw-counts-v1"
AXIS_SIZE = 18164
FAM153B_INDEX = 5091
FAM153B_HGNC = "HGNC:27323"
FAM153B_TENV_ENSEMBL = "ENSG00000182230"
FAM153B_TAHOE_ENSEMBL = "ENSG00000289731"
FAM153B_TAHOE_TOKEN = "61482"

TENV_ITEM = "pert_data/GSE262419"
TAHOE_ITEM = "pert_data/Tahoe-100M"
AXIS_PATH = "stage2/01_bases/output/gene-identity/tables/genes.csv"
CROSSWALK_PATH = "stage2/01_bases/output/gene-identity/tables/gene_crosswalk.csv"
TENV_ACCEPTANCE = "stage2/01_bases/log/records/native_tenv_verify.json"
TAHOE_ACCEPTANCE = "stage2/01_bases/log/records/tahoe_counts_source_acceptance.json"

# Expected real-data populations. They are asserted, never inferred.
TENV_COUNTS = {"profiles": 4510, "genes": 19613, "treated": 4206,
               "vehicle": 272, "blank": 32}
TAHOE_COUNTS = {"profiles": 67018, "genes": 18164, "treated": 65618,
                "vehicle": 1400, "cells": 95624334, "plates": 14, "shards": 3388}

RECIPE = {
    "tenv": ("GSE262419 TempO-Seq plate CSVs summed by the resolved native HGNC "
             "probe identity; int32; no normalization, no delta, no axis "
             "subsetting. The probe evidence is retained and the frozen-axis "
             "bridge is sealed as tenv_axis_mapping.csv."),
    "tahoe": ("Tahoe-100M cells summed per stable (plate, well, cell_line_id) "
              "profile as int32 UMI counts, projected to the frozen 18,164 axis "
              "before normalization; marker tokens dropped; the original "
              "Ensembl/token identity is retained in token_mapping.csv. No "
              "normalization, no delta."),
    "provenance": ("Default mode rehomes the source-accepted payloads and re-seals "
                   "them after re-checking their sha256 against the recorded "
                   "source-replay acceptance record and re-deriving the Tahoe "
                   "token mapping from the immutable source metadata; "
                   "--full-rebuild rebuilds both payloads from source/."),
}

TENV_PAYLOAD = ("tenv.h5ad", "tenv_probe_mapping.csv")
TAHOE_PAYLOAD = ("tahoe.h5ad", "token_mapping.csv")
PRODUCT_FILES = ("tenv.h5ad", "tenv_probe_mapping.csv", "tenv_axis_mapping.csv",
                 "tahoe.h5ad", "token_mapping.csv")


class RawCountsError(ReleaseError):
    """A raw-count release violates its contract."""


# --------------------------------------------------------------------------- #
# reading helpers
# --------------------------------------------------------------------------- #
def _axis_frame(axis_path, *, fixture=False):
    axis = pd.read_csv(axis_path, dtype={"symbol": str, "hgnc_id": str,
                                         "ensembl_id": str}, keep_default_na=False)
    required = {"symbol", "hgnc_id", "ensembl_id"}
    if not required <= set(axis.columns):
        raise RawCountsError(f"axis requires columns {sorted(required)}")
    if not fixture:
        if len(axis) != AXIS_SIZE:
            raise RawCountsError(f"expected {AXIS_SIZE} frozen genes, got {len(axis)}")
        if sha256(axis_path) != AXIS_SHA256:
            raise RawCountsError("frozen gene axis digest mismatch")
    for field in ("symbol", "hgnc_id", "ensembl_id"):
        values = axis[field]
        if values.eq("").any() or values.duplicated().any():
            raise RawCountsError(f"axis {field} missing, empty or duplicated")
    axis = axis.reset_index(drop=True)
    axis.insert(0, "gene_index", np.arange(len(axis), dtype=np.int64))
    return axis


def tenv_axis_mapping(axis, var):
    """One native column per frozen gene, with the exact-HGNC resolution evidence."""
    axis = pd.DataFrame(axis)
    var = pd.DataFrame(var)
    if "hgnc_id" not in var or "ensembl_id" not in var:
        raise RawCountsError("native tenv var lacks hgnc_id/ensembl_id")
    native_index = {}
    for column, (hgnc, ensembl) in enumerate(
            zip(var.hgnc_id.astype(str), var.ensembl_id.astype(str), strict=True)):
        if hgnc in ("", "nan"):
            continue
        if hgnc in native_index:
            raise RawCountsError(f"ambiguous native tenv HGNC identity: {hgnc}")
        native_index[hgnc] = (column, ensembl)
    missing = [hgnc for hgnc in axis.hgnc_id if hgnc not in native_index]
    if missing:
        raise RawCountsError(f"tenv missing frozen HGNC identities: {missing[:5]}")
    rows = []
    for row in axis.itertuples(index=False):
        column, ensembl = native_index[row.hgnc_id]
        rows.append({"gene_index": int(row.gene_index), "symbol": row.symbol,
                     "hgnc_id": row.hgnc_id, "ensembl_id": ensembl,
                     "native_column": column, "native_symbol": row.symbol,
                     "resolution": "exact_hgnc"})
    mapping = pd.DataFrame(rows)
    columns = mapping.native_column.to_numpy()
    if len(set(columns.tolist())) != len(columns):
        raise RawCountsError("multiple frozen genes resolve to one native column")
    if mapping.gene_index.tolist() != list(range(len(axis))):
        raise RawCountsError("bridge must cover every frozen gene_index in order")
    return mapping


def _check_fam153b(axis, bridge, mapping):
    if len(axis) <= FAM153B_INDEX:
        return
    if str(axis.iloc[FAM153B_INDEX].hgnc_id) != FAM153B_HGNC:
        raise RawCountsError("frozen FAM153B gene_index moved")
    if bridge is not None:
        if str(bridge.iloc[FAM153B_INDEX].ensembl_id) != FAM153B_TENV_ENSEMBL:
            raise RawCountsError("tenv FAM153B native identity differs")
    if mapping is not None:
        row = mapping.iloc[FAM153B_INDEX]
        if (str(row.native_ensembl_id), str(row.native_token_id)) != (
                FAM153B_TAHOE_ENSEMBL, FAM153B_TAHOE_TOKEN):
            raise RawCountsError("Tahoe FAM153B native identity differs")


def _h5_header(path):
    with h5py.File(path, "r") as handle:
        if "obs" not in handle or "var" not in handle or "X" not in handle:
            raise RawCountsError(f"not a counts payload: {path}")
        return read_elem(handle["obs"]), read_elem(handle["var"])


def _x_kind(matrix):
    if isinstance(matrix, h5py.Group):
        if matrix.attrs.get("encoding-type") != "csr_matrix":
            raise RawCountsError("sparse raw counts must be CSR")
        return "csr"
    return "dense"


def _x_dtype(matrix):
    return matrix["data"].dtype if isinstance(matrix, h5py.Group) else matrix.dtype


def _counts_rows(matrix, start, stop):
    """Dense int32 window of one raw payload's X matrix."""
    if isinstance(matrix, h5py.Group):
        pointer = matrix["indptr"]
        lo, hi = int(pointer[start]), int(pointer[stop])
        data = matrix["data"][lo:hi].astype(np.int32)
        indices = matrix["indices"][lo:hi].astype(np.int32)
        local = pointer[start:stop + 1].astype(np.int64) - lo
        block = csr_matrix((data, indices, local),
                           shape=(stop - start, int(matrix.attrs["shape"][1])))
        return block.toarray().astype(np.int32, copy=False)
    return np.asarray(matrix[start:stop], dtype=np.int32)


# --------------------------------------------------------------------------- #
# publication
# --------------------------------------------------------------------------- #
def _record(path, repo_root):
    """An acceptance record: repo-relative pointer plus its digest and content."""
    path = Path(path)
    return {"path": _relative_path(path, repo_root), "absolute": str(path),
            "sha256": sha256(path), "value": read_json(path)}


def _copy_payload(candidate, payload_dir, names, expected_digests, label):
    candidate, payload_dir = Path(candidate), Path(payload_dir)
    copied = {}
    for name in names:
        source = payload_dir / name
        if not source.is_file():
            raise RawCountsError(f"{label}: missing accepted payload: {source}")
        actual = sha256(source)
        expected = expected_digests.get(name)
        if expected is not None and actual != expected:
            raise RawCountsError(
                f"{label}: payload {name} digest {actual} does not match the "
                f"recorded source-replay acceptance digest {expected}")
        destination = candidate / name
        if source.resolve() != destination.resolve():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, destination)
        copied[name] = {"sha256": actual, "bytes": destination.stat().st_size}
    return copied


def _seal_tenv(candidate, payload_dir, manifest_name, axis, *, fixture,
               payload_digests=None, accepted_manifest_digest=None):
    payload_dir = Path(payload_dir)
    files = _copy_payload(candidate, payload_dir, TENV_PAYLOAD,
                          payload_digests or {}, "tenv")
    build_manifest = read_json(payload_dir / manifest_name)
    if not fixture:
        if build_manifest.get("schema") != "stage2-bases-native-tenv-v1":
            raise RawCountsError("unexpected native tenv build schema")
        if build_manifest.get("shape") != [TENV_COUNTS["profiles"], TENV_COUNTS["genes"]]:
            raise RawCountsError("unexpected native tenv payload shape")
        if build_manifest.get("populations") != {
                "blank": TENV_COUNTS["blank"], "control": TENV_COUNTS["vehicle"],
                "treated": TENV_COUNTS["treated"]}:
            raise RawCountsError("unexpected native tenv population")
        if build_manifest.get("files") != {
                name: files[name]["sha256"] for name in TENV_PAYLOAD}:
            raise RawCountsError("native tenv manifest payload digests differ")
        if (accepted_manifest_digest is not None
                and sha256(payload_dir / manifest_name) != accepted_manifest_digest):
            raise RawCountsError("native tenv manifest is not the accepted one")
    obs, var = _h5_header(Path(candidate) / "tenv.h5ad")
    if not fixture and len(obs) != TENV_COUNTS["profiles"]:
        raise RawCountsError("native tenv profile count differs")
    bridge = tenv_axis_mapping(axis, var)
    bridge_path = Path(candidate) / "tenv_axis_mapping.csv"
    bridge.to_csv(bridge_path, index=False, lineterminator="\n")
    files["tenv_axis_mapping.csv"] = {"sha256": sha256(bridge_path),
                                      "bytes": bridge_path.stat().st_size}
    return obs, var, files, build_manifest, bridge


def _seal_tahoe(candidate, payload_dir, manifest_name, axis, *, fixture,
                payload_digests=None, gene_metadata=None, crosswalk=None):
    payload_dir = Path(payload_dir)
    files = _copy_payload(candidate, payload_dir, TAHOE_PAYLOAD,
                          payload_digests or {}, "tahoe")
    build_manifest = read_json(payload_dir / manifest_name)
    if not fixture:
        if build_manifest.get("schema") != "bases-tahoe-shared-counts-v1":
            raise RawCountsError("unexpected Tahoe build schema")
        counts = {"profiles": build_manifest["profiles"],
                  "cells": build_manifest["cells"], "controls": build_manifest["controls"],
                  "treated": build_manifest["treated"], "genes": build_manifest["genes"]}
        if counts != {"profiles": TAHOE_COUNTS["profiles"],
                      "cells": TAHOE_COUNTS["cells"], "controls": TAHOE_COUNTS["vehicle"],
                      "treated": TAHOE_COUNTS["treated"], "genes": TAHOE_COUNTS["genes"]}:
            raise RawCountsError(f"unexpected Tahoe build counts: {counts}")
        if build_manifest.get("h5ad_sha256") != files["tahoe.h5ad"]["sha256"]:
            raise RawCountsError("Tahoe payload digest differs from build manifest")
        if build_manifest.get("mapping_sha256") != files["token_mapping.csv"]["sha256"]:
            raise RawCountsError("Tahoe token mapping digest differs from build manifest")
    obs, var = _h5_header(Path(candidate) / "tahoe.h5ad")
    mapping = pd.read_csv(Path(candidate) / "token_mapping.csv",
                          dtype={"native_token_id": str, "native_ensembl_id": str,
                                 "hgnc_id": str, "ensembl_id": str},
                          keep_default_na=False)
    if not fixture:
        if len(obs) != TAHOE_COUNTS["profiles"] or len(var) != TAHOE_COUNTS["genes"]:
            raise RawCountsError("unexpected Tahoe payload shape")
        if int(obs.is_control.sum()) != TAHOE_COUNTS["vehicle"]:
            raise RawCountsError("unexpected Tahoe vehicle control count")
        if int(obs.n_cells.sum()) != TAHOE_COUNTS["cells"]:
            raise RawCountsError("Tahoe cell count is not conserved")
        if obs.plate.nunique() != TAHOE_COUNTS["plates"]:
            raise RawCountsError("unexpected Tahoe plate count")
        if mapping.gene_index.tolist() != list(range(AXIS_SIZE)):
            raise RawCountsError("Tahoe token mapping must cover the frozen axis in order")
        if gene_metadata is None or crosswalk is None:
            raise RawCountsError("token-mapping replay requires source gene metadata")
        expected = raw.tahoe_axis_mapping(axis, crosswalk, gene_metadata,
                                          fixture=fixture)
        columns = ["gene_index", "symbol", "hgnc_id", "ensembl_id",
                   "native_column", "native_symbol", "native_ensembl_id",
                   "native_token_id", "resolution"]
        actual = mapping[columns].astype(str).reset_index(drop=True)
        reference = expected[columns].astype(str).reset_index(drop=True)
        if not actual.equals(reference):
            raise RawCountsError("Tahoe token mapping differs from source gene metadata")
    return obs, var, files, build_manifest, mapping


def _provenance(mode, *, acceptance=None, payloads=None, counts=None,
                source_manifest=None):
    provenance = {"mode": mode}
    if acceptance is not None:
        value = acceptance["value"]
        provenance.update({"acceptance_record": acceptance["path"],
                           "acceptance_record_sha256": acceptance["sha256"],
                           "recorded_result": value.get("result")})
        for key in ("manifest_sha256", "source_rows_replayed", "scope",
                    "immutable_source_manifest_sha256", "raw_token_mapping_sha256",
                    "source_shards"):
            if key in value:
                provenance[f"recorded_{key}"] = value[key]
        replay_record = value.get("independent_replay")
        if replay_record:
            provenance["recorded_source_replay"] = {
                key: replay_record.get(key) for key in
                ("schema", "input_sha256", "h5ad_sha256", "verified_profiles",
                 "verified_cells")}
    if source_manifest is not None:
        provenance.update(source_manifest)
    if payloads is not None:
        provenance["payload_sha256"] = {name: record["sha256"]
                                        for name, record in payloads.items()}
    if counts is not None:
        provenance["counts"] = counts
    return provenance


def _relative_path(path, repo_root):
    """Manifest keys are repo-relative; fixture keys are bare file names.

    ``absolute()`` (not ``resolve()``) keeps a declared read-only input reachable
    through a symlink without leaking the symlink target into the manifest.
    """
    path = Path(path)
    if repo_root is None:
        return path.name
    root = Path(repo_root).absolute()
    try:
        return str(path.absolute().relative_to(root))
    except ValueError as error:
        raise RawCountsError(f"release input outside the repository: {path}") from error


def _input_entry(path, inputs, repo_root):
    path = Path(path)
    key = _relative_path(path, repo_root)
    if key in inputs:
        raise RawCountsError(f"duplicate release input path: {key}")
    inputs[key] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    return key


def _outputs(candidate, names):
    outputs = {}
    for name in names:
        path = Path(candidate) / name
        record = {"sha256": sha256(path), "bytes": path.stat().st_size}
        if name.endswith(".csv"):
            frame = pd.read_csv(path, keep_default_na=False)
            record.update({"rows": len(frame), "columns": frame.columns.tolist()})
        outputs[name] = record
    return outputs


def _prepare_acceptance(mode, fixture, data_root, tenv_payload, tahoe_payload,
                        tenv_acceptance, tahoe_acceptance, repo_root):
    """Validate the recorded source replay and return (records, payload digests)."""
    if mode != "source-accepted-replay" or fixture:
        return None, None, {}
    tenv_record = _record(tenv_acceptance, repo_root)
    tahoe_record = _record(tahoe_acceptance, repo_root)
    if tenv_record["value"].get("result") != "PASS":
        raise RawCountsError("the recorded tenv source replay is not PASS")
    tenv_build = read_json(Path(tenv_payload) / "tenv_manifest.json")
    if tenv_record["value"].get("manifest_sha256") != sha256(
            Path(tenv_payload) / "tenv_manifest.json"):
        raise RawCountsError("tenv payload is not the source-accepted payload")
    if tenv_record["value"].get("source_rows_replayed") != TENV_COUNTS["profiles"]:
        raise RawCountsError("tenv recorded source replay is incomplete")
    tahoe_value = tahoe_record["value"]
    if tahoe_value.get("result") != "PASS":
        raise RawCountsError("the recorded Tahoe source replay is not PASS")
    if tahoe_value.get("raw_release_manifest", {}).get("sha256") != sha256(
            Path(tahoe_payload) / "manifest.json"):
        raise RawCountsError("Tahoe payload is not the source-accepted payload")
    recorded = tahoe_value.get("raw_release", {})
    expected_tahoe = {"tahoe.h5ad": recorded.get("h5ad_sha256"),
                      "token_mapping.csv": tahoe_value.get("raw_token_mapping_sha256")}
    recorded_tenv = {name: tenv_build["files"][name] for name in TENV_PAYLOAD}
    if tahoe_value.get("source_shards") != TAHOE_COUNTS["shards"]:
        raise RawCountsError("Tahoe recorded source replay did not cover every shard")
    if tahoe_value.get("independent_replay", {}).get("verified_profiles") != (
            TAHOE_COUNTS["profiles"]):
        raise RawCountsError("Tahoe recorded independent replay profile count differs")
    if tahoe_value.get("immutable_source_manifest_sha256") != sha256(
            Path(data_root) / "Tahoe-100M/manifest.sha256"):
        raise RawCountsError("Tahoe immutable source manifest digest differs")
    return (tenv_record, tahoe_record), (recorded_tenv, expected_tahoe), {}


def build(candidate, *, mode="source-accepted-replay", data_root, axis_path,
          crosswalk_path, tenv_payload, tahoe_payload, tenv_acceptance=None,
          tahoe_acceptance=None, gene_metadata=None, fixture=False,
          checkpoints_root=None, repo_root=None, tahoe_expected_shards=None,
          tahoe_expected_genes=None):
    """Write a complete raw-count candidate outside ``output/``."""
    if mode not in ("source-accepted-replay", "source-rebuild"):
        raise RawCountsError(f"unsupported publication mode: {mode}")
    if fixture and mode == "source-rebuild" and tahoe_expected_shards is None:
        raise RawCountsError("a fixture source-rebuild needs explicit shard counts")
    candidate = Path(candidate)
    candidate.mkdir(parents=True, exist_ok=True)
    data_root = Path(data_root)
    axis = _axis_frame(axis_path, fixture=fixture)
    crosswalk = pd.read_csv(crosswalk_path, dtype=str, keep_default_na=False)
    if not fixture and repo_root is None:
        raise RawCountsError("a real build requires repo_root for relative manifest paths")
    inputs = {}
    for path in (axis_path, crosswalk_path):
        _input_entry(path, inputs, repo_root)
    records, recorded, _ = _prepare_acceptance(
        mode, fixture, data_root, tenv_payload, tahoe_payload,
        tenv_acceptance, tahoe_acceptance, repo_root)
    if mode == "source-accepted-replay":
        for directory in (tenv_payload, tahoe_payload):
            for name in TENV_PAYLOAD + TAHOE_PAYLOAD + (
                    "tenv_manifest.json", "manifest.json"):
                path = Path(directory) / name
                if path.is_file() and str(path) not in inputs:
                    _input_entry(path, inputs, repo_root)
        tenv_digests, tahoe_digests = recorded if recorded else ({}, {})
        obs_t, var_t, files_t, build_t, bridge = _seal_tenv(
            candidate, tenv_payload, "tenv_manifest.json", axis, fixture=fixture,
            payload_digests=tenv_digests,
            accepted_manifest_digest=(records[0]["value"]["manifest_sha256"]
                                      if records else None))
        obs_h, var_h, files_h, build_h, mapping = _seal_tahoe(
            candidate, tahoe_payload, "manifest.json", axis, fixture=fixture,
            payload_digests=tahoe_digests,
            gene_metadata=gene_metadata, crosswalk=crosswalk)
    else:
        # Stronger mode: rebuild both payloads from the immutable source trees.
        if gene_metadata is None:
            raise RawCountsError("a full rebuild needs the source gene metadata")
        native_build = native.build_tenv(
            data_root / "GSE262419/source",
            data_root / "Annotation/source/hgnc_complete_set.txt",
            data_root / "Annotation/local/cm_compound_annotation.csv", candidate,
            manifest="tenv_build.json")
        mapping_frame = raw.tahoe_axis_mapping(axis, crosswalk, gene_metadata,
                                               fixture=fixture)
        options = {}
        if tahoe_expected_shards is not None:
            options["expected_shards"] = tahoe_expected_shards
        if tahoe_expected_genes is not None:
            options["expected_genes"] = tahoe_expected_genes
        tahoe_build = tahoe.build_tahoe(
            data_root / "Tahoe-100M/source", mapping_frame,
            data_root / "Tahoe-100M/manifest.sha256", candidate,
            manifest="tahoe_build.json",
            checkpoints=(None if checkpoints_root is None else Path(checkpoints_root)),
            **options)
        obs_t, var_t, files_t, build_t, bridge = _seal_tenv(
            candidate, candidate, "tenv_build.json", axis, fixture=fixture,
            payload_digests={name: native_build["files"][name] for name in TENV_PAYLOAD})
        obs_h, var_h, files_h, build_h, mapping = _seal_tahoe(
            candidate, candidate, "tahoe_build.json", axis, fixture=fixture,
            payload_digests={"tahoe.h5ad": tahoe_build["h5ad_sha256"],
                             "token_mapping.csv": tahoe_build["mapping_sha256"]},
            gene_metadata=gene_metadata, crosswalk=crosswalk)
    if not fixture:
        _check_fam153b(axis, bridge, mapping)
    tenv_counts = {"profiles": int(len(obs_t)), "genes": int(len(var_t)),
                   "treated": int((~obs_t.is_control & ~obs_t.is_blank).sum()),
                   "vehicle": int(obs_t.is_control.sum()),
                   "blank": int(obs_t.is_blank.sum())}
    tahoe_counts = {"profiles": int(len(obs_h)), "genes": int(len(var_h)),
                    "treated": int((~obs_h.is_control).sum()),
                    "vehicle": int(obs_h.is_control.sum()),
                    "cells": int(obs_h.n_cells.sum())}
    if not fixture:
        if tenv_counts != TENV_COUNTS:
            raise RawCountsError(f"unexpected tenv population: {tenv_counts}")
        if tahoe_counts != {key: TAHOE_COUNTS[key] for key in tahoe_counts}:
            raise RawCountsError(f"unexpected Tahoe population: {tahoe_counts}")
        for path in (data_root / "GSE262419/manifest.sha256",
                     data_root / "Tahoe-100M/manifest.sha256", gene_metadata):
            _input_entry(path, inputs, repo_root)
        for record in records or ():
            _input_entry(record["absolute"], inputs, repo_root)
    names = list(PRODUCT_FILES)
    for extra in ("tenv_build.json", "tahoe_build.json"):
        if (candidate / extra).is_file():
            names.append(extra)
    outputs = _outputs(candidate, names)
    provenance_mode = "fixture" if fixture else mode
    item_manifests = {}
    if not fixture:
        item_manifests = {
            "tenv": {"path": f"{TENV_ITEM}/manifest.sha256",
                     "sha256": sha256(data_root / "GSE262419/manifest.sha256")},
            "tahoe": {"path": f"{TAHOE_ITEM}/manifest.sha256",
                      "sha256": sha256(data_root / "Tahoe-100M/manifest.sha256")}}
    manifest = {
        "schema": SCHEMA,
        "recipe": RECIPE,
        "axis": {"size": int(len(axis)),
                 "path": _relative_path(axis_path, repo_root),
                 "sha256": inputs[_relative_path(axis_path, repo_root)]["sha256"]},
        "inputs": inputs,
        "outputs": outputs,
        "sources": {
            "tenv": {
                "counts": tenv_counts,
                "bridge": {"axis_genes_mapped": int(len(bridge)),
                           "native_genes": int(len(var_t)),
                           "native_genes_off_axis": int(len(var_t) - len(bridge)),
                           "resolution": {"exact_hgnc": int(len(bridge))}},
                "item_manifest": item_manifests.get("tenv"),
                "provenance": _provenance(
                    provenance_mode, acceptance=(records[0] if records else None),
                    payloads=files_t, counts=tenv_counts),
            },
            "tahoe": {
                "counts": tahoe_counts,
                "token_mapping": {
                    "axis_genes_mapped": int(len(mapping)),
                    "native_tokens": int(len(mapping)),
                    "resolution_counts": ({str(k): int(v) for k, v in
                                           mapping.resolution.value_counts().items()}
                                          if "resolution" in mapping
                                          else {})},
                "item_manifest": item_manifests.get("tahoe"),
                "provenance": _provenance(
                    provenance_mode, acceptance=(records[1] if records else None),
                    payloads=files_h, counts=tahoe_counts),
            },
        },
        "counts": {"axis_size": int(len(axis)),
                   "tenv_profiles": tenv_counts["profiles"],
                   "tahoe_profiles": tahoe_counts["profiles"]},
    }
    write_json(candidate / "manifest.json", manifest)
    return manifest


def verify(candidate):
    """Re-check a candidate's manifest, payload digests, schema and index space."""
    candidate = Path(candidate)
    manifest_path = candidate / "manifest.json"
    if not manifest_path.is_file():
        raise RawCountsError("candidate has no manifest.json")
    manifest = read_json(manifest_path)
    if manifest.get("schema") != SCHEMA:
        raise RawCountsError("unexpected raw-count release schema")
    check_manifest(candidate, manifest)
    axis_size = int(manifest["counts"]["axis_size"])
    tenv_obs, tenv_var = _h5_header(candidate / "tenv.h5ad")
    tahoe_obs, tahoe_var = _h5_header(candidate / "tahoe.h5ad")
    with h5py.File(candidate / "tenv.h5ad") as handle:
        if _x_kind(handle["X"]) != "dense":
            raise RawCountsError("tenv raw counts must be dense")
        if _x_dtype(handle["X"]).kind not in "iu":
            raise RawCountsError("tenv raw counts must be integer")
        if tuple(handle["X"].shape) != (len(tenv_obs), len(tenv_var)):
            raise RawCountsError("tenv count shape differs from metadata")
        if handle["X"][:min(8, len(tenv_obs))].min() < 0:
            raise RawCountsError("negative tenv raw counts")
        if "layers" in handle and len(handle["layers"]):
            raise RawCountsError("raw counts must not carry normalized layers")
    with h5py.File(candidate / "tahoe.h5ad") as handle:
        if _x_kind(handle["X"]) != "csr":
            raise RawCountsError("Tahoe raw counts must be CSR")
        if _x_dtype(handle["X"]).kind not in "iu":
            raise RawCountsError("Tahoe raw counts must be integer")
        if tuple(int(v) for v in handle["X"].attrs["shape"]) != (
                len(tahoe_obs), len(tahoe_var)):
            raise RawCountsError("Tahoe count shape differs from metadata")
    bridge = pd.read_csv(candidate / "tenv_axis_mapping.csv",
                         dtype={"hgnc_id": str, "ensembl_id": str},
                         keep_default_na=False)
    if len(bridge) != axis_size or bridge.gene_index.tolist() != list(range(axis_size)):
        raise RawCountsError("tenv bridge must cover the frozen axis in order")
    if not bridge.resolution.eq("exact_hgnc").all():
        raise RawCountsError("tenv bridge resolution is not exact_hgnc")
    if bridge.native_column.duplicated().any():
        raise RawCountsError("tenv bridge maps two frozen genes to one native column")
    mapping = pd.read_csv(candidate / "token_mapping.csv",
                          dtype={"native_token_id": str, "native_ensembl_id": str},
                          keep_default_na=False)
    if len(mapping) != axis_size or mapping.gene_index.tolist() != list(range(axis_size)):
        raise RawCountsError("token mapping must cover the frozen axis in order")
    if mapping.native_token_id.duplicated().any():
        raise RawCountsError("token mapping reuses one source token")
    for source, obs in (("tenv", tenv_obs), ("tahoe", tahoe_obs)):
        if obs.index.has_duplicates or obs.index.hasnans:
            raise RawCountsError(f"{source}: duplicate or missing profile identities")
    if tenv_obs.is_control.dtype != np.dtype(bool):
        raise RawCountsError("tenv vehicle flags are not boolean")
    if int(tenv_obs.is_control.sum()) != int(
            manifest["sources"]["tenv"]["counts"]["vehicle"]):
        raise RawCountsError("tenv vehicle population differs from the manifest")
    if int((~tahoe_obs.is_control).sum()) != int(
            manifest["sources"]["tahoe"]["counts"]["treated"]):
        raise RawCountsError("Tahoe treated population differs from the manifest")
    for source in ("tenv", "tahoe"):
        provenance = manifest["sources"][source]["provenance"]
        if provenance.get("mode") not in ("fixture", "source-accepted-replay",
                                          "source-rebuild"):
            raise RawCountsError(f"{source}: unknown provenance mode")
        declared = provenance.get("payload_sha256", {})
        for name, record in declared.items():
            if sha256(candidate / name) != record:
                raise RawCountsError(f"{source}: payload digest differs: {name}")
    return {"result": "PASS", "manifest_sha256": sha256(manifest_path),
            "outputs": sorted(output_digests(manifest))}


def publish_release(output, *, root=None, **kwargs):
    """Build, verify and replace the accepted raw-count release (candidate-then-replace)."""
    return publish(output, lambda stage: build(stage, **kwargs),
                   lambda stage: verify(stage), root=root)


def replay_release(output, *, root=None, **kwargs):
    """Determinism gate: rebuild from the pinned inputs and compare every digest."""
    return replay(output, lambda stage: build(stage, **kwargs),
                  lambda stage: verify(stage), root=root)


class RawCounts:
    """Bounded read-only access to one raw-count payload in a sealed release.

    ``open_tenv``/``open_tahoe`` are the public entry points. Iterate with
    :meth:`batches`; a dense observation x 18,164 matrix is never materialised.
    """

    def __init__(self, root, source):
        self.root = Path(root)
        self.source = source
        self.manifest = read_json(self.root / "manifest.json")
        if self.manifest.get("schema") != SCHEMA:
            raise RawCountsError("unsupported raw-count release schema")
        self.payload = self.root / ("tenv.h5ad" if source == "tenv" else "tahoe.h5ad")
        self._handle = h5py.File(self.payload, "r")
        self.obs = read_elem(self._handle["obs"])
        self.var = read_elem(self._handle["var"])
        self._matrix = self._handle["X"]
        self.shape = (len(self.obs), len(self.var))

    def close(self):
        self._handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def bridge(self):
        """The source->frozen-axis mapping evidence for this payload."""
        name = "tenv_axis_mapping.csv" if self.source == "tenv" else "token_mapping.csv"
        return pd.read_csv(self.root / name, keep_default_na=False, dtype={
            "native_token_id": str, "native_ensembl_id": str,
            "hgnc_id": str, "ensembl_id": str, "symbol": str})

    def counts(self, start=0, stop=None):
        stop = self.shape[0] if stop is None else int(stop)
        return _counts_rows(self._matrix, int(start), stop)

    def batches(self, batch_size=256):
        if batch_size <= 0:
            raise RawCountsError("batch_size must be positive")
        for start in range(0, self.shape[0], batch_size):
            stop = min(start + batch_size, self.shape[0])
            yield start, self.counts(start, stop)


def open_tenv(root):
    return RawCounts(root, "tenv")


def open_tahoe(root):
    return RawCounts(root, "tahoe")
