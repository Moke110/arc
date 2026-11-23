"""Map fixed gene identities to source-native count columns."""
from collections import defaultdict

import numpy as np
import pandas as pd

from stage2_util.derivation.release import ReleaseError, sha256
from stage2_util.identity import AXIS_SHA256

AXIS_SIZE = 18164
FAM153B = 'HGNC:27323'

class RawCountsError(ReleaseError):
    """A raw-count release violates its contract."""


# --------------------------------------------------------------------------- #
# reading helpers
# --------------------------------------------------------------------------- #
def axis_frame(axis_path, *, fixture=False):
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


def tahoe_axis_mapping(axis, crosswalk, gene_metadata, fixture=False):
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
    mapping = gene_mapping('tahoe', axis, var, crosswalk, fixture)
    mapping.insert(0, 'gene_index', np.arange(len(mapping), dtype='int32'))
    return mapping


def read_csv(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


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
