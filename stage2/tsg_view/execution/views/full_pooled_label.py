"""The `full-pooled-label` view: the `full` inputs, only the label pooled.

Every row, expression value and control value is the `full` view's; nothing is copied. The
single difference is the training label: each treated profile's Δ becomes its condition's
`merge` Δ, so a condition's repeats share one label and the label noise can move without
moving the profile count (ADR 0082, plan C7.3/C7.4).

Provenance is read from the canonical products, never re-derived here:

* ``full`` product manifest -> the accepted 01 payload/manifest digests the ``X``/``control``
  alias points at; row identity is asserted by the ``profile_ids`` digest.
* ``merge`` product manifest + delta layer -> the pooled Δ transport. The merge payload is
  re-hashed, its rows must be exactly the canonical condition rows, and its member table
  must reconcile to the conditions the label is projected onto.

Transport is a row-aligned float32 ``<source>.delta.npy``: a treated row carries its
condition's merge Δ, a vehicle row keeps the full view's exact ``X - control``. The payload
is therefore deliberately not ``X - control`` for treated rows; the manifest declares
``delta_kind: pooled`` and names the formula the freeze uses.
"""
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import (BATCH, CONTROL_RULE, RECIPE, SOURCES, HarnessError, condition_table,
               merge_profile_id, relative, sha256, write_json)

# Released-population expectations (reconnaissance notes-bases-views §0/§5); asserted only
# for a non-fixture build so a corpus/release drift fails closed instead of publishing a
# differently shaped label transport.
EXPECTED_ROWS = {"tenv": 4478, "tahoe_bulk": 66718}
EXPECTED_TREATED = {"tenv": 4206, "tahoe_bulk": 65318}
EXPECTED_VEHICLE = {"tenv": 272, "tahoe_bulk": 1400}
EXPECTED_CONDITIONS = {"tenv": 1377, "tahoe_bulk": 52777}
EXPECTED_MEMBER_HISTOGRAM = {
    "tenv": {"1": 0, "2": 0, "3": 1320, "4": 48, "5": 0, "6": 9},
    "tahoe_bulk": {"1": 43036, "2": 7941, "3": 1400, "4": 350, "5": 0, "6": 0, "7": 0,
                   "8": 0, "9": 0, "10": 0, "11": 0, "12": 0, "13": 0, "14": 0, "15": 0,
                   "16": 50},
}
VEHICLE_RULE = ("vehicle row: the full view's exact X - control (a vehicle row has no "
                "condition)")
TRANSFER = ("treated row delta = the delta layer of the merge row whose "
            "(cell_line, drug_parent_inchikey, dose_um) equals the row's condition")
TRANSPORT = ("row-aligned float32 .npy; treated rows carry their condition's merge delta; "
             "X/control aliased to the accepted full payload digests; delta is not "
             "X - control for treated rows by design")


# --------------------------------------------------------------------------- #
# provenance products
# --------------------------------------------------------------------------- #
def _product_manifest(context, view):
    """Locate a previously built product manifest (same views root, else the canonical repo)."""
    name = f"{context.source}.json"
    candidates = [Path(context.output_dir) / view / name]
    if not context.fixture:
        candidates.append(Path(context.root) / "stage2/tsg_view/output/views" / view / name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise HarnessError(
        f"the {view} product is unavailable: {candidates[0]}; build `{view}` before "
        "`full-pooled-label`")


def _read_manifest(context, view):
    path = _product_manifest(context, view)
    manifest = json.loads(path.read_text())
    if (manifest.get("schema") != "tsg-view-product-v1" or manifest.get("view") != view
            or manifest.get("source") != context.source):
        raise HarnessError(f"the {view} product manifest identity is wrong: {path}")
    return path, manifest


def _accepted_digests(context, full):
    """Re-hash the accepted release this view aliases; assert the recorded provenance."""
    manifest_path = context.accepted_dir / "manifest.json"
    if not manifest_path.is_file():
        raise HarnessError(f"the accepted release manifest is missing: {manifest_path}")
    accepted = json.loads(manifest_path.read_text())
    accepted_sha = sha256(manifest_path)
    alias = full.get("alias") or {}
    if accepted_sha != alias.get("manifest_sha256"):
        raise HarnessError(
            "the accepted release manifest digest differs from the full product alias")
    payload = f"{context.source}.h5ad"
    record = (accepted.get("outputs") or {}).get(payload)
    if not record:
        raise HarnessError(f"the accepted release does not declare {payload}")
    payload_sha = sha256(context.accepted_dir / payload)
    if payload_sha != record["sha256"]:
        raise HarnessError("the accepted payload digest differs from its release manifest")
    if payload_sha != (alias.get("payload_sha256") or {}).get(context.source):
        raise HarnessError(
            "the accepted payload digest differs from the full product alias")
    return accepted, accepted_sha, payload_sha


def _merge_payload(context, merge_path, merge):
    """Re-hash the merge payload the pooled Δ is transported from."""
    payload = merge_path.parent / merge["payload"]["file"]
    if not payload.is_file():
        raise HarnessError(f"the merge payload is missing: {payload}")
    digest = sha256(payload)
    if digest != merge["payload"]["sha256"]:
        raise HarnessError(
            f"the merge payload digest {digest} differs from its manifest "
            f"{merge['payload']['sha256']}")
    if int(merge["payload"]["rows"]) != int(merge["rows"]):
        raise HarnessError("the merge payload declares unexpected rows")
    if int(merge["payload"]["cols"]) != int(context.size):
        raise HarnessError("the merge payload declares an unexpected gene axis")
    return payload, digest


def _merge_obs(payload, merge):
    """The merge product's row space, in its declared payload order."""
    import h5py
    from anndata.io import read_elem

    with h5py.File(payload, "r") as handle:
        obs = pd.DataFrame(read_elem(handle["obs"]))
        if "layers/delta" not in handle:
            raise HarnessError("the merge payload carries no delta layer")
    obs = obs.reset_index(drop=True)
    obs["profile_id"] = obs["profile_id"].astype(str)
    if len(obs) != int(merge["rows"]):
        raise HarnessError(
            f"the merge payload holds {len(obs)} rows, manifest says {merge['rows']}")
    return obs


def _merge_members(merge_path, merge):
    record = merge.get("members")
    if not record:
        raise HarnessError("the merge product does not declare its member table")
    table = pd.read_csv(merge_path.parent / record["file"],
                        dtype={"profile_id": str, "member_profile_id": str})
    return list(zip(table["profile_id"], table["member_profile_id"]))


def _condition_identities(context, conditions):
    return _condition_identities_by_source(conditions, context.source)


def _route_conditions(context, conditions, merge_obs):
    """Prove the merge rows are exactly the canonical condition rows, then route them."""
    identities = _condition_identities(context, conditions)
    observed = list(merge_obs["profile_id"])
    if observed != identities:
        raise HarnessError(
            "the merge payload rows are not the canonical condition rows this view labels")
    index = {identity: position for position, identity in enumerate(observed)}
    return np.asarray([index[identity] for identity in identities], dtype=np.int64)


def _condition_membership(context, conditions, merge_path, merge):
    """Reconcile the label projection's membership to the merge member table."""
    expected = []
    for identity, condition in zip(_condition_identities(context, conditions),
                                   conditions.itertuples()):
        for member in condition.profile_ids:
            expected.append((identity, str(member)))
    observed = _merge_members(merge_path, merge)
    if sorted(expected) != sorted(observed):
        raise HarnessError(
            "the conditions this view labels do not reconcile to the merge member table")


def _check_released(source, rows, conditions, merge):
    """Fail closed when the released row/condition population moves under the label."""
    treated = int((rows["role"].astype(str) == "treated").sum())
    observed = {"rows": int(len(rows)), "treated": treated,
                "vehicle": int(len(rows)) - treated, "conditions": int(len(conditions))}
    expected = {"rows": EXPECTED_ROWS[source], "treated": EXPECTED_TREATED[source],
                "vehicle": EXPECTED_VEHICLE[source],
                "conditions": EXPECTED_CONDITIONS[source]}
    if observed != expected:
        raise HarnessError(
            f"{source}: the released row space {observed} differs from {expected}")
    histogram = {}
    for count in range(1, int(conditions["member_count"].max()) + 1):
        histogram[str(count)] = int((conditions["member_count"] == count).sum())
    if histogram != EXPECTED_MEMBER_HISTOGRAM[source]:
        raise HarnessError(
            f"{source}: the member-count histogram {histogram} differs from "
            f"{EXPECTED_MEMBER_HISTOGRAM[source]}")
    merged = {key: value for key, value in
              merge["special_conditions"]["member_count_histogram"].items() if value}
    if {key: value for key, value in histogram.items() if value} != merged:
        raise HarnessError(
            f"{source}: the conditions differ from the merge product's member counts")


# --------------------------------------------------------------------------- #
# transport
# --------------------------------------------------------------------------- #
def _vehicle_delta(context, rows):
    """The full view's exact ``X - control`` for the vehicle rows (no condition to pool)."""
    import h5py

    payload = context.accepted_dir / f"{context.source}.h5ad"
    positions = np.flatnonzero(rows["role"].astype(str).to_numpy() != "treated")
    columns = int(context.size)
    delta = np.empty((len(positions), columns), dtype=np.float32)
    order = np.argsort(positions, kind="stable")
    selected = positions[order]
    with h5py.File(payload, "r") as handle:
        x_layer, control_layer = handle["X"], handle["layers/control"]
        if x_layer.shape[1] != columns or control_layer.shape[1] != columns:
            raise HarnessError("the accepted payload does not carry the frozen gene axis")
        for start in range(0, len(selected), BATCH):
            chunk = [int(row) for row in selected[start:start + BATCH]]
            x = np.asarray(x_layer[chunk], dtype=np.float32)
            control = np.asarray(control_layer[chunk], dtype=np.float32)
            delta[order[start:start + BATCH]] = x - control
    return positions, delta


def _merge_delta_matrix(payload, rows, columns):
    """Extract the merge delta layer into a row-indexed memmap (canonical merge order)."""
    import h5py

    descriptor, name = tempfile.mkstemp(prefix=".merge-delta-", suffix=".npy")
    os.close(descriptor)
    temporary = Path(name)
    matrix = np.lib.format.open_memmap(temporary, mode="w+", dtype=np.float32,
                                       shape=(rows, columns))
    with h5py.File(payload, "r") as handle:
        layer = handle["layers/delta"]
        if tuple(layer.shape) != (rows, columns):
            raise HarnessError("the merge delta layer shape differs from the payload")
        for start in range(0, rows, BATCH):
            stop = min(start + BATCH, rows)
            matrix[start:stop] = np.asarray(layer[start:stop], dtype=np.float32)
    matrix.flush()
    return temporary, matrix


def _delta_blocks(targets, merge_matrix, vehicle_positions, vehicle_delta, columns):
    """Stream product-order ``{"delta": ...}`` blocks for the shared product writer."""
    rows = len(targets)
    vehicle_slot = np.full(rows, -1, dtype=np.int64)
    vehicle_slot[vehicle_positions] = np.arange(len(vehicle_positions), dtype=np.int64)
    for start in range(0, rows, BATCH):
        stop = min(start + BATCH, rows)
        block = np.empty((stop - start, columns), dtype=np.float32)
        target = targets[start:stop]
        labelled = target >= 0
        if labelled.any():
            block[labelled] = merge_matrix[target[labelled]]
        if not labelled.all():
            slots = vehicle_slot[start:stop][~labelled]
            if (slots < 0).any():
                raise HarnessError(
                    "an unlabelled product row is not a declared vehicle row")
            block[~labelled] = vehicle_delta[slots]
        yield start, {"delta": block}


def _members_table(conditions, rows):
    """Treated product row -> every member profile id of its condition."""
    ids = rows["profile_id"].astype(str).to_numpy()
    entries = []
    for condition in conditions.itertuples():
        members = [str(member) for member in condition.profile_ids]
        for row in condition.members:
            identity = str(ids[int(row)])
            for member in members:
                entries.append({"profile_id": identity, "member_profile_id": member})
    return pd.DataFrame(entries, columns=["profile_id", "member_profile_id"])


def _conditions_table(conditions, rows, source):
    """Per-condition member counts with both identities, digest-pinned as a companion."""
    ids = rows["profile_id"].astype(str).to_numpy()
    identities = _condition_identities_by_source(conditions, source)
    return pd.DataFrame({
        "condition_key": [str(condition.condition_key)
                          for condition in conditions.itertuples()],
        "merge_profile_id": identities,
        "cell_line": conditions["cell_line"].astype(str).to_numpy(),
        "drug_parent_inchikey": conditions["drug_parent_inchikey"].astype(str).to_numpy(),
        "dose_um": conditions["dose_um"].astype(float).to_numpy(),
        "member_count": conditions["member_count"].astype(np.int64).to_numpy(),
        "member_profile_ids": [";".join(str(member) for member in condition.profile_ids)
                               for condition in conditions.itertuples()],
        "profile_ids": [";".join(str(ids[int(row)]) for row in condition.members)
                        for condition in conditions.itertuples()],
    })


def _condition_identities_by_source(conditions, source):
    return [merge_profile_id(source, condition.cell_line,
                             condition.drug_parent_inchikey, condition.dose_um)
            for condition in conditions.itertuples()]


def _histogram(conditions):
    counts = conditions["member_count"].astype(np.int64)
    return {str(count): int((counts == count).sum())
            for count in range(1, int(counts.max()) + 1)}


def build(context):
    """Write the pooled-label product for ``context.source`` and return its manifest."""
    source, root = context.source, Path(context.root)
    columns = int(context.size)
    full_path, full = _read_manifest(context, "full")
    merge_path, merge = _read_manifest(context, "merge")
    if full.get("payload") is not None:
        raise HarnessError("the full product must publish manifests only")
    accepted, accepted_sha, accepted_payload_sha = _accepted_digests(context, full)
    merge_payload, merge_payload_sha = _merge_payload(context, merge_path, merge)
    merge_manifest_sha = sha256(merge_path)
    merge_obs = _merge_obs(merge_payload, merge)

    # The row space is the full view's rows verbatim: identity comes from the accepted obs.
    with context.accepted_source() as accepted_source:
        rows = accepted_source.obs.copy()
    rows["profile_id"] = rows.index.astype(str)
    roles = set(rows["role"].astype(str))
    if roles - {"treated", "vehicle"}:
        raise HarnessError(f"unexpected full-view roles: {sorted(roles)}")
    if len(rows) != int(full["rows"]) or rows["profile_id"].astype(str).duplicated().any():
        raise HarnessError("the full row space differs from the full product manifest")

    conditions, row_to_condition = condition_table(rows, source)
    route = _route_conditions(context, conditions, merge_obs)
    _condition_membership(context, conditions, merge_path, merge)
    if not context.fixture:
        _check_released(source, rows, conditions, merge)

    # Product row -> merge row (treated) or vehicle slot (no condition to pool).
    treated = np.flatnonzero(rows["role"].astype(str).to_numpy() == "treated")
    targets = np.full(len(rows), -1, dtype=np.int64)
    targets[treated] = route[row_to_condition[treated]]
    if (targets[treated] < 0).any():
        raise HarnessError("a treated row has no merge label")
    vehicle_positions, vehicle_delta = _vehicle_delta(context, rows)

    table = _conditions_table(conditions, rows, source)
    table_path = context.staging_dir / f"{source}.conditions.csv"
    table.to_csv(table_path, index=False, lineterminator="\n")
    members = _members_table(conditions, rows)

    temporary, merge_matrix = _merge_delta_matrix(merge_payload, len(merge_obs), columns)
    try:
        manifest = context.write_product(
            rows=rows,
            blocks=_delta_blocks(targets, merge_matrix, vehicle_positions, vehicle_delta,
                                 columns),
            mode="delta", recipe=RECIPE["full-pooled-label"],
            control_rule=CONTROL_RULE["full-pooled-label"],
            alias={
                "release": relative(context.accepted_dir, root),
                "manifest_sha256": accepted_sha,
                "payload_sha256": {name: accepted["outputs"][f"{name}.h5ad"]["sha256"]
                                   for name in SOURCES
                                   if f"{name}.h5ad" in accepted["outputs"]},
                "outputs": {name: record["sha256"]
                            for name, record in sorted(accepted["outputs"].items())},
                "recipe": accepted["recipe"],
                "normalization": accepted["normalization"],
                "control": accepted["control"],
            },
            members=members, control_sources=None, delta_kind="pooled", seed=None,
            selection=None)
    finally:
        del merge_matrix
        Path(temporary).unlink(missing_ok=True)

    if manifest["profile_ids"]["sha256"] != full["profile_ids"]["sha256"]:
        raise HarnessError(
            "the pooled-label row identities differ from the full product's rows")
    histogram = _histogram(conditions)
    pooled_label = {
        "delta_kind": "pooled",
        "label_source": {
            "view": "merge",
            "product": merge["product"],
            "rows": int(merge["rows"]),
            "conditions": int(merge["conditions"]),
            "manifest": {"file": relative(merge_path, root),
                         "sha256": merge_manifest_sha},
            "payload": {"file": merge["payload"]["file"],
                        "sha256": merge_payload_sha,
                        "bytes": int(merge["payload"]["bytes"]),
                        "rows": int(merge["payload"]["rows"]),
                        "cols": int(merge["payload"]["cols"]),
                        "dtype": merge["payload"]["dtype"],
                        "layers": list(merge["payload"]["layers"])},
            "layer": "delta",
            "transfer": TRANSFER,
            "row_key": ("cell_line|drug_parent_inchikey|dose_um -> "
                        "merge::<source>::<cell_line>::<parent>::<dose>"),
            "recipe": merge["recipe"],
            "control_rule": merge["control_rule"],
        },
        "transport": TRANSPORT,
        "vehicle_rule": VEHICLE_RULE,
        "conditions": int(len(conditions)),
        "treated_rows": int(len(treated)),
        "vehicle_rows": int(len(rows) - len(treated)),
        "member_count_histogram": histogram,
        "member_count": {
            "min": int(conditions["member_count"].min()),
            "median": float(conditions["member_count"].median()),
            "max": int(conditions["member_count"].max()),
            "sum": int(conditions["member_count"].sum()),
        },
        "conditions_table": {
            "file": table_path.name, "sha256": sha256(table_path),
            "bytes": table_path.stat().st_size, "rows": int(len(table)),
            "columns": list(table.columns),
        },
    }
    manifest["full"] = {
        "view": "full",
        "product": full["product"],
        "rows": int(full["rows"]),
        "treated": int(full["treated"]),
        "manifest": {"file": relative(full_path, root),
                     "sha256": sha256(full_path)},
        "profile_ids": full["profile_ids"],
        "payload": {"file": f"{source}.h5ad",
                    "sha256": accepted_payload_sha,
                    "bytes": int(accepted["outputs"][f"{source}.h5ad"]["bytes"]),
                    "rows": int(full["rows"]), "cols": columns, "dtype": "float32",
                    "layers": ["X", "control"], "aliased": True},
        "alias": full.get("alias"),
    }
    manifest["pooled_label"] = pooled_label
    write_json(context.staging_dir / f"{source}.json", manifest)
    context.evidence.update({
        "rows": int(len(rows)),
        "treated": int(len(treated)),
        "conditions": int(len(conditions)),
        "delta_payload": {"file": manifest["payload"]["file"],
                          "sha256": manifest["payload"]["sha256"],
                          "bytes": int(manifest["payload"]["bytes"])},
        "full_payload_sha256": accepted_payload_sha,
        "full_manifest_sha256": manifest["full"]["manifest"]["sha256"],
        "merge_payload_sha256": merge_payload_sha,
        "merge_manifest_sha256": merge_manifest_sha,
        "member_count_histogram": histogram,
        "member_count": pooled_label["member_count"],
        "vehicle_rows": pooled_label["vehicle_rows"],
        "conditions_table": pooled_label["conditions_table"],
        "members_rows": int(len(members)),
    })
    return manifest
