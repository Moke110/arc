"""The `merge` view: one replicate-pooled profile per condition (ADR 0082, plan C7.1-C7.2).

Treated expression is the log-normalized **sum** of each condition's member profiles'
retained G_5000 integer counts: the integer counts are added first and the pooled row is
normalized on its own selected-axis total, so a pooled profile is never the mean of the
members' normalized values. The control is pooled by the same rule over the union of the
members' matched control groups' vehicle wells. A repeat that spans two plates therefore
becomes one profile with one deepened control.

The recipe itself is delegated to the shared harness primitives; this module only orders
them, asserts the released-population invariants, names the tenv trap conditions in the
product manifest and records the pooling-identity probe that proves the sum-based rule.
"""
import numpy as np
import pandas as pd

from . import (BATCH, CONTROL_RULE, RECIPE, HarnessError, condition_table,
               control_source_table, merge_profile_id, merge_rows, member_table,
               pooled_condition_blocks, pooled_condition_counts, pooled_control,
               vehicle_group_sums, write_json)

# Released-replica expectations from the reconnaissance (notes-bases-views §0/§5). They are
# asserted only for a non-fixture build: a released corpus change must fail closed here
# rather than silently publish a differently shaped view.
EXPECTED_CONDITIONS = {"tenv": 1377, "tahoe_bulk": 52777}
EXPECTED_MEMBER_HISTOGRAM = {
    "tenv": {"3": 1320, "4": 48, "6": 9},
    "tahoe_bulk": {"1": 43036, "2": 7941, "3": 1400, "4": 350, "16": 50},
}
EXPECTED_CROSS_PLATE = {"tenv": 6, "tahoe_bulk": 8741}
EXPECTED_TRAPS = {"tenv": {"four_well": 48, "six_well": 9}}
IDENTITY_SAMPLE = 256
CONTROL_TOKEN = "pooled-sum"


def _entry(row, source):
    """One named trap condition, compact enough for a tracked manifest."""
    return {
        "profile_id": merge_profile_id(source, row["cell_line"],
                                       row["drug_parent_inchikey"], row["dose_um"]),
        "cell_line": str(row["cell_line"]),
        "drug_parent_inchikey": str(row["drug_parent_inchikey"]),
        "drug_parent_smiles": str(row["drug_parent_smiles"]),
        "dose_um": float(row["dose_um"]),
        "compound_names": [str(name) for name in row["compound_names"]],
        "plates": [str(plate) for plate in row["plates"]],
        "control_groups": [str(group) for group in row["control_groups"]],
        "member_count": int(row["member_count"]),
        "member_profile_ids": [str(member) for member in row["profile_ids"]],
    }


def special_conditions(conditions, source):
    """Name the repeats a reader must be able to find in the manifest.

    Tenv carries two trap families the spec calls out: the 48 four-well conditions (each
    plate's single reference compound at one dose) and the 9 six-well conditions (two test
    articles of one parent — Tamoxifen/Tamoxifen citrate on one plate, the PFOA and PFOS
    free-acid/salt pairs across two plates). A member-count histogram and the cross-plate
    condition count are recorded for every source so the pooled row space reconciles to the
    members it replaced.
    """
    frame = pd.DataFrame(conditions)
    if not len(frame):
        raise HarnessError(f"{source}: no condition to describe")
    member_counts = frame["member_count"].to_numpy(dtype=np.int64)
    histogram = {}
    for count in range(1, int(member_counts.max()) + 1):
        histogram[str(count)] = int((member_counts == count).sum())
    cross_plate = int((frame["plates"].apply(len) > 1).sum())
    tenv = source == "tenv"
    four = frame[member_counts == 4] if tenv else frame.iloc[0:0]
    six = frame[member_counts == 6] if tenv else frame.iloc[0:0]
    return {
        "source_has_named_traps": bool(tenv and (len(four) or len(six))),
        "member_count_histogram": histogram,
        "cross_plate_conditions": cross_plate,
        "four_well_reference_compounds": [_entry(row, source) for _, row in four.iterrows()],
        "six_well_salt_pairs": [_entry(row, source) for _, row in six.iterrows()],
    }


def _check_released(source, conditions, special):
    """Fail closed if the released corpus no longer matches the locked view population."""
    observed = int(len(conditions))
    if observed != EXPECTED_CONDITIONS[source]:
        raise HarnessError(
            f"{source}: the released corpus yields {observed} conditions, expected "
            f"{EXPECTED_CONDITIONS[source]}")
    histogram = {key: value for key, value in
                 special["member_count_histogram"].items() if value}
    if histogram != EXPECTED_MEMBER_HISTOGRAM[source]:
        raise HarnessError(
            f"{source}: the member-count histogram {histogram} differs from "
            f"{EXPECTED_MEMBER_HISTOGRAM[source]}")
    if special["cross_plate_conditions"] != EXPECTED_CROSS_PLATE[source]:
        raise HarnessError(
            f"{source}: {special['cross_plate_conditions']} cross-plate conditions, "
            f"expected {EXPECTED_CROSS_PLATE[source]}")
    if source in EXPECTED_TRAPS:
        traps = EXPECTED_TRAPS[source]
        four = special["four_well_reference_compounds"]
        six = special["six_well_salt_pairs"]
        if len(four) != traps["four_well"] or len(six) != traps["six_well"]:
            raise HarnessError(
                f"{source}: the named traps do not reconcile ({len(four)} four-well, "
                f"{len(six)} six-well)")
        if any(len(entry["compound_names"]) != 1 for entry in four):
            raise HarnessError("a four-well reference-compound condition names several articles")
        if any(len(entry["compound_names"]) != 2 for entry in six):
            raise HarnessError("a six-well salt-pair condition does not name two articles")


def _pooling_identity(reader, genes, row_to_condition, pooled, conditions, view_api,
                      limit=IDENTITY_SAMPLE):
    """Prove on a bounded sample that a pooled row is not the members' normalized mean.

    The pooled accumulator is already in memory; this probe streams the corpus once more and
    accumulates ``mean(normalized(member))`` for at most ``limit`` repeated conditions. It
    fails closed when the two coincide, which is the one way the normalized-sum rule could
    silently degrade to a mean-of-normalized-values rule.
    """
    repeated = np.flatnonzero(conditions["member_count"].to_numpy(dtype=np.int64) > 1)
    if not len(repeated):
        return {"checked_conditions": 0, "checked_members": 0, "differing_rows": 0,
                "max_abs_difference": 0.0, "rule": "normalize-of-sum"}
    sample = repeated[:limit]
    slot = {int(condition): position for position, condition in enumerate(sample)}
    accumulator = np.zeros((len(sample), len(genes)), dtype=np.float64)
    members = np.zeros(len(sample), dtype=np.int64)
    library_totals = np.zeros(len(reader.obs), dtype=np.int64)
    for condition in conditions.itertuples():
        library_totals[np.asarray(condition.members, dtype=np.int64)] = (
            condition.library_totals.astype(np.int64))
    for start, block in reader.batches(BATCH):
        stop = start + len(block)
        local = np.flatnonzero(np.isin(row_to_condition[start:stop], sample))
        for offset in local:
            row = start + int(offset)
            position = slot[int(row_to_condition[row])]
            counts = block[offset][genes].astype(np.int64).reshape(1, -1)
            total = np.asarray([library_totals[row]], dtype=np.int64)
            if int(total[0]) <= 0:
                raise HarnessError(
                    f"{reader.source}: profile {reader.obs.index[row]} has no retained total")
            accumulator[position] += view_api.log1p_cp10k(counts, total)[0]
            members[position] += 1
    mean = accumulator / members[:, None]
    pooled_counts = np.asarray(pooled[sample], dtype=np.int64)
    pooled_x = view_api.log1p_cp10k(pooled_counts, view_api.library_totals(pooled_counts))
    difference = np.abs(pooled_x - mean)
    differing = int((difference.max(axis=1) > 1e-6).sum())
    if not differing:
        raise HarnessError(
            "the pooled expression equals the mean of the members' normalized values on "
            f"{len(sample)} sampled conditions; the normalize-of-sum rule is not in effect")
    return {"checked_conditions": int(len(sample)), "checked_members": int(members.sum()),
            "differing_rows": differing, "max_abs_difference": float(difference.max()),
            "rule": ("sum the members' retained integer counts, then normalize the sum on "
                     "its own selected-axis total; never a mean of normalized values")}


def build(context):
    """Write the pooled product for ``context.source`` and return its manifest."""
    source = context.source
    view_api = context.modules["view"]
    with context.accepted_source() as accepted:
        obs = accepted.obs.copy()
    obs["profile_id"] = obs.index.astype(str)
    conditions, row_to_condition = condition_table(obs, source)
    genes = context.genes

    with context.corpus(source) as reader:
        pooled = pooled_condition_counts(reader, genes, row_to_condition, conditions)
        group_sums = vehicle_group_sums(reader, genes)
        identity = _pooling_identity(reader, genes, row_to_condition, pooled, conditions,
                                     view_api)

    # The pooled treated total is the sum of the members' retained totals; the pooled
    # control is normalized on the sum of the union of the members' vehicle wells' counts.
    control, control_totals = pooled_control(conditions, group_sums)
    expected_totals = np.asarray(
        [int(condition.library_totals.sum(dtype=np.int64))
         for condition in conditions.itertuples()], dtype=np.int64)
    treated_totals = np.asarray(pooled.sum(axis=1), dtype=np.int64)
    if not np.array_equal(treated_totals, expected_totals):
        raise HarnessError(
            f"{source}: a pooled retained total differs from the sum of its members' totals")

    source_sha = context.source_manifests()["accepted_view_sha256"]
    rows = merge_rows(conditions, source, treated_totals, source_sha)
    members = member_table(
        conditions, lambda condition: merge_profile_id(
            source, condition.cell_line, condition.drug_parent_inchikey, condition.dose_um))
    controls = control_source_table(conditions, rows["profile_id"].tolist(),
                                    CONTROL_TOKEN, control_totals)
    manifest = context.write_product(
        rows=rows, blocks=pooled_condition_blocks(pooled, control), mode="payload",
        recipe=RECIPE["merge"], control_rule=CONTROL_RULE["merge"], members=members,
        control_sources=controls, seed=None, selection=None)

    special = special_conditions(conditions, source)
    if not context.fixture:
        _check_released(source, conditions, special)
    manifest["special_conditions"] = special
    manifest["pooling_identity"] = identity
    manifest["pooled_retained_total"] = int(control_totals.sum(dtype=np.int64))
    write_json(context.staging_dir / f"{source}.json", manifest)
    context.evidence.update({
        "pooling_rule": CONTROL_RULE["merge"],
        "conditions": int(len(conditions)),
        "members": int(rows["member_count"].sum()),
        "cross_plate_conditions": special["cross_plate_conditions"],
        "member_count_histogram": special["member_count_histogram"],
        "named_traps": {"four_well_reference_compounds": len(
                            special["four_well_reference_compounds"]),
                        "six_well_salt_pairs": len(special["six_well_salt_pairs"]),
                        "source_has_named_traps": special["source_has_named_traps"]},
        "pooling_identity": identity,
        "pooled_treated_retained_total": int(treated_totals.sum(dtype=np.int64)),
        "pooled_control_retained_total": int(control_totals.sum(dtype=np.int64)),
    })
    return manifest
