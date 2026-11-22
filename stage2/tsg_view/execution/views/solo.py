"""The member-choice views: `solo` (draw) and `solo-deepest` (largest retained total).

Both views keep exactly one treated profile per condition and reuse that member's
**full-view** row verbatim: its retained total, its ``X`` and its matching-control value are
the accepted release's values at that row, so the view is a pure row subset of ``full`` and
the vehicle tail is unchanged. The product rows are ordered by the chosen ``profile_id``,
and the choice itself consumes only the canonical condition table plus the frozen seed, so
folds, corpus row order and training seeds cannot change a single chosen member.

``solo`` draws with :func:`views.solo_choice` from the frozen ``SOLO_SEED``; ``solo-deepest``
is the sensitivity arm that takes the member with the largest retained selected-axis total
(ties break to the smallest ``profile_id``) and therefore equals ``solo`` on singletons.
"""
import hashlib

import numpy as np

from . import (CONTROL_RULE, RECIPE, SOLO_SEED, condition_table, control_source_table,
               member_table, selected_blocks, selected_groups, selected_rows, solo_choice)

SELECTION = {
    "draw": ("one member per condition drawn from sha256(SOLO_SEED|source|condition_key) "
             "modulo the member count (members ordered by profile_id); independent of "
             "folds, corpus row order and training seeds"),
    "deepest": ("one member per condition: the largest retained selected-axis library_total, "
                "ties breaking to the smallest profile_id"),
}


def chosen_profile_ids(conditions, chosen):
    """The chosen member's own ``profile_id`` per condition, in condition order."""
    ids = []
    for position, row in enumerate(np.asarray(chosen, dtype=np.int64)):
        condition = conditions.iloc[position]
        members = [int(member) for member in condition["members"]]
        profile_ids = [str(profile) for profile in condition["profile_ids"]]
        ids.append(profile_ids[members.index(int(row))])
    return np.asarray(ids, dtype=object)


def _ordered_conditions(conditions, order, groups):
    """Conditions in row order, with each row's own matched control group.

    ``solo`` keeps the chosen member's full-view control, which averages that member's own
    ``matched_control_group`` — not the union of every member's group. Recording the
    condition union would overstate the vehicles behind a cross-plate row, so the
    ``control_sources.csv`` provenance is narrowed to the chosen member's own group.
    """
    ordered = conditions.iloc[order].reset_index(drop=True).copy()
    ordered["control_groups"] = [[str(group)] if str(group) else [] for group in groups]
    return ordered


def materialize(context, *, rule, selection, seed, recipe, control_rule):
    """Write one member-choice product; the harness owns the schema and the payload."""
    source = context.source
    with context.accepted_source() as reader:
        obs = reader.obs.copy()
        conditions, _ = condition_table(obs, source)
        chosen = solo_choice(conditions, rule=rule, seed=SOLO_SEED, source=source)
        profile_ids = chosen_profile_ids(conditions, chosen)
        order = np.argsort(profile_ids, kind="stable")
        positions = np.asarray(chosen, dtype=np.int64)[order]
        identities = profile_ids[order]
        groups = selected_groups(reader, positions)
        ordered = _ordered_conditions(conditions, order, groups)
        rows = selected_rows(reader, positions)
        members = member_table(ordered, identities)
        controls = control_source_table(ordered, identities, control_rule)
        if not len(rows) == len(identities) == len(conditions):
            raise AssertionError("one product row per condition is required")
        manifest = context.write_product(
            rows=rows, blocks=selected_blocks(reader, positions), mode="payload",
            recipe=recipe, control_rule=control_rule, members=members,
            control_sources=controls, seed=seed, selection=selection)
    discarded = int(conditions["member_count"].sum()) - len(conditions)
    context.evidence.update({
        "selection_rule": rule,
        "selection": selection,
        "seed": seed,
        "source_conditions": int(len(conditions)),
        "discarded_members": discarded,
        "multi_member_conditions": int((conditions["member_count"] > 1).sum()),
        "chosen_profile_id_sha256": hashlib.sha256(
            "\n".join(str(identity) for identity in identities).encode()).hexdigest(),
        "fold_and_training_seed_independent": True,
    })
    return manifest


def build(context):
    return materialize(context, rule="draw", selection=SELECTION["draw"],
                       seed=SOLO_SEED, recipe=RECIPE["solo"],
                       control_rule=CONTROL_RULE["solo"])
