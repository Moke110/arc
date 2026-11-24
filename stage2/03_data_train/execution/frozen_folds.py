"""Project Step 07 frozen member labels to pooled Step 03 rows."""
from training_data import STRIPS

PARTITIONS = {"train": "train", "validation": "valid", "test": "test"}


def group_labels(folds, observations, groups, *, strip_label="excluded"):
    member_labels = folds.partitions.reindex(observations.index)
    identities = folds.identities.reindex(observations.index)
    if (member_labels.isna().any().any()
            or not identities.drug_parent.eq(observations.drug_parent_inchikey).all()
            or not identities.dose_um.eq(observations.dose_um).all()):
        raise ValueError("missing frozen fold assignment or mismatched member identity")
    per_group = member_labels.groupby(groups, sort=False)
    if not per_group.nunique().eq(1).all().all():
        raise ValueError("frozen fold assignments disagree between members of a pooled row")
    return [
        [strip_label if label in STRIPS else PARTITIONS[label] for label in original] + [
            STRIPS.get(label, "") for label in original[3:]
        ]
        for original in per_group.first().itertuples(index=False, name=None)
    ]
