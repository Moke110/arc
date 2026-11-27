"""Paired TranSiGen architecture scores on the same folds."""

from tsg_axis.metrics import summarize_folds


def pair_shape_scores(fold_scores, *, control):
    baseline = fold_scores.loc[fold_scores["shape"] == control, ["fold", "metric", "value"]]
    pairs = fold_scores.merge(baseline.rename(columns={"value": "control_score"}),
                              on=["fold", "metric"], how="left", validate="many_to_one", indicator=True)
    if not pairs["_merge"].eq("both").all():
        raise ValueError("every score needs the original architecture on the same fold and metric")
    pairs = pairs.drop(columns="_merge")
    pairs["score_difference"] = pairs.value - pairs.control_score
    summary = summarize_folds(pairs, ["shape", "metric"])
    difference = pairs.groupby(["shape", "metric"], sort=False).score_difference.agg(["mean", "count"])
    summary = summary.merge(difference.rename(columns={"mean": "score_difference", "count": "n_paired_folds"}),
                            on=["shape", "metric"], validate="one_to_one")
    return pairs, summary
