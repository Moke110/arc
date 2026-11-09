"""Plot frozen ARC group sizes using only CSVs and Tahoe observation metadata."""

import json
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.ticker import PercentFormatter

OUTPUT = Path(__file__).resolve().parent
ROOT = next(parent for parent in OUTPUT.parents if (parent / ".git").exists())
STAGE = ROOT / "stage2"
CORPUS = STAGE / "02_dataset_corpus/output/corpus/v2"
STEP4 = STAGE / "04_pdiff_hyperparameter"
RELEASES = {
    "Tahoe": STEP4 / "execution/output/tahoe-smoke-ticket68/release",
    "tenv": STAGE / "03_modes_tenv/output/input",
}


def read_groups(dataset):
    release = RELEASES[dataset]
    wells = pd.read_csv(release / "wells.csv")
    controls = pd.read_csv(release / "controls.csv").set_index("control_index")
    test = pd.concat(
        [pd.read_csv(release / f"fold_{fold}/test_cases.csv") for fold in range(3)],
        ignore_index=True,
    )
    keys = ["case_id", "group_id", "control_index"]
    pd.testing.assert_frame_equal(
        wells[keys].sort_values("case_id").reset_index(drop=True),
        test[keys].sort_values("case_id").reset_index(drop=True),
    )
    assert wells.case_id.is_unique
    members = wells[["group_id", "control_index"]].drop_duplicates()
    members = members.join(controls[["n_wells"]], on="control_index")
    groups = wells.groupby("group_id").size().rename("treated_profiles").to_frame()
    groups["control_profiles"] = members.groupby("group_id").n_wells.sum()
    groups["baseline_vectors"] = members.groupby("group_id").size()

    plate_groups = None
    if dataset == "Tahoe":
        with h5py.File(CORPUS / "tahoe.h5ad") as source:
            obs = source["obs"]

            def column(key):
                entry = obs[key]
                if isinstance(entry, h5py.Group):
                    return entry["categories"].asstr()[:][entry["codes"][:]]
                return entry.asstr()[:] if entry.dtype.kind in "OS" else entry[:]

            metadata = pd.DataFrame({
                key: column(key) for key in
                ("profile_id", "n_cells", "drug", "is_control", "plate", "cell_line_id")
            })
        treated = wells.join(
            metadata.set_index("profile_id")[["n_cells"]], on="case_id", validate="one_to_one"
        )
        dmso = metadata.loc[metadata.is_control & metadata.drug.eq("DMSO_TF")]
        pools = dmso.groupby(["plate", "cell_line_id"]).n_cells.agg(["size", "sum"])
        controls = controls.join(pools, on=["plate", "cell_line_id"], validate="one_to_one")
        assert controls.n_wells.eq(controls["size"]).all() and treated.n_cells.notna().all()
        groups["treated_cells"] = treated.groupby("group_id").n_cells.sum()
        groups["control_cells"] = (
            members.join(controls[["sum"]], on="control_index").groupby("group_id")["sum"].sum()
        )
        plate_groups = (
            treated.groupby(["group_id", "control_index"]).n_cells.sum()
            .rename("treated_cells").reset_index()
            .join(controls[["sum"]], on="control_index")
            .rename(columns={"sum": "control_cells"})
        )
    return groups, plate_groups


def stats(values):
    return {
        "min": int(values.min()), "p05": float(values.quantile(0.05)),
        "median": float(values.median()), "p95": float(values.quantile(0.95)),
        "max": int(values.max()),
    }


def draw_bars(axis, groups, columns, edges, title, xlabel):
    labels = ["Treated", "Matched DMSO pool"]
    colors = ["#236894", "#c26331"]
    slots = np.arange(len(edges) - 1)
    rows = []
    for offset, field, label, color in zip((-0.19, 0.19), columns, labels, colors):
        counts = np.histogram(groups[field], bins=edges)[0]
        assert counts.sum() == len(groups)
        fraction = counts / len(groups)
        axis.bar(slots + offset, fraction, width=0.37, color=color, label=label, zorder=3)
        rows.extend(zip([field] * len(slots), edges[:-1], edges[1:], counts, fraction))
    axis.set(xlabel=xlabel, ylabel="Share of groups")
    axis.set_title(title, loc="left", fontsize=13, fontweight="bold")
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    axis.set_xticks(slots)
    axis.grid(axis="y", color="#e4e9ee", zorder=0)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, loc="upper right", fontsize=9)
    return rows


def main():
    tahoe, plate_groups = read_groups("Tahoe")
    tenv, _ = read_groups("tenv")
    datasets = {"Tahoe": tahoe, "tenv": tenv}
    summary = {}
    for dataset, groups in datasets.items():
        widths = groups.treated_profiles
        summary[dataset] = {
            "groups": len(groups),
            "treated_wells": int(widths.sum()),
            "groups_by_treated_wells": widths.value_counts().sort_index().to_dict(),
            "groups_sampled_with_replacement": int(widths.lt(3).sum()),
            "groups_with_equal_treated_and_control_wells": int(widths.eq(groups.control_profiles).sum()),
            **{field: stats(groups[field]) for field in groups},
        }
        if dataset == "Tahoe":
            summary[dataset]["groups_with_equal_treated_and_control_cells"] = int(
                groups.treated_cells.eq(groups.control_cells).sum()
            )
    summary["Tahoe_plate_stratified"] = {
        "groups": len(plate_groups),
        "treated_cells": stats(plate_groups.treated_cells),
        "control_cells": stats(plate_groups.control_cells),
        "groups_with_fewer_than_256_treated_cells": int(plate_groups.treated_cells.lt(256).sum()),
        "groups_with_fewer_than_256_control_cells": int(plate_groups.control_cells.lt(256).sum()),
    }
    (OUTPUT / "perturbdiff_population_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    fig = plt.figure(figsize=(13, 8), layout="constrained")
    grid = fig.add_gridspec(2, 2)
    fig.suptitle("PerturbDiff population sizes in the ARC benchmark", fontsize=18, fontweight="bold")
    rows = []
    for pos, (dataset, groups) in enumerate(datasets.items()):
        axis = fig.add_subplot(grid[0, pos])
        edges = [1, 2, 3, 4, 5, 6, 11, 16, 21, 26, 31]
        rows.extend((dataset, *row) for row in draw_bars(
            axis, groups, ("treated_profiles", "control_profiles"), edges,
            f"{'AB'[pos]}   {dataset} | {len(groups):,} groups",
            "Source well profiles per group (before sampling or averaging)",
        ))
        axis.set_xticklabels(["1", "2", "3", "4", "5", "6–10", "11–15", "16–20", "21–25", "26–30"])
        axis.set_ylim(0, 1)
    edges = [1, 100, 250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 160000]
    for pos, (name, groups, title) in enumerate((
        ("Tahoe", tahoe, "C   Tahoe | Cells per ARC condition"),
        ("Tahoe_plate_stratified", plate_groups, "D   Tahoe | Cells per condition within plate"),
    )):
        axis = fig.add_subplot(grid[1, pos])
        rows.extend((name, *row) for row in draw_bars(
            axis, groups, ("treated_cells", "control_cells"), edges, title,
            "Single cells underlying the well profiles",
        ))
        axis.set_xticklabels(
            ["1–99", "100–249", "250–499", "500–999", "1–2k",
             "2–4k", "4–8k", "8–16k", "16–32k", "32k+"],
            fontsize=8,
        )
        axis.set_ylim(0, 0.55)
    fig.savefig(OUTPUT / "perturbdiff_population_sizes.png", dpi=180, facecolor="white")
    pd.DataFrame(rows, columns=[
        "dataset", "field", "bin_start", "bin_stop", "groups", "fraction"
    ]).to_csv(OUTPUT / "perturbdiff_population_histogram.csv", index=False)


if __name__ == "__main__":
    main()
