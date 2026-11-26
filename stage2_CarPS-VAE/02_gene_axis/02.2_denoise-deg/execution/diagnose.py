"""Measure tenv detection and DMSO variability without applying gene exclusions."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def measure(source, genes):
    columns = source.native_columns(genes.gene_index)
    obs = source.obs
    dmso = (obs.role == "vehicle") & (obs.vehicle == "DMSO")
    biological = obs.role.isin(["treated", "vehicle"])
    retained = np.zeros(len(obs), dtype=bool)
    reads = np.zeros(len(genes), dtype=np.int64)
    detected = reads.copy()
    normalized = np.zeros(len(genes))
    control_blocks = []
    control_rows = []
    library_sizes = np.zeros(len(obs), dtype=np.int64)
    for start, counts in source.batches(128):
        rows = slice(start, start + len(counts))
        depth = counts.sum(axis=1)
        library_sizes[rows] = depth
        eligible = biological.to_numpy()[rows] & (depth > 0)
        retained[rows] = eligible
        group = dmso.to_numpy()[rows] & eligible
        common = counts[:, columns]
        values = common[eligible]
        reads += values.sum(axis=0)
        detected += (values > 0).sum(axis=0)
        normalized += (values * (1e6 / depth[eligible, None])).sum(axis=0)
        control_rows.extend(np.flatnonzero(group) + start)
        control_blocks.append(np.log2(1 + common[group] * (1e6 / depth[group, None])))
    controls = np.concatenate(control_blocks)
    metadata = obs.iloc[control_rows][["plate", "well"]].reset_index(drop=True)
    metadata["native_library_size"] = library_sizes[control_rows]
    metadata["median_peer_correlation"] = np.nan

    means, within_ss, sizes = [], [], []
    for _, positions in metadata.groupby("plate", sort=False).indices.items():
        group = controls[positions]
        means.append(group.mean(axis=0))
        within_ss.append(((group - means[-1]) ** 2).sum(axis=0))
        sizes.append(len(group))
        correlation = np.corrcoef(group)
        np.fill_diagonal(correlation, np.nan)
        metadata.loc[positions, "median_peer_correlation"] = np.nanmedian(correlation, axis=1)
    sizes = np.asarray(sizes)
    means = np.asarray(means)
    n = sizes.sum()
    plates = len(sizes)
    if plates < 2 or (sizes < 2).any():
        raise ValueError("plate variance decomposition needs replicated controls on multiple plates")
    grand_mean = np.average(means, axis=0, weights=sizes)
    within_var = np.sum(within_ss, axis=0) / (n - plates)
    between_ms = ((means - grand_mean) ** 2 * sizes[:, None]).sum(axis=0) / (plates - 1)
    effective_n = (n - (sizes ** 2).sum() / n) / (plates - 1)
    plate_var = np.maximum((between_ms - within_var) / effective_n, 0)
    total = plate_var + within_var

    evidence = genes[["gene_index", "symbol"]].copy()
    evidence["detection_fraction"] = detected / retained.sum()
    evidence["mean_raw_count"] = reads / retained.sum()
    evidence["mean_cpm"] = normalized / retained.sum()
    # Exact per-gene quantiles, with only one small gene block resident.
    q95 = np.empty(len(genes))
    for lo in range(0, len(genes), 2048):
        hi = min(lo + 2048, len(genes))
        abundance = np.empty((retained.sum(), hi - lo), dtype=np.float32)
        offset = 0
        for start, counts in source.batches(128):
            rows = slice(start, start + len(counts))
            keep = retained[rows]
            depth = library_sizes[rows][keep]
            next_offset = offset + len(depth)
            abundance[offset:next_offset] = (
                counts[keep][:, columns[lo:hi]] * (1e6 / depth[:, None])
            )
            offset = next_offset
        q95[lo:hi] = np.quantile(abundance, .95, axis=0, overwrite_input=True)
    evidence["p95_cpm"] = q95
    evidence["dmso_mean_log2_cpm1"] = grand_mean
    evidence["plate_sd"] = np.sqrt(plate_var)
    evidence["within_plate_sd"] = np.sqrt(within_var)
    evidence["plate_variance_fraction"] = np.divide(
        plate_var, total, out=np.full_like(total, np.nan), where=total > 0
    )
    return evidence, metadata, {
        "observed_roles": obs.role.value_counts().to_dict(),
        "expression_wells": int(retained.sum()),
        "expression_wells_by_role": obs.loc[retained, "role"].value_counts().to_dict(),
        "zero_depth_nonblank_wells": int((biological.to_numpy() & (library_sizes == 0)).sum()),
        "dmso_wells": int(n),
        "plate_count": plates,
        "dmso_wells_per_plate": metadata.plate.value_counts().to_dict(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--step", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo / "stage2_CarPS-VAE/01_data_raw/output/api"))
    from native_corpus import open_corpus

    genes = pd.read_csv(args.repo / "stage2_CarPS-VAE/util/output/gene-identity/genes.csv")
    with open_corpus(args.repo / "stage2_CarPS-VAE/01_data_raw/output") as corpus:
        evidence, controls, summary = measure(corpus.source("tenv"), genes)
        summary["corpus_manifest_sha256"] = corpus.manifest_sha256
    output = args.step / "figure/data"
    output.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(output / "diagnostics_genes.csv.gz", index=False)
    controls.to_csv(output / "diagnostics_dmso_wells.csv", index=False)
    x = evidence.dmso_mean_log2_cpm1
    measured = x > 0
    cases = []
    exclusions = evidence[["gene_index", "symbol"]].assign(
        dmso_measured=measured, low_detection=evidence.detection_fraction < .10
    )
    for target in (.05, .075, .10):
        excluded = []
        for metric in ("plate_sd", "within_plate_sd"):
            y = evidence[metric]
            k = (y[measured] / x[measured]).quantile(1 - target)
            above = measured & (y > k * x)
            excluded.append(above)
            cases.append({"target": target, "metric": metric, "k": k,
                          "dmso_measured_genes": int(measured.sum()), "above": int(above.sum())})
            if target == .10:
                exclusions[metric + "_above_boundary"] = above
        cases[-1]["union_above"] = int((excluded[0] | excluded[1]).sum())
        cases[-1]["intersection_above"] = int((excluded[0] & excluded[1]).sum())
    pd.DataFrame(cases).to_csv(output / "diagnostics_variance_boundaries.csv", index=False)
    exclusions["excluded"] = (
        exclusions.low_detection | exclusions.plate_sd_above_boundary
        | exclusions.within_plate_sd_above_boundary
    )
    exclusions.to_csv(output / "diagnostics_gene_exclusions.csv", index=False)
    (output / "diagnostics_population.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
