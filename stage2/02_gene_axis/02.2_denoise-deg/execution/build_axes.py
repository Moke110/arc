"""Publish the five nested Denoise-DEG axes from the recorded selection evidence."""

import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
import pandas as pd

from diagnose_scores import disease_evidence

SIZES = (1000, 2000, 3000, 4000, 5000)
RULES = ("low_detection", "plate_sd_above_boundary", "within_plate_sd_above_boundary")


def sha256(path):
    with path.open("rb") as handle:
        return file_digest(handle, "sha256").hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def axis_score_histograms(ranking):
    for metric in ("p_tenv", "p_disease", "score"):
        values = ranking[metric]
        width = max(1, int(np.ceil((values.max() + 1) / 50)))
        edges = np.arange(-0.5, values.max() + width + 0.5, width)
        total, _ = np.histogram(values, bins=edges)
        for size in SIZES:
            selected, _ = np.histogram(values.iloc[:size], bins=edges)
            for left, right, count, in_axis in zip(edges[:-1], edges[1:], total, selected, strict=True):
                yield {"axis_size": size, "metric": metric, "bin_left": left, "bin_right": right,
                       "retained_genes": count, "axis_genes": in_axis}


def select(repo, step):
    data = step / "figure/data"
    corpus_manifest = repo / "stage2/01_data_raw/output/manifest.json"
    disease_root = repo / "stage2/02_gene_axis/02.1_hvg-deg/output/release"
    disease_manifest = disease_root / "manifest.json"
    previous = read_json(disease_manifest)
    if sha256(corpus_manifest) != read_json(data / "diagnostics_population.json")["corpus_manifest_sha256"]:
        raise ValueError("DMSO and abundance evidence came from a different raw corpus")
    if read_json(corpus_manifest)["common_axis"] != previous["common_axis"]:
        raise ValueError("raw corpus and disease evidence use different common gene axes")
    for filename, entry in previous["outputs"].items():
        if filename.startswith("evidence/disease-fusion/") and sha256(disease_root / filename) != entry["sha256"]:
            raise ValueError(f"02.1 disease evidence checksum differs: {filename}")
    disease, disease_pins = disease_evidence(repo)
    if disease_pins["disease_catalog_sha256"] != read_json(data / "diagnostics_score_summary.json")["disease_catalog_sha256"]:
        raise ValueError("score evidence came from a different disease release")
    if (data / "diagnostics_disease_top200.csv").read_text() != disease.to_csv(index=False):
        raise ValueError("DEG200 evidence differs from released disease ranks")
    if disease.duplicated(["disease", "gene_index"]).any() or not disease.groupby("disease").size().eq(200).all():
        raise ValueError("each eligible disease must contain exactly 200 distinct genes")

    identities = ["gene_index", "symbol"]
    genes_path = repo / "stage2/util/output/gene-identity/genes.csv"
    genes = pd.read_csv(genes_path, usecols=identities)
    flags = pd.read_csv(data / "diagnostics_gene_exclusions.csv")
    scores = pd.read_csv(data / "diagnostics_scores.csv")
    measured = pd.read_csv(data / "diagnostics_genes.csv.gz")
    table = genes.merge(flags, on=identities, validate="one_to_one").merge(
        scores, on=identities, validate="one_to_one"
    ).merge(measured, on=identities, validate="one_to_one")
    if len(table) != previous["common_axis"]["size"] or any(
        len(frame) != len(table) for frame in (genes, flags, scores, measured)
    ):
        raise ValueError("selection evidence does not exactly cover the common gene axis")
    if not table.low_detection.equals(table.detection_fraction.lt(.10)):
        raise ValueError("nonzero-well exclusion flags differ from the 0.10 threshold")
    x = table.dmso_mean_log2_cpm1
    if not table.dmso_measured.equals(x.gt(0)):
        raise ValueError("DMSO calibration population differs from measured expression")
    for metric in ("plate_sd", "within_plate_sd"):
        measured_ratio = table.loc[table.dmso_measured, metric] / x[table.dmso_measured]
        k = measured_ratio.quantile(.90)
        if not table[metric + "_above_boundary"].equals(x.gt(0) & (table[metric] > k * x)):
            raise ValueError(f"{metric} flags differ from the calibrated threshold")
    if not table[list(RULES)].any(axis=1).equals(table.excluded):
        raise ValueError("exclusion flags do not match the union of screening rules")
    p_disease = table.gene_index.map(disease.gene_index.value_counts()).fillna(0).astype(int)
    if not p_disease.equals(table.p_disease) or not table.score.eq(table.p_tenv * (2 ** p_disease)).all():
        raise ValueError("scores differ from the selected disease evidence or scoring formula")

    columns = [*identities, "p_tenv", "p_disease", "score"]
    ranking = table.loc[~table.excluded, columns].sort_values(
        ["score", "gene_index"], ascending=[False, True], ignore_index=True
    )
    if ranking.score.gt(0).sum() < max(SIZES):
        raise ValueError("fewer than 5,000 retained genes have positive scores")
    ranking.insert(0, "rank", range(1, len(ranking) + 1))
    output = step / "output/release"
    output.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(output / "gene_ranking.csv", index=False)
    axes = {}
    coverage = []
    for size in SIZES:
        selected = ranking.iloc[:size]
        axis_path = output / f"gene_axis_{size}.csv"
        selected[identities].sort_values("gene_index").to_csv(axis_path, index=False)
        axes[str(size)] = sha256(axis_path)
        for name, top in disease.groupby("disease", sort=False):
            hits = int(top.gene_index.isin(selected.gene_index).sum())
            coverage.append({"disease": name, "axis_size": size, "hits": hits,
                             "total": len(top), "coverage": hits / len(top)})
    coverage_path = data / "disease-DEG200-coverage.csv"
    pd.DataFrame(coverage).to_csv(coverage_path, index=False)
    histogram_path = data / "score-in-axis.csv"
    pd.DataFrame(axis_score_histograms(ranking)).to_csv(histogram_path, index=False)
    evidence_names = (
        "diagnostics_genes.csv.gz", "diagnostics_gene_exclusions.csv",
        "diagnostics_scores.csv", "diagnostics_disease_top200.csv",
        "diagnostics_population.json", "diagnostics_score_summary.json",
        "diagnostics_variance_boundaries.csv", coverage_path.name, histogram_path.name,
    )
    manifest = {
        "schema": "stage2-denoise-deg-axes-v1",
        "common_axis": previous["common_axis"],
        "inputs": {str(path.relative_to(repo)): sha256(path)
                   for path in (corpus_manifest, disease_manifest, genes_path)},
        "screening": {"nonzero_fraction_min": .10, "dmso_variance_upper_fraction": .10},
        "score": "p_tenv * 2 ** p_disease",
        "tie_break": "ascending frozen gene_index",
        "axis_order": "ascending frozen gene_index",
        "axes": axes,
        "ranked_genes": len(ranking),
        "ranking_sha256": sha256(output / "gene_ranking.csv"),
        "evidence": {f"figure/data/{name}": sha256(data / name) for name in evidence_names},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--step", type=Path, required=True)
    args = parser.parse_args()
    select(args.repo.resolve(), args.step.resolve())
