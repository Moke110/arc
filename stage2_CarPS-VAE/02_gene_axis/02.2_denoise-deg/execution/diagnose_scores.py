"""Measure preliminary score distributions before applying gene exclusions."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd


def disease_evidence(repo):
    root = repo / "stage2_CarPS-VAE/02_gene_axis/02.1_hvg-deg/output/release/evidence/disease-fusion"
    catalog = json.loads((root / "catalog.json").read_text())
    tables = []
    for row in catalog["diseases"]:
        if row["eligible"]:
            frame = pd.read_csv(root / row["deg1000_file"])
            if len(frame) != row["deg1000"] or frame.gene_index.duplicated().any():
                raise ValueError(f"invalid published DEG1000 list for {row['category']}")
            tables.append(frame.assign(disease=row["category"]))
    return pd.concat(tables, ignore_index=True), {
        "disease_catalog_sha256": hashlib.sha256((root / "catalog.json").read_bytes()).hexdigest(),
    }


def plate_evidence(source, plate, positions, columns, command, tmp_path):
    positions = np.sort(positions)
    obs = source.obs.iloc[positions]
    if not np.all(np.diff(positions) == 1):
        raise ValueError(f"{plate} source rows are not contiguous")
    depth, matrix = [], []
    for start in range(positions[0], positions[-1] + 1, 128):
        counts = source.counts(start, min(start + 128, positions[-1] + 1))
        depth.extend(counts.sum(axis=1))
        matrix.append(counts[:, columns])
    depth = np.asarray(depth)
    dmso = (obs.role == "vehicle") & (obs.vehicle == "DMSO")
    treated = (obs.role == "treated") & (obs.match_status == "matched")
    retained = (dmso | treated).to_numpy() & (depth > 0)
    dmso_rows = positions[dmso & retained]
    for index in positions[treated & retained]:
        if not np.isin(dmso_rows, source.matched_control_rows(index)).all():
            raise ValueError(f"unexpected matched controls on {plate}")
    if len(dmso_rows) < 2:
        raise ValueError(f"insufficient DMSO controls on {plate}")

    condition = obs[["drug_parent_inchikey", "dose_um"]].copy()
    condition["condition"] = (
        condition.drug_parent_inchikey + "@" + condition.dose_um.astype(str)
    )
    metadata = condition.assign(
        group=condition.condition.mask(dmso, "DMSO"), depth=depth
    ).loc[retained]
    counts = np.concatenate(matrix)[retained]
    np.asarray(counts, dtype="<i4").tofile(tmp_path / "plate.bin")
    metadata.to_csv(tmp_path / "metadata.csv", index=False)
    run = subprocess.run(
        [*command, "plate.bin", "metadata.csv", str(len(columns)), "deg.csv"],
        cwd=tmp_path, capture_output=True, text=True, check=True,
    )
    if run.stderr:
        print(run.stderr, file=sys.stderr)
    frame = pd.read_csv(tmp_path / "deg.csv").rename(columns={"group": "condition"})
    return frame.assign(plate=plate), metadata.loc[metadata.group != "DMSO", "group"].unique()


def score_genes(genes, hits, disease):
    scores = genes[["gene_index", "symbol"]].copy()
    for key, votes in (("p_tenv", hits), ("p_disease", disease)):
        scores[key] = scores.gene_index.map(votes.gene_index.value_counts()).fillna(0).astype(int)
        assert scores[key].sum() == len(votes)
    scores["score"] = scores.p_tenv * (2 ** scores.p_disease)
    return scores


def score_summary(scores):
    return {key: int((scores[key] > 0).sum()) for key in ("p_tenv", "p_disease", "score")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--step", type=Path, required=True)
    args = parser.parse_args()
    repo, step = args.repo.resolve(), args.step.resolve()
    sys.path[:0] = [str(repo / "stage2_CarPS-VAE/01_data_raw/output/api"),
                    str(repo / "stage2_CarPS-VAE/util/output/api")]
    from native_corpus import open_corpus
    from stage2_util.derivation.deg_rank import deg_effect
    from stage2_util.derivation.r_env import entry_point

    r, home = entry_point()
    os.environ["R_HOME"] = str(home)
    command = [str(r), "--vanilla", "--slave", "-f",
               str(step / "execution/plate_deg.R"), "--args"]
    genes = pd.read_csv(repo / "stage2_CarPS-VAE/util/output/gene-identity/genes.csv")
    disease, summary = disease_evidence(repo)
    all_frames, plate_counts = [], []
    with TemporaryDirectory(prefix="plate-deg-", dir=step) as scratch:
        with open_corpus(repo / "stage2_CarPS-VAE/01_data_raw/output") as corpus:
            source = corpus.source("tenv")
            columns = source.native_columns(genes.gene_index)
            for plate, positions in source.obs.groupby("plate").indices.items():
                frame, conditions = plate_evidence(
                    source, plate, positions, columns, command, Path(scratch)
                )
                all_frames.append(frame)
                plate_counts.append(pd.DataFrame({"plate": plate, "condition": conditions}))
                print(f"{plate}: {len(conditions)} contrasts, {len(frame)} qualifying rows", flush=True)
    population = pd.concat(plate_counts, ignore_index=True)
    hits = pd.concat(all_frames, ignore_index=True)
    hits["effect"] = deg_effect(hits.log2fc, hits.padj)
    hits["direction"] = np.where(hits.log2fc > 0, "up", "down")
    hits = hits.sort_values(
        ["condition", "effect", "gene_index", "plate"], ascending=[True, False, True, True]
    ).drop_duplicates(["condition", "gene_index"])
    hits = hits.groupby(["condition", "direction"], sort=False).head(100)
    scores = score_genes(genes, hits, disease)

    output = step / "figure/data"
    scores.to_csv(output / "diagnostics_scores.csv", index=False)
    hits.to_csv(output / "diagnostics_tenv_top_deg.csv.gz", index=False)
    disease.to_csv(output / "diagnostics_disease_deg1000.csv", index=False)
    summary.update(
        plate_contrasts=len(population),
        unique_conditions=population.condition.nunique(),
        conditions_without_hits=population.condition.nunique() - hits.condition.nunique(),
        tenv_votes=len(hits),
        disease_votes=len(disease),
        diseases=disease.groupby("disease").size().to_dict(),
        positive_scores=score_summary(scores),
    )
    (output / "diagnostics_score_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
