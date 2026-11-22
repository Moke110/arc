#!/usr/bin/env python3
"""Build and independently verify the source-native Stage 2 corpus."""
from __future__ import annotations

import argparse
import fcntl
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem
import pyarrow.parquet as pq

STEP = Path(__file__).resolve().parents[1]
REPO = next(parent for parent in STEP.parents if (parent / ".git").exists())
sys.path[:0] = [str(path) for path in (
    STEP / "output/api", REPO / "stage2/01_bases/output/api", REPO / "stage2/util/output/api",
)]

from bases import corpus as legacy_corpus, native, raw, raw_counts, release
from cardiac_sources import build_cardiac_sources
from stage2_util.publication import publication_lock, publish_candidate
import native_tahoe as tahoe
from native_corpus import CorpusError, open_corpus, seal_corpus, sha256

GSE_TENV = "GSE262419"
TAHOE = "Tahoe-100M"
EXCLUSIONS = {
    "XPert-L1000": "Landmark-gene expression profiles lack verified native raw integer counts.",
    "Annotation": "Reference metadata only, with no measured perturbation profiles.",
}


def pin_source(item, path):
    """Verify one consumed file against the external source item's immutable manifest."""
    relative = path.relative_to(item / "source").as_posix()
    for entry in (item / "manifest.sha256").read_text().splitlines():
        expected, filename = entry.split(maxsplit=1)
        if filename.removeprefix("*").removeprefix("./") == relative:
            if sha256(path) != expected:
                raise CorpusError(f"source digest differs: {relative}")
            return expected
    raise CorpusError(f"unsealed source file: {relative}")


def write_metadata(payload, source, item, drug_metadata=None):
    with h5py.File(payload) as handle:
        raw_obs = read_elem(handle["obs"])
    if source == "tenv":
        records, excluded = legacy_corpus._tenv_rows(raw_obs, item)
        role = np.select([raw_obs.is_blank, raw_obs.is_control], ["blank", "vehicle"], "treated")
        compatible = ~raw_obs.is_blank & (raw_obs.is_control | raw_obs.vehicle.eq("DMSO"))
        groups = ("tenv::" + raw_obs.plate.astype(str) + "::" + raw_obs.vehicle.astype(str)).where(
            compatible, ""
        )
    else:
        structures = legacy_corpus.tahoe_structures(drug_metadata)
        records, excluded = legacy_corpus._tahoe_rows(raw_obs, item, structures)
        role = np.where(raw_obs.is_control, "vehicle", "treated")
        groups = "tahoe::" + raw_obs.plate.astype(str) + "::" + raw_obs.cell_line_id.astype(str)
    canonical = pd.DataFrame(records).set_index("profile_id", drop=False)
    obs = canonical.reindex(raw_obs.index)
    obs["profile_id"] = raw_obs.index
    obs["role"] = role
    obs["matched_control_group"] = groups
    obs["source_item"] = item["name"]
    obs["source_item_manifest_sha256"] = item["sha256"]
    obs["source_row"] = np.arange(len(obs))
    if source == "tenv":
        fallback_compound = raw_obs.compound_name
        fallback_cell = "hiPSC-CM"
    else:
        fallback_compound = raw_obs.drug
        fallback_cell = raw_obs.cell_name
    obs["compound_name"] = obs.compound_name.where(obs.compound_name.notna(), fallback_compound)
    obs["cell_line"] = obs.cell_line.fillna(fallback_cell)
    for field in ("plate", "dose_um"):
        obs[field] = raw_obs[field]
    obs["view_exclusion_reason"] = ""
    for row in excluded:
        obs.loc[row["profile_id"], "view_exclusion_reason"] = row["reason"]
    for field in obs.select_dtypes(include="object"):
        obs[field] = obs[field].fillna("")
    with h5py.File(payload.with_suffix(".obs.h5"), "w") as handle:
        write_elem(handle, "obs", obs)


def axis_inputs():
    stem = REPO / "stage2/util/output/gene-identity/tables"
    axis_path = stem / "genes.csv"
    axis = raw_counts._axis_frame(axis_path)
    crosswalk = pd.read_csv(stem / "gene_crosswalk.csv", dtype=str, keep_default_na=False)
    return axis, crosswalk


def tahoe_mappings(data, axis, crosswalk):
    gene_file = data / TAHOE / "source/metadata/gene_metadata.parquet"
    pin_source(data / TAHOE, gene_file)
    common = raw.tahoe_axis_mapping(axis, crosswalk, gene_file)
    genes = pq.read_table(gene_file).to_pandas().sort_values("token_id").reset_index(drop=True)
    mapping = pd.DataFrame({
        "gene_index": np.arange(len(genes)),
        "symbol": genes.gene_symbol.fillna("").astype(str),
        "native_ensembl_id": genes.ensembl_id,
        "native_token_id": genes.token_id,
    })
    positions = dict(zip(genes.token_id, mapping.gene_index, strict=True))
    indices = {int(row.gene_index): int(positions[int(row.native_token_id)])
               for row in common.itertuples(index=False)}
    return mapping, indices


def tenv_args(data):
    return (
        data / GSE_TENV / "source",
        data / "Annotation/source/hgnc_complete_set.txt",
        data / "Annotation/local/cm_compound_annotation.csv",
    )


def verify_product(root):
    with open_corpus(root) as corpus:
        release.check_manifest(root, corpus.manifest)
        for name in corpus.manifest["sources"]:
            data = corpus.source(name)
            if name in {"tenv", "tahoe_bulk"} and not data.shared_axis_complete:
                raise CorpusError(f"{name}: common axis coverage is incomplete")
        return corpus.manifest_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--materialized-root", type=Path, default=REPO)
    parser.add_argument("--output", type=Path, default=STEP / "output/native")
    parser.add_argument("--tahoe-built", type=Path,
                        help="reuse a native Tahoe build, independently replayed before use")
    parser.add_argument("--staging", type=Path, default=STEP / "execution/checkpoints")
    args = parser.parse_args()
    source_root = args.materialized_root.resolve()
    data = source_root / "pert_data"
    axis, crosswalk = axis_inputs()
    mapping, tahoe_columns = tahoe_mappings(data, axis, crosswalk)
    drug_metadata = data / TAHOE / "source/metadata/drug_metadata.parquet"
    pin_source(data / TAHOE, drug_metadata)

    def build(candidate):
        inputs = {
            "common_axis_sha256": raw_counts.AXIS_SHA256,
            "crosswalk_sha256": sha256(
                REPO / "stage2/util/output/gene-identity/tables/gene_crosswalk.csv"
            ),
            "source_manifests": {
                item.name: sha256(item / "manifest.sha256")
                for item in sorted(data.iterdir()) if (item / "manifest.sha256").is_file()
            },
        }

        pin_source(data / "Annotation", tenv_args(data)[1])
        tenv_record = native.build_tenv(*tenv_args(data), candidate)
        for filename, record in tenv_record["inputs"].items():
            if filename in {path.name for path in tenv_args(data)[1:]}:
                continue
            if pin_source(data / GSE_TENV, data / GSE_TENV / "source" / filename) != record["sha256"]:
                raise CorpusError(f"tenv source replay differs: {filename}")
        native.verify_tenv(*tenv_args(data), candidate)
        with h5py.File(candidate / "tenv.h5ad") as handle:
            var = read_elem(handle["var"])
        tenv_bridge = raw_counts.tenv_axis_mapping(axis, var)
        tenv_columns = {int(row.gene_index): int(row.native_column)
                        for row in tenv_bridge.itertuples(index=False)}
        inputs["tenv_source_replay"] = tenv_record

        tahoe_source = data / TAHOE / "source"
        tahoe_manifest = data / TAHOE / "manifest.sha256"
        options = {"expected_genes": len(mapping)}
        tahoe_candidate = args.tahoe_built
        if tahoe_candidate is None:
            tahoe_candidate = candidate / "tahoe-build"
            tahoe.build_tahoe(
                tahoe_source, mapping, tahoe_manifest, tahoe_candidate,
                checkpoints=args.staging / "tahoe-native-plates", **options,
            )
        inputs["tahoe_source_replay"] = tahoe.verify_tahoe(
            tahoe_source, mapping, tahoe_manifest, tahoe_candidate, **options
        )
        payload = candidate / "tahoe_bulk.h5ad"
        if args.tahoe_built is None:
            (tahoe_candidate / "tahoe.h5ad").rename(payload)
            shutil.rmtree(tahoe_candidate)
        else:
            shutil.copyfile(tahoe_candidate / "tahoe.h5ad", payload)
        if sha256(payload) != inputs["tahoe_source_replay"]["h5ad_sha256"]:
            raise CorpusError("Tahoe corpus payload differs from independently replayed build")

        for source, item in (("tenv", GSE_TENV), ("tahoe_bulk", TAHOE)):
            write_metadata(candidate / f"{source}.h5ad", source,
                           {"name": item, "sha256": inputs["source_manifests"][item]},
                           drug_metadata)
        columns = build_cardiac_sources(candidate, data, axis, pin_source, inputs)
        columns.update(tenv=tenv_columns, tahoe_bulk=tahoe_columns)
        seal_corpus(candidate, columns, len(axis), raw_counts.AXIS_SHA256, inputs, EXCLUSIONS,
                    artifacts=("tenv_probe_mapping.csv", "tenv_manifest.json"))

    args.staging.mkdir(parents=True, exist_ok=True)
    with publication_lock(args.output) as output, (args.staging / "native-corpus.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.action == "publish":
            publish_candidate(release, output, build, verify_product, staging=args.staging)
        else:
            release.replay(output, build, verify_product, root=args.staging)
        print(verify_product(output))


if __name__ == "__main__":
    main()
