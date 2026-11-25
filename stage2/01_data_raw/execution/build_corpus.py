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
    STEP / "output/api", REPO / "stage2/util/output/api",
)]

from corpus_sources import identities, mapping as gene_mapping, native
from stage2_util.derivation import release
from cardiac_sources import PRIMARY_MEMBERSHIP, build_cardiac_sources
from corpus_migration import predecessor_record, stage_predecessor, verify_equivalence
from stage2_util.publication import publication_lock, publish_candidate
import native_tahoe as tahoe
from native_corpus import DATASETS, SCHEMA, CorpusError, open_corpus, seal_corpus, sha256

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
        records, excluded = identities._tenv_rows(raw_obs, item)
        role = np.select([raw_obs.is_blank, raw_obs.is_control], ["blank", "vehicle"], "treated")
        compatible = ~raw_obs.is_blank & (raw_obs.is_control | raw_obs.vehicle.eq("DMSO"))
        groups = ("tenv::" + raw_obs.plate.astype(str) + "::" + raw_obs.vehicle.astype(str)).where(
            compatible, ""
        )
    else:
        structures = identities.tahoe_structures(drug_metadata)
        records, excluded = identities._tahoe_rows(raw_obs, item, structures)
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
    axis = gene_mapping.axis_frame(axis_path)
    crosswalk = pd.read_csv(stem / "gene_crosswalk.csv", dtype=str, keep_default_na=False)
    return axis, crosswalk


def tahoe_mappings(data, axis, crosswalk):
    gene_file = data / TAHOE / "source/metadata/gene_metadata.parquet"
    pin_source(data / TAHOE, gene_file)
    common = gene_mapping.tahoe_axis_mapping(axis, crosswalk, gene_file)
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


SUPPLEMENTS = {
    "tenv": ("tenv_probe_mapping.csv", "tenv_manifest.json"),
    "dtoxs_primary": (PRIMARY_MEMBERSHIP,),
}
DEFAULT = STEP / "output"
HISTORY = STEP.parent / "history"
AXES = Path("02_gene_axis/02.1_hvg-deg")
TRAINING = Path("03_data_train")
TAHOE_VIEW = TRAINING / "03.2_tahoe_hvg-5k"


def _seal(candidate, columns, common_axis, inputs, exclusions, predecessor):
    seal_corpus(candidate, columns, common_axis["size"], common_axis["sha256"], inputs,
                exclusions, supplements=SUPPLEMENTS,
                predecessor=predecessor_record(predecessor) if predecessor is not None else None)


def _preserve_api(candidate, output):
    if (output / "api").is_dir():
        shutil.copytree(output / "api", candidate / "api", ignore=shutil.ignore_patterns("__pycache__"))


def _verify_primary_membership(root, source):
    ledger = pd.read_csv(root / "dtoxs_primary" / PRIMARY_MEMBERSHIP, dtype=str,
                         keep_default_na=False)
    if (len(source.obs) != 842 or len(ledger) != 956 or ledger.source_table.nunique() != 8
            or ledger[["source_table", "source_row"]].duplicated().any()
            or set(ledger.profile_id) != set(source.obs.index)):
        raise CorpusError("DToxS-primary source membership or observation count differs")


def migrate(predecessor, candidate):
    """Build a candidate from the independently accepted flat native corpus."""
    release.check_manifest(predecessor.root, predecessor.manifest)
    axis, _ = axis_inputs()
    columns = stage_predecessor(predecessor, candidate, axis)
    _seal(candidate, columns, predecessor.manifest["common_axis"],
          predecessor.manifest["inputs"], predecessor.manifest["exclusions"], predecessor)


def verify_product(root):
    with open_corpus(root) as corpus:
        if corpus.manifest["schema"] != SCHEMA or set(corpus.manifest["sources"]) != DATASETS:
            raise CorpusError("unexpected corpus dataset inventory")
        expected = {
            f"{name}/{filename}" for name in DATASETS
            for filename in (f"{name}.h5ad", "manifest.json", *SUPPLEMENTS.get(name, ()))
        }
        if set(corpus.manifest["outputs"]) != expected:
            raise CorpusError("unexpected corpus file inventory")
        files = {path.relative_to(root).as_posix() for path in root.glob("*/*")
                 if path.parent.name != "api"}
        if files != expected or {path.name for path in root.iterdir()} - {"api", "manifest.json"} != DATASETS:
            raise CorpusError("corpus payload directory differs from release inventory")
        release.check_manifest(root, corpus.manifest)
        for name in corpus.manifest["sources"]:
            dataset = corpus.source(name)
            if name in {"tenv", "tahoe_bulk"} and not dataset.shared_axis_complete:
                raise CorpusError(f"{name}: common axis coverage is incomplete")
            if name == "dtoxs_primary":
                _verify_primary_membership(root, dataset)
        return corpus.manifest_sha256


def _select_predecessor(output, explicit):
    active = output / "native"
    if (active / "manifest.json").is_file():
        if explicit is not None and sha256(explicit / "manifest.json") != sha256(active / "manifest.json"):
            raise CorpusError("selected predecessor differs from the installed legacy corpus")
        return active
    return explicit


def verify_archived_legacy(legacy):
    """Require a verified local archive handoff before replacing a legacy directory."""
    inventory = HISTORY / "archive-manifest.json"
    error = "replacing the installed legacy corpus requires a verified archive handoff"
    if not legacy.is_relative_to(REPO) or not inventory.is_file() or inventory.is_symlink():
        raise CorpusError(error)
    original = legacy.relative_to(REPO)
    handoff = release.read_json(inventory)
    destination = Path(handoff.get("path_map", {}).get(original.as_posix(), "."))
    archive = REPO / destination
    if (handoff.get("schema") != "stage2-local-history-v1" or not archive.is_relative_to(HISTORY)
            or not archive.resolve().is_relative_to(HISTORY.resolve())
            or not (archive / "manifest.json").is_file() or not (legacy / "manifest.json").is_file()
            or sha256(archive / "manifest.json") != sha256(legacy / "manifest.json")):
        raise CorpusError(error)
    record = release.read_json(legacy / "manifest.json")
    declared = set(release.output_digests(record)) | {"manifest.json"}
    if {path.relative_to(legacy).as_posix() for path in legacy.rglob("*")
            if path.is_file() or path.is_symlink()} != declared:
        raise CorpusError(f"{error}: unexpected legacy files")
    entries = handoff.get("entries", {})
    for name in declared:
        entry = entries.get((original / name).as_posix(), {})
        if (entry.get("destination") != (destination / name).as_posix()
                or (archive / name).is_symlink()):
            raise CorpusError(f"{error}: missing archived entry {name}")
    release.check_manifest(archive, record)


def verify_downstream(candidate, output):
    """Check the active training releases before replacing the accepted corpus."""
    if not ((output / "manifest.json").is_file() or (output / "native/manifest.json").is_file()):
        return
    stage2 = output.parent.parent
    sys.path[:0] = [str(STEP.parent / directory / "output/api")
                    for directory in (AXES, TRAINING, TAHOE_VIEW)]
    from prediction_axis import open_axes
    from tahoe_expression import open_tahoe
    from training_data import TENV_AXIS_STEPS, open_tenv_axis_training

    with open_corpus(candidate) as corpus:
        for axis_id in TENV_AXIS_STEPS:
            with open_tenv_axis_training(axis_id, stage2 / TRAINING) as training:
                if training.manifest["inputs"]["raw_sha256"] != corpus.manifest_sha256:
                    raise CorpusError(f"{axis_id}: accepted tenv training pins a different native corpus")
        axis = open_axes(stage2 / AXES / "output/release").axis(5000)
        projection = axis.project(corpus, "tahoe_bulk")
        with open_tahoe(stage2 / TAHOE_VIEW / "output/release") as view:
            if tuple(view.var.gene_index) != axis.gene_indices or not view.shape[0]:
                raise CorpusError("Tahoe: accepted well-level axis or population differs")
            for index in (0, view.shape[0] - 1):
                identity = view.obs.index[index]
                native_row = projection.source.obs.index.get_loc(identity)
                counts = projection.counts(native_row, native_row + 1)
                totals = counts.sum(axis=1, dtype=np.int64)
                if (totals <= 0).any():
                    raise CorpusError(f"Tahoe: empty selected-axis library for {identity}")
                expected = np.log1p(counts.astype(np.float64) * (10000.0 / totals[:, None])).astype(np.float32)
                expression, control, delta = (view.read(index, index + 1, layer)
                                               for layer in ("X", "control", "delta"))
                if (not np.array_equal(expression, expected) or not np.isfinite(control).all()
                        or not np.array_equal(delta, expression - control)):
                    raise CorpusError(f"Tahoe: accepted well-level view differs from candidate at {identity}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify", "migrate"))
    parser.add_argument("--materialized-root", type=Path, default=REPO)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--predecessor", type=Path,
                        help="sealed v1 release used for migration or equivalence checks")
    parser.add_argument("--tahoe-built", type=Path,
                        help="reuse a native Tahoe build, independently replayed before use")
    parser.add_argument("--staging", type=Path, default=STEP / "execution/checkpoints")
    args = parser.parse_args()
    def build(candidate, predecessor):
        data = args.materialized_root.resolve() / "pert_data"
        axis, crosswalk = axis_inputs()
        tahoe_mapping, tahoe_columns = tahoe_mappings(data, axis, crosswalk)
        drug_metadata = data / TAHOE / "source/metadata/drug_metadata.parquet"
        pin_source(data / TAHOE, drug_metadata)
        inputs = {
            "common_axis_sha256": gene_mapping.AXIS_SHA256,
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
        tenv_bridge = gene_mapping.tenv_axis_mapping(axis, var)
        tenv_columns = {int(row.gene_index): int(row.native_column)
                        for row in tenv_bridge.itertuples(index=False)}
        inputs["tenv_source_replay"] = tenv_record

        tahoe_source = data / TAHOE / "source"
        tahoe_manifest = data / TAHOE / "manifest.sha256"
        options = {"expected_genes": len(tahoe_mapping)}
        tahoe_candidate = args.tahoe_built
        if tahoe_candidate is None:
            tahoe_candidate = candidate / "tahoe-build"
            tahoe.build_tahoe(
                tahoe_source, tahoe_mapping, tahoe_manifest, tahoe_candidate,
                checkpoints=args.staging / "tahoe-native-plates", **options,
            )
        inputs["tahoe_source_replay"] = tahoe.verify_tahoe(
            tahoe_source, tahoe_mapping, tahoe_manifest, tahoe_candidate, **options
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
        if predecessor is not None:
            inputs["publication_method"] = "external-source-replay"
        _seal(candidate, columns, {"size": len(axis), "sha256": gene_mapping.AXIS_SHA256},
              inputs, EXCLUSIONS, predecessor)

    with publication_lock(args.output) as output:
        stable = output == DEFAULT.resolve()
        if stable and (output / "manifest.json").is_file() and args.action != "verify":
            if (args.action != "publish" or args.predecessor is not None or args.tahoe_built is not None
                    or args.materialized_root.resolve() != REPO):
                raise CorpusError("the installed corpus generation is pinned; build with a scratch output")
            installed = verify_product(output)
            verify_downstream(output, output)
            print(installed)
            return

        previous = _select_predecessor(output, args.predecessor)
        legacy = output / "native"
        if args.action != "verify" and legacy.exists():
            verify_archived_legacy(legacy)
        if args.action == "migrate" and previous is None:
            parser.error("migration requires a sealed v1 release")

        args.staging.mkdir(parents=True, exist_ok=True)
        with (args.staging / "native-corpus.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if previous is None:
                verified = _run(args, output, None, build)
            else:
                with open_corpus(previous) as predecessor:
                    verified = _run(args, output, predecessor, build)
        print(verified)


def _run(args, output, predecessor, build):
    migrated = args.action == "migrate"
    if args.action == "verify":
        accepted = release.read_json(output / "manifest.json")
        migrated = ("predecessor" in accepted
                    and accepted["inputs"].get("publication_method") != "external-source-replay")
        if migrated and predecessor is None:
            return verify_product(output)

    def prepare(candidate):
        if migrated:
            migrate(predecessor, candidate)
        else:
            build(candidate, predecessor)
        _preserve_api(candidate, output)

    def verify(candidate):
        verify_product(candidate)
        if predecessor is not None:
            verify_equivalence(predecessor, candidate)

    if args.action == "verify":
        return release.replay(output, prepare, verify, root=args.staging)["manifest_sha256"]

    def verify_before_install(candidate):
        verify(candidate)
        if output == DEFAULT.resolve():
            verify_downstream(candidate, output)

    installed = publish_candidate(release, output, prepare, verify_before_install, staging=args.staging)
    return sha256(installed / "manifest.json")


if __name__ == "__main__":
    main()
