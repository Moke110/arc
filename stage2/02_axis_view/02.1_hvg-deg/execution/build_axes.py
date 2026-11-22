"""Publish and replay the four accepted HVG-DEG axes without expression data."""
import argparse
import json
import shutil
import sys
from dataclasses import replace
from pathlib import Path

STEP = Path(__file__).resolve().parents[1]
STAGE2 = STEP.parents[1]
sys.path.insert(0, str(STEP / "output/api"))
sys.path.insert(0, str(STAGE2 / "01_bases/execution"))
sys.path.insert(0, str(STAGE2 / "01_corpus/output/api"))

import build_gene_axis as legacy
from native_corpus import open_corpus
from prediction_axis import COMMON_SHA256, COMMON_SIZE, SCHEMA, open_axes, sha256

OUTPUT = STEP / "output/release"
METHOD_FILES = (
    "stage2/01_bases/execution/build_gene_axis.py",
    "stage2/01_bases/output/api/bases/gene_axis.py",
)


def _copy(source, target, expected):
    if sha256(source) != expected:
        raise ValueError(f"HVG-DEG evidence digest differs: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def _candidate(source_root, output):
    config = legacy.default_config(source_root)
    accepted_path = config.output / "manifest.json"
    accepted = legacy.release.read_json(accepted_path)
    if accepted["schema"] != legacy.SCHEMA or accepted["release"]["algorithm"] != legacy.ALGORITHM:
        raise ValueError("unrecognized historical axis method")
    for path, record in accepted["inputs"].items():
        if sha256(source_root / path) != record["sha256"]:
            raise ValueError(f"pinned historical axis input differs: {path}")
    common_table = config.genes.parent / "tables/genes.csv"
    if sha256(common_table) != COMMON_SHA256:
        raise ValueError("fixed Stage 2 common gene identity differs")
    products = {
        "evidence/gene-axis/accepted_manifest.json": (accepted_path, sha256(accepted_path)),
        "evidence/tenv-hvg/manifest.json": (
            config.tenv_hvg.with_name("manifest.json"), sha256(config.tenv_hvg.with_name("manifest.json"))
        ),
    }
    axes = {}
    for name, record in accepted["outputs"].items():
        source = source_root / name
        if source.parent != config.output:
            continue
        if source.name.startswith("gene_axis_"):
            size = int(source.stem.rsplit("_", 1)[1])
            products[source.name] = (source, record["sha256"])
            axes[str(size)] = record["sha256"]
        else:
            products[f"evidence/gene-axis/{source.name}"] = (source, record["sha256"])
    for name, record in accepted["inputs"].items():
        source = source_root / name
        if source == config.tenv_hvg:
            products[f"evidence/tenv-hvg/{source.name}"] = (source, record["sha256"])
        elif source.parent == config.fusion_root:
            products[f"evidence/disease-fusion/{source.name}"] = (source, record["sha256"])
    for target, (source, digest) in products.items():
        _copy(source, output / target, digest)
    manifest = {
        "schema": SCHEMA,
        "common_axis": {"size": COMMON_SIZE, "sha256": COMMON_SHA256},
        "axes": axes,
        "algorithm": legacy.ALGORITHM,
        "method": {name: sha256(STAGE2.parent / name) for name in METHOD_FILES},
        "outputs": {name: {"sha256": digest} for name, (_, digest) in products.items()},
    }
    legacy.release.write_json(output / "manifest.json", manifest)
    return verify_release(output)


def verify_release(output):
    published = open_axes(output)
    manifest = published.manifest
    if manifest["algorithm"] != legacy.ALGORITHM or manifest["method"] != {
        name: sha256(STAGE2.parent / name) for name in METHOD_FILES
    }:
        raise ValueError("published HVG-DEG selection method differs")
    legacy.release.check_manifest(output, manifest)
    evidence = output / "evidence"
    config = replace(legacy.default_config(STAGE2.parent),
                     tenv_hvg=evidence / "tenv-hvg/descriptive.csv",
                     fusion_root=evidence / "disease-fusion")
    if sha256(config.genes.parent / "tables/genes.csv") != COMMON_SHA256:
        raise ValueError("fixed Stage 2 common gene identity differs")
    axis = legacy.read_axis(config.genes)
    if len(axis) != COMMON_SIZE:
        raise ValueError("selection input does not cover the fixed common genes")
    tenv_hvg = legacy.read_tenv(config.tenv_hvg, COMMON_SIZE)
    fusion = legacy.read_fusion(config, COMMON_SIZE)
    selected = legacy.compute(axis, tenv_hvg, fusion, config.axis_sizes)["axes"]
    prior = set()
    for size in published.sizes:
        prediction = published.axis(size)
        expected = tuple(selected[size].tolist())
        if prediction.gene_indices != expected or not prior.issubset(expected):
            raise ValueError(f"recomputed ordered HVG-DEG axis differs: {size}")
        if prediction.symbols != tuple(axis.symbol.iloc[list(expected)]):
            raise ValueError(f"HVG-DEG symbols differ from fixed gene identity: {size}")
        prior = set(expected)
    return published


def verify_coverage(release, corpus_path, required_sources):
    report = {}
    with open_corpus(corpus_path) as corpus:
        entries = corpus.manifest["sources"]
        missing = set(required_sources).difference(entries)
        if missing:
            raise ValueError(f"corpus is missing required sources: {sorted(missing)}")
        for name, metadata in entries.items():
            if metadata["shared_axis_complete"] or name in required_sources:
                projection = release.axis(release.sizes[-1]).project(corpus, name)
                observations = projection.source.shape[0]
                sample = projection.counts(0, min(1, observations))
                report[name] = {"shared_axis_complete": True, "observations": observations,
                                "projected_sample_shape": list(sample.shape),
                                "projected_sample_dtype": str(sample.dtype)}
            else:
                report[name] = {"shared_axis_complete": False}
        return {"manifest_sha256": corpus.manifest_sha256, "sources": report}


def run(action, source_root, output, *, corpus=None, required_sources=("tenv", "tahoe_bulk")):
    candidate = legacy.release.candidate(output) if action == "publish" else None
    try:
        released = (_candidate(source_root.resolve(), candidate) if candidate is not None
                    else verify_release(output))
        result = {"axes": list(released.sizes), "manifest_sha256": released.manifest_sha256}
        if corpus is not None:
            result["corpus"] = verify_coverage(released, corpus, required_sources)
        if candidate is not None:
            legacy.release.replace(candidate, output)
        return result
    finally:
        if candidate is not None and candidate.exists():
            shutil.rmtree(candidate)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--source-root", type=Path, default=STAGE2.parent)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--corpus", type=Path)
    parser.add_argument("--required-source", action="append", dest="required_sources")
    args = parser.parse_args()
    result = run(args.action, args.source_root, args.output, corpus=args.corpus,
                 required_sources=args.required_sources or ("tenv", "tahoe_bulk"))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
