"""Publish and replay the four accepted HVG-DEG axes without expression data."""
import argparse
import json
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

STEP = Path(__file__).resolve().parents[1]
REPO = next(parent for parent in STEP.parents if (parent / ".git").exists())
STAGE2 = REPO / "stage2"
NATIVE_CORPUS = STAGE2 / "01_data_raw/output"
sys.path.insert(0, str(STEP / "output/api"))
sys.path.insert(0, str(STAGE2 / "01_data_raw/output/api"))
sys.path.insert(0, str(STAGE2 / "util/output/api"))

import axis_selection
from axis_selection_kernel import ALGORITHM, AXIS_SIZES
from stage2_util.derivation import release
from native_corpus import open_corpus
from stage2_util.publication import publication_lock, replace_release
from prediction_axis import COMMON_SHA256, COMMON_SIZE, open_axes, sha256

OUTPUT = STEP / "output/release"

@dataclass(frozen=True)
class SelectionInputs:
    genes: Path
    tenv_hvg: Path
    fusion_root: Path


def _candidate(source, output):
    """Copy the verified owner release with the original provenance intact."""
    accepted = verify_release(source)
    for name, record in accepted.manifest["outputs"].items():
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, target)
        if sha256(target) != record["sha256"]:
            raise ValueError(f"axis evidence changed during publication: {name}")
    shutil.copyfile(source / "manifest.json", output / "manifest.json")
    if sha256(output / "manifest.json") != accepted.manifest_sha256:
        raise ValueError("axis manifest changed during publication")
    return verify_release(output)


def verify_release(output):
    published = open_axes(output)
    manifest = published.manifest
    if manifest["algorithm"] != ALGORITHM:
        raise ValueError("unrecognized published HVG-DEG selection algorithm")
    release.check_manifest(output, manifest)
    evidence = output / "evidence"
    identity = STAGE2 / "util/output/gene-identity"
    config = SelectionInputs(
        genes=identity / "genes.csv",
        tenv_hvg=evidence / "tenv-hvg/descriptive.csv",
        fusion_root=evidence / "disease-fusion",
    )
    if sha256(identity / "tables/genes.csv") != COMMON_SHA256:
        raise ValueError("fixed Stage 2 common gene identity differs")
    axis = axis_selection.read_axis(config.genes)
    if len(axis) != COMMON_SIZE:
        raise ValueError("selection input does not cover the fixed common genes")
    tenv_hvg = axis_selection.read_tenv(config.tenv_hvg, COMMON_SIZE)
    fusion = axis_selection.read_fusion(config, COMMON_SIZE)
    selected = axis_selection.compute(axis, tenv_hvg, fusion, AXIS_SIZES)["axes"]
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
    if action == "publish":
        corpus = NATIVE_CORPUS if corpus is None else Path(corpus)
        if not (corpus / "manifest.json").is_file():
            raise ValueError(f"axis publication requires a native corpus manifest: {corpus}")
        if Path(output).resolve() == OUTPUT.resolve():
            required_sources = set(required_sources) | {"tenv", "tahoe_bulk"}
    with publication_lock(output) as destination:
        candidate = release.candidate(destination) if action == "publish" else None
        try:
            source = source_root.resolve() / STEP.relative_to(REPO) / "output/release"
            released = (_candidate(source, candidate) if candidate is not None
                        else verify_release(destination))
            result = {"axes": list(released.sizes), "manifest_sha256": released.manifest_sha256}
            if corpus is not None:
                result["corpus"] = verify_coverage(released, corpus, required_sources)
            if candidate is not None:
                replace_release(candidate, destination, release.replace,
                                verify_existing=verify_release)
            return result
        finally:
            if candidate is not None and candidate.exists():
                shutil.rmtree(candidate)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--source-root", type=Path, default=REPO)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--corpus", type=Path)
    parser.add_argument("--required-source", action="append", dest="required_sources")
    args = parser.parse_args()
    result = run(args.action, args.source_root, args.output, corpus=args.corpus,
                 required_sources=args.required_sources or ("tenv", "tahoe_bulk"))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
