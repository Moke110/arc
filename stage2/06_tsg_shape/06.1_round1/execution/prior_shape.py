"""Load the accepted #347 fit provenance without rebinding its historical runs."""

from itertools import product

from shape_model import STEP, axis_run
from shape_roster import PRIOR_NAMES

PRIOR_RECEIPT = STEP / "log/347-formal.json"
PRIOR_RECEIPT_SHA256 = "1b4e65c6a764fb984201d8270b3d8b7fefcc6b510552fecdcd61c26c3098031f"
PRIOR_RELEASE_REVISION = "9dee5fd90ea6d2bed775655b0205875c8812938b"
PRIOR_RUN_SCHEMA = "stage2-tenv-transigen-shape-round1-run-v1"
PRIOR_PROTOCOL = "matched-plate-control-transigen-shape-round1-v1"
PRIOR_MATRIX_SCHEMA = "stage2-tenv-transigen-shape-round1-matrix-v1"
PRIOR_MATRIX = STEP / "output/formal/prior-347/matrix"


def verify_matrix(matrix, coordinates, schema):
    receipt = axis_run.read_json(matrix / "manifest.json")
    if (receipt["schema"] != schema
            or receipt["preflight_sha256"] != axis_run.digest(matrix / "preflight.json")
            or [(fit["shape"], fit["fold"]) for fit in receipt["fits"]] != list(coordinates)):
        raise ValueError("formal matrix receipt has an invalid preflight or coordinate roster")
    paths, fits = {}, {}
    for (shape, fold), fit in zip(coordinates, receipt["fits"], strict=True):
        relative = f"runs/{shape}/fold-{fold}"
        path = matrix / relative
        if (fit["path"] != relative
                or fit["run_manifest_sha256"] != axis_run.digest(path / "manifest.json")
                or fit["registration_sha256"] != axis_run.digest(path / "registration.json")):
            raise ValueError(f"sealed matrix run differs: {shape}, fold {fold}")
        paths[shape, fold], fits[shape, fold] = path, fit
    return paths, fits


def load_prior(matrix, receipt_path=PRIOR_RECEIPT):
    if axis_run.digest(receipt_path) != PRIOR_RECEIPT_SHA256:
        raise ValueError("prior evidence differs from the sealed #347 receipt")
    receipt = axis_run.read_json(receipt_path)
    report = matrix.parent / "report"
    for path, field in ((matrix / "manifest.json", "matrix_manifest_sha256"),
                        (matrix / "preflight.json", "preflight_sha256"),
                        (report / "manifest.json", "report_manifest_sha256")):
        if axis_run.digest(path) != receipt[field]:
            raise ValueError(f"prior evidence differs from the sealed #347 source: {path}")
    if axis_run.read_json(report / "manifest.json")["outputs"] != receipt["report_outputs_sha256"]:
        raise ValueError("sealed #347 report output identities differ")
    axis_run.check_outputs(report, {"outputs": receipt["report_outputs_sha256"]})

    source_paths, source_fits = verify_matrix(
        matrix, tuple(product(PRIOR_NAMES, range(3))), PRIOR_MATRIX_SCHEMA)
    paths, fits = {}, {}
    for (old, fold), path in source_paths.items():
        canonical = PRIOR_NAMES[old]
        paths[canonical, fold] = path
        fits[canonical, fold] = {**source_fits[old, fold], "shape": canonical, "source_shape": old}
    provenance = {
        "ticket": receipt["ticket"], "release_revision": PRIOR_RELEASE_REVISION,
        "matrix_code_revision": receipt["matrix_code_revision"],
        "release_receipt_sha256": PRIOR_RECEIPT_SHA256,
        "matrix_manifest_sha256": receipt["matrix_manifest_sha256"],
        "preflight_sha256": receipt["preflight_sha256"],
        "report_manifest_sha256": receipt["report_manifest_sha256"],
        "report_outputs_sha256": receipt["report_outputs_sha256"],
        "coordinates": list(fits.values()),
        "step_05_2_reference_parity": receipt["reference_parity"],
    }
    return paths, fits, provenance
