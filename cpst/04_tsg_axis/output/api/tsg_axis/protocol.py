"""Published identifiers and stopping budget for the matched-control comparison."""

from cpst_util.tenv_folds import COLD_DRUG, open_tenv_folds

PROTOCOL = "matched-plate-control-native-loss-v2"
RUN_SCHEMA = "stage2-tenv-transigen-axis-run-v3"
REPORT_SCHEMA = "stage2-tenv-transigen-axis-report-v3"
MAX_EPOCHS = 1000
PATIENCE = 30


def shared_protocol():
    """Read the versioned fold and seed roster used for the rerun."""
    return open_tenv_folds().protocol(COLD_DRUG)
