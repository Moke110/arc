"""The registered Step 06.4 arms, training seeds and bounded gate."""

from itertools import product

CONTROL = "XX-XX"
DIRECT_DELTA = "XX-XXΔ"
CONTROL_ONLY = "X-XΔ"
MEASURED_SHARED = "XΔ-XΔ"
MEASURED_SEPARATE = "XΔ-XΔΔ"
SHAPES = (CONTROL, DIRECT_DELTA, CONTROL_ONLY, MEASURED_SHARED, MEASURED_SEPARATE)
SOURCE_SHAPES = {CONTROL: "0_original", DIRECT_DELTA: "1_output-delta"}
MEASURED_DELTA = (MEASURED_SHARED, MEASURED_SEPARATE)
SEEDS = (17, 42, 2026)
SMOKE_FOLDS = (0, 1)
SMOKE_SEEDS = (17, 42)


def coordinates(smoke):
    folds, seeds = (SMOKE_FOLDS, SMOKE_SEEDS) if smoke else (range(3), SEEDS)
    return tuple(product(SHAPES, folds, seeds))


def fit_coordinates(smoke):
    return tuple(key for key in coordinates(smoke) if smoke or key[0] not in SOURCE_SHAPES)


def objective_label(shape):
    if shape == CONTROL_ONLY:
        return "Control KL and reconstruction, double Delta MSE"
    if shape in MEASURED_DELTA:
        return "Control and measured Delta reconstruction and Gaussian alignment"
    return "Step 05.2"
