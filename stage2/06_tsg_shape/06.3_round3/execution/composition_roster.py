"""The registered Step 06.3 comparison roster and bounded CUDA smoke."""

from itertools import product

CONTROL = "0_original"
DIRECT_DELTA = "1_output-delta"
NO_X2 = "2_output-delta-no-x2"
FILM = "3_z-interact_FiLM"
DIRECT_DELTA_FILM = "4_1+3"
NO_X2_FILM = "5_2+3"
SHAPES = (CONTROL, DIRECT_DELTA, NO_X2, FILM, DIRECT_DELTA_FILM, NO_X2_FILM)
ROUND_TWO_SHAPES = dict(zip(SHAPES, (
    "original", "output-delta", "output-delta-nox2", "z-interact_FiLM",
    "output-delta", "output-delta-nox2",
), strict=True))
COMBINATIONS = (DIRECT_DELTA_FILM, NO_X2_FILM)
SEEDS = (17, 42, 2026, 43)
SMOKE_FOLDS = (0, 1)
SMOKE_SEEDS = (17, 42)
REUSED_SEED = 42
NO_X2_SHAPES = (NO_X2, NO_X2_FILM)
HEAD_PARAMETERS = {shape: 106488 if shape in (FILM, *COMBINATIONS) else 106600 for shape in SHAPES}


def coordinates(smoke):
    folds, seeds = (SMOKE_FOLDS, SMOKE_SEEDS) if smoke else (range(3), SEEDS)
    return tuple(product(SHAPES, folds, seeds))


def fit_coordinates(smoke):
    return tuple(key for key in coordinates(smoke)
                 if smoke or key[0] in COMBINATIONS or key[2] != REUSED_SEED)


def objective_label(shape):
    return "No X2 MSE, double Delta MSE" if shape in NO_X2_SHAPES else "Step 05.2"
