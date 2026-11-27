"""The registered order and parameters of the eight Step 06.2 coordinates."""

from itertools import product

CONTROL = "original"
DIRECT_DELTA = "output-delta"
NO_X2 = "output-delta-nox2"
MLP = "z-interact_MLP"
FILM = "z-interact_FiLM"
SHARED = "shared-z"
RESIDUAL = "z-res"
BOTH = "shared-z-res"
SHAPES = (CONTROL, DIRECT_DELTA, NO_X2, MLP, FILM, SHARED, RESIDUAL, BOTH)
HEAD_PARAMETERS = {shape: {MLP: 106485, FILM: 106488}.get(shape, 106600) for shape in SHAPES}
BASE_OBJECTIVE = "Step 05.2"
NO_X2_OBJECTIVE = "No X2 MSE, double Delta MSE"


def coordinates(smoke=False, *, include_original=True):
    shapes = SHAPES if include_original else SHAPES[1:]
    return tuple(product(shapes, (0,) if smoke else range(3)))


def objective_label(shape):
    return NO_X2_OBJECTIVE if shape == NO_X2 else BASE_OBJECTIVE
