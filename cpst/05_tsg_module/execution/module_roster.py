"""Register all five modules on the shared cold-drug roster."""

from itertools import product

from util.training.measured_roster import SHAPES


def protocol_roster():
    from tsg_axis.protocol import shared_protocol

    protocol = shared_protocol()
    return tuple(fold.number for fold in protocol.folds), protocol.seeds


def fit_coordinates(smoke):
    folds, seeds = protocol_roster()
    if smoke:
        return tuple(product(SHAPES, folds[:1], seeds[:1]))
    return tuple(product(SHAPES, folds, seeds))


coordinates = fit_coordinates
