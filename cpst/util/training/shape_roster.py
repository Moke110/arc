"""Current and historical Step 06.1 coordinate identities."""

from itertools import product

SHAPES = {
    "original": ((1200,), (800,), 100),
    "widen_z-110": ((1200,), (800,), 110),
    "widen_z-120": ((1200,), (800,), 120),
    "widen_z-130": ((1200,), (800,), 130),
    "widen_enc-2k": ((2000,), (800,), 100),
    "widen_enc-4k": ((4000,), (800,), 100),
    "widen_enc-6k": ((6000,), (800,), 100),
    "deepen_enc": ((1200, 1200), (800,), 100),
    "deepen_dec": ((1200,), (800, 800), 100),
    "deepen_both": ((1200, 1200), (800, 800), 100),
}
CONTROL = "original"
NEW_WIDTHS = ("widen_enc-2k", "widen_enc-4k")
SMOKE_SHAPES = (CONTROL, *NEW_WIDTHS)
PRIOR_NAMES = {
    "widen_z-100": CONTROL,
    "widen_z-110": "widen_z-110",
    "widen_z-120": "widen_z-120",
    "widen_z-130": "widen_z-130",
    "widen_enc": "widen_enc-6k",
    "deepen_enc": "deepen_enc",
    "deepen_dec": "deepen_dec",
    "deepen-both": "deepen_both",
}


def coordinates(smoke=False):
    return tuple(product(SMOKE_SHAPES if smoke else SHAPES, (0,) if smoke else range(3)))


def dimensions(shape):
    encoder, decoder, latent = SHAPES[shape]
    return dict(encoder=list(encoder), decoder=list(decoder), n_latent=latent,
                predictive_head_dim=latent + 400 + 32)
