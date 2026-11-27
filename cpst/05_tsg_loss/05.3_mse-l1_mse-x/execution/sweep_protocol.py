"""Step 05.3 coordinate and objective constants."""

from fractions import Fraction

ALPHAS = ("0", "0.5", "0.75", "0.9", "1", "1.1", "1.25", "1.5", "2")
ALPHA_LABELS = {Fraction(label): label for label in ALPHAS}
REFERENCE_ALPHA = "1"
SMOKE_ALPHAS = ("0", REFERENCE_ALPHA)
DELTA_L1_WEIGHT = Fraction(7, 4)
