"""Decide Step 06 promotion from the paired overall validation AP@100."""
from decimal import Decimal

SCREEN_SEED = 17
REPLICATION_SEEDS = (42, 2026)
THRESHOLD = Decimal('0.010')


def finite_decimal(value):
    result = Decimal(str(value))
    if not result.is_finite():
        raise ValueError('expected a finite validation effect')
    return result


def paired_gain(intervention, original):
    """Subtract the published metric values without binary threshold error."""
    if intervention is None or original is None:
        return None
    return finite_decimal(intervention) - finite_decimal(original)


def screen_decision(gains):
    """Omit pending seeds and mark terminal unavailable effects with None."""
    def decision(status, reason, next_seeds=()):
        return {'status': status, 'reason': reason, 'next_seeds': list(next_seeds)}

    if set(gains) - {SCREEN_SEED, *REPLICATION_SEEDS}:
        raise ValueError('unexpected screening seed')
    if SCREEN_SEED not in gains:
        return decision('pending_screen', 'seed 17 overall paired effect is pending', (SCREEN_SEED,))
    primary = gains[SCREEN_SEED]
    if primary is None:
        return decision('incomplete', 'seed 17 overall paired effect is unavailable')
    if finite_decimal(primary) < THRESHOLD:
        return decision('below_threshold', 'seed 17 overall paired effect is below +0.010 AP@100')
    pending = [seed for seed in REPLICATION_SEEDS if seed not in gains]
    if pending:
        return decision('pending_replication', 'seed 17 cleared the +0.010 AP@100 gate', pending)
    unavailable = [seed for seed in REPLICATION_SEEDS if gains[seed] is None]
    if unavailable:
        return decision('incomplete', f'overall paired effect unavailable for seeds {unavailable}')
    nonpositive = [seed for seed in REPLICATION_SEEDS if finite_decimal(gains[seed]) <= 0]
    if nonpositive:
        return decision('inconsistent', f'overall paired effect is nonpositive for seeds {nonpositive}')
    return decision('promoted', 'seed 17 cleared the gate and both added effects are positive')
