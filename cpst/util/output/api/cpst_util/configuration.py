"""Stable hashes of ordered registered configurations."""

from hashlib import sha256
import json


def ordered_configurations_sha256(configurations):
    values = [json.dumps(config, sort_keys=True) for config in configurations]
    return sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()
