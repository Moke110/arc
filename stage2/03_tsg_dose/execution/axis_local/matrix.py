"""Print the registered dual-source formal origin matrix."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'output/api'))
from tsg_round1.matrix import plan
from tsg_round1.protocol import SCENARIOS


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['plan'])
    parser.add_argument('--scenario', choices=SCENARIOS)
    args = parser.parse_args()
    print(json.dumps(plan(args.scenario), indent=2))
