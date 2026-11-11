#!/usr/bin/env python3
"""Replay a standalone tenv PR panel from its sealed numeric source."""
import argparse
from pathlib import Path
from bases.plots import replay_retrieval

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    replay_retrieval(args.report,args.output)
