#!/usr/bin/env python3
"""How many station-days of a given kind are still unpicked.

The wrapper needs an exact remaining count to know when to stop; scraping it out
of the log is fragile (and once got it wrong, ending a run at 2,184 of 13,948).

Usage:  remaining.py ../full/done_land.txt EDL,SAM
"""
import os, sys

from mxio import load_index

prog, kinds = sys.argv[1], set(sys.argv[2].split(","))
dmin = int(sys.argv[3]) if len(sys.argv) > 3 else 127
dmax = int(sys.argv[4]) if len(sys.argv) > 4 else 282
min_bytes = int(sys.argv[5]) if len(sys.argv) > 5 else 5_000_000

done = set(open(prog).read().split()) if os.path.exists(prog) else set()
work = [r for r in load_index()
        if r["kind"] in kinds and dmin <= r["day"] <= dmax and r["nbytes"] >= min_bytes]
print(sum(1 for r in work if f"{r['station']}/{r['day']}" not in done))
