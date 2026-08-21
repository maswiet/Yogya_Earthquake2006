#!/usr/bin/env python3
"""Parse a NLLoc summary .hyp file into a catalog CSV with quality flags.

Usage:
  parse_nll.py --hyp ../nll/loc/pilot.sum.grid0.loc.hyp --out ../pilot/catalog_nll.csv
"""
import argparse, os, re
from datetime import datetime, timezone

import numpy as np
import pandas as pd

GEO = re.compile(r"^GEOGRAPHIC\s+OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                 r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
QUAL = re.compile(r"RMS\s+([\d.eE+-]+)\s+Nphs\s+(\d+)\s+Gap\s+([\d.]+)\s+Dist\s+([\d.]+)")
STAT = re.compile(r"CovXX\s+(-?[\d.eE+-]+).*?YY\s+(-?[\d.eE+-]+)\s+YZ\s+-?[\d.eE+-]+\s+"
                  r"ZZ\s+(-?[\d.eE+-]+)")


def parse(path):
    rows, cur = [], {}
    for line in open(path):
        if line.startswith("GEOGRAPHIC"):
            m = GEO.search(line)
            if m:
                y, mo, d, h, mi = map(int, m.groups()[:5])
                sec = float(m.group(6))
                cur = {"time": datetime(y, mo, d, h, mi, 0, tzinfo=timezone.utc).timestamp() + sec,
                       "latitude": float(m.group(7)), "longitude": float(m.group(8)),
                       "depth": float(m.group(9))}
        elif line.startswith("QUALITY"):
            m = QUAL.search(line)
            if m:
                cur.update(rms=float(m.group(1)), nphs=int(m.group(2)),
                           gap=float(m.group(3)), dist=float(m.group(4)))
        elif line.startswith("STATISTICS"):
            m = STAT.search(line)
            if m:
                cxx, cyy, czz = map(float, m.groups())
                cur["errh_km"] = round(float(np.sqrt(max(cxx, 0) + max(cyy, 0))), 2)
                cur["errz_km"] = round(float(np.sqrt(max(czz, 0))), 2)
            if cur:
                rows.append(cur)
            cur = {}
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hyp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-gap", type=float, default=180.0)
    ap.add_argument("--max-errh", type=float, default=5.0)
    ap.add_argument("--max-rms", type=float, default=0.5)
    ap.add_argument("--min-nphs", type=int, default=8)
    a = ap.parse_args()

    df = parse(a.hyp)
    if df.empty:
        print("no hypocentres parsed from", a.hyp)
        return
    df["time_utc"] = pd.to_datetime(df["time"], unit="s", utc=True)
    df = df.sort_values("time").reset_index(drop=True)
    good = df[(df.gap <= a.max_gap) & (df.errh_km <= a.max_errh) &
              (df.rms <= a.max_rms) & (df.nphs >= a.min_nphs)].copy()
    df["quality_pass"] = df.index.isin(good.index)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    df.to_csv(a.out, index=False)
    good.to_csv(a.out.replace(".csv", "_good.csv"), index=False)
    print(f"parsed {len(df)} hypocentres -> {a.out}")
    print(f"  RMS median {df.rms.median():.3f} s | gap median {df.gap.median():.0f} deg "
          f"| errH median {df.errh_km.median():.1f} km | errZ median {df.errz_km.median():.1f} km")
    print(f"  quality pass (gap<={a.max_gap}, errH<={a.max_errh}, RMS<={a.max_rms}, "
          f"Nphs>={a.min_nphs}): {len(good)}")
    if len(good):
        print(f"  depth: median {good.depth.median():.1f} km, "
              f"5-95% {good.depth.quantile(0.05):.1f}-{good.depth.quantile(0.95):.1f} km")


if __name__ == "__main__":
    main()
