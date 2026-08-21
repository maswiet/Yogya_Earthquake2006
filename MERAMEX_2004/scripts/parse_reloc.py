#!/usr/bin/env python3
"""Turn HypoDD `hypoDD.reloc` and GrowClust `out.growclust_cat` into catalog
CSVs with the same column names as `catalog_nll.csv`, so every downstream plot
can take any of the three.

Usage:
  parse_reloc.py --hypodd ../hypodd/hypoDD.reloc --out ../full/catalog_hypodd.csv
  parse_reloc.py --growclust ../growclust/OUT/out.growclust_cat \
      --out ../full/catalog_growclust.csv
"""
import argparse, os

import pandas as pd
from obspy import UTCDateTime

HYPODD_COLS = ["id", "latitude", "longitude", "depth", "x", "y", "z",
               "ex", "ey", "ez", "year", "month", "day", "hour", "minute",
               "second", "mag", "nccp", "nccs", "nctp", "ncts", "rcc", "rct", "cid"]
GC_COLS = ["year", "month", "day", "hour", "minute", "second", "id",
           "latitude", "longitude", "depth", "mag", "qid", "cid", "nbranch",
           "qnpair", "qndiffP", "qndiffS", "rmsP", "rmsS", "eh", "ez", "et",
           "latC", "lonC", "depC"]


def stamp(df):
    t = []
    for r in df.itertuples():
        sec = float(r.second)
        base = UTCDateTime(int(r.year), int(r.month), int(r.day),
                           int(r.hour), int(r.minute), 0)
        t.append(str(base + sec))
    df.insert(0, "time_utc", t)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hypodd")
    ap.add_argument("--growclust")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.hypodd:
        df = pd.read_csv(a.hypodd, sep=r"\s+", header=None, names=HYPODD_COLS)
        df = stamp(df)
        keep = ["time_utc", "latitude", "longitude", "depth", "ex", "ey", "ez", "cid", "id"]
        df = df[keep].rename(columns={"ex": "errx_m", "ey": "erry_m", "ez": "errz_m"})
        # HypoDD reports errors in metres when ISOLV=2 (LSQR); make that explicit
        df["errh_km"] = ((df.errx_m ** 2 + df.erry_m ** 2) ** 0.5) / 1000.0
        df["errz_km"] = df.errz_m / 1000.0
    elif a.growclust:
        df = pd.read_csv(a.growclust, sep=r"\s+", header=None, names=GC_COLS)
        df = stamp(df)
        df = df[["time_utc", "latitude", "longitude", "depth", "eh", "ez",
                 "nbranch", "cid", "id", "rmsP", "rmsS"]]
        df = df.rename(columns={"eh": "errh_km", "ez": "errz_km"})
    else:
        raise SystemExit("give --hypodd or --growclust")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"{len(df)} relocated events -> {a.out}")
    if "cid" in df:
        multi = df.groupby("cid").size()
        print(f"  {len(multi)} clusters; largest {multi.max()} events; "
              f"{int((multi > 1).sum())} clusters with >1 event")
    print(f"  median errH {df.errh_km.median():.2f} km, errZ {df.errz_km.median():.2f} km")


if __name__ == "__main__":
    main()
