#!/usr/bin/env python3
"""Convert ph2dt's catalog differential times (`dt.ct`) into the GrowClust
`xcordata` format.

GrowClust was written for waveform cross-correlation differential times, but it
accepts any differential-time file in the same layout. Feeding it the catalog
times (with the ph2dt weight standing in for the correlation coefficient) gives
a catalog-only relocation - useful as the first pass, and as the baseline that
a later cross-correlation run has to beat.

  dt.ct    : "# ID1 ID2"        then  "STA TT1 TT2 WEIGHT PHA"
  xcordata : "# ID1 ID2 OTC"    then  "STA TDIF RXCOR PHA"   (tdif_fmt 12)

Usage:  ct2xcor.py --ct ../hypodd/dt.ct --out ../growclust/IN/xcordata.txt
"""
import argparse, os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ct", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-weight", type=float, default=0.0)
    a = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    npair = nobs = 0
    with open(a.ct) as src, open(a.out, "w") as dst:
        for line in src:
            if line.startswith("#"):
                _, id1, id2 = line.split()[:3]
                dst.write(f"# {id1} {id2} 0.000\n")
                npair += 1
                continue
            p = line.split()
            if len(p) < 5:
                continue
            sta, tt1, tt2, wt, pha = p[0], float(p[1]), float(p[2]), float(p[3]), p[4]
            if wt < a.min_weight:
                continue
            dst.write(f"{sta:<7}{tt1 - tt2:9.4f} {wt:6.3f} {pha}\n")
            nobs += 1
    print(f"{npair} pairs, {nobs} differential times -> {a.out}")


if __name__ == "__main__":
    main()
