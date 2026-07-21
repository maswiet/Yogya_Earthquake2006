#!/usr/bin/env python3
"""Quality-screen the located catalogue, and estimate the false-detection rate.

RMS alone does not separate real events from spurious associations: NLLoc
down-weights badly-fitting phases, so an event built partly from noise picks
can still report a healthy RMS. evid 3319 is the type example -- RMS 0.178 s,
but with S residuals of -2.48 s (TF09b) and -3.76 s (TF11b).

This script parses every hypocentre file once and records, per event, the
metrics that DO separate them: the largest absolute phase residual, how many
phases are badly fitted, the P/S counts, and the usual gap/RMS/error terms.
It then reports how the failure rate varies with ML, which is the number a
reviewer will ask for when the catalogue claims sub-completeness detections.

Writes full/catalog_quality.csv (one row per event, all metrics + pass flag).
"""
import os, re, glob, argparse
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HYP = sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))

QUAL = re.compile(r"RMS\s+([\d.eE+-]+)\s+Nphs\s+(\d+)\s+Gap\s+([\d.]+)\s+Dist\s+([\d.]+)")
STAT = re.compile(r"Hyp\s+.*?ErrX\s+([\d.eE+-]+)\s+.*?ErrY\s+([\d.eE+-]+)\s+.*?ErrZ\s+([\d.eE+-]+)")

# A phase is "badly fitted" beyond this; ~5x the catalogue's typical RMS.
RES_BAD = 0.5
# Screening thresholds. Deliberately mild -- the point is to flag events that
# cannot be defended, not to sculpt the catalogue down to the cleanest core.
#
# max|res| is recorded but NOT used as a criterion: it scales with the number
# of phases (median nphs 10 -> 21 from ML<0 to ML>1, median max|res| 0.61 ->
# 1.41 over the same range), so thresholding it rejects precisely the
# best-recorded events. Screening on it inverted the pass rate with magnitude
# (58%/65%/61%); the scale-free criteria below give the physically expected
# 65%/78%/73%, rising with magnitude.
MAX_FRAC_BAD  = 0.25    # at most a quarter of phases badly fitted
MAX_RMS       = 0.5     # s
MIN_NPHS      = 8
MIN_NS        = 2       # S picks, needed for depth and for ML
MAX_GAP       = 180.


def parse_all():
    rows = []
    for evid, path in enumerate(HYP):
        rms = gap = dist = np.nan; nphs = 0
        res = []; nP = nS = 0
        for line in open(path):
            if line.startswith("QUALITY"):
                m = QUAL.search(line)
                if m:
                    rms = float(m.group(1)); nphs = int(m.group(2))
                    gap = float(m.group(3)); dist = float(m.group(4))
            elif " > " in line and "GAU" in line:
                left, right = line.split(" > "); lf = left.split(); rf = right.split()
                ph = lf[4]
                if ph not in ("P", "S"):
                    continue
                nP += ph == "P"; nS += ph == "S"
                res.append(abs(float(rf[1])))
        res = np.array(res) if res else np.array([np.nan])
        nbad = int(np.sum(res > RES_BAD))
        rows.append(dict(evid=evid, rms=rms, nphs=nphs, gap=gap, dist_min=dist,
                         nP=nP, nS=nS, max_res=float(np.nanmax(res)),
                         med_res=float(np.nanmedian(res)), n_bad=nbad,
                         frac_bad=nbad/max(len(res), 1)))
        if (evid+1) % 2000 == 0:
            print(f"  parsed {evid+1}/{len(HYP)}", flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{ROOT}/full/catalog_quality.csv")
    a = ap.parse_args()

    q = parse_all()
    q["pass"] = ((q.frac_bad <= MAX_FRAC_BAD) & (q.rms <= MAX_RMS) &
                 (q.nphs >= MIN_NPHS) & (q.nS >= MIN_NS) & (q.gap <= MAX_GAP))

    mag = f"{ROOT}/full/catalog_magnitude.csv"
    if os.path.exists(mag):
        q = q.merge(pd.read_csv(mag)[["evid", "ML", "ML_std"]], on="evid", how="left")

    q.to_csv(a.out, index=False)

    n = len(q); npass = int(q["pass"].sum())
    print(f"\nevents            : {n}")
    print(f"pass quality screen: {npass}  ({100*npass/n:.1f}%)")
    print(f"rejected          : {n-npass}  ({100*(n-npass)/n:.1f}%)")
    print("\nrejection reasons (not exclusive):")
    for lab, m in [(f">{MAX_FRAC_BAD:.0%} phases badly fitted", q.frac_bad > MAX_FRAC_BAD),
                   (f"rms > {MAX_RMS} s", q.rms > MAX_RMS),
                   (f"nphs < {MIN_NPHS}", q.nphs < MIN_NPHS),
                   (f"nS < {MIN_NS}", q.nS < MIN_NS),
                   (f"gap > {MAX_GAP:g}", q.gap > MAX_GAP)]:
        print(f"  {lab:32s} {int(m.sum()):6d}")

    if "ML" in q.columns:
        print("\nfailure rate vs magnitude:")
        print(f"  {'band':>14} {'n':>7} {'pass':>7} {'reject %':>9} {'med max|res|':>13}")
        for lo, hi, lab in [(-9, -1, "ML < -1"), (-1, -0.5, "-1..-0.5"),
                            (-0.5, 0, "-0.5..0"), (0, 0.5, "0..0.5"),
                            (0.5, 1, "0.5..1"), (1, 9, "ML > 1")]:
            m = (q.ML >= lo) & (q.ML < hi)
            if not m.sum():
                continue
            print(f"  {lab:>14} {int(m.sum()):7d} {int(q['pass'][m].sum()):7d} "
                  f"{100*(1-q['pass'][m].mean()):8.1f}% {q.max_res[m].median():13.2f}")
    print(f"\nwrote {os.path.relpath(a.out, ROOT)}")


if __name__ == "__main__":
    main()
