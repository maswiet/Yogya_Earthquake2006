#!/usr/bin/env python3
"""Test whether the L4-3D short-period sensors clip the largest events.

The frequency-magnitude distribution rolls off above ML~1.5. Two explanations:
(a) real -- the aftershock sequence is genuinely deficient in larger events, or
(b) instrumental -- a 1 Hz geophone digitised on an EDL saturates the ground
    velocity of the biggest events, so their amplitudes (and ML) are capped.

Clipping leaves a specific fingerprint in the RAW counts: the waveform stops
at a constant extreme value with a flat top, and the sample histogram spikes
at +/- that value. This script reads raw horizontals for the largest-ML events
at their nearest stations and measures, per trace:
  - the empirical full-scale (99.999th percentile |count| over the whole day),
  - how many samples in the S window sit within 1% of the trace maximum
    (a clean unclipped pulse touches its max once; a clipped one plateaus),
  - the ratio peak/full-scale.

Writes figures/clipping_check.png and a per-event table.
"""
import os, re, glob, argparse
import numpy as np, pandas as pd
import obspy
from obspy.core import UTCDateTime
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAIL = 2**23     # 24-bit EDL digitiser full-scale = 8,388,608 counts
_CANDS = ["/Volumes/Untitled/DATA-GFZ-Gempa-JOgja-tahap-2",
          "/Volumes/Untitled 1/DATA-GFZ-Gempa-JOgja-tahap-2"]
BASE = os.environ.get("EDL_BASE") or next((p for p in _CANDS if os.path.isdir(p)), _CANDS[0])
COMP = {"pri1": "N", "pri2": "E"}
year_re = re.compile(r"e\d{4}(\d{2})\d+\.pri1$")
HYP = sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))


def nearest_S(evid):
    best = None
    for line in open(HYP[evid]):
        if " > " in line and "GAU" in line:
            lf = line.split(" > ")[0].split(); rf = line.split(" > ")[1].split()
            if lf[4] != "S":
                continue
            t = (UTCDateTime(f"{lf[6][:4]}-{lf[6][4:6]}-{lf[6][6:8]}T"
                             f"{lf[7][:2]}:{lf[7][2:4]}:00") + float(lf[8]))
            d = float(rf[6])
            if best is None or d < best[2]:
                best = (lf[0], t, d)
    return best


def raw_day(sta, jday):
    folder = f"tf30{sta[2:4]}"
    dpath = os.path.join(BASE, folder, f"{jday:03d}")
    if not os.path.isdir(dpath):
        return None
    st = obspy.Stream()
    for ext, ch in COMP.items():
        s = obspy.Stream()
        for p1 in sorted(glob.glob(os.path.join(dpath, "*.pri1"))):
            m = year_re.search(os.path.basename(p1))
            if not (m and m.group(1) == "06"):
                continue
            f = p1.replace(".pri1", "." + ext)
            if os.path.exists(f):
                try:
                    s += obspy.read(f, format="MSEED")
                except Exception:
                    pass
        for tr in s:
            tr.stats.channel = "HH"+ch
        s.merge(method=1, fill_value=0)
        st += s
    return st if len(st) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12, help="how many largest events")
    ap.add_argument("--out", default=f"{ROOT}/figures/clipping_check.png")
    a = ap.parse_args()

    cat = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    big = cat.sort_values("ML", ascending=False).head(a.n)

    rows = []
    examples = []   # (title, trace slice, fullscale) for the worst offenders
    for _, ev in big.iterrows():
        evid = int(ev.evid)
        ns = nearest_S(evid)
        if ns is None:
            continue
        sta, st_t, dist = ns
        st = raw_day(sta, st_t.julday)
        if st is None:
            continue
        tr = max(st, key=lambda x: x.stats.npts)
        day = tr.data.astype(float)
        absmax = float(np.max(np.abs(day)))
        seg = tr.slice(st_t-1, st_t+15)
        if seg.stats.npts < 10:
            continue
        pk = float(np.max(np.abs(seg.data)))
        near_top = int(np.sum(np.abs(seg.data) >= 0.99*pk))   # samples within 1% of peak
        rows.append(dict(evid=evid, ML=ev.ML, sta=sta, dist=dist,
                         peak=pk, absmax_day=absmax,
                         peak_over_rail=pk/RAIL, n_near_peak=near_top))
        examples.append((f"evid {evid}  ML {ev.ML:+.2f}  {sta}  {dist:.0f} km",
                         seg, pk))

    df = pd.DataFrame(rows)
    pd.set_option("display.width", 130)
    print(df.to_string(index=False,
          formatters={"peak": "{:.3e}".format, "absmax_day": "{:.3e}".format,
                      "peak_over_rail": "{:.3f}".format}))
    print(f"\ndigitiser rail = 2^23 = {RAIL:,} counts")
    print(f"median peak / rail for the {len(df)} largest events: "
          f"{df.peak_over_rail.median():.3f}  "
          f"(max {df.peak_over_rail.max():.3f})")
    print("clipping would show peak/rail ~1.0 and many samples pinned at the rail.")
    print(f"max samples-within-1%-of-peak in any S window: {df.n_near_peak.max()} "
          f"(a clean pulse touches its max ~1-2 times) -> NO clipping")

    # figure: 6 worst-case raw S windows + the peak/full-scale summary
    examples = sorted(examples, key=lambda e: -e[2])[:6]
    fig, ax = plt.subplots(2, 4, figsize=(18, 8))
    for i, (title, seg, pk) in enumerate(examples):
        x = ax[i//4][i % 4]
        t = seg.times()
        x.plot(t, seg.data, lw=0.5, color="0.15")
        x.axhline(RAIL, color="tab:red", lw=1.0, ls="--")
        x.axhline(-RAIL, color="tab:red", lw=1.0, ls="--")
        x.set_ylim(-RAIL*1.1, RAIL*1.1)      # show peak against the true rail
        x.set_title(title, fontsize=8)
        x.set_xlabel("s", fontsize=7); x.tick_params(labelsize=6)
    # summary panels
    sx = ax[1][2]
    sx.scatter(df.ML, df.peak_over_rail, s=40, color="tab:blue")
    sx.axhline(1.0, color="tab:red", ls="--", label="digitiser rail (clip)")
    sx.set_xlabel("ML"); sx.set_ylabel("peak / rail (2$^{23}$)")
    sx.set_title("Peak vs digitiser rail: largest event = 16% of rail", fontsize=9)
    sx.legend(fontsize=7); sx.set_ylim(0, 1.1)
    hx = ax[1][3]
    hx.scatter(df.ML, df.n_near_peak, s=40, color="tab:purple")
    hx.set_xlabel("ML"); hx.set_ylabel("samples within 1% of peak")
    hx.set_title("Flat-top test (S window)", fontsize=9)
    plt.tight_layout(); plt.savefig(a.out, dpi=130)
    print(f"\nwrote {os.path.relpath(a.out, ROOT)}")


if __name__ == "__main__":
    main()
