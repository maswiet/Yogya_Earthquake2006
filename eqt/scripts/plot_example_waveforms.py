#!/usr/bin/env python3
"""Waveform evidence for what the small-magnitude tail of the catalogue is.

Plots Wood-Anderson horizontals at the nearest stations for three representative
events -- a well-recorded ML~1.5 event, an event at the completeness magnitude,
and one of the ML<-1 events -- with the P/S picks, the amplitude measurement
window, and the pre-P noise window used to compute SNR.

Motivation: 39% of the ML<-1 events fall within +/-10 min of midnight (a 28x
enrichment confined entirely to that magnitude band), which points at a
day-boundary artifact in the per-julian-day amplitude measurement rather than
at genuine microseismicity. These panels show what that looks like on the trace.

Runs in the `eqt` env (obspy). Needs the raw EDL volume mounted.
"""
import os, re, glob, math, argparse
import numpy as np, pandas as pd
import obspy
from obspy.core import UTCDateTime
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE = os.environ.get("EDL_BASE", "/Volumes/Untitled/DATA-GFZ-Gempa-JOgja-tahap-2")
COMP = {"pri1": "N", "pri2": "E"}
year_re = re.compile(r"e\d{4}(\d{2})\d+\.pri1$")

w0 = 2*math.pi*1.0; h = 0.707
PAZ_L4 = {"poles": [complex(-h*w0,  w0*math.sqrt(1-h*h)),
                    complex(-h*w0, -w0*math.sqrt(1-h*h))],
          "zeros": [0j, 0j], "gain": 1.0, "sensitivity": 1.7e8}
PAZ_WA = {"poles": [complex(-6.2832, -4.7124), complex(-6.2832, 4.7124)],
          "zeros": [0j, 0j], "gain": 1.0, "sensitivity": 2080.0}
BP_LO, BP_HI, BP_CORNERS = 1.0, 20.0, 4      # must match build_amplitudes.py

HYP = sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))


def picks_of(evid):
    """P and S picks for one event from its NLLoc hypocentre file."""
    rows = []
    for line in open(HYP[evid]):
        if " > " in line and "GAU" in line:
            left, right = line.split(" > "); lf = left.split(); rf = right.split()
            sta, ph = lf[0], lf[4]
            if ph not in ("P", "S"):
                continue
            t = (UTCDateTime(f"{lf[6][:4]}-{lf[6][4:6]}-{lf[6][6:8]}T"
                             f"{lf[7][:2]}:{lf[7][2:4]}:00") + float(lf[8]))
            rows.append(dict(sta=sta, phase=ph, t=t, dist=float(rf[6])))
    return pd.DataFrame(rows)


def read_wa(sta, t0, t1):
    """Wood-Anderson horizontals for one station spanning [t0,t1].

    Reads every julian-day folder the window touches, so a window crossing
    midnight is assembled from both days rather than truncated at the boundary.
    """
    folder = f"tf30{sta[2:4]}"
    st = obspy.Stream()
    for jd in sorted({t0.julday, t1.julday}):
        dpath = os.path.join(BASE, folder, f"{jd:03d}")
        if not os.path.isdir(dpath):
            continue
        for ext, ch in COMP.items():
            s = obspy.Stream()
            for p1 in sorted(glob.glob(os.path.join(dpath, "*.pri1"))):
                m = year_re.search(os.path.basename(p1))
                if not (m and m.group(1) == "06"):
                    continue
                f = p1.replace(".pri1", "." + ext)
                if not os.path.exists(f):
                    continue
                try:
                    tr = obspy.read(f, format="MSEED", headonly=True)[0]
                    if tr.stats.endtime < t0-5 or tr.stats.starttime > t1+5:
                        continue
                    s += obspy.read(f, format="MSEED")
                except Exception:
                    pass
            for tr in s:
                tr.stats.channel = "HH"+ch
            st += s
    if len(st) == 0:
        return None
    st.merge(method=1, fill_value=None)
    st.trim(t0-5, t1+5)
    if len(st) == 0 or all(tr.stats.npts < 10 for tr in st):
        return None
    for tr in st:
        if hasattr(tr.data, "filled"):
            tr.data = tr.data.filled(0.0)
    st.detrend("demean")
    try:
        st.simulate(paz_remove=PAZ_L4, paz_simulate=PAZ_WA, water_level=10)
        # Same band as build_amplitudes.py, so the traces shown are the ones
        # the magnitudes were actually measured on (TF10b carries 50 Hz mains).
        st.filter("bandpass", freqmin=BP_LO, freqmax=BP_HI,
                  corners=BP_CORNERS, zerophase=True)
    except Exception:
        return None
    return st


def snr(tr, ptime, stime):
    """Peak WA amplitude in the S window over RMS of the pre-P noise window."""
    sig = tr.slice(stime-1.0, stime+15.0)
    noi = tr.slice(ptime-12.0, ptime-2.0)
    if sig.stats.npts < 10 or noi.stats.npts < 10:
        return np.nan, np.nan, np.nan
    pk = float(np.max(np.abs(sig.data)))
    rms = float(np.sqrt(np.mean(noi.data**2)))
    return pk, rms, (pk/rms if rms > 0 else np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default=None,
                    help="comma-separated evids; default = one per ML band")
    ap.add_argument("--nsta", type=int, default=3)
    ap.add_argument("--out", default=f"{ROOT}/figures/example_waveforms.png")
    a = ap.parse_args()

    cat = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    cat["t"] = pd.to_datetime(cat.time)

    if a.events:
        evids = [int(x) for x in a.events.split(",")]
    else:
        def pick(lo, hi):
            s = cat[(cat.ML >= lo) & (cat.ML < hi) & (cat.gap < 160)]
            return int(s.sort_values("n_sta", ascending=False).iloc[0].evid)
        evids = [pick(1.4, 2.2), pick(-0.6, -0.4), pick(-9, -1.0)]

    labels = ["well recorded", "at completeness Mc", "small-magnitude tail"]
    fig, axes = plt.subplots(len(evids), a.nsta,
                             figsize=(4.6*a.nsta, 3.1*len(evids)), squeeze=False)

    for row, evid in enumerate(evids):
        ev = cat[cat.evid == evid].iloc[0]
        pk = picks_of(evid)
        s_pk = pk[pk.phase == "S"].sort_values("dist")
        used = 0
        for _, p in s_pk.iterrows():
            if used >= a.nsta:
                break
            ptt = pk[(pk.sta == p.sta) & (pk.phase == "P")]
            ptime = UTCDateTime(ptt.iloc[0].t) if len(ptt) else UTCDateTime(p.t)-2
            st = read_wa(p.sta, ptime-15, UTCDateTime(p.t)+20)
            if st is None:
                continue
            tr = max(st, key=lambda x: x.stats.npts)
            ax = axes[row][used]
            rel = tr.times() + (tr.stats.starttime - UTCDateTime(p.t))
            ax.plot(rel, tr.data, lw=0.55, color="0.15")
            pkA, rms, s2n = snr(tr, ptime, UTCDateTime(p.t))
            ax.axvspan(-1, 15, color="tab:orange", alpha=0.16, zorder=0)
            ax.axvspan(float(ptime-UTCDateTime(p.t))-12,
                       float(ptime-UTCDateTime(p.t))-2,
                       color="tab:blue", alpha=0.13, zorder=0)
            ax.axvline(float(ptime-UTCDateTime(p.t)), color="tab:blue", lw=1.3)
            ax.axvline(0, color="tab:red", lw=1.3)
            ax.text(0.02, 0.95, f"{p.sta} {tr.stats.channel}  {p.dist:.0f} km",
                    transform=ax.transAxes, va="top", fontsize=8, weight="bold")
            ax.text(0.02, 0.06,
                    "SNR " + ("n/a" if np.isnan(s2n) else f"{s2n:.1f}"),
                    transform=ax.transAxes, fontsize=8,
                    color="tab:red" if (np.isnan(s2n) or s2n < 3) else "tab:green")
            if tr.stats.npts and np.count_nonzero(tr.data) < 0.8*tr.stats.npts:
                ax.text(0.98, 0.95, "DATA GAP", transform=ax.transAxes,
                        ha="right", va="top", fontsize=8, color="tab:red",
                        weight="bold")
            ax.set_xlim(rel.min(), rel.max()); ax.tick_params(labelsize=7)
            if used == 0:
                ax.set_ylabel(f"ML {ev.ML:+.2f}\n{labels[row]}\nWA disp. (mm)",
                              fontsize=8)
            if row == len(evids)-1:
                ax.set_xlabel("time relative to S pick (s)", fontsize=8)
            used += 1
        axes[row][0].set_title(
            f"evid {evid}   {ev.time[:19]} UTC   ML {ev.ML:+.2f}   "
            f"{int(ev.n_sta)} sta   gap {ev.gap:.0f}$\\degree$",
            fontsize=9, loc="left")
        for k in range(used, a.nsta):
            axes[row][k].axis("off")

    fig.suptitle("Wood-Anderson horizontals: P (blue) / S (red) picks, "
                 "S-amplitude window (orange), pre-P noise window (blue)",
                 fontsize=10)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(a.out, dpi=140)
    print("events plotted:", evids)
    print("wrote", os.path.relpath(a.out, ROOT))


if __name__ == "__main__":
    main()
