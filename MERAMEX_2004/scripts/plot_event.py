#!/usr/bin/env python3
"""Record section for one associated event: waveforms straight from the archive
with the EQTransformer picks overlaid, sorted by epicentral distance.

This is the visual check that the associations are real earthquakes and that the
p0/p1/p2 -> Z/N/E convention gives sensible P-on-vertical arrivals.

Usage:
  plot_event.py --event-idx 42 --assignments ../pilot/events_assignments.csv \
      --catalog ../pilot/catalog_nll.csv --out ../figures/event_42.png
"""
import argparse, json, math, os, sys, warnings

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy import UTCDateTime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mxio import build_stream, load_index

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event-idx", type=int, required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--catalog", default=None,
                    help="NLL catalog; used for the hypocentre in the title")
    ap.add_argument("--periods", default=os.path.join(HERE, "..", "config",
                                                      "stations_periods.json"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--pre", type=float, default=10.0)
    ap.add_argument("--post", type=float, default=60.0)
    ap.add_argument("--component", default="Z")
    ap.add_argument("--freqmin", type=float, default=2.0)
    ap.add_argument("--freqmax", type=float, default=20.0)
    ap.add_argument("--max-traces", type=int, default=24)
    a = ap.parse_args()

    asg = pd.read_csv(a.assignments)
    g = asg[asg.event_idx == a.event_idx]
    if g.empty:
        sys.exit(f"event_idx {a.event_idx} not in {a.assignments}")
    t0 = UTCDateTime(float(g.time.min()))

    periods = json.load(open(a.periods))

    def site(pid):
        base = pid[:-1] if pid[-1].islower() and pid[:-1] in periods else pid
        sites = periods[base]["sites"]
        k = 0 if len(sites) == 1 else ord(pid[-1]) - ord("a")
        return base, sites[min(k, len(sites) - 1)]

    hyp = None
    if a.catalog and os.path.exists(a.catalog):
        cat = pd.read_csv(a.catalog)
        cat["t"] = pd.to_datetime(cat.time_utc, utc=True).astype("int64") / 1e9
        near = cat.iloc[(cat.t - t0.timestamp).abs().argmin()]
        if abs(near.t - t0.timestamp) < 60:
            hyp = near

    index = {(r["station"], r["day"]): r for r in load_index()}

    rows = []
    for pid, gg in g.groupby("station"):
        base, s = site(pid)
        if hyp is not None:
            d = math.hypot((s["lat"] - hyp.latitude) * 111.19,
                           (s["lon"] - hyp.longitude) * 111.19 *
                           math.cos(math.radians(s["lat"])))
        else:
            d = float(gg.time.min() - g.time.min())
        rows.append((d, base, pid, gg))
    rows.sort()
    rows = rows[:a.max_traces]

    fig, ax = plt.subplots(figsize=(11, 0.42 * len(rows) + 2.6))
    plotted = 0
    for i, (d, base, pid, gg) in enumerate(rows):
        rec = index.get((base, t0.julday))
        if rec is None:
            continue
        st = build_stream(rec["path"], base, rec["kind"])
        if st is None:
            continue
        st = st.select(component=a.component)
        if not len(st):
            continue
        st = st.slice(t0 - a.pre, t0 + a.post).copy()
        if not len(st) or st[0].stats.npts < 10:
            continue
        tr = st[0]
        tr.detrend("demean")
        tr.filter("bandpass", freqmin=a.freqmin, freqmax=a.freqmax,
                  corners=4, zerophase=True)
        y = tr.data.astype(float)
        m = np.max(np.abs(y)) or 1.0
        tt = tr.stats.starttime - t0 + np.arange(tr.stats.npts) / tr.stats.sampling_rate
        ax.plot(tt, y / m * 0.45 + i, lw=0.45, color="0.2")
        for _, r in gg.iterrows():
            dt = float(r.time) - t0.timestamp
            ax.plot([dt, dt], [i - 0.45, i + 0.45], lw=1.4,
                    color="#d62728" if r.phase == "P" else "#1f77b4")
        ax.text(-a.pre - 1.0, i, f"{pid}  {d:5.0f} km", ha="right", va="center",
                fontsize=7.5, family="monospace")
        plotted += 1

    ax.set_yticks([])
    ax.set_xlim(-a.pre - 14, a.post)
    ax.set_xlabel(f"Seconds relative to first pick ({t0.strftime('%Y-%m-%d %H:%M:%S')} UTC)")
    title = f"event {a.event_idx} — {plotted} stations, {a.component} comp, " \
            f"{a.freqmin:g}–{a.freqmax:g} Hz bandpass"
    if hyp is not None:
        title += (f"\nNLL: {hyp.latitude:.3f}, {hyp.longitude:.3f}, "
                  f"{hyp.depth:.1f} km — RMS {hyp.rms:.2f} s, gap {hyp.gap:.0f}°, "
                  f"Nphs {int(hyp.nphs)}")
    ax.set_title(title, fontsize=11)
    ax.plot([], [], color="#d62728", lw=1.4, label="P pick")
    ax.plot([], [], color="#1f77b4", lw=1.4, label="S pick")
    ax.legend(loc="lower right", fontsize=8)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=150, bbox_inches="tight")
    print("wrote", os.path.abspath(a.out), f"({plotted} traces)")


if __name__ == "__main__":
    main()
