#!/usr/bin/env python3
"""Four-panel summary of a MERAMEX pilot catalog: map, depth section,
depth histogram and daily rate.

Usage:
  plot_pilot.py --catalog ../pilot/catalog_nll.csv --stations ../pilot/events_stations.csv \
      --out ../figures/pilot_summary.png
"""
import argparse, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
COAST = os.path.join(HERE, "..", "data", "coastline.xy")


def coast_segments(path):
    segs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if len(cur) > 1:
                segs.append(np.array(cur))
            cur = []
            continue
        p = line.split()
        if len(p) >= 2:
            try:
                cur.append((float(p[0]), float(p[1])))
            except ValueError:
                pass
    if len(cur) > 1:
        segs.append(np.array(cur))
    return segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--stations", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="MERAMEX 2004 pilot — EQTransformer + PyOcto + NonLinLoc")
    ap.add_argument("--good-only", action="store_true")
    a = ap.parse_args()

    cat = pd.read_csv(a.catalog)
    if a.good_only and "quality_pass" in cat:
        cat = cat[cat.quality_pass]
    sta = pd.read_csv(a.stations)
    cat["time_utc"] = pd.to_datetime(cat["time_utc"], utc=True)

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.4, 1], hspace=0.28, wspace=0.24)

    # --- map ---
    ax = fig.add_subplot(gs[0, 0])
    if os.path.exists(COAST):
        for s in coast_segments(COAST):
            ax.plot(s[:, 0], s[:, 1], color="0.55", lw=0.7, zorder=1)
    sc = ax.scatter(cat.longitude, cat.latitude, c=cat.depth, s=9, cmap="plasma_r",
                    vmin=0, vmax=max(40, np.nanpercentile(cat.depth, 95)),
                    edgecolors="none", alpha=0.85, zorder=3)
    ax.scatter(sta.longitude, sta.latitude, marker="^", s=34, c="none",
               edgecolors="k", linewidths=0.7, zorder=4, label=f"{len(sta)} stations")
    # frame both the array and the (mostly offshore) seismicity
    lon_lo = min(sta.longitude.min(), cat.longitude.quantile(0.01))
    lon_hi = max(sta.longitude.max(), cat.longitude.quantile(0.99))
    lat_lo = min(sta.latitude.min(), cat.latitude.quantile(0.01))
    lat_hi = max(sta.latitude.max(), cat.latitude.quantile(0.99))
    ax.set_xlim(lon_lo - 0.25, lon_hi + 0.25)
    ax.set_ylim(lat_lo - 0.25, lat_hi + 0.25)
    ax.set_aspect(1 / np.cos(np.radians(float(sta.latitude.mean()))))
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
    ax.set_title(f"{len(cat)} located events")
    ax.legend(loc="lower left", fontsize=8)
    plt.colorbar(sc, ax=ax, label="Depth (km)", shrink=0.85)

    # --- N-S cross-section ---
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(cat.latitude, cat.depth, s=8, c=cat.depth, cmap="plasma_r",
               vmin=0, vmax=max(40, np.nanpercentile(cat.depth, 95)), alpha=0.8)
    ax.scatter(sta.latitude, np.zeros(len(sta)) - 2, marker="^", s=26,
               c="none", edgecolors="k", linewidths=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("Latitude (°N)"); ax.set_ylabel("Depth (km)")
    ax.set_title("N–S section (depth vs latitude)")
    ax.grid(alpha=0.25)

    # --- depth histogram ---
    ax = fig.add_subplot(gs[1, 0])
    ax.hist(cat.depth.dropna(), bins=np.arange(0, max(60, cat.depth.max()) + 2, 2),
            color="#3b6ea5", edgecolor="white", linewidth=0.4)
    ax.set_xlabel("Depth (km)"); ax.set_ylabel("Events")
    ax.set_title(f"Depth distribution (median {cat.depth.median():.1f} km)")
    ax.grid(alpha=0.25)

    # --- daily rate ---
    ax = fig.add_subplot(gs[1, 1])
    daily = cat.set_index("time_utc").resample("1D").size()
    ax.bar(daily.index, daily.values, width=0.85, color="#c1582f")
    ax.set_ylabel("Events / day")
    ax.set_title(f"Daily rate (mean {daily.mean():.0f}/day)")
    ax.tick_params(axis="x", rotation=30, labelsize=8)
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle(a.title, fontsize=13, y=0.975)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=160, bbox_inches="tight")
    print("wrote", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
