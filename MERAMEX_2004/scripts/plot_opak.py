#!/usr/bin/env python3
"""2004 MERAMEX seismicity vs the 2006 Yogyakarta rupture zone.

Left  : Central Java overview — where the 2004 earthquakes actually are.
Right : zoom on the future 2006 rupture, with the MERAMEX station geometry and
        a detection-capability proxy (distance to the 4th-nearest station, which
        is what sets whether a small event can be located at all).

Usage:
  plot_opak.py --catalog ../wide11/catalog_nll.csv \
      --stations ../wide11/events_stations.csv \
      --aftershocks ../../eqt/full/catalog_growclust.csv --out ../figures/opak_2004_vs_2006.png
"""
import argparse, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_pilot import coast_segments

# 2006 aftershock-zone centroid and box, from the repo's own GrowClust catalog
BOX = dict(lat=(-8.01, -7.87), lon=(110.35, 110.53))
CENTRE = (-7.940, 110.460)


def km_grid(sta, lat_lo, lat_hi, lon_lo, lon_hi, n=180, k=4):
    """Distance (km) to the k-th nearest station on a lat/lon grid."""
    la = np.linspace(lat_lo, lat_hi, n)
    lo = np.linspace(lon_lo, lon_hi, n)
    LO, LA = np.meshgrid(lo, la)
    d = np.stack([np.hypot((LA - r.latitude) * 111.19,
                           (LO - r.longitude) * 111.19 * np.cos(np.radians(LA)))
                  for r in sta.itertuples()])
    return LO, LA, np.sort(d, axis=0)[k - 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--stations", required=True)
    ap.add_argument("--aftershocks", required=True)
    ap.add_argument("--coast", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "data", "coastline.xy"))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cat = pd.read_csv(a.catalog)
    sta = pd.read_csv(a.stations).drop_duplicates("id")
    aft = pd.read_csv(a.aftershocks)
    coast = coast_segments(a.coast) if a.coast and os.path.exists(a.coast) else []

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(15, 6.6))

    # ---------- left: regional overview ----------
    for seg in coast:
        axL.plot(seg[:, 0], seg[:, 1], color="0.55", lw=0.7, zorder=1)
    axL.scatter(sta.longitude, sta.latitude, marker="^", s=26, facecolor="none",
                edgecolor="0.35", lw=0.6, label=f"{len(sta)} MERAMEX sites", zorder=3)
    sc = axL.scatter(cat.longitude, cat.latitude, c=cat.depth, s=34, cmap="plasma_r",
                     vmin=0, vmax=150, edgecolor="k", lw=0.3, zorder=4,
                     label=f"{len(cat)} events, 3–13 Jun 2004")
    plt.colorbar(sc, ax=axL, label="Depth (km)", shrink=0.85)
    axL.add_patch(plt.Rectangle((BOX["lon"][0], BOX["lat"][0]),
                                BOX["lon"][1] - BOX["lon"][0],
                                BOX["lat"][1] - BOX["lat"][0],
                                fill=False, ec="#d62728", lw=2, zorder=5))
    axL.annotate("2006 rupture zone", CENTRE[::-1], xytext=(110.0, -7.55),
                 color="#d62728", fontsize=10, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.3), zorder=6)
    axL.set_xlim(108.7, 112.8); axL.set_ylim(-9.8, -5.6)
    axL.set_xlabel("Longitude (°E)"); axL.set_ylabel("Latitude (°N)")
    axL.set_title("Central Java, 11 days of 2004 MERAMEX recording\n"
                  "seismicity concentrates offshore; almost none on land", fontsize=11)
    axL.legend(loc="lower left", fontsize=8.5)
    axL.set_aspect(1 / np.cos(np.radians(7.7)))

    # ---------- right: zoom on the 2006 rupture ----------
    lat_lo, lat_hi, lon_lo, lon_hi = -8.25, -7.60, 110.10, 110.85
    LO, LA, D4 = km_grid(sta, lat_lo, lat_hi, lon_lo, lon_hi)
    cf = axR.contourf(LO, LA, D4, levels=[0, 10, 15, 20, 30, 50], cmap="Greens_r", alpha=0.45)
    plt.colorbar(cf, ax=axR, label="Distance to 4th-nearest station (km)", shrink=0.85)

    axR.scatter(aft.lon, aft.lat, s=1.2, color="#d62728", alpha=0.10, lw=0,
                label=f"{len(aft):,} aftershocks of the 2006 Mw 6.3 (this repo)", zorder=2)
    axR.scatter(sta.longitude, sta.latitude, marker="^", s=95, facecolor="w",
                edgecolor="k", lw=1.0, zorder=4, label="MERAMEX 2004 stations")
    for r in sta.itertuples():
        if lat_lo < r.latitude < lat_hi and lon_lo < r.longitude < lon_hi:
            axR.annotate(r.id, (r.longitude, r.latitude), xytext=(4, 4),
                         textcoords="offset points", fontsize=7, zorder=5)
    inzoom = cat[cat.latitude.between(lat_lo, lat_hi) & cat.longitude.between(lon_lo, lon_hi)]
    axR.scatter(inzoom.longitude, inzoom.latitude, s=70, marker="*", color="#1f77b4",
                edgecolor="k", lw=0.4, zorder=6,
                label=f"2004 events in this window: {len(inzoom)}")
    axR.set_xlim(lon_lo, lon_hi); axR.set_ylim(lat_lo, lat_hi)
    axR.set_xlabel("Longitude (°E)"); axR.set_ylabel("Latitude (°N)")
    axR.set_title("The future 2006 rupture had the densest coverage in the network\n"
                  "— and produced no located event in these 11 days", fontsize=11)
    axR.legend(loc="lower left", fontsize=8.5, framealpha=0.92)
    axR.set_aspect(1 / np.cos(np.radians(7.94)))

    fig.tight_layout()
    fig.savefig(a.out, dpi=150, bbox_inches="tight")
    print("wrote", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
