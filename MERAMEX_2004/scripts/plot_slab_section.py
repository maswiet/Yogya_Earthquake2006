#!/usr/bin/env python3
"""Trench-perpendicular depth sections through the MERAMEX 2004 catalogue,
drawn with PyGMT (env `gmt`).

Each panel is a north–south swath at a fixed longitude: the topography and
bathymetry along the profile on top, the hypocentres below. Stacked west to
east they show how the Wadati–Benioff zone bends beneath Java and whether the
dip changes along strike.

Usage:
  plot_slab_section.py --catalog ../full/catalog_nll.csv \
      --profiles 109.6,110.5,111.4 --halfwidth 0.5 --out ../figures/meramex_slab.png
"""
import argparse, os

import numpy as np
import pandas as pd
import pygmt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

LAT0, LAT1 = -10.9, -5.4
W_MAP, H_TOPO = 15.0, 2.2                 # cm
KM_PER_DEG = 111.19


def topo_profile(grid, lon, lat0, lat1, n=400):
    """Elevation (m) sampled along a meridian of the pre-loaded relief grid."""
    track = pd.DataFrame({"lon": np.full(n, lon),
                          "lat": np.linspace(lat0, lat1, n)})
    return pygmt.grdtrack(points=track, grid=grid, newcolname="elev")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default=f"{ROOT}/wide11/catalog_nll.csv")
    ap.add_argument("--stations", default=f"{ROOT}/config/stations_info.csv")
    ap.add_argument("--profiles", default="109.6,110.5,111.4")
    ap.add_argument("--halfwidth", type=float, default=0.5, help="degrees of longitude")
    ap.add_argument("--zmax", type=float, default=250.0)
    ap.add_argument("--exaggeration", type=float, default=1.0,
                    help="vertical exaggeration; 1.0 = true scale, the honest "
                         "choice when the point of the figure is slab dip")
    ap.add_argument("--out", default=f"{ROOT}/figures/meramex_slab.png")
    a = ap.parse_args()

    e = pd.read_csv(a.catalog)
    sites = pd.read_csv(a.stations).drop_duplicates(subset="sta")
    lons = [float(v) for v in a.profiles.split(",")]
    pad = a.halfwidth + 0.2
    relief = pygmt.datasets.load_earth_relief(
        resolution="30s", registration="pixel",
        region=[min(lons) - pad, max(lons) + pad, LAT0, LAT1])

    # true-scale section height: same km per cm horizontally and vertically
    km_per_cm = (LAT1 - LAT0) * KM_PER_DEG / W_MAP
    h_sec = a.zmax / km_per_cm * a.exaggeration

    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL="10p", FONT_ANNOT_PRIMARY="9p",
                 MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="1p")

    for i, lon in enumerate(lons):
        band = e[(e.longitude - lon).abs() <= a.halfwidth]
        sband = sites[(sites.lon - lon).abs() <= a.halfwidth]

        # ---- topography / bathymetry strip ----
        tp = topo_profile(relief, lon, LAT0, LAT1)
        rel_km = tp.elev / 1000.0
        fig.basemap(region=[LAT0, LAT1, -7.5, 3.5], projection=f"X{W_MAP}c/{H_TOPO}c",
                    frame=["Wsne", "ya5f1+lkm"])
        fig.plot(x=[LAT0, LAT1], y=[0, 0], pen="0.5p,gray50,-")
        fig.plot(x=np.r_[LAT0, tp.lat, LAT1], y=np.r_[-7.5, rel_km, -7.5],
                 fill="lightbrown", pen="0.8p,black")
        fig.text(x=LAT0 + 0.08, y=2.4, text=f"{lon:g}°E  ±{a.halfwidth:g}°",
                 font="10p,Helvetica-Bold,black", justify="TL",
                 fill="white@25", pen="0.4p,black")
        # stations in the swath
        if len(sband):
            fig.plot(x=sband.lat, y=np.full(len(sband), 2.6),
                     style="t0.22c", fill="white", pen="0.5p,black")

        # ---- hypocentre section ----
        fig.shift_origin(yshift=f"-{h_sec + 0.15}c")
        frame = ["WSne" if i == len(lons) - 1 else "Wsne",
                 "ya50f10+ldepth (km)"]
        if i == len(lons) - 1:
            frame.append("xa1f0.5+llatitude (°N)")
        fig.basemap(region=[LAT0, LAT1, 0, a.zmax], projection=f"X{W_MAP}c/-{h_sec}c",
                    frame=frame)
        if len(band):
            pygmt.makecpt(cmap="hot", series=[0, 40], reverse=True)
            fig.plot(x=band.latitude, y=band.depth, fill=band.nphs, cmap=True,
                     style="c0.17c", pen="0.25p,black")
        ve = "true scale" if a.exaggeration == 1.0 else f"VE x{a.exaggeration:g}"
        fig.text(x=LAT1 - 0.08, y=a.zmax * 0.94, text=f"n = {len(band)}   ({ve})",
                 font="9p,Helvetica,black", justify="BR", fill="white@30")

        if i < len(lons) - 1:
            fig.shift_origin(yshift=f"-{H_TOPO + 1.5}c")

    fig.colorbar(cmap=True, frame=["x+lphases used in the location", "y+lN"],
                 position=f"JBC+o0c/1.5c+w{W_MAP * 0.7}c/0.35c+h")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=250)
    for lon in lons:
        n = ((e.longitude - lon).abs() <= a.halfwidth).sum()
        print(f"  {lon:g}°E ±{a.halfwidth:g}°: {n} events")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
