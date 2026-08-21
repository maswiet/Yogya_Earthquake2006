#!/usr/bin/env python3
"""Regional PyGMT map of the MERAMEX 2004 catalogue over shaded relief and
Java-margin bathymetry (env `gmt`).

The sea is the point here — the trench, the forearc high and the Java Sea shelf
are what the OBS were deployed to constrain — so the ocean is shown with real
bathymetry rather than a flat blue fill, unlike the onshore Yogyakarta figures.

Usage:
  plot_map_gmt.py --catalog ../full/catalog_nll.csv --out ../figures/meramex_map.png
"""
import argparse, os

import numpy as np
import pandas as pd
import pygmt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
REGION = [108.4, 113.0, -10.9, -5.4]


def load_sites(path):
    """Station table -> (land, obs) frames with one row per site."""
    s = pd.read_csv(path)
    s = s.drop_duplicates(subset="sta")
    obs = s[s.kind.isin(["OBH", "OBS"])]
    land = s[~s.kind.isin(["OBH", "OBS"])]
    return land, obs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default=f"{ROOT}/wide11/catalog_nll.csv")
    ap.add_argument("--stations", default=f"{ROOT}/config/stations_info.csv")
    ap.add_argument("--out", default=f"{ROOT}/figures/meramex_map.png")
    ap.add_argument("--region", default=None, help="lon0/lon1/lat0/lat1")
    ap.add_argument("--resolution", default="30s",
                    help="earth_relief resolution; 30s pixel-registered is the "
                         "one whose tile cache covers the whole Java margin")
    ap.add_argument("--registration", default="pixel", choices=["pixel", "gridline"])
    ap.add_argument("--dmax", type=float, default=250.0)
    ap.add_argument("--profiles", default=None,
                    help="draw profile lines, e.g. '109.5,110.5,111.5'")
    ap.add_argument("--labels", action="store_true", help="label land stations")
    ap.add_argument("--title", default="MERAMEX 2004 — EQTransformer catalogue")
    a = ap.parse_args()

    region = [float(v) for v in a.region.split("/")] if a.region else REGION
    e = pd.read_csv(a.catalog)
    e = e[e.latitude.between(region[2], region[3]) &
          e.longitude.between(region[0], region[1])]
    e = e.sort_values("depth", ascending=False)      # shallow drawn on top
    land, obs = load_sites(a.stations)

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="13p,Helvetica-Bold", MAP_FRAME_TYPE="fancy")

    grid = pygmt.datasets.load_earth_relief(resolution=a.resolution, region=region,
                                            registration=a.registration)
    shade = pygmt.grdgradient(grid=grid, radiance=[315, 30], normalize="t0.6")
    # "geo" spans bathymetry and topography in one continuous scale
    pygmt.makecpt(cmap="geo", series=[-7000, 3000], continuous=True)
    fig.grdimage(grid=grid, region=region, projection="M17c", cmap=True,
                 shading=shade, frame=["WSne", "xa1f0.5", "ya1f0.5"])
    fig.coast(region=region, projection="M17c", shorelines="0.5p,black",
              resolution="i")

    # trench trace: the 6 km isobath is a fair proxy south of Java
    fig.grdcontour(grid=grid, levels=[-6000], pen="0.5p,black,--", annotation=None)

    if a.profiles:
        for lon in (float(v) for v in a.profiles.split(",")):
            fig.plot(x=[lon, lon], y=[region[2], region[3]], pen="1.2p,red,-")
            fig.text(x=lon, y=region[2] + 0.12, text=f"{lon:g}°E",
                     font="9p,Helvetica-Bold,red", justify="BC", fill="white@30")

    pygmt.makecpt(cmap="turbo", series=[0, a.dmax], reverse=False)
    fig.plot(x=e.longitude, y=e.latitude, fill=e.depth, cmap=True,
             style="c0.16c", pen="0.2p,black", transparency=10)

    fig.plot(x=land.lon, y=land.lat, style="t0.28c", fill="white", pen="0.6p,black")
    fig.plot(x=obs.lon, y=obs.lat, style="i0.40c", fill="yellow", pen="0.9p,black")
    for r in obs.itertuples():
        fig.text(x=r.lon, y=r.lat, text=r.sta, font="7p,Helvetica-Bold,black",
                 justify="CT", offset="0c/-0.30c", fill="white@30")
    if a.labels:
        for r in land.itertuples():
            fig.text(x=r.lon, y=r.lat, text=r.sta, font="5p,Helvetica,black",
                     justify="CB", offset="0c/0.18c")

    fig.colorbar(cmap=True, frame=["x+lhypocentre depth", "y+lkm"],
                 position="JMR+o0.9c/0c+w10c/0.4c")
    pygmt.makecpt(cmap="geo", series=[-7000, 3000])
    fig.colorbar(cmap=True, frame=["x+ltopography / bathymetry", "y+lm"],
                 position="JBC+o0c/1.3c+w11c/0.35c+h")

    fig.text(x=region[0] + 0.06, y=region[3] - 0.06,
             text=f"{a.title}  |  n = {len(e):,}",
             font="12p,Helvetica-Bold,black", justify="TL",
             fill="white@25", pen="0.5p,black", offset="0.1c/-0.1c")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=250)
    print(f"n={len(e)} events, {len(land)} land sites, {len(obs)} OBS/OBH; wrote {a.out}")


if __name__ == "__main__":
    main()
