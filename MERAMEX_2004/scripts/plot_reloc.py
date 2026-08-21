#!/usr/bin/env python3
"""Before/after relocation maps and a depth section, in PyGMT (env `gmt`).

Absolute locations scatter by their own error ellipse; relative relocation
collapses the common part of that error, so structure that was a cloud becomes a
line or a plane. This figure puts the two side by side at identical scale so the
change is legible rather than asserted.

Usage:
  plot_reloc.py --before ../full/catalog_nll.csv \
      --after ../full/catalog_growclust.csv --label GrowClust \
      --out ../figures/full_reloc.png
"""
import argparse, os

import numpy as np
import pandas as pd
import pygmt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))


def panel(fig, df, region, title, width, dmax, relief, shade, stations=None):
    pygmt.makecpt(cmap="geo", series=[-7000, 3000], continuous=True)
    fig.grdimage(grid=relief, region=region, projection=f"M{width}c", cmap=True,
                 shading=shade, frame=["WSne", "xa0.5f0.1", "ya0.5f0.1"])
    fig.coast(region=region, projection=f"M{width}c", shorelines="0.5p,black",
              resolution="i")
    if stations is not None and len(stations):
        fig.plot(x=stations.lon, y=stations.lat, style="t0.22c", fill="white",
                 pen="0.5p,black")
    pygmt.makecpt(cmap="turbo", series=[0, dmax])
    if len(df):
        d = df.sort_values("depth", ascending=False)
        fig.plot(x=d.longitude, y=d.latitude, fill=d.depth, cmap=True,
                 style="c0.13c", pen="0.15p,black")
    fig.text(x=region[0] + 0.02, y=region[3] - 0.02, text=f"{title}  (n={len(df)})",
             font="11p,Helvetica-Bold,black", justify="TL", offset="0.1c/-0.1c",
             fill="white@25", pen="0.5p,black")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True)
    ap.add_argument("--after", required=True)
    ap.add_argument("--label", default="relocated")
    ap.add_argument("--stations", default=f"{ROOT}/config/stations_info.csv")
    ap.add_argument("--region", default=None, help="lon0/lon1/lat0/lat1")
    ap.add_argument("--dmax", type=float, default=150.0)
    ap.add_argument("--width", type=float, default=11.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    before = pd.read_csv(a.before)
    after = pd.read_csv(a.after)
    # keep only the events that actually relocated, and the same set before
    if "nbranch" in after:
        after = after[after.nbranch >= 2]
    if a.region:
        region = [float(v) for v in a.region.split("/")]
    else:
        lon = pd.concat([before.longitude, after.longitude])
        lat = pd.concat([before.latitude, after.latitude])
        pad = 0.15
        region = [lon.quantile(0.01) - pad, lon.quantile(0.99) + pad,
                  lat.quantile(0.01) - pad, lat.quantile(0.99) + pad]
    for df in (before, after):
        df.drop(df[~(df.longitude.between(region[0], region[1]) &
                     df.latitude.between(region[2], region[3]))].index, inplace=True)

    sites = pd.read_csv(a.stations).drop_duplicates(subset="sta")
    sites = sites[sites.lon.between(region[0], region[1]) &
                  sites.lat.between(region[2], region[3])]

    # 15s is pixel-registered only, and its tile cache stops at 10 deg S; fall
    # back to 30s when the frame reaches further south.
    try:
        relief = pygmt.datasets.load_earth_relief(resolution="15s", region=region,
                                                  registration="pixel")
    except Exception:
        relief = pygmt.datasets.load_earth_relief(resolution="30s", region=region,
                                                  registration="pixel")
    shade = pygmt.grdgradient(grid=relief, radiance=[315, 30], normalize="t0.6")

    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL="10p", FONT_ANNOT_PRIMARY="9p")
    panel(fig, before, region, "NonLinLoc (absolute)", a.width, a.dmax,
          relief, shade, sites)
    fig.shift_origin(xshift=f"{a.width + 1.0}c")
    panel(fig, after, region, a.label, a.width, a.dmax, relief, shade, sites)
    fig.colorbar(cmap=True, frame=["x+ldepth", "y+lkm"],
                 position="JMR+o0.8c/0c+w8c/0.35c")
    fig.savefig(a.out, dpi=250)

    print(f"wrote {a.out}")
    for name, df in (("before", before), (a.label, after)):
        if "errh_km" in df and len(df):
            print(f"  {name:<12} n={len(df):5d}  median errH "
                  f"{df.errh_km.median():6.2f} km  errZ {df.errz_km.median():6.2f} km")


if __name__ == "__main__":
    main()
