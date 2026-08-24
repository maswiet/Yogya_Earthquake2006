#!/usr/bin/env python3
"""Ray-coverage maps and sections for the MERAMEX tomography assessment.

Hit count per cell at a series of depths, drawn over the same relief/bathymetry
base as the catalogue map, with the stations and the resolved-area contour.

Usage:
  plot_tomo_coverage.py --tomo ../tomo5 --depths 5,15,25,40,60,100 \
      --out ../figures/tomo_coverage.png
"""
import argparse, os, sys

import numpy as np
import pygmt

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tomo", default=os.path.join(HERE, "..", "tomo5"))
    ap.add_argument("--stations", default=os.path.join(HERE, "..", "full",
                                                       "events_land_stations.csv"))
    ap.add_argument("--obs", default=os.path.join(HERE, "..", "config", "stations_info.csv"))
    ap.add_argument("--depths", default="5,15,25,40,60,100")
    ap.add_argument("--minhit", type=int, default=5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import pandas as pd
    g = np.load(os.path.join(a.tomo, "coverage.npz"))
    nlat, nlon, nz = int(g["nlat"]), int(g["nlon"]), int(g["nz"])
    hits = g["hits"].reshape(nlat, nlon, nz)
    lat = g["lat0"] + (np.arange(nlat) + 0.5) * g["dlat"]
    lon = g["lon0"] + (np.arange(nlon) + 0.5) * g["dlon"]
    dz = float(g["dz"])

    sta = pd.read_csv(a.stations).drop_duplicates("id")
    obs = pd.read_csv(a.obs)
    obs = obs[obs.kind.isin(["OBH", "OBS"])].drop_duplicates("sta")

    depths = [float(x) for x in a.depths.split(",")]
    region = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="11p", FONT_LABEL="9p", FONT_ANNOT_PRIMARY="8p",
                 MAP_FRAME_TYPE="plain")
    pygmt.makecpt(cmap="lajolla", series=[0, 200, 5], reverse=True, background=True)

    # Mercator aspect ratio: height/width ~ 1.189 at this latitude
    # panels need to be taller to keep proper shape
    with fig.subplot(nrows=2, ncols=3, figsize=("24c", "18c"), margins=["0.4c", "0.7c"],
                     frame="lrtb"):
        for panel, zc in enumerate(depths):
            k = int(np.clip(zc / dz, 0, nz - 1))
            layer = hits[:, :, k]
            with fig.set_panel(panel=panel):
                grd = pygmt.xyz2grd(
                    x=np.repeat(lon, nlat), y=np.tile(lat, nlon),
                    z=layer.T.ravel(), region=region,
                    spacing=(float(g["dlon"]), float(g["dlat"])))
                fig.grdimage(grd, region=region, projection="M8.0c",
                             cmap=True, nan_transparent=True)
                fig.grdcontour(grd, levels=[a.minhit], pen="0.9p,black",
                               annotation=None, region=region, projection="M8.0c")
                fig.coast(region=region, projection="M8.0c", shorelines="0.4p,gray30",
                          frame=["WSne", "xa1f0.5", "ya1f0.5"])
                fig.plot(x=sta.longitude, y=sta.latitude, style="t0.16c",
                         fill="white", pen="0.3p,black", region=region, projection="M8.0c")
                fig.plot(x=obs.lon, y=obs.lat, style="i0.20c",
                         fill="yellow", pen="0.4p,black", region=region, projection="M8.0c")
                nres = int((layer >= a.minhit).sum())
                fig.text(position="TL", offset="0.15c/-0.15c", justify="TL",
                         text=f"{zc:.0f} km   {nres} cells >= {a.minhit} rays",
                         font="8.5p,Helvetica-Bold", fill="white@25", pen="0.3p,black")
    fig.colorbar(position="JBC+w12c/0.35c+h+o0c/1.1c",
                 frame=['x+l"rays per cell"'])
    fig.savefig(a.out, dpi=200)
    print("wrote", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
