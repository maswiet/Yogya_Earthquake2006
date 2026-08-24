#!/usr/bin/env python3
"""Input vs recovered checkerboard at a series of depths.

Usage:
  plot_checker.py --npz ../tomo5/checker_20km.npz --depths 5,15,25,40 \
      --out ../figures/checker_20km.png
"""
import argparse, os

import numpy as np
import pygmt

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--stations", default=os.path.join(HERE, "..", "full",
                                                       "events_land_stations.csv"))
    ap.add_argument("--depths", default="5,15,25,40")
    ap.add_argument("--phase", default="P")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default=None)
    a = ap.parse_args()

    import pandas as pd
    z = np.load(a.npz)
    nlat, nlon, nz = int(z["nlat"]), int(z["nlon"]), int(z["nz"])
    dz = float(z["dz"])
    recflat = z["recP" if a.phase == "P" else "recS"]
    # show the input only where the inversion was allowed to solve, so the two
    # rows are compared over the same domain
    trueflat = np.where(np.isfinite(recflat), z["true_v"], np.nan)
    true = trueflat.reshape(nlat, nlon, nz) * 100
    rec = recflat.reshape(nlat, nlon, nz) * 100
    lat = z["lat0"] + (np.arange(nlat) + 0.5) * z["dlat"]
    lon = z["lon0"] + (np.arange(nlon) + 0.5) * z["dlon"]
    sta = pd.read_csv(a.stations).drop_duplicates("id")
    depths = [float(x) for x in a.depths.split(",")]
    region = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]
    amp = float(z["amp"])

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="10p", FONT_ANNOT_PRIMARY="7p", MAP_FRAME_TYPE="plain")
    pygmt.makecpt(cmap="polar", series=[-amp, amp, amp / 20], background=True)

    n = len(depths)
    with fig.subplot(nrows=2, ncols=n, figsize=(f"{5.6*n}c", "11c"),
                     margins=["0.35c", "0.7c"], frame="lrtb"):
        for row, (field, label) in enumerate(((true, "input"), (rec, "recovered"))):
            for col, zc in enumerate(depths):
                k = int(np.clip(zc / dz, 0, nz - 1))
                layer = field[:, :, k]
                with fig.set_panel(panel=row * n + col):
                    grd = pygmt.xyz2grd(
                        x=np.repeat(lon, nlat), y=np.tile(lat, nlon),
                        z=layer.T.ravel(), region=region,
                        spacing=(float(z["dlon"]), float(z["dlat"])))
                    fig.grdimage(grd, region=region, projection="M5.2c",
                                 cmap=True, nan_transparent=True)
                    fig.coast(region=region, projection="M5.2c",
                              shorelines="0.4p,gray20",
                              frame=["WSne", "xa1f0.5", "ya1f0.5"])
                    fig.plot(x=sta.longitude, y=sta.latitude, style="t0.11c",
                             fill="white", pen="0.2p,black")
                    fig.text(position="TL", offset="0.12c/-0.12c", justify="TL",
                             text=f"{label}  {zc:.0f} km",
                             font="8p,Helvetica-Bold", fill="white@25",
                             pen="0.3p,black")
    fig.colorbar(position="JBC+w10c/0.32c+h+o0c/1.0c",
                 frame=[f'x+l"{a.phase}-velocity anomaly (%)"'])
    fig.savefig(a.out, dpi=200)
    print("wrote", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
