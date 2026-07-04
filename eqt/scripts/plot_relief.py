#!/usr/bin/env python3
"""Plot the latest (VELEST-model) relocation over a grayscale shaded-relief map of
Yogyakarta. Relief is neutral gray so the depth-coloured hypocentres stand out.
Runs in the `gmt` env (pygmt). Saves eqt/figures/aftershock_relief_map.png
"""
import os, json, argparse
import numpy as np, pandas as pd
import pygmt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default=f"{ROOT}/full/catalog_velest.csv")
    ap.add_argument("--out", default=f"{ROOT}/figures/aftershock_relief_map.png")
    ap.add_argument("--wellconstrained", action="store_true", default=True)
    a = ap.parse_args()

    e = pd.read_csv(a.catalog)
    if a.wellconstrained and "gap" in e:
        e = e[(e.gap < 180) & (e.get("errh_km", 0) < 5) & (e.rms < 0.5)]
    region = [110.15, 110.72, -8.12, -7.70]

    # station sites
    periods = json.load(open(f"{ROOT}/config/stations_periods.json"))
    slon=[]; slat=[]
    for info in periods.values():
        for s in info.get("sites", []): slon.append(s["lon"]); slat.append(s["lat"])

    fig = pygmt.Figure()
    grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)
    # neutral GRAYSCALE relief + hillshade -> contrasts with coloured points
    pygmt.makecpt(cmap="gray", series=[-200, 3000], continuous=True)
    fig.grdimage(grid, region=region, projection="M16c", cmap=True,
                 shading="+a315+nt0.8", frame=["WSne", "xa0.1f0.05", "ya0.1f0.05"])
    fig.coast(region=region, projection="M16c", shorelines="0.5p,black",
              rivers="a/0.5p,steelblue", resolution="f")
    # hypocentres coloured by depth with a vivid cpt
    pygmt.makecpt(cmap="turbo", series=[0, 20], reverse=False)
    fig.plot(x=e.longitude, y=e.latitude, fill=e.depth, cmap=True,
             style="c0.07c", pen="0.1p,black", transparency=15)
    # stations
    fig.plot(x=slon, y=slat, style="t0.45c", fill="white", pen="1.2p,black")
    fig.colorbar(cmap=True, frame=["x+lhypocentre depth", "y+lkm"],
                 position="JMR+o0.8c/0c+w9c/0.4c")
    # topo colourbar (gray)
    pygmt.makecpt(cmap="gray", series=[-200, 3000])
    fig.colorbar(cmap=True, frame=["x+ltopography", "y+lm"],
                 position="JBC+o0c/1.2c+w8c/0.35c+h")
    fig.text(x=region[0]+0.02, y=region[3]-0.02, text=f"Yogyakarta 2006 aftershocks (VELEST reloc, n={len(e)})",
             font="12p,Helvetica-Bold,black", justify="TL", offset="0.1c/-0.1c", fill="white@30")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=200)
    print(f"n={len(e)} plotted; wrote {a.out}")

if __name__ == "__main__":
    main()
