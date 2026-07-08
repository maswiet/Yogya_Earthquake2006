#!/usr/bin/env python3
"""Plot the latest (VELEST-model) relocation over a LIGHT grayscale shaded-relief
map of Yogyakarta, with a moderate-blue sea and XN station-code labels.
Runs in the `gmt` env (pygmt). Saves eqt/figures/aftershock_relief_map.png
"""
import os, argparse
import numpy as np, pandas as pd
import pygmt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def load_stations():
    rows = []
    for l in open(f"{ROOT}/config/stations_xn_meta.txt"):
        p = l.strip().split("|")
        if len(p) < 6: continue
        rows.append((p[1], float(p[4]), float(p[5])))   # name, lat, lon
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default=f"{ROOT}/full/catalog_velest.csv")
    ap.add_argument("--out", default=f"{ROOT}/figures/aftershock_relief_map.png")
    ap.add_argument("--title", default="Yogyakarta 2006 aftershocks (VELEST reloc")
    a = ap.parse_args()

    e = pd.read_csv(a.catalog)
    e = e.rename(columns={"lat": "latitude", "lon": "longitude"})   # accept HypoDD cols
    if "gap" in e:
        e = e[(e.gap < 180) & (e.get("errh_km", 0) < 5) & (e.rms < 0.5)]
    region = [110.15, 110.72, -8.12, -7.70]
    stations = load_stations()

    fig = pygmt.Figure()
    grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)
    # LIGHT grayscale relief: wide series pushes land into the bright half; gentle shade
    pygmt.makecpt(cmap="gray", series=[-1600, 1600], continuous=True)
    fig.grdimage(grid, region=region, projection="M16c", cmap=True,
                 shading="+a315+nt0.35", frame=["WSne", "xa0.1f0.05", "ya0.1f0.05"])
    # moderate-blue sea (covers the dark relief over water), coastline + rivers
    fig.coast(region=region, projection="M16c", water="151/193/226",
              shorelines="0.6p,black", rivers="a/0.4p,steelblue", resolution="f")
    # hypocentres coloured by depth (vivid -> contrasts with light gray land)
    pygmt.makecpt(cmap="turbo", series=[0, 20])
    fig.plot(x=e.longitude, y=e.latitude, fill=e.depth, cmap=True,
             style="c0.07c", pen="0.1p,black", transparency=20)
    # stations + code labels
    slon = [s[2] for s in stations]; slat = [s[1] for s in stations]
    fig.plot(x=slon, y=slat, style="t0.45c", fill="white", pen="1.2p,black")
    for name, la, lo in stations:
        fig.text(x=lo, y=la, text=name, font="9p,Helvetica-Bold,black",
                 justify="CB", offset="0c/0.28c", fill="white@25", pen="0.25p,black")
    fig.colorbar(cmap=True, frame=["x+lhypocentre depth", "y+lkm"],
                 position="JMR+o0.8c/0c+w9c/0.4c")
    pygmt.makecpt(cmap="gray", series=[-1600, 1600])
    fig.colorbar(cmap=True, frame=["x+ltopography", "y+lm"],
                 position="JBC+o0c/1.2c+w8c/0.35c+h")
    fig.text(x=region[0]+0.015, y=region[3]-0.015,
             text=f"{a.title}, n={len(e)})",
             font="12p,Helvetica-Bold,black", justify="TL", offset="0.1c/-0.1c",
             fill="white@20", pen="0.5p,black")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=200)
    print(f"n={len(e)} events, {len(stations)} stations; wrote {a.out}")

if __name__ == "__main__":
    main()
