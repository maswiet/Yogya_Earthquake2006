#!/usr/bin/env python3
"""Aftershock hypocentres over Yogyakarta relief: marker SIZE ~ ML, COLOR ~ depth.
Runs in the `gmt` env. Saves eqt/figures/aftershock_magnitude_map.png
"""
import os
import numpy as np, pandas as pd
import pygmt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def load_stations():
    rows=[]
    for l in open(f"{ROOT}/config/stations_xn_meta.txt"):
        p=l.strip().split("|")
        if len(p)>=6: rows.append((p[1],float(p[4]),float(p[5])))
    return rows

def main():
    e=pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    e=e[(e.gap<180) & (e.n_sta>=4)].copy()
    e=e[e.ML>=-0.5]                       # drop the unreliable low tail
    # marker size (cm) scaled with ML
    mlc=e.ML.clip(-0.5,3.5)
    e["size"]=0.04 + 0.085*(mlc-(-0.5))
    region=[110.15,110.72,-8.12,-7.70]
    stations=load_stations()

    fig=pygmt.Figure()
    grid=pygmt.datasets.load_earth_relief(resolution="03s", region=region)
    pygmt.makecpt(cmap="gray", series=[-1600,1600], continuous=True)
    fig.grdimage(grid, region=region, projection="M16c", cmap=True,
                 shading="+a315+nt0.35", frame=["WSne","xa0.1f0.05","ya0.1f0.05"])
    fig.coast(region=region, projection="M16c", water="151/193/226",
              shorelines="0.6p,black", rivers="a/0.4p,steelblue", resolution="f")
    pygmt.makecpt(cmap="turbo", series=[0,20])
    fig.plot(x=e.longitude, y=e.latitude, size=e["size"], fill=e.depth, cmap=True,
             style="cc", pen="0.2p,black", transparency=15)
    slon=[s[2] for s in stations]; slat=[s[1] for s in stations]
    fig.plot(x=slon, y=slat, style="t0.45c", fill="white", pen="1.2p,black")
    for name,la,lo in stations:
        fig.text(x=lo,y=la,text=name,font="9p,Helvetica-Bold,black",justify="CB",
                 offset="0c/0.28c",fill="white@25",pen="0.25p,black")
    fig.colorbar(cmap=True, frame=["x+lhypocentre depth","y+lkm"],
                 position="JMR+o0.8c/0c+w9c/0.4c")
    # ML size legend (reference circles near bottom-left)
    lx = region[0]+0.05
    for i,ml in enumerate([1.0,2.0,3.0]):
        ly = region[2]+0.06+i*0.05
        fig.plot(x=[lx], y=[ly], size=[0.04+0.085*(ml-(-0.5))], style="cc",
                 fill="gray40", pen="0.4p,black")
        fig.text(x=lx+0.03, y=ly, text=f"ML {ml:.0f}", font="8p,Helvetica,black", justify="LM")
    fig.text(x=region[0]+0.015,y=region[3]-0.015,
             text=f"Yogyakarta 2006 aftershocks — size~ML, colour~depth (n={len(e)})",
             font="12p,Helvetica-Bold,black",justify="TL",offset="0.1c/-0.1c",
             fill="white@20",pen="0.5p,black")
    out=f"{ROOT}/figures/aftershock_magnitude_map.png"
    fig.savefig(out,dpi=200); print(f"n={len(e)}; wrote {out}")

if __name__=="__main__":
    main()
