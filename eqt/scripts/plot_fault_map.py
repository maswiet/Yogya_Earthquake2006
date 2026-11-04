#!/usr/bin/env python3
"""GrowClust aftershock catalog on hillshade + USER-DIGITIZED fault traces
(Opak, Mataram, Dengkeng — eqt/data/faults_digitized/faults.csv). Runs in `gmt`.
Saves eqt/figures/growclust_fault_map.png
"""
import os
import numpy as np, pandas as pd
import pygmt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FCSV = os.path.join(ROOT, "data/faults_digitized/faults.csv")
REGION = [110.15, 110.75, -8.15, -7.63]; PROJ = "M17c"

# fault -> (colour, label, label lon, lat, angle)
FINFO = {
    "Opak":    ("red",      "Opak Fault",    110.345, -7.918, 48),
    "Mataram": ("purple",   "Mataram Fault", 110.388, -7.758,  3),
    "Dengkeng":("darkgreen","Dengkeng Fault",110.590, -7.772,-13),
}
CITIES = [("Yogyakarta",110.366,-7.797),("Bantul",110.328,-7.888),("Sleman",110.353,-7.716),
          ("Klaten",110.606,-7.703),("Kulon Progo",110.175,-7.840),("Gunungkidul",110.601,-7.966)]
MECA=[("GFZ", 110.420,-8.005, 110.470,-8.055, "red",       6.3, (229,85,-9)),
      ("USGS",110.468,-7.960, 110.580,-8.055, "gold",      6.5, (231,87,3)),
      ("GCMT",110.555,-8.010, 110.688,-8.055, "steelblue", 6.4, (232,86,-13))]

def load_xn():
    r=[]
    for l in open(f"{ROOT}/config/stations_xn_meta.txt"):
        p=l.strip().split("|")
        if len(p)>=6: r.append((p[1],float(p[4]),float(p[5])))
    return r

def main():
    f = pd.read_csv(FCSV)
    e = pd.read_csv(f"{ROOT}/full/catalog_growclust.csv").rename(
            columns={"lat":"latitude","lon":"longitude","dep":"depth"})

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="15p,Helvetica-Bold")
    grid = pygmt.datasets.load_earth_relief(resolution="03s", region=REGION)
    pygmt.makecpt(cmap="gray", series=[-3000,2400], continuous=True)
    fig.grdimage(grid, region=REGION, projection=PROJ, cmap=True, shading="+a315+nt0.7",
                 frame=["WSne+tYogyakarta 2006 aftershocks (GrowClust) + digitized faults "
                        f"(n={len(e)})", "xa0.25f0.05","ya0.25f0.05"])
    fig.coast(region=REGION, projection=PROJ, water="170/200/224",
              shorelines="0.7p,black", resolution="f")

    # aftershocks coloured by focal depth
    pygmt.makecpt(cmap="jet", series=[0,25], reverse=True)
    fig.plot(x=e.longitude, y=e.latitude, fill=e.depth, cmap=True, style="c0.05c",
             pen="0.04p,gray30", transparency=12)

    # digitized fault traces (bold, white dash down the centre for contrast)
    for name,(col,lab,lo,la,ang) in FINFO.items():
        g = f[f.fault==name].sort_values("vertex")
        fig.plot(x=g.longitude, y=g.latitude, pen=f"2.8p,{col}")
        fig.plot(x=g.longitude, y=g.latitude, pen="0.6p,white,-")

    xn=load_xn()
    fig.plot(x=[s[2] for s in xn], y=[s[1] for s in xn], style="i0.30c", fill="white", pen="0.8p,black")

    for nm,slo,sla,plo,pla,col,mw,(st,dp,rk) in MECA:
        fig.meca(spec={"strike":st,"dip":dp,"rake":rk,"magnitude":mw}, scale="1.0c",
                 longitude=slo, latitude=sla, depth=12, plot_longitude=plo, plot_latitude=pla,
                 offset="0.6p,gray20", compressionfill=col, extensionfill="white", pen="0.6p,black")
        fig.text(x=plo, y=pla-0.052, text=nm, font="7.5p,Helvetica-Bold,black", justify="MC", fill="white@35")

    for nm,lo,la in CITIES:
        fig.plot(x=[lo],y=[la],style="s0.26c",fill="firebrick",pen="0.6p,black")
        fig.text(x=lo,y=la,text=nm,font="9p,Helvetica-Bold,20/20/90",justify="LM",offset="0.2c/0c",fill="white@45")

    for name,(col,lab,lo,la,ang) in FINFO.items():
        fig.text(x=lo,y=la,text=lab,angle=ang,font=f"10p,Helvetica-BoldOblique,{col}",
                 fill="white@25",pen=f"0.3p,{col}")

    fig.colorbar(cmap=True, frame=["x+lFocal depth","y+lkm"], position="JMR+o0.7c/2.6c+w7c/0.35c")
    fig.basemap(map_scale="jTL+w20k+o0.6c/0.6c+f+u")
    fig.text(position="BC",
             text="fault traces: author digitization (Opak / Mataram / Dengkeng)  \\267  "
                  "relocation: GrowClust (Trugman & Shearer 2017)",
             font="8p,Helvetica-Oblique,gray20", offset="0c/-1.0c", no_clip=True)
    out=f"{ROOT}/figures/growclust_fault_map.png"
    fig.savefig(out,dpi=230); print(f"n={len(e)}; faults={f.fault.nunique()}; wrote {out}")

if __name__=="__main__":
    main()
