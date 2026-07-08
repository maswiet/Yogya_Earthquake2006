#!/usr/bin/env python3
"""GrowClust aftershock catalog on a hillshade + faults/structures map,
styled after Ramdhan et al. (2025) Fig. 1. Runs in the `gmt` env.
Faults/lineament digitized approximately from that figure; geology shown as
generalized domains (limestone Southern Mts vs Bantul-basin sediment) — the full
Rahardjo et al. (1995) lithology needs the source shapefile.
Saves eqt/figures/growclust_tectonic_map.png
"""
import os
import numpy as np, pandas as pd
import pygmt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REGION = [110.15, 110.75, -8.15, -7.63]
PROJ = "M17c"

# --- fault traces (lon,lat) approx from Ramdhan et al. 2025 Fig.1 ---
FAULTS = {
    "Opak Fault":    [(110.400,-8.030),(110.430,-7.900),(110.455,-7.780),(110.470,-7.700)],
    "Oyo Fault":     [(110.428,-7.985),(110.470,-7.920),(110.495,-7.875)],
    "Ngalang Fault": [(110.460,-7.990),(110.500,-7.855),(110.535,-7.715)],
    "Nglipar Fault": [(110.550,-7.895),(110.615,-7.800),(110.690,-7.730)],
}
FAULT_LABEL = {  # (lon,lat,angle)
    "Opak Fault":    (110.447,-7.815, 68),
    "Oyo Fault":     (110.463,-7.930, 52),
    "Ngalang Fault": (110.512,-7.835, 66),
    "Nglipar Fault": (110.628,-7.795, 47),
}
MURIA_PROGO = [(110.190,-7.660),(110.245,-7.820),(110.310,-8.060)]

# generalized geology domains (approx; after Rahardjo et al. 1995)
LIMESTONE = [(110.44,-8.15),(110.75,-8.15),(110.75,-7.80),(110.60,-7.83),
             (110.52,-7.92),(110.46,-8.00),(110.44,-8.05)]           # Southern Mts limestone (SE)
BASIN     = [(110.30,-8.02),(110.44,-8.00),(110.44,-7.74),(110.30,-7.74)]  # Bantul basin sediment (W)

CITIES = [("Yogyakarta",110.366,-7.797),("Bantul",110.328,-7.888),
          ("Sleman",110.353,-7.716),("Klaten",110.606,-7.703),
          ("Kulon Progo",110.175,-7.840),("Gunungkidul",110.601,-7.966)]
# 2006 mainshock, four agencies (approx positions from the figure)
MAINSHOCK = [("BMKG",110.330,-8.050,"yellow"),("GFZ",110.420,-8.005,"red"),
             ("USGS",110.468,-7.960,"gold"),("GCMT",110.555,-8.010,"steelblue")]

def load_xn():
    r=[]
    for l in open(f"{ROOT}/config/stations_xn_meta.txt"):
        p=l.strip().split("|")
        if len(p)>=6: r.append((p[1],float(p[4]),float(p[5])))
    return r

def main():
    e = pd.read_csv(f"{ROOT}/full/catalog_growclust.csv").rename(columns={"lat":"latitude","lon":"longitude","dep":"depth"})
    fig = pygmt.Figure()
    grid = pygmt.datasets.load_earth_relief(resolution="03s", region=REGION)
    pygmt.makecpt(cmap="gray", series=[-2000,1900], continuous=True)
    fig.grdimage(grid, region=REGION, projection=PROJ, cmap=True,
                 shading="+a315+nt0.6", frame=["WSne","xa0.25f0.05","ya0.25f0.05"])
    fig.coast(region=REGION, projection=PROJ, water="176/206/228",
              shorelines="0.7p,black", resolution="f")
    # generalized geology domains (semi-transparent)
    fig.plot(data=[[x,y] for x,y in LIMESTONE], projection=PROJ, region=REGION,
             fill="135/206/235@70", close=True, pen="0.3p,70/130/180@40")
    # faults
    for nm,pts in FAULTS.items():
        fig.plot(data=[[x,y] for x,y in pts], pen="2.2p,black")
    fig.plot(data=[[x,y] for x,y in MURIA_PROGO], pen="1.6p,black,--")
    # aftershocks coloured by focal depth (red shallow -> blue deep, like the paper)
    pygmt.makecpt(cmap="jet", series=[0,25], reverse=True)
    fig.plot(x=e.longitude, y=e.latitude, fill=e.depth, cmap=True,
             style="c0.055c", pen="0.05p,gray30", transparency=10)
    # XN stations
    xn=load_xn()
    fig.plot(x=[s[2] for s in xn], y=[s[1] for s in xn], style="i0.32c",
             fill="yellow", pen="0.8p,black")
    # mainshock stars
    for nm,lo,la,c in MAINSHOCK:
        fig.plot(x=[lo], y=[la], style="a0.6c", fill=c, pen="1p,black")
        fig.text(x=lo, y=la, text=nm, font="8p,Helvetica-Bold,black",
                 justify="TL", offset="0.18c/-0.05c", fill="white@40")
    # fault labels
    for nm,(lo,la,ang) in FAULT_LABEL.items():
        fig.text(x=lo, y=la, text=nm, angle=ang, font="11p,Helvetica-BoldOblique,black",
                 fill="white@50", pen="0.2p")
    fig.text(x=110.44,y=-7.66,text="Muria–Progo Lineament", angle=-60,
             font="9p,Helvetica-Oblique,black", justify="ML")
    # cities
    for nm,lo,la in CITIES:
        fig.plot(x=[lo],y=[la],style="s0.28c",fill="firebrick",pen="0.6p,black")
        fig.text(x=lo,y=la,text=nm,font="9p,Helvetica-Bold,20/20/90",justify="LM",
                 offset="0.22c/0c",fill="white@50")
    # domain labels
    fig.text(x=110.63,y=-8.06,text="Southern Mts (limestone)",font="8p,Helvetica-Oblique,30/60/110",justify="MC")
    # colorbar + scale
    fig.colorbar(cmap=True, frame=["x+lFocal depth","y+lkm"], position="JMR+o0.7c/0c+w9c/0.4c")
    fig.basemap(map_scale="jBL+w20k+o0.6c/0.6c+f+u")
    fig.text(x=REGION[0]+0.01,y=REGION[3]-0.01,
             text=f"Yogyakarta 2006 aftershocks (GrowClust, n={len(e)}) | faults after Ramdhan et al. 2025",
             font="10p,Helvetica-Bold,black",justify="TL",offset="0.15c/-0.15c",fill="white@25",pen="0.5p")
    out=f"{ROOT}/figures/growclust_tectonic_map.png"
    fig.savefig(out,dpi=220); print(f"n={len(e)}; wrote {out}")

if __name__=="__main__":
    main()
