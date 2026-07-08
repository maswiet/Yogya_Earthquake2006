#!/usr/bin/env python3
"""GrowClust aftershock catalog on hillshade + REAL geology (Rahardjo et al. 1995
shapefile) + fault structures (after Ramdhan et al. 2025). Runs in the `gmt` env.
Saves eqt/figures/growclust_geology_map.png
"""
import os
import numpy as np, pandas as pd, geopandas as gpd
import pygmt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SHP = os.path.join(ROOT,"data/geology_yogyakarta/Geology Yogyakarta.shp")
if not os.path.exists(SHP): SHP = "/Users/maswiet/Downloads/14. Geology Yogyakarta/Geology Yogyakarta.shp"
REGION = [110.15, 110.75, -8.15, -7.63]; PROJ = "M17c"

# map symbol -> paper-style age/lithology group
SYM2GRP = {"Qa":"alluvium","Qc":"alluvium","Qt3":"alluvium",
    "Qmi":"volccover","Qvu3":"volccover","Tpdi":"diorite","a1":"diorite",
    "Tmwl":"limestone","Tmj":"limestone","Tmps":"limestone","Tmpk2":"limestone",
    "Tmng":"lowmioc","Tms3":"lowmioc","Tmse":"lowmioc","Tmss":"lowmioc","Tmo":"lowmioc",
    "Tomk":"oligomioc","Tomm3":"oligomioc","Teon":"eocene","KTm1":"metamorph"}
GRP = [("alluvium","200/200/200","Quaternary alluvium"),
       ("volccover","224/208/140","Quaternary volcanic cover"),
       ("oligomioc","216/116/70","Oligo-Miocene volcanics"),
       ("lowmioc","178/142/96","Lower-Mid Miocene volcaniclastics"),
       ("limestone","150/200/226","Mio-Pliocene limestone/marl"),
       ("diorite","220/70/55","Diorite / andesite intrusion"),
       ("eocene","176/210/140","Eocene siliciclastics"),
       ("metamorph","182/112/192","Pre-Tertiary metamorphics")]
GCOLOR = {k:c for k,c,_ in GRP}

FAULTS = {
    "Opak Fault":    [(110.400,-8.030),(110.430,-7.900),(110.455,-7.780),(110.470,-7.700)],
    "Oyo Fault":     [(110.428,-7.985),(110.470,-7.920),(110.495,-7.875)],
    "Ngalang Fault": [(110.460,-7.990),(110.500,-7.855),(110.535,-7.715)],
    "Nglipar Fault": [(110.550,-7.895),(110.615,-7.800),(110.690,-7.730)]}
FLABEL = {"Opak Fault":(110.447,-7.815,68),"Oyo Fault":(110.463,-7.930,52),
          "Ngalang Fault":(110.512,-7.835,66),"Nglipar Fault":(110.628,-7.795,47)}
MURIA = [(110.190,-7.660),(110.245,-7.820),(110.310,-8.060)]
CITIES=[("Yogyakarta",110.366,-7.797),("Bantul",110.328,-7.888),("Sleman",110.353,-7.716),
        ("Klaten",110.606,-7.703),("Kulon Progo",110.175,-7.840),("Gunungkidul",110.601,-7.966)]
MAIN=[("BMKG",110.330,-8.050,"yellow"),("GFZ",110.420,-8.005,"red"),
      ("USGS",110.468,-7.960,"gold"),("GCMT",110.555,-8.010,"steelblue")]

def load_xn():
    r=[]
    for l in open(f"{ROOT}/config/stations_xn_meta.txt"):
        p=l.strip().split("|")
        if len(p)>=6: r.append((p[1],float(p[4]),float(p[5])))
    return r

def main():
    g = gpd.read_file(SHP)
    import shapely.geometry as sg
    g = gpd.clip(g, sg.box(*[REGION[0],REGION[2],REGION[1],REGION[3]]))
    g["grp"] = g["SYMBOLS"].map(SYM2GRP)
    e = pd.read_csv(f"{ROOT}/full/catalog_growclust.csv").rename(columns={"lat":"latitude","lon":"longitude","dep":"depth"})

    fig = pygmt.Figure()
    grid = pygmt.datasets.load_earth_relief(resolution="03s", region=REGION)
    pygmt.makecpt(cmap="gray", series=[-3000,2400], continuous=True)
    fig.grdimage(grid, region=REGION, projection=PROJ, cmap=True,
                 shading="+a315+nt0.7", frame=["WSne","xa0.25f0.05","ya0.25f0.05"])
    # geology polygons (semi-transparent so hillshade + events show)
    for k,color,_ in GRP:
        sub = g[g.grp==k]
        if len(sub): fig.plot(data=sub, fill=color, transparency=42, pen="0.2p,80/80/80@50")
    fig.coast(region=REGION, projection=PROJ, water="170/200/224",
              shorelines="0.7p,black", resolution="f")
    # faults
    for nm,pts in FAULTS.items(): fig.plot(data=[[x,y] for x,y in pts], pen="2.2p,black")
    fig.plot(data=[[x,y] for x,y in MURIA], pen="1.6p,black,--")
    # aftershocks
    pygmt.makecpt(cmap="jet", series=[0,25], reverse=True)
    fig.plot(x=e.longitude, y=e.latitude, fill=e.depth, cmap=True, style="c0.05c",
             pen="0.04p,gray30", transparency=8)
    xn=load_xn()
    fig.plot(x=[s[2] for s in xn], y=[s[1] for s in xn], style="i0.30c", fill="white", pen="0.8p,black")
    for nm,lo,la,c in MAIN:
        fig.plot(x=[lo],y=[la],style="a0.55c",fill=c,pen="1p,black")
        fig.text(x=lo,y=la,text=nm,font="8p,Helvetica-Bold,black",justify="TL",offset="0.16c/-0.04c",fill="white@40")
    for nm,(lo,la,ang) in FLABEL.items():
        fig.text(x=lo,y=la,text=nm,angle=ang,font="10p,Helvetica-BoldOblique,black",fill="white@55",pen="0.2p")
    for nm,lo,la in CITIES:
        fig.plot(x=[lo],y=[la],style="s0.26c",fill="firebrick",pen="0.6p,black")
        fig.text(x=lo,y=la,text=nm,font="9p,Helvetica-Bold,20/20/90",justify="LM",offset="0.2c/0c",fill="white@45")
    # depth colorbar
    fig.colorbar(cmap=True, frame=["x+lFocal depth","y+lkm"], position="JMR+o0.7c/2.6c+w7c/0.35c")
    # geology legend (panel over the SW sea)
    fig.plot(data=[[110.156,-7.858],[110.388,-7.858],[110.388,-8.075],[110.156,-8.075]],
             fill="white@12", close=True, pen="0.6p,black")
    x0,ytop=110.172,-7.885; dy=0.0262
    for i,(k,color,lab) in enumerate(GRP):
        yy=ytop-i*dy
        fig.plot(x=[x0],y=[yy],style="s0.26c",fill=color,pen="0.4p,black")
        fig.text(x=x0+0.013,y=yy,text=lab,font="7.3p,Helvetica,black",justify="LM")
    fig.basemap(map_scale="jBR+w20k+o0.6c/0.6c+f+u")
    fig.text(x=REGION[0]+0.008,y=REGION[3]-0.008,
             text=f"GrowClust aftershocks (n={len(e)}) \\267 geology: Rahardjo 1995 \\267 faults: Ramdhan 2025",
             font="9.5p,Helvetica-Bold,black",justify="TL",offset="0.12c/-0.12c",fill="white@20",pen="0.5p")
    out=f"{ROOT}/figures/growclust_geology_map.png"
    fig.savefig(out,dpi=230); print(f"n={len(e)}, geology polys={len(g)}; wrote {out}")

if __name__=="__main__":
    main()
