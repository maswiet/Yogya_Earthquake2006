#!/usr/bin/env python3
"""Plot good-quality BayesISOLA focal mechanisms over the GrowClust catalog +
digitized faults. Runs in `gmt` env. Quality filter: VR>=VRMIN, DC%>=DCMIN, nsta>=6.
Saves eqt/figures/mt_mechanisms_map.png
"""
import os, numpy as np, pandas as pd, pygmt
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
REGION=[110.15,110.75,-8.15,-7.63]; PROJ="M17c"
VRMIN=float(os.environ.get("VRMIN","0.6")); DCMIN=float(os.environ.get("DCMIN","50"))

FINFO={"Opak":("red","Opak Fault",110.345,-7.918,48),
       "Mataram":("purple","Mataram Fault",110.388,-7.758,3),
       "Dengkeng":("darkgreen","Dengkeng Fault",110.590,-7.772,-13)}

def main():
    e=pd.read_csv(f"{ROOT}/full/catalog_growclust.csv").rename(columns={"lat":"latitude","lon":"longitude","dep":"depth"})
    mt=pd.read_csv(f"{ROOT}/mt/results/mt_catalog.csv")
    mt=mt[(mt.status=="ok")&(mt.get('VR',0)>=VRMIN)&(mt.get('dc_perc',0)>=DCMIN)&(mt.nsta>=6)].copy()
    f=pd.read_csv(f"{ROOT}/data/faults_digitized/faults.csv")
    fig=pygmt.Figure(); pygmt.config(FONT_TITLE="15p,Helvetica-Bold")
    grid=pygmt.datasets.load_earth_relief(resolution="03s",region=REGION)
    pygmt.makecpt(cmap="gray",series=[-3000,2400],continuous=True)
    fig.grdimage(grid,region=REGION,projection=PROJ,cmap=True,shading="+a315+nt0.7",
                 frame=["WSne+tYogyakarta 2006 aftershock focal mechanisms (BayesISOLA MT, VR>="+f"{VRMIN})","xa0.25f0.05","ya0.25f0.05"])
    fig.coast(region=REGION,projection=PROJ,water="170/200/224",shorelines="0.7p,black",resolution="f")
    pygmt.makecpt(cmap="jet",series=[0,25],reverse=True)
    fig.plot(x=e.longitude,y=e.latitude,fill="gray70",style="c0.03c",transparency=55)
    for name,(col,lab,lo,la,ang) in FINFO.items():
        g=f[f.fault==name].sort_values("vertex")
        fig.plot(x=g.longitude,y=g.latitude,pen=f"2.4p,{col}")
        fig.text(x=lo,y=la,text=lab,angle=ang,font=f"9p,Helvetica-BoldOblique,{col}",fill="white@30")
    # beachballs coloured by centroid depth
    pygmt.makecpt(cmap="jet",series=[0,25],reverse=True)
    for _,r in mt.iterrows():
        try:
            spec=dict(strike=r.s1,dip=r.d1,rake=r.r1,magnitude=r.Mw)
        except Exception: continue
        fig.meca(spec=spec, scale="0.9c", longitude=r.lon, latitude=r.lat, depth=r.get('z',10),
                 compressionfill="black", extensionfill="white", pen="0.5p,black")
    fig.colorbar(cmap=True,frame=["x+lDepth","y+lkm"],position="JMR+o0.7c/2.6c+w7c/0.35c")
    fig.basemap(map_scale="jTL+w20k+o0.6c/0.6c+f+u")
    fig.text(position="BC",text=f"{len(mt)} well-resolved mechanisms (VR>={VRMIN}, DC>={DCMIN:.0f}%, nsta>=6)  \\267  BayesISOLA (Vackar et al. 2017)",
             font="8p,Helvetica-Oblique,gray20",offset="0c/-0.9c",no_clip=True)
    out=f"{ROOT}/figures/mt_mechanisms_map.png"; fig.savefig(out,dpi=230)
    print(f"{len(mt)} mechanisms plotted -> {out}")

if __name__=="__main__": main()
