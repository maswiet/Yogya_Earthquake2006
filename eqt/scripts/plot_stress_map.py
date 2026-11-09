#!/usr/bin/env python3
"""Statistical stress direction (SHmax) from the aftershock moment-tensor catalog,
overlaid on topography. Runs in `gmt` env (pygmt).

Reads  eqt/mt/results/stress_result.json   (from stress_inversion.py)
       eqt/mt/results/mt_catalog.csv        (quality beachballs at locations)
       eqt/data/faults_digitized/faults.csv (Opak/Mataram/Dengkeng traces)
Saves  eqt/figures/stress_map.png

SHmax drawn World-Stress-Map style: oriented bars on a lattice (regime-coloured),
with a 68%-confidence fan; quality beachballs at their epicentres.
"""
import os, json
import numpy as np, pandas as pd
import pygmt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES  = f"{ROOT}/mt/results"
REGION=[110.15,110.75,-8.15,-7.63]; PROJ="M18c"
VRMIN=float(os.environ.get("VRMIN","0.3")); DCMIN=float(os.environ.get("DCMIN","50"))
NUSEDMIN=int(os.environ.get("NUSEDMIN","5"))   # gate on stations actually inverted

FINFO={"Opak":("red",110.345,-7.918,48),"Mataram":("purple",110.388,-7.758,3),
       "Dengkeng":("darkgreen",110.590,-7.772,-13)}

def regime(s1p, s3p):
    """WSM regime from sigma1/sigma3 plunge -> (name, colour)."""
    if s1p >= 52: return "Normal","255/70/70"
    if s3p >= 52: return "Thrust","70/110/255"
    return "Strike-slip","40/170/70"

def load_xn():
    r=[]
    p=f"{ROOT}/config/stations_xn_meta.txt"
    if os.path.exists(p):
        for l in open(p):
            q=l.strip().split("|")
            if len(q)>=6: r.append((float(q[4]),float(q[5])))
    return r

def main():
    S=json.load(open(f"{RES}/stress_result.json"))
    pt=S["ptaxes"]; R=S["R"]
    # robust, plane-independent: SHmax perpendicular to the well-clustered horizontal
    # T-axis (extension). Uncertainty fan = T-axis axial SD, NOT the degenerate formal cone.
    az=pt["SHmax_from_T"]; sd=pt["T_trend_sd"]; azlo,azhi=az-sd, az+sd
    text_az=pt["T_trend"]        # extension azimuth (Shmin)
    s1=S["sigma1"]; s3=S["sigma3"]; reg,regcol=regime(s1["plunge"],s3["plunge"])

    d=pd.read_csv(f"{RES}/mt_catalog.csv"); d=d[d.status=="ok"]
    q=d[(d.VR>=VRMIN)&(d.dc_perc>=DCMIN)&(d.nsta_used>=NUSEDMIN)].copy()
    e=pd.read_csv(f"{ROOT}/full/catalog_growclust.csv").rename(
        columns={"lat":"latitude","lon":"longitude","dep":"depth"})
    f=pd.read_csv(f"{ROOT}/data/faults_digitized/faults.csv")

    fig=pygmt.Figure(); pygmt.config(FONT_TITLE="15p,Helvetica-Bold")
    grid=pygmt.datasets.load_earth_relief(resolution="03s",region=REGION)
    pygmt.makecpt(cmap="gray",series=[-3000,2400],continuous=True)
    fig.grdimage(grid,region=REGION,projection=PROJ,cmap=True,shading="+a315+nt0.7",
        frame=["WSne+tYogyakarta 2006 — statistical SHmax from aftershock moment tensors",
               "xa0.25f0.05","ya0.25f0.05"])
    fig.coast(region=REGION,projection=PROJ,water="170/200/224",shorelines="0.6p,black",resolution="f")

    # faint aftershock cloud for context
    fig.plot(x=e.longitude,y=e.latitude,style="c0.03c",fill="gray35",transparency=55)
    # fault traces
    for nm,(col,lo,la,ang) in FINFO.items():
        g=f[f.fault==nm].sort_values("vertex")
        fig.plot(x=g.longitude,y=g.latitude,pen=f"2.4p,{col}")
        fig.plot(x=g.longitude,y=g.latitude,pen="0.5p,white,-")

    # SHmax bars on a lattice (WSM style), all at regional SHmax with CI fan
    lons=np.linspace(REGION[0]+0.06,REGION[1]-0.06,7)
    lats=np.linspace(REGION[2]+0.05,REGION[3]-0.05,6)
    L=0.028   # half-length in degrees (approx; drawn in geographic)
    def seg(lon,lat,azimuth,half):
        a=np.radians(azimuth); dlat=half*np.cos(a); dlon=half*np.sin(a)/np.cos(np.radians(lat))
        return [lon-dlon,lon+dlon],[lat-dlat,lat+dlat]
    for lo in lons:
        for la in lats:
            # CI fan (light)
            for aa in (azlo,azhi):
                xs,ys=seg(lo,la,aa,L*0.92); fig.plot(x=xs,y=ys,pen=f"1.2p,{regcol}@70")
            xs,ys=seg(lo,la,az,L); fig.plot(x=xs,y=ys,pen=f"3.2p,{regcol}")

    # quality beachballs at epicentres
    for _,r in q.iterrows():
        fig.meca(spec={"strike":r.s1,"dip":r.d1,"rake":r.r1,"magnitude":r.get("Mw",3.0)},
                 scale="0.55c",longitude=r.lon,latitude=r.lat,depth=float(r.z)/1000 if r.get("z") else 10,
                 compressionfill="black",extensionfill="white",pen="0.4p,black")

    xn=load_xn()
    if xn: fig.plot(x=[s[1] for s in xn],y=[s[0] for s in xn],style="i0.28c",fill="white",pen="0.7p,black")

    for nm,(col,lo,la,ang) in FINFO.items():
        fig.text(x=lo,y=la,text=f"{nm} Fault",angle=ang,
                 font=f"9p,Helvetica-BoldOblique,{col}",fill="white@25",pen=f"0.3p,{col}")

    # info box
    txt=(f"SHmax = N{az:.0f}\\260E  (\\261{sd:.0f}\\260)    extension (Shmin) = N{text_az:.0f}\\260E    "
         f"regime: {reg}-dominant transtension    R = {R:.2f}")
    fig.text(position="TC",text=txt,font="10p,Helvetica-Bold,black",
             fill="white@20",pen="0.6p,black",offset="0c/-0.55c",no_clip=True)
    fig.text(position="BL",offset="0.3c/0.3c",fill="white@25",pen="0.4p,black",
             text=f"@;{regcol};|@;; SHmax bar (\\136 horizontal T-axis)   "
                  f"n={len(q)} MT  (>={NUSEDMIN} sta used, VR>={VRMIN}, DC>={DCMIN}%)",
             font="8p,Helvetica,black")
    fig.basemap(map_scale="jBR+w20k+o0.6c/0.6c+f+u")
    out=f"{ROOT}/figures/stress_map.png"
    fig.savefig(out,dpi=230); print(f"n_MT={len(q)} SHmax={az:.0f} R={R:.2f} regime={reg}; wrote {out}")

if __name__=="__main__":
    main()
