#!/usr/bin/env python3
"""Parse the VELEST-model relocation, compare to the Central-Java NLLoc catalog,
and write catalog_velest.csv + full/velest_reloc_compare.png.
"""
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SUM = f"{ROOT}/nll/loc_v/yogya_v.sum.grid0.loc.hyp"

geo = re.compile(r"^GEOGRAPHIC\s+OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                 r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
qual = re.compile(r"RMS\s+([\d.eE+-]+)\s+Nphs\s+(\d+)\s+Gap\s+([\d.]+)\s+Dist\s+([\d.]+)")
stat = re.compile(r"CovXX\s+(-?[\d.eE+-]+).*?YY\s+(-?[\d.eE+-]+)\s+YZ\s+-?[\d.eE+-]+\s+ZZ\s+(-?[\d.eE+-]+)")
from datetime import datetime, timezone
rows=[]; cur={}
for line in open(SUM):
    if line.startswith("GEOGRAPHIC"):
        m=geo.search(line)
        if m:
            y,mo,d,h,mi=map(int,m.groups()[:5]); sec=float(m.group(6))
            cur={"time":datetime(y,mo,d,h,mi,0,tzinfo=timezone.utc).timestamp()+sec,
                 "latitude":float(m.group(7)),"longitude":float(m.group(8)),"depth":float(m.group(9))}
    elif line.startswith("QUALITY"):
        m=qual.search(line)
        if m: cur.update(rms=float(m.group(1)),nphs=int(m.group(2)),gap=float(m.group(3)))
    elif line.startswith("STATISTICS"):
        m=stat.search(line)
        if m:
            cxx,cyy,czz=map(float,m.groups())
            cur["errh_km"]=round(float(np.sqrt(max(cxx,0)+max(cyy,0))),2)
            cur["errz_km"]=round(float(np.sqrt(max(czz,0))),2)
        rows.append(cur); cur={}
v=pd.DataFrame(rows); v["time_utc"]=pd.to_datetime(v["time"],unit="s",utc=True)
v.to_csv(f"{ROOT}/full/catalog_velest.csv",index=False)
good=v[(v.gap<180)&(v.errh_km<5)&(v.rms<0.5)]
print(f"VELEST reloc: {len(v)} events, {len(good)} well-constrained")
print(f"RMS median {v.rms.median():.3f}s | depth median {v.depth.median():.1f} km "
      f"(good {good.depth.median():.1f}) | errH median {v.errh_km.median():.1f} km")

# compare to Central-Java NLLoc
c=pd.read_csv(f"{ROOT}/full/catalog_nll.csv")
EP=pd.Timestamp("1970-01-01",tz="UTC")
c["time"]=(pd.to_datetime(c["time_utc"],utc=True)-EP).dt.total_seconds()
m=pd.merge_asof(v.sort_values("time"),
                c[["time","latitude","longitude","depth"]].sort_values("time").rename(
                    columns={"latitude":"c_lat","longitude":"c_lon","depth":"c_dep"}),
                on="time",direction="nearest",tolerance=5).dropna(subset=["c_lat"])
dkm=np.sqrt(((m.latitude-m.c_lat)*111)**2+((m.longitude-m.c_lon)*111*np.cos(np.radians(m.latitude)))**2)
ddep=m.depth-m.c_dep
print(f"matched {len(m)} | epicenter shift median {dkm.median():.1f} km | "
      f"depth shift median {ddep.median():+.1f} km (VELEST - CentralJava)")

fig,ax=plt.subplots(2,2,figsize=(14,11))
sc=ax[0,0].scatter(good.longitude,good.latitude,s=5,c=good.depth,cmap="turbo",vmin=0,vmax=20,alpha=0.5)
plt.colorbar(sc,ax=ax[0,0],label="depth km"); ax[0,0].set_aspect("equal","box")
ax[0,0].set_title(f"VELEST-model reloc (n={len(good)})"); ax[0,0].set_xlabel("Lon"); ax[0,0].set_ylabel("Lat")
ax[0,1].hist(v.depth,bins=np.arange(0,26,1),color="darkgreen",alpha=0.6,label=f"VELEST (med {v.depth.median():.1f})")
ax[0,1].hist(c.depth,bins=np.arange(0,26,1),color="gray",alpha=0.5,label=f"CentralJava (med {c.depth.median():.1f})")
ax[0,1].set_xlabel("depth km"); ax[0,1].set_title("Depth: VELEST vs Central-Java model"); ax[0,1].legend()
ax[1,0].hist(v.rms,bins=np.arange(0,0.5,0.02),color="seagreen",edgecolor="k")
ax[1,0].axvline(v.rms.median(),color="k",ls="--"); ax[1,0].set_xlabel("RMS s")
ax[1,0].set_title(f"RMS (median {v.rms.median():.3f} s)")
ax[1,1].hist(ddep,bins=np.arange(-8,8,0.5),color="indianred",edgecolor="k")
ax[1,1].axvline(ddep.median(),color="k",ls="--")
ax[1,1].set_xlabel("depth shift VELEST - CentralJava (km)"); ax[1,1].set_title(f"Depth change (median {ddep.median():+.1f} km)")
plt.tight_layout(); plt.savefig(f"{ROOT}/full/velest_reloc_compare.png",dpi=110)
print("wrote full/catalog_velest.csv and full/velest_reloc_compare.png")
