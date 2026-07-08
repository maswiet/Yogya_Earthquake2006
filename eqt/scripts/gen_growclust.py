#!/usr/bin/env python3
"""Build GrowClust inputs from the full refined (loc_v) catalog:
  growclust/IN/evlist.txt   - event list (yr mo dy hr mi sec lat lon dep mag eh ez rms id)
  growclust/IN/stlist.txt   - station list (code lat lon)
  growclust/IN/vzmodel.txt  - 1-D model (refined VELEST), piecewise constant
  growclust/hypoDD.pha      - phase file for ph2dt (-> dt.ct -> xcordata)
  growclust/ph2dt.inp       - ph2dt control
  growclust/growclust.inp   - GrowClust control
Selects all well-associated events (gap<gapmax, nphs>=minphs).
"""
import os, re, glob, json, string, math
import numpy as np, pandas as pd
from obspy import UTCDateTime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GC = f"{ROOT}/growclust"
for d in ("IN", "OUT", "TT"): os.makedirs(f"{GC}/{d}", exist_ok=True)
VPVS = 1.735
# refined VELEST model layers (top depth, Vp) -> piecewise constant to 26 km
LAYERS = [(0.0,2.90),(0.7,4.30),(2.0,4.65),(4.0,5.49),(10.0,6.30),
          (13.0,6.39),(16.0,6.55),(22.0,6.80),(26.0,7.20),(30.0,7.20),(40.0,8.00)]
GAPMAX=float(os.environ.get("GC_GAP","200")); NPHMIN=int(os.environ.get("GC_NPH","6"))
NMAX=int(os.environ.get("GC_NMAX","16500"))

def parse_events():
    geo=re.compile(r"OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                   r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
    qual=re.compile(r"RMS\s+([\d.eE+-]+)\s+Nphs\s+(\d+)\s+Gap\s+([\d.]+)")
    stat=re.compile(r"CovXX\s+(-?[\d.eE+-]+).*?YY\s+(-?[\d.eE+-]+)\s+YZ\s+-?[\d.eE+-]+\s+ZZ\s+(-?[\d.eE+-]+)")
    evs=[]
    for path in sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp")):
        ev={"phs":[]}; ot=None
        for line in open(path):
            if line.startswith("GEOGRAPHIC"):
                m=geo.search(line)
                if m:
                    y,mo,d,h,mi=map(int,m.groups()[:5]); s=float(m.group(6))
                    ot=UTCDateTime(f"{y:04d}-{mo:02d}-{d:02d}T{h:02d}:{mi:02d}:00")+s
                    ev.update(ot=ot,lat=float(m.group(7)),lon=float(m.group(8)),dep=float(m.group(9)))
            elif line.startswith("QUALITY"):
                mm=qual.search(line)
                if mm: ev["rms"]=float(mm.group(1)); ev["nphs"]=int(mm.group(2)); ev["gap"]=float(mm.group(3))
            elif line.startswith("STATISTICS"):
                ms=stat.search(line)
                if ms:
                    cxx,cyy,czz=map(float,ms.groups())
                    ev["eh"]=round(math.sqrt(max(cxx,0)+max(cyy,0)),2); ev["ez"]=round(math.sqrt(max(czz,0)),2)
            elif " > " in line and "GAU" in line:
                lf=line.split(" > ")[0].split(); rf=line.split(" > ")[1].split()
                tt=float(rf[0])+float(rf[1])
                if tt>0: ev["phs"].append((lf[0],lf[4],tt))
        if ot and ev["phs"]: evs.append(ev)
    return evs

def main():
    evs=parse_events()
    mag=pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    mt=pd.to_datetime(mag.time,utc=True).astype("int64").to_numpy()/1e9; mlv=mag.ML.to_numpy()
    def ml_of(ot):
        i=np.argmin(np.abs(mt-ot.timestamp)); return mlv[i] if abs(mt[i]-ot.timestamp)<2 else 0.0
    sel=[e for e in evs if e.get("gap",999)<GAPMAX and e.get("nphs",0)>=NPHMIN and e["dep"]<26][:NMAX]
    print(f"{len(evs)} located -> {len(sel)} for GrowClust")

    # evlist + pha (shared sequential IDs)
    fev=open(f"{GC}/IN/evlist.txt","w"); fph=open(f"{GC}/hypoDD.pha","w")
    for eid,e in enumerate(sel,1):
        ot=e["ot"]; ml=ml_of(ot)
        fev.write(f"{ot.year:4d} {ot.month:2d} {ot.day:2d} {ot.hour:2d} {ot.minute:2d} "
                  f"{ot.second+ot.microsecond*1e-6:6.3f} {e['lat']:9.5f} {e['lon']:10.5f} "
                  f"{e['dep']:8.3f} {ml:6.3f} {e.get('eh',0):6.3f} {e.get('ez',0):6.3f} "
                  f"{e.get('rms',0):6.3f} {eid:9d}\n")
        fph.write(f"# {ot.year} {ot.month} {ot.day} {ot.hour} {ot.minute} "
                  f"{ot.second+ot.microsecond*1e-6:5.2f} {e['lat']:.4f} {e['lon']:.4f} "
                  f"{e['dep']:.2f} {ml:.1f} {e.get('eh',0):.2f} {e.get('ez',0):.2f} {e.get('rms',0):.2f} {eid}\n")
        for sta,ph,tt in e["phs"]:
            fph.write(f"{sta:<7}{tt:7.3f} {1.0 if ph=='P' else 0.5:.2f} {ph}\n")
    fev.close(); fph.close()

    # station list
    per=json.load(open(f"{ROOT}/config/stations_periods.json"))
    with open(f"{GC}/IN/stlist.txt","w") as f:
        for code,info in per.items():
            sites=info.get("sites",[])
            if len(sites)<=1:
                if sites: f.write(f"{code:<6}{sites[0]['lat']:9.4f} {sites[0]['lon']:10.4f}\n")
            else:
                for suf,s in zip(string.ascii_lowercase,sites):
                    f.write(f"{code+suf:<6}{s['lat']:9.4f} {s['lon']:10.4f}\n")

    # vzmodel (piecewise constant depth Vp Vs=0 -> use vpvs)
    with open(f"{GC}/IN/vzmodel.txt","w") as f:
        for i in range(len(LAYERS)-1):
            d0,v0=LAYERS[i]; d1,_=LAYERS[i+1]
            f.write(f"{d0:6.2f} {v0:5.2f} 0.00\n{d1:6.2f} {v0:5.2f} 0.00\n")

    # ph2dt.inp
    open(f"{GC}/ph2dt.inp","w").write(
        "* ph2dt.inp\nIN/stlist.txt\nhypoDD.pha\n"
        "*MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS\n"
        "   0      120     10     40     8      8     80\n")

    # growclust.inp
    open(f"{GC}/growclust.inp","w").write(f"""* GrowClust control - EQT Yogya 2006 (refined VELEST model)
* evlist_fmt
1
* fin_evlist
IN/evlist.txt
* stlist_fmt
1
* fin_stlist
IN/stlist.txt
* xcordat_fmt  tdif_fmt
1  12
* fin_xcordat
IN/xcordata.txt
* fin_vzmdl
IN/vzmodel.txt
* fout_vzfine
TT/vzfine.txt
* fout_pTT
TT/tt.pg
* fout_sTT
TT/tt.sg
* vpvs_factor  rayparam_min
  {VPVS}       0.0
* tt_dep0 tt_dep1 tt_ddep
  0.  40.  1.
* tt_del0 tt_del1 tt_ddel
  0.  200.  2.
* rmin delmax rmsmax
  0.0  80    0.30
* rpsavgmin rmincut ngoodmin iponly
  0  0  0  0
* nboot nbranch_min
  100  1
* fout_cat
OUT/out.growclust_cat
* fout_clust
OUT/out.growclust_clust
* fout_log
OUT/out.growclust_log
* fout_boot
OUT/out.growclust_boot
""")
    print(f"wrote GrowClust inputs in {GC}")

if __name__=="__main__":
    main()
