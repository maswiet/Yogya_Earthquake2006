#!/usr/bin/env python3
"""Build HypoDD inputs from the refined NLLoc (loc_v) catalog:
  hypodd/hypoDD.pha  - phase file (event headers + P/S travel times) for ph2dt
  hypodd/station.dat - station list
  hypodd/ph2dt.inp, hypodd/hypoDD.inp - control files (refined VELEST 1-D model)

Selects a high-quality subset (gap<160, rms<0.30, nphs>=10, main cluster) up to
MAXEVE=6500.
"""
import os, re, glob, json, string, math
import numpy as np, pandas as pd
from obspy import UTCDateTime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HD = f"{ROOT}/hypodd"; os.makedirs(HD, exist_ok=True)
CLAT, CLON = -7.92, 110.44
# refined VELEST P model (top depth km, Vp) and Vp/Vs
MTOP = [0.0, 0.7, 2.0, 4.0, 7.0, 10.0, 13.0, 16.0, 22.0, 30.0]
MVEL = [2.90, 4.30, 4.65, 5.49, 5.49, 6.30, 6.39, 6.55, 6.80, 7.20]
VPVS = 1.735

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
                sta,ph=lf[0],lf[4]; tt=float(rf[0])+float(rf[1])
                if tt>0: ev["phs"].append((sta,ph,tt))
        if ot and ev["phs"]: evs.append(ev)
    return evs

def main():
    evs=parse_events()
    # ML by time
    mag=pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    mt=pd.to_datetime(mag.time,utc=True).astype("int64").to_numpy()/1e9
    mlv=mag.ML.to_numpy()
    def ml_of(ot):
        i=np.argmin(np.abs(mt-ot.timestamp))
        return mlv[i] if abs(mt[i]-ot.timestamp)<2 else 0.0
    # select high quality, main cluster
    sel=[]
    for e in evs:
        d=math.hypot(e["lat"]-CLAT, e["lon"]-CLON)
        if e.get("gap",999)<150 and e.get("rms",9)<0.22 and e.get("nphs",0)>=12 and d<0.30 and e["dep"]<25:
            sel.append(e)
    sel=sel[:3500]
    print(f"{len(evs)} events -> {len(sel)} selected for HypoDD")

    # station.dat (per-period codes)
    per=json.load(open(f"{ROOT}/config/stations_periods.json"))
    with open(f"{HD}/station.dat","w") as f:
        for code,info in per.items():
            sites=info.get("sites",[])
            if len(sites)<=1:
                if sites: f.write(f"{code:<7}{sites[0]['lat']:.4f} {sites[0]['lon']:.4f}\n")
            else:
                for suf,s in zip(string.ascii_lowercase,sites):
                    f.write(f"{code+suf:<7}{s['lat']:.4f} {s['lon']:.4f}\n")

    # hypoDD.pha
    with open(f"{HD}/hypoDD.pha","w") as f:
        for eid,e in enumerate(sel,1):
            ot=e["ot"]; ml=ml_of(ot)
            f.write(f"# {ot.year} {ot.month} {ot.day} {ot.hour} {ot.minute} "
                    f"{ot.second+ot.microsecond*1e-6:5.2f} {e['lat']:.4f} {e['lon']:.4f} "
                    f"{e['dep']:.2f} {ml:.1f} {e.get('eh',0):.2f} {e.get('ez',0):.2f} "
                    f"{e.get('rms',0):.2f} {eid}\n")
            for sta,ph,tt in e["phs"]:
                wt = 1.0 if ph=="P" else 0.5
                f.write(f"{sta:<7}{tt:7.3f} {wt:.2f} {ph}\n")

    # ph2dt.inp
    open(f"{HD}/ph2dt.inp","w").write(
        "* ph2dt.inp\nstation.dat\nhypoDD.pha\n"
        "*MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS\n"
        "   0      120     10     50     8      8     100\n")

    # hypoDD.inp with refined 1-D model (exact format: ONLY the dt.cc line is blank)
    top=" ".join(f"{d:.2f}" for d in MTOP); vel=" ".join(f"{v:.2f}" for v in MVEL)
    open(f"{HD}/hypoDD.inp","w").write(f"""* hypoDD.inp (EQT Yogya 2006, refined VELEST model)
*--- input file selection
* cross correlation diff times:

*
* catalog P diff times:
dt.ct
*
* event file:
event.sel
*
* station file:
station.dat
*
*--- output file selection
* original locations:
hypoDD.loc
* relocations:
hypoDD.reloc
* station information:
hypoDD.sta
* residual information:
hypoDD.res
* source parameter information:
hypoDD.src
*
*--- data type selection:
* IDAT IPHA DIST
   2     3    120
*
*--- event clustering:
* OBSCC OBSCT
   0     8
*
*--- solution control:
* ISTART ISOLV NSET
   2       2     4
*
*--- data weighting and re-weighting:
* NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP
   5     -9    -9    -9   -9    1     1     8    -9   60
   5     -9    -9    -9   -9    1     1     6    5    60
   5     -9    -9    -9   -9    1    0.8    4    3    60
   5     -9    -9    -9   -9    1    0.8    3    2    60
*
*--- 1D model:
* NLAY RATIO
  {len(MTOP)}   {VPVS}
* TOP
{top}
* VEL
{vel}
*
*--- event selection:
* CID
   0
* ID
""")
    print(f"wrote {HD}/hypoDD.pha, station.dat, ph2dt.inp, hypoDD.inp")

if __name__=="__main__":
    main()
