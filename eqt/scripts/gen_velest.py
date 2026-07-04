#!/usr/bin/env python3
"""Generate VELEST inputs from NLLoc-located events for a minimum-1D inversion.

Writes into eqt/velest/:
  yogya.sta  - station file (6-char codes, per-period sites)
  yogya.pha  - CNV phase file (selected high-quality events)
  yogya.mod  - starting P & S 1-D model (data-informed Central Java)
  velest.cmn - control file (simultaneous inversion, station corrections)

Selection: gap<gapmax, rms<rmsmax, nphs>=minphs; subsample evenly in time to
max_events (VELEST minimum-1D uses a well-distributed subset).
"""
import os, re, glob, json, argparse, string
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OLAT, OLON = -7.92, 110.44
VPVS = 1.735
# Geologically-constrained model: depth(km, top of layer), Vp(km/s), vdamp.
#  - SHALLOW sediment FIXED a priori (Merapi volcaniclastics, >700 m base): damp 999.
#    (lateral W/E sediment-vs-limestone contrast is carried by STATION CORRECTIONS.)
#  - DEEP crust FIXED to realistic values (removes the earlier 7.16 km/s artifact): damp 999.
#  - Only the well-sampled MID-CRUST (2-13 km) is inverted (damp 1).
MODEL = [(-2.0,2.50,999.0),   # near-surface sediment (fixed)
         ( 0.0,2.90,999.0),   # Merapi sediment (fixed)
         ( 0.7,4.30,100.0),   # base of ~700 m sediment -> bedrock (lightly damped)
         ( 2.0,5.20,1.0),     # --- invert mid-crust ---
         ( 4.0,5.70,1.0),
         ( 7.0,5.95,1.0),
         (10.0,6.20,1.0),
         (13.0,6.40,20.0),    # transition (semi-fixed)
         (16.0,6.55,999.0),   # --- fixed deep crust ---
         (22.0,6.80,999.0),
         (30.0,7.20,999.0),
         (40.0,8.00,999.0)]

def parse_hyps():
    events = []
    for path in sorted(glob.glob(f"{ROOT}/nll/loc/yogya.2*.grid0.loc.hyp")):
        ot=lat=lon=dep=None; gap=9999.; rms=9.; phs=[]
        for line in open(path):
            if line.startswith("GEOGRAPHIC"):
                m=re.search(r"OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                            r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)",line)
                if m:
                    ot=tuple(int(x) for x in m.groups()[:5])+(float(m.group(6)),)
                    lat=float(m.group(7)); lon=float(m.group(8)); dep=float(m.group(9))
            elif line.startswith("QUALITY"):
                mg=re.search(r"RMS\s+([\d.eE+-]+).*Nphs\s+(\d+)\s+Gap\s+([\d.]+)",line)
                if mg: rms=float(mg.group(1)); nph=int(mg.group(2)); gap=float(mg.group(3))
            elif " > " in line and "GAU" in line:
                left,right=line.split(" > "); lf=left.split(); rf=right.split()
                tt=float(rf[0])+float(rf[1])
                phs.append((lf[0],lf[4],tt))
        if ot and phs:
            events.append(dict(ot=ot,lat=lat,lon=lon,dep=dep,gap=gap,rms=rms,phs=phs))
    return events

# station elevations (m) by our code (from XN metadata; unmapped early sites approx)
ELEV = {"TF18":44,"TF13":130,"TF14":50,"TF16":230,"TF19":175,"TF12":120,
        "TF15a":260,"TF10a":172,"TF10b":73,"TF17":200,"TF11b":200,"TF07b":120,
        "TF09b":55,"TF07a":110,"TF09a":60,"TF11a":200,"TF15b":172}

def write_stations(path):
    periods=json.load(open(f"{ROOT}/config/stations_periods.json"))
    lines=["(a6,f7.4,a1,1x,f8.4,a1,1x,i4,1x,i1,1x,i3,1x,f5.2,2x,f5.2)"]
    idx=1; codes=[]
    def add(code,s):
        nonlocal idx
        la=s["lat"]; lo=s["lon"]; ns="N"; ew="E"
        if la<0: ns="S"; la=-la
        if lo<0: ew="W"; lo=-lo
        elev=ELEV.get(code,120)
        lines.append(f"{code:<6}{la:7.4f}{ns} {lo:8.4f}{ew} {elev:4d} {1:1d} {idx:3d} {0.0:5.2f}  {0.0:5.2f}")
        idx+=1; codes.append(code)
    for code,info in periods.items():
        sites=info.get("sites",[])
        if len(sites)<=1:
            if sites: add(code,sites[0])
        else:
            for suf,s in zip(string.ascii_lowercase,sites): add(f"{code}{suf}",s)
    lines.append("")
    open(path,"w").write("\n".join(lines)+"\n")
    return set(codes)

def write_model(path):
    L=["Yogya 2006 minimum-1D start (data-informed Central Java)     Ref TF12"]
    L.append(f"{len(MODEL):3d}        vel,depth,vdamp,phase (f5.2,5x,f7.2,2x,f7.3,3x,a1)")
    for i,(d,vp,vd) in enumerate(MODEL):
        tag="\t\t\tP-VELOCITY MODEL" if i==0 else ""
        L.append(f"{vp:5.2f}     {d:7.2f}  {vd:7.3f}{tag}")
    L.append(f"{len(MODEL):3d}")
    for i,(d,vp,vd) in enumerate(MODEL):
        vs=vp/VPVS
        tag="\t\t\tS-VELOCITY MODEL" if i==0 else ""
        L.append(f"{vs:5.2f}     {d:7.2f}  {vd:7.3f}{tag}")
    open(path,"w").write("\n".join(L)+"\n")

def write_phases(path, events, known):
    out=[]
    for k,e in enumerate(events):
        yy,mo,dd,hh,mi,sec=e["ot"]; yy%=100
        la=e["lat"]; lo=e["lon"]; ns="N"; ew="E"
        if la<0: ns="S"; la=-la
        if lo<0: ew="W"; lo=-lo
        out.append("")
        out.append(f"{yy:2d}{mo:02d}{dd:02d} {hh:2d}{mi:02d} {sec:5.2f} "
                   f"{la:7.4f}{ns} {lo:8.4f}{ew} {e['dep']:7.2f}  {k/100.0:5.2f}")
        for sta,ph,tt in e["phs"]:
            if sta not in known: continue
            iwt = 0 if ph=="P" else 1
            out.append(f"  {sta:<6}  {ph:<1}   {iwt:1d}   {tt:6.2f}")
    out.append("")
    open(path,"w").write("\n".join(out)+"\n")

def write_cmn(path, neqs):
    cmn=f"""******* CONTROL-FILE FOR PROGRAM  V E L E S T *******
*** ( all lines starting with  *  are ignored! )
*** next line contains a title:
Yogya 2006 EQT minimum-1D inversion
***  olat       olon   icoordsystem      zshift   itrial ztrial    ised
  {OLAT:9.4f}  {-OLON:9.4f}      0            0.0      0     0.00       1
*** neqs   nshot   rotate
  {neqs:4d}      0      0.0
*** isingle   iresolcalc
       0          0
*** dmax    itopo    zmin     veladj    zadj   lowveloclay
   200.0      0      -2.00       0.20    5.00       0
*** nsp    swtfac   vpvs       nmod
     2      0.75    {VPVS:6.3f}        1
***   othet   xythet    zthet    vthet   stathet
      0.01    0.01      0.01     1.00     1.00
*** nsinv   nshcor   nshfix     iuseelev    iusestacorr
       1       0       0           1            1
*** iturbo    icnvout   istaout   ismpout
       1         1         1        0
*** irayout   idrvout   ialeout   idspout   irflout   irfrout   iresout
       0         0         0         0         0         0         0
*** delmin   ittmax   invertratio
    0.010      99          1
*** Modelfile:
yogya.mod
*** Stationfile:
yogya.sta
*** Seismofile:

*** File with region names:
regionsnamen.dat
*** File with region coordinates:
regionskoord.dat
*** File #1 with topo data:

*** File #2 with topo data:

*** DATA INPUT files:
*** File with Earthquake data:
yogya.pha
*** File with Shot data:

*** OUTPUT files:
*** Main print output file:
yogya.OUT
*** File with single event locations:
yogya.CNV
*** File with final hypocenters in *.cnv format:
yogya.finalcnv
*** File with new station corrections:
yogya.STA
"""
    open(path,"w").write(cmn)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--max_events", type=int, default=1200)
    ap.add_argument("--gapmax", type=float, default=160)
    ap.add_argument("--rmsmax", type=float, default=0.30)
    ap.add_argument("--minphs", type=int, default=10)
    a=ap.parse_args()
    vdir=f"{ROOT}/velest"; os.makedirs(vdir, exist_ok=True)
    ev=parse_hyps()
    print(f"parsed {len(ev)} located events")
    sel=[e for e in ev if e["gap"]<a.gapmax and e["rms"]<a.rmsmax and len(e["phs"])>=a.minphs]
    print(f"quality events (gap<{a.gapmax}, rms<{a.rmsmax}, nphs>={a.minphs}): {len(sel)}")
    if len(sel)>a.max_events:
        step=len(sel)/a.max_events
        sel=[sel[int(i*step)] for i in range(a.max_events)]
    print(f"selected {len(sel)} events for inversion")
    known=write_stations(f"{vdir}/yogya.sta")
    write_model(f"{vdir}/yogya.mod")
    write_phases(f"{vdir}/yogya.pha", sel, known)
    write_cmn(f"{vdir}/velest.cmn", len(sel))
    # copy region files (required by cmn)
    import shutil
    src=f"{ROOT}/tools/REAL/demo_real/VELEST"
    for f in ("regionsnamen.dat","regionskoord.dat"):
        shutil.copy(f"{src}/{f}", f"{vdir}/{f}")
    print(f"wrote VELEST inputs in {vdir}")

if __name__=="__main__":
    main()
