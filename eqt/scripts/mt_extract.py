#!/usr/bin/env python3
"""Extract per-event 3-component windows from raw EDL for BayesISOLA MT.
For an event, writes <outdir>/sac/YK.TFxx..HH[ZNE].sac and <outdir>/network.stn
with date-correct site coordinates. p0=Z, p1=N, p2=E; 200 Hz -> decimate to 100.
"""
import os, re, glob, json
import numpy as np, obspy
from obspy import UTCDateTime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW  = "/Volumes/Untitled 1/DATA-GFZ-Gempa-JOgja-tahap-2"
PER  = json.load(open(f"{ROOT}/config/stations_periods.json"))
PRE, POST = 90.0, 90.0           # window seconds around origin (long pre for noise slice)
COMP = {"pri0":("HHZ",0.0,-90.0),"pri1":("HHN",0.0,0.0),"pri2":("HHE",90.0,0.0)}
fn_re = re.compile(r"e\d{4}(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\.pri0$")

def seg_start(path):
    m=fn_re.search(os.path.basename(path))
    yy,mo,dd,hh,mi,ss=map(int,m.groups())
    return UTCDateTime(2000+yy,mo,dd,hh,mi,ss)

def site_coords(code, jd):
    for s in PER[code]["sites"]:
        if s["day_start"]<=jd<=s["day_end"]:
            return s["lat"], s["lon"]
    return None

def extract_event(t0, outdir):
    """t0: UTCDateTime origin. Returns list of (code,lat,lon) with written data."""
    sacdir=f"{outdir}/sac"; os.makedirs(sacdir, exist_ok=True)
    w0, w1 = t0-PRE, t0+POST
    jd = int(t0.julday)
    got=[]
    for code,info in PER.items():
        folder=info["folder"]; cc=site_coords(code,jd)
        if cc is None: continue
        lat,lon=cc
        # candidate day folders (handle window crossing midnight)
        days=sorted(set([int((w0).julday), int((w1).julday)]))
        p0s=[]
        for d in days:
            p0s += glob.glob(f"{RAW}/{folder}/{d}/*.pri0")
        # pick segments overlapping [w0,w1] (segments are 1800 s)
        sel=[p for p in p0s if seg_start(p)<=w1 and seg_start(p)+1800.0>=w0]
        if not sel: continue
        comp_ok=0
        for ext,(chan,az,dip) in COMP.items():
            st=obspy.Stream()
            for p0 in sorted(sel):
                f=p0.replace(".pri0","."+ext)
                if os.path.exists(f):
                    try: st+=obspy.read(f,format="MSEED")
                    except Exception: pass
            if len(st)==0: continue
            st.merge(method=1, fill_value=0)
            try: st.decimate(2, no_filter=False)   # 200->100 Hz
            except Exception: pass
            st.trim(w0, w1, pad=True, fill_value=0)
            tr=st[0]
            tr.stats.network="YK"; tr.stats.station=code; tr.stats.location=""; tr.stats.channel=chan
            tr.stats.sac=obspy.core.AttribDict(dict(stla=lat,stlo=lon,stel=0.0,cmpaz=az,
                         cmpinc=dip+90.0, knetwk="YK", kstnm=code, kcmpnm=chan))
            tr.write(f"{sacdir}/YK.{code}..{chan}.sac", format="SAC")
            comp_ok+=1
        if comp_ok==3: got.append((code,lat,lon))
    with open(f"{outdir}/network.stn","w") as f:
        for code,lat,lon in got:
            f.write(f"YK:{code}::HH   {lat:.5f}   {lon:.5f}\n")
    return got

if __name__=="__main__":
    import sys
    t0=UTCDateTime(sys.argv[1]); outdir=sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    got=extract_event(t0, outdir)
    print(f"event {t0}: {len(got)} stations with 3-comp data:", [g[0] for g in got])
