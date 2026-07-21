#!/usr/bin/env python3
"""Measure Wood-Anderson amplitudes for every S pick, to compute ML.

Step A: parse NLLoc (VELEST-reloc) hyp files -> picks_for_amp.csv
        (event, station TF-code, phase, arrival epoch, epicentral dist, depth).
Step B: per (logger folder, julian day): read .pri horizontals, remove L4-3D
        response, simulate Wood-Anderson, measure peak |A| in an S-window.
        -> amplitudes.csv (event, station, amp_mm, hypo_dist_km). Resumable.

Runs in the `eqt` env (obspy).
"""
import os, re, glob, sys, csv, argparse, math
import numpy as np
import obspy
from obspy.core import UTCDateTime

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# The EDL volume has mounted under different labels across sessions; take the
# first one that exists, or an explicit EDL_BASE override.
_CANDIDATES = ["/Volumes/Untitled/DATA-GFZ-Gempa-JOgja-tahap-2",
               "/Volumes/Untitled 1/DATA-GFZ-Gempa-JOgja-tahap-2"]
BASE = os.environ.get("EDL_BASE") or next(
    (p for p in _CANDIDATES if os.path.isdir(p)), _CANDIDATES[0])
COMP = {"pri1": "N", "pri2": "E"}          # horizontals (p1=N, p2=E)
year_re = re.compile(r"e\d{4}(\d{2})\d+\.pri1$")

# L4-3D geophone: f0=1 Hz, damping h=0.707, velocity sensitivity 1.7e8 counts/(m/s)
w0 = 2*math.pi*1.0; h = 0.707
PAZ_L4 = {"poles": [complex(-h*w0,  w0*math.sqrt(1-h*h)),
                    complex(-h*w0, -w0*math.sqrt(1-h*h))],
          "zeros": [0j, 0j], "gain": 1.0, "sensitivity": 1.7e8}
# Wood-Anderson (IASPEI): f0=1.25 Hz, h=0.8, static magnification 2080
PAZ_WA = {"poles": [complex(-6.2832,-4.7124), complex(-6.2832,4.7124)],
          "zeros": [0j, 0j], "gain": 1.0, "sensitivity": 2080.0}

def parse_picks():
    out = f"{ROOT}/full/picks_for_amp.csv"
    if os.path.exists(out):
        return out
    rows = []
    geo = re.compile(r"OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                     r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
    for evid, path in enumerate(sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))):
        depth=None; ot=None
        for line in open(path):
            if line.startswith("GEOGRAPHIC"):
                m=geo.search(line)
                if m:
                    depth=float(m.group(9))
            elif " > " in line and "GAU" in line:
                left,right=line.split(" > "); lf=left.split(); rf=right.split()
                sta,ph=lf[0],lf[4]
                if ph!="S":            # ML from S-wave (horizontals)
                    continue
                yyyymmdd=lf[6]; hhmm=lf[7]; sec=float(lf[8])
                t=UTCDateTime(f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:8]}T"
                              f"{hhmm[:2]}:{hhmm[2:4]}:00")+sec
                dist=float(rf[6])
                rows.append((evid, sta, round(t.timestamp,3), round(dist,2),
                             round(depth,2), t.julday))
    import pandas as pd
    pd.DataFrame(rows, columns=["evid","sta","arr","dist","depth","jday"]).to_csv(out,index=False)
    print(f"parsed {len(rows)} S picks -> {out}")
    return out

def folder_for(sta):            # TF07a -> tf3007
    return f"tf30{sta[2:4]}"

# Half-length of the window deconvolved around each pick. The instrument
# simulation must NOT be run on the merged day-long trace: obspy's simulate()
# applies a 5% cosine taper, which over 24 h suppresses everything within
# ~72 min of the day boundary (up to ~200x, i.e. >2 magnitude units, for events
# in the first few minutes of a day). Deconvolving a short window instead keeps
# the taper far away from the measured S phase.
PAD_S = 60.0

# Amplitudes are measured on a band-passed WA trace. Without this the peak is
# taken from the raw simulation, and stations sitting near power lines (TF10b
# carries a strong 50 Hz mains tone, sampled at 200 Hz) have their small-event
# amplitudes set by the hum rather than by the S phase -- a +0.27 ML bias at
# ML<0 for TF10b. 1-20 Hz spans the S-wave corner frequencies of local
# microearthquakes at 10-30 km while putting 50 Hz far into the stop band.
BP_LO, BP_HI, BP_CORNERS = 1.0, 20.0, 4


def build_day_raw(folder, jday):
    """Read + merge the raw horizontals for one logger-day (no deconvolution)."""
    dpath=os.path.join(BASE, folder, f"{jday:03d}")
    if not os.path.isdir(dpath): return None
    st=obspy.Stream()
    for comp_ext,ch in COMP.items():
        s=obspy.Stream()
        for p1 in sorted(glob.glob(os.path.join(dpath,"*.pri1"))):
            m=year_re.search(os.path.basename(p1))
            if not (m and m.group(1)=="06"): continue
            f=p1.replace(".pri1","."+comp_ext)
            if os.path.exists(f):
                try: s+=obspy.read(f,format="MSEED")
                except Exception: pass
        if len(s)==0: continue
        for tr in s: tr.stats.channel="HH"+ch
        s.merge(method=1, fill_value=0)
        st+=s
    if len(st)==0: return None
    st.detrend("demean")
    return st


def wa_window(raw, t0, t1, folder="", jday=0):
    """Wood-Anderson displacement (mm) for [t0,t1], deconvolved on a padded
    window so the simulate() taper never reaches the measurement interval."""
    seg=raw.slice(t0-PAD_S, t1+PAD_S)
    seg=obspy.Stream([tr for tr in seg if tr.stats.npts > 100])
    if len(seg)==0: return None
    seg=seg.copy(); seg.detrend("demean")
    try:
        seg.simulate(paz_remove=PAZ_L4, paz_simulate=PAZ_WA, water_level=10)
        seg.filter("bandpass", freqmin=BP_LO, freqmax=BP_HI,
                   corners=BP_CORNERS, zerophase=True)
    except Exception as ex:
        print(f"  sim fail {folder} {jday}: {ex}", file=sys.stderr); return None
    return seg

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{ROOT}/full/amplitudes.csv")
    ap.add_argument("--progress", default=f"{ROOT}/full/amp_done.txt")
    ap.add_argument("--only", default=None, help="limit to folder:jday for testing")
    a=ap.parse_args()
    import pandas as pd
    picks=pd.read_csv(parse_picks())
    picks["folder"]=picks["sta"].apply(folder_for)
    done=set(open(a.progress).read().split()) if os.path.exists(a.progress) else set()
    groups=picks.groupby(["folder","jday"])
    keys=list(groups.groups.keys())
    if a.only:
        f,j=a.only.split(":"); keys=[(f,int(j))]
    new = not os.path.exists(a.out)
    fout=open(a.out,"a",newline=""); w=csv.writer(fout)
    if new: w.writerow(["evid","sta","amp_mm","hypo_km"])
    fp=open(a.progress,"a")
    for i,(folder,jday) in enumerate(keys,1):
        key=f"{folder}:{jday}"
        if key in done and not a.only: continue
        g=groups.get_group((folder,jday))
        raw=build_day_raw(folder,jday)
        if raw is None:
            fp.write(key+"\n"); fp.flush(); continue
        n=0
        for _,r in g.iterrows():
            t0=UTCDateTime(r["arr"])-1.0; t1=UTCDateTime(r["arr"])+15.0
            wa=wa_window(raw, t0, t1, folder, jday)
            if wa is None: continue
            peak=0.0
            for tr in wa:
                seg=tr.slice(t0,t1)
                if seg.stats.npts>0:
                    peak=max(peak, float(np.max(np.abs(seg.data))))
            if peak>0:
                hypo=math.sqrt(r["dist"]**2+r["depth"]**2)
                w.writerow([int(r["evid"]), r["sta"], f"{peak:.5e}", round(hypo,2)]); n+=1
        fout.flush(); fp.write(key+"\n"); fp.flush()
        print(f"[{i}/{len(keys)}] {key}: {n} amps", flush=True)
    fout.close(); fp.close()

if __name__=="__main__":
    main()
