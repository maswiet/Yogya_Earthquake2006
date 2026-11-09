#!/usr/bin/env python3
"""Batch BayesISOLA MT inversion for Yogya 2006 aftershocks (ML>=MLMIN).
Run FROM eqt/mt (needs ./green with compiled gr_xyz, elemse).
Writes eqt/mt/results/mt_catalog.csv (checkpointed per event) + per-event outputs.
"""
import sys, os, glob, shutil, time, traceback
CLONE="/Users/maswiet/Work/Yogya_Earthquake2006/eqt/tools/BayesISOLA"
sys.path.insert(0, CLONE)
sys.path.insert(0, "/Users/maswiet/Work/Yogya_Earthquake2006/eqt/scripts")
import numpy as np, pandas as pd, BayesISOLA
from obspy import UTCDateTime
from mt_extract import extract_event

ROOT="/Users/maswiet/Work/Yogya_Earthquake2006/eqt"
MT=f"{ROOT}/mt"; RES=f"{MT}/results"; EVROOT=f"{MT}/events"
os.makedirs(RES, exist_ok=True); os.makedirs(EVROOT, exist_ok=True)
MLMIN=float(os.environ.get("MLMIN","1.5"))
NSTA_MIN=int(os.environ.get("NSTA_MIN","6"))

def band(mag):
    if mag>=2.5: return 0.8,2.0
    if mag>=2.0: return 1.0,2.5
    return 1.2,3.0

def stf_dur(mag):   # triangle half-duration ~ source duration (s)
    return max(0.2, 0.5*10**(0.5*(mag-3.0)))   # ~0.5 s at M3, shorter for small

def clear_green():
    for p in glob.glob(f"{MT}/green/*.dat")+glob.glob(f"{MT}/green/elemse*.txt")+[f"{MT}/green/grdat.hed"]:
        try: os.remove(p)
        except OSError: pass

def run_event(evid, t0, lat, lon, depth, mag):
    evdir=f"{EVROOT}/{evid}"; outdir=f"{evdir}/out"
    if os.path.isdir(evdir): shutil.rmtree(evdir)
    os.makedirs(evdir)
    got=extract_event(t0, evdir)
    if len(got)<NSTA_MIN: return dict(evid=evid, status=f"skip_nsta={len(got)}")
    fmin,fmax=band(mag)
    clear_green()
    inp=BayesISOLA.load_data(outdir=outdir)
    inp.set_event_info(lat=lat,lon=lon,depth=max(depth,1.0),mag=mag,t=t0,agency='EQT')
    inp.set_source_time_function('triangle', stf_dur(mag))
    inp.read_network_coordinates(f"{evdir}/network.stn")
    inp.read_crust(f"{MT}/crustal.dat")
    inp.load_files(f"{evdir}/sac/", xml_dir=f"{MT}/stations.xml")
    grid=BayesISOLA.grid(inp, location_unc=1500, depth_unc=12000, time_unc=1.5,
                         step_x=1000, step_z=1000, max_points=300, circle_shape=True, rupture_velocity=1000)
    data=BayesISOLA.process_data(inp, grid, threads=8, use_precalculated_Green='auto', fmax=fmax, fmin=fmin)
    cova=BayesISOLA.covariance_matrix(data)
    cova.covariance_matrix_noise(crosscovariance=True, save_non_inverted=True)
    sol=BayesISOLA.resolve_MT(data, cova, deviatoric=True)
    c=sol.centroid
    g=lambda k: c.get(k) if isinstance(c,dict) else getattr(c,k,None)
    r=dict(evid=evid, time=str(t0), lat=lat, lon=lon, dep_cat=depth, ML=mag, nsta=len(got),
           fmin=fmin, fmax=fmax, status="ok")
    for k in ('VR','dc_perc','clvd_perc','iso_perc','Mw','mom','s1','d1','r1','s2','d2','r2',
              'shift','x','y','z','det_Ca','CN'):
        v=g(k)
        if v is not None: r[k]=v
    return r

def main():
    cat=pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    cat=cat[(cat.gap<180)&(cat.ML>=MLMIN)].copy().sort_values("ML",ascending=False).reset_index(drop=True)
    # better absolute loc: merge GrowClust if time-matched (optional)
    print(f"{len(cat)} events ML>={MLMIN}")
    outcsv=f"{RES}/mt_catalog.csv"; rows=[]
    done=set()
    if os.path.exists(outcsv):
        prev=pd.read_csv(outcsv); rows=prev.to_dict('records'); done=set(prev.evid)
    for i,e in cat.iterrows():
        t0=UTCDateTime(e['time']); evid=t0.strftime("%Y%m%dT%H%M%S")
        if evid in done: continue
        t=time.time()
        try:
            r=run_event(evid, t0, float(e.latitude), float(e.longitude), float(e.depth), float(e.ML))
        except Exception as ex:
            r=dict(evid=evid, time=str(t0), ML=float(e.ML), status="ERR:"+str(ex)[:80])
            traceback.print_exc()
        r['sec']=round(time.time()-t,1)
        rows.append(r); pd.DataFrame(rows).to_csv(outcsv, index=False)
        vr=r.get('VR'); dc=r.get('dc_perc')
        print(f"[{i+1}/{len(cat)}] {evid} ML{e.ML:.2f} {r['status']} "
              f"VR={vr if vr is None else round(vr,3)} DC%={dc if dc is None else round(dc)} "
              f"s/d/r={r.get('s1')},{r.get('d1')},{r.get('r1')} ({r['sec']}s)", flush=True)
    print("DONE batch")

if __name__=="__main__":
    main()
