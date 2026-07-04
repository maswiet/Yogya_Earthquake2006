#!/usr/bin/env python3
"""Compute ML per event from Wood-Anderson amplitudes, merge into a magnitude
catalog, and estimate the Gutenberg-Richter b-value and completeness Mc.

Per-pick ML = log10(A_mm) + 1.110*log10(R/100) + 0.00189*(R-100) + 3.0   (Hutton & Boore 1987)
Per-event ML = median over stations (>= min_sta amplitudes).
"""
import os, re, glob, argparse
import numpy as np, pandas as pd
from obspy import UTCDateTime
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def ml_pick(A_mm, R_km):
    return np.log10(A_mm) + 1.110*np.log10(R_km/100.) + 0.00189*(R_km-100.) + 3.0

def event_meta():
    """evid -> (time,lat,lon,depth,gap,rms) using the SAME sorted order as amplitudes."""
    geo=re.compile(r"OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                   r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
    qual=re.compile(r"RMS\s+([\d.eE+-]+)\s+Nphs\s+\d+\s+Gap\s+([\d.]+)")
    rows={}
    for evid,path in enumerate(sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))):
        t=lat=lon=dep=None; gap=999.; rms=9.
        for line in open(path):
            if line.startswith("GEOGRAPHIC"):
                m=geo.search(line)
                if m:
                    y,mo,d,h,mi=map(int,m.groups()[:5]); s=float(m.group(6))
                    t=UTCDateTime(f"{y:04d}-{mo:02d}-{d:02d}T{h:02d}:{mi:02d}:00")+s
                    lat=float(m.group(7)); lon=float(m.group(8)); dep=float(m.group(9))
            elif line.startswith("QUALITY"):
                mm=qual.search(line)
                if mm: rms=float(mm.group(1)); gap=float(mm.group(2))
        rows[evid]=(t.isoformat(),lat,lon,dep,gap,rms)
    return rows

def maxc_mc(mags, dM=0.1):
    """Maximum-curvature Mc + 0.2 correction."""
    bins=np.arange(np.floor(mags.min()/dM)*dM, mags.max()+dM, dM)
    h,_=np.histogram(mags,bins=bins)
    return bins[np.argmax(h)]+0.2

def b_value(mags, mc, dM=0.1):
    m=mags[mags>=mc-1e-9]
    if len(m)<20: return np.nan, np.nan, len(m)
    b=np.log10(np.e)/(m.mean()-(mc-dM/2))
    b_err=b/np.sqrt(len(m))           # Aki/Shi-Bolt approx
    return b, b_err, len(m)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--amp", default=f"{ROOT}/full/amplitudes.csv")
    ap.add_argument("--out", default=f"{ROOT}/full/catalog_magnitude.csv")
    ap.add_argument("--fig", default=f"{ROOT}/figures/magnitude_gutenberg_richter.png")
    ap.add_argument("--min_sta", type=int, default=3)
    a=ap.parse_args()
    amp=pd.read_csv(a.amp)
    # use REFINED-catalog distances: merge measured amp_mm with refined dist+depth
    pk=pd.read_csv(f"{ROOT}/full/picks_for_amp.csv")
    pk["hypo_ref"]=np.sqrt(pk["dist"]**2 + pk["depth"]**2)
    amp=amp.merge(pk[["evid","sta","hypo_ref"]], on=["evid","sta"], how="left")
    amp["hypo_km"]=amp["hypo_ref"].fillna(amp["hypo_km"])   # fallback to stored if unmatched
    amp["ml"]=ml_pick(amp.amp_mm, amp.hypo_km)
    g=amp.groupby("evid")
    ev=g.agg(ML=("ml","median"), ML_std=("ml","std"), n_sta=("ml","size")).reset_index()
    ev=ev[ev.n_sta>=a.min_sta]
    meta=event_meta()
    ev["time"]=ev.evid.map(lambda i: meta[i][0])
    ev["latitude"]=ev.evid.map(lambda i: meta[i][1])
    ev["longitude"]=ev.evid.map(lambda i: meta[i][2])
    ev["depth"]=ev.evid.map(lambda i: meta[i][3])
    ev["gap"]=ev.evid.map(lambda i: meta[i][4])
    ev["ML"]=ev.ML.round(2)
    ev.to_csv(a.out, index=False)

    M=ev.ML.values
    mc=maxc_mc(M); b,berr,nb=b_value(M,mc)
    print(f"events with ML: {len(ev)}  (>= {a.min_sta} stations)")
    print(f"ML range {M.min():.2f}..{M.max():.2f}  median {np.median(M):.2f}")
    print(f"Mc = {mc:.2f} | b-value = {b:.2f} +/- {berr:.2f}  (N>=Mc = {nb})")

    # figure: FMD + histogram + ML-time + ML-map
    fig,ax=plt.subplots(2,2,figsize=(14,11))
    dM=0.1; bins=np.arange(np.floor(M.min()),M.max()+dM,dM)
    h,edges=np.histogram(M,bins=bins); centers=(edges[:-1]+edges[1:])/2
    cum=np.cumsum(h[::-1])[::-1]
    ax[0,0].semilogy(centers,cum,"s",ms=4,color="navy",label="cumulative")
    ax[0,0].semilogy(centers,h,"o",ms=3,color="orange",label="non-cumulative")
    mm=np.linspace(mc,M.max(),20); a_val=np.log10(cum[np.argmin(abs(centers-mc))])+b*mc
    ax[0,0].semilogy(mm,10**(a_val-b*mm),"r-",lw=2,label=f"b={b:.2f}±{berr:.2f}")
    ax[0,0].axvline(mc,color="gray",ls="--",label=f"Mc={mc:.2f}")
    ax[0,0].set_xlabel("ML"); ax[0,0].set_ylabel("N"); ax[0,0].legend()
    ax[0,0].set_title("Frequency-magnitude distribution")
    ax[0,1].hist(M,bins=bins,color="steelblue",edgecolor="k")
    ax[0,1].axvline(mc,color="red",ls="--"); ax[0,1].set_xlabel("ML"); ax[0,1].set_ylabel("events")
    ax[0,1].set_title(f"ML histogram (median {np.median(M):.2f}, max {M.max():.2f})")
    t=pd.to_datetime(ev.time)
    ax[1,0].scatter(t,ev.ML,s=4,alpha=0.3,color="darkred")
    ax[1,0].set_ylabel("ML"); ax[1,0].set_title("ML vs time")
    plt.setp(ax[1,0].get_xticklabels(),rotation=45,ha="right",fontsize=7)
    sc=ax[1,1].scatter(ev.longitude,ev.latitude,s=6,c=ev.ML,cmap="hot_r",vmin=mc,vmax=3)
    plt.colorbar(sc,ax=ax[1,1],label="ML"); ax[1,1].set_aspect("equal","box")
    ax[1,1].set_title("Epicentres coloured by ML"); ax[1,1].set_xlabel("Lon"); ax[1,1].set_ylabel("Lat")
    plt.tight_layout(); plt.savefig(a.fig,dpi=120)
    print("wrote",a.out,"and",a.fig)

if __name__=="__main__":
    main()
