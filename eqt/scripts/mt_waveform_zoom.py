#!/usr/bin/env python3
"""Zoomed observed-vs-synthetic waveform fit for one event, window [0, TMAX] s,
y-axis auto-scaled per station to the peak inside the window. Reconstructs the
synthetic exactly as BayesISOLA.plot.plot_seismo does.

Run from eqt/mt:  python ../scripts/mt_waveform_zoom.py [EVID] [TMAX]
"""
import sys, os
CLONE="/Users/maswiet/Work/Yogya_Earthquake2006/eqt/tools/BayesISOLA"
sys.path.insert(0, CLONE)
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import BayesISOLA
from obspy import UTCDateTime
from BayesISOLA.fileformats import read_elemse, read_elemse_from_files
from BayesISOLA.helpers import my_filter

ROOT="/Users/maswiet/Work/Yogya_Earthquake2006/eqt"; MT=f"{ROOT}/mt"; FIG=f"{ROOT}/figures"

def main():
    evid = sys.argv[1] if len(sys.argv)>1 else "20060609T083059"
    tmax = float(sys.argv[2]) if len(sys.argv)>2 else 10.0
    r = pd.read_csv(f"{MT}/results/mt_catalog.csv").query("evid==@evid").iloc[0]
    evdir=f"{MT}/events/{evid}"; outdir=f"{evdir}/out"
    t0=UTCDateTime(r['time']); fmin,fmax=float(r['fmin']),float(r['fmax'])
    import glob
    for p in glob.glob(f"{MT}/green/*.dat")+glob.glob(f"{MT}/green/elemse*.txt")+[f"{MT}/green/grdat.hed"]:
        try: os.remove(p)
        except OSError: pass

    inp=BayesISOLA.load_data(outdir=outdir)
    inp.set_event_info(lat=float(r['lat']),lon=float(r['lon']),depth=max(float(r['dep_cat']),1.0),
                       mag=float(r['ML']),t=t0,agency='EQT')
    inp.set_source_time_function('triangle', 0.2)
    inp.read_network_coordinates(f"{evdir}/network.stn")
    inp.read_crust(f"{MT}/crustal.dat")
    inp.load_files(f"{evdir}/sac/", xml_dir=f"{MT}/stations.xml")
    grid=BayesISOLA.grid(inp, location_unc=1500, depth_unc=12000, time_unc=1.5,
                         step_x=1000, step_z=1000, max_points=300, circle_shape=True, rupture_velocity=1000)
    data=BayesISOLA.process_data(inp, grid, threads=8, use_precalculated_Green='auto', fmax=fmax, fmin=fmin)
    cova=BayesISOLA.covariance_matrix(data); cova.covariance_matrix_noise(crosscovariance=True, save_non_inverted=True)
    sol=BayesISOLA.resolve_MT(data, cova, deviatoric=True)
    c=sol.centroid
    print(f"{evid} VR={c['VR']:.3f} DC%={c['dc_perc']:.0f} Mw={c['Mw']:.2f} "
          f"s/d/r={c['s1']:.0f}/{c['d1']:.0f}/{c['r1']:.0f} depth={c['z']/1000:.1f}km shift={c['shift']:.2f}s")

    # reconstruct observed + synthetic (mirrors plot_seismo)
    obs = data.data_shifts[c['shift_idx']]
    npts = data.npts_slice; sr = data.samprate
    if c['path']:
        elemse = read_elemse_from_files(inp.nr, c['path'], inp.stations, inp.event['t'], sr, data.npts_elemse, data.invert_displacement)
    else:
        elemse = read_elemse(inp.nr, data.npts_elemse, 'green/elemse'+c['id']+'.dat', inp.stations, data.invert_displacement)
    for rr in range(inp.nr):
        for e in range(6):
            my_filter(elemse[rr][e], inp.stations[rr]['fmin'], inp.stations[rr]['fmax'])
            elemse[rr][e].trim(UTCDateTime(0)+data.elemse_start_origin)
    t=np.arange(0,(npts-0.5)/sr,1./sr)
    win=t<=tmax

    used=[s for s in range(inp.nr) if (inp.stations[s]['useZ'] or inp.stations[s]['useN'] or inp.stations[s]['useE'])]
    comps=[0,1,2]; cname=['Z','N','E']
    fig,ax=plt.subplots(len(used),3,figsize=(13,2.4*len(used)+0.6),squeeze=False,sharex=True)
    for ri,s in enumerate(used):
        SY={}
        for comp in range(3):
            SY[comp]=np.zeros(npts)
            for e in range(6):
                SY[comp]+=elemse[s][e][comp].data[0:npts]*c['a'][e,0]
        pk=0.0
        for comp in comps:
            d=obs[s][comp].data[0:len(t)]
            use=inp.stations[s][{0:'useZ',1:'useN',2:'useE'}[comp]]
            if use: pk=max(pk, np.abs(d[win]).max(), np.abs(SY[comp][win]).max())
        pk=pk or 1.0
        for comp in comps:
            a=ax[ri][comp]
            d=obs[s][comp].data[0:len(t)]; use=inp.stations[s][{0:'useZ',1:'useN',2:'useE'}[comp]]
            a.plot(t[win],d[win],'k',lw=1.6,label='observed')
            if use:
                a.plot(t[win],SY[comp][win],'r',lw=1.6,label='synthetic')
            else:
                a.plot(t[win],SY[comp][win],color='0.6',lw=1.2)
            a.set_xlim(0,tmax); a.set_ylim(-1.15*pk,1.15*pk)
            a.axhline(0,color='0.8',lw=0.5)
            if ri==0: a.set_title(cname[comp],fontsize=15)
            a.tick_params(labelsize=10)
            if comp==0:
                a.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
            else:
                a.set_yticklabels([])
            if not use and comp==1:
                a.text(0.5,0.9,'(not inverted)',transform=a.transAxes,ha='center',fontsize=9,color='0.5')
        ax[ri][0].set_ylabel(f"{inp.stations[s]['code']}\n{inp.stations[s]['dist']/1000:.0f} km, {inp.stations[s]['az']:.0f}°",fontsize=11)
    for comp in range(3): ax[-1][comp].set_xlabel('time [s]',fontsize=12)
    h,l=ax[0][1].get_legend_handles_labels()
    fig.legend(h,l,loc='upper center',ncol=2,fontsize=12,frameon=True,bbox_to_anchor=(0.5,1.0))
    fig.suptitle(f"{evid}  ML{r['ML']}  Mw{c['Mw']:.2f}  VR={c['VR']:.2f}  DC={c['dc_perc']:.0f}%  "
                 f"({len(used)} sta used)  —  first {tmax:.0f} s",y=1.0 if len(used)>1 else 1.06,fontsize=12)
    fig.tight_layout(rect=[0,0,1,0.97 if len(used)>1 else 0.8])
    out=f"{FIG}/mt_waveform_zoom_{evid}.png"
    fig.savefig(out,dpi=200,bbox_inches='tight'); print("wrote",out)

if __name__=="__main__":
    import multiprocessing as _mp; _mp.freeze_support(); main()
