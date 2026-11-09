#!/usr/bin/env python3
"""Re-run one BayesISOLA MT inversion and plot observed-vs-synthetic seismograms
(waveform fit) + the best-fitting MT beachball, for illustration.

Run from eqt/mt:  python ../scripts/mt_waveform_plot.py [EVID]
Default EVID = best-VR event. Reuses existing events/<EVID>/sac + network.stn.
Figures land in events/<EVID>/out/ ; the key ones are copied to eqt/figures/.
"""
import sys, os, shutil
CLONE="/Users/maswiet/Work/Yogya_Earthquake2006/eqt/tools/BayesISOLA"
sys.path.insert(0, CLONE)
import pandas as pd
import BayesISOLA
from obspy import UTCDateTime

ROOT="/Users/maswiet/Work/Yogya_Earthquake2006/eqt"
MT=f"{ROOT}/mt"; FIG=f"{ROOT}/figures"

def main():
    evid = sys.argv[1] if len(sys.argv)>1 else "20060609T083059"
    cat = pd.read_csv(f"{MT}/results/mt_catalog.csv")
    row = cat[cat.evid==evid]
    if row.empty:
        print("evid not in catalog:", evid); return 1
    r = row.iloc[0]
    evdir=f"{MT}/events/{evid}"; outdir=f"{evdir}/out"
    if not os.path.isdir(f"{evdir}/sac"):
        print("no SAC for", evid, "- re-run batch extraction first"); return 1
    t0=UTCDateTime(r['time']); lat,lon=float(r['lat']),float(r['lon'])
    depth=max(float(r['dep_cat']),1.0); mag=float(r['ML'])
    fmin,fmax=float(r['fmin']),float(r['fmax'])
    print(f"{evid}  ML{mag}  VR(catalog)={r['VR']:.3f}  band {fmin}-{fmax} Hz  nsta={int(r['nsta'])}")

    # clear old Green's functions
    import glob
    for p in glob.glob(f"{MT}/green/*.dat")+glob.glob(f"{MT}/green/elemse*.txt")+[f"{MT}/green/grdat.hed"]:
        try: os.remove(p)
        except OSError: pass

    inp=BayesISOLA.load_data(outdir=outdir)
    inp.set_event_info(lat=lat,lon=lon,depth=depth,mag=mag,t=t0,agency='EQT')
    inp.set_source_time_function('triangle', 0.2)
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
    print(f"re-run: VR={c['VR']:.3f}  DC%={c['dc_perc']:.0f}  Mw={c['Mw']:.2f}  "
          f"s/d/r={c['s1']:.0f}/{c['d1']:.0f}/{c['r1']:.0f}  depth={c['z']/1000:.1f} km")

    # plot: waveform fit (obs vs synthetic) + beachball, skip the slow map/PPDF plots
    BayesISOLA.plot(sol, maps=False, slices=False, maps_sum=False, MT=True,
                    uncertainty=None, seismo=True, seismo_sharey=True, seismo_cova=True,
                    noise=False, spectra=False, stations=False, covariance_matrix=False)

    os.makedirs(FIG, exist_ok=True)
    produced=[]
    for src,dst in (("seismo.png",       f"mt_waveform_fit_{evid}.png"),
                    ("seismo_sharey.png", f"mt_waveform_fit_sharey_{evid}.png"),
                    ("seismo_cova.png",   f"mt_waveform_fit_whitened_{evid}.png"),
                    ("centroid.png",      f"mt_beachball_{evid}.png")):
        s=f"{outdir}/{src}"
        if os.path.exists(s):
            shutil.copy(s, f"{FIG}/{dst}"); produced.append(f"{FIG}/{dst}")
    print("PRODUCED:")
    for p in produced: print("  ", p)

if __name__=="__main__":
    import multiprocessing as _mp; _mp.freeze_support(); main()
