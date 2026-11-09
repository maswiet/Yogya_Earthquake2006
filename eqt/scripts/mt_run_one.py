#!/usr/bin/env python3
"""Run one BayesISOLA MT inversion. Run from eqt/mt (needs ./green with binaries)."""
import sys, os
CLONE="/Users/maswiet/Work/Yogya_Earthquake2006/eqt/tools/BayesISOLA"
sys.path.insert(0, CLONE)   # use patched clone
import numpy as np, BayesISOLA
from obspy import UTCDateTime

def main():
    evdir, outdir = sys.argv[1], sys.argv[2]
    lat,lon,depth,mag = map(float, sys.argv[3:7]); t0=UTCDateTime(sys.argv[7])
    fmin,fmax = (float(sys.argv[8]), float(sys.argv[9])) if len(sys.argv)>9 else (0.7,2.5)
    MT="/Users/maswiet/Work/Yogya_Earthquake2006/eqt/mt"

    inp=BayesISOLA.load_data(outdir=outdir)
    inp.set_event_info(lat=lat,lon=lon,depth=depth,mag=mag,t=t0,agency='EQT')
    inp.set_source_time_function('triangle', 0.5)
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
    print("CENTROID KEYS:", sorted(c.keys()))
    for k in ('VR','dc_perc','clvd_perc','iso_perc','Mw','mom','s1','d1','r1','s2','d2','r2','shift','z','x','y'):
        if k in c: print(f"  {k} = {c[k]}")
    print("max_VR:", getattr(sol,'max_VR',None))


if __name__=='__main__':
    import multiprocessing as _mp
    _mp.freeze_support()
    main()
