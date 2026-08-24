#!/usr/bin/env python3
"""What the ocean-bottom stations would add to the tomography.

The OBS arrival times are not usable yet - no GPS underwater, no drift file on
the archive - so they are absent from the catalogue. This builds the ray set
they *would* contribute if the clocks were corrected: every catalogue event is
paired with each OBS site within a distance cutoff, and those rays are appended
to the land ray set. Comparing coverage and checkerboard recovery with and
without them says whether the clock-correction work is worth doing.

The OH sites are hydrophones and get P only; OS12 is left out (see OBS_QC.md).

Usage:
  tomo_obs_gain.py --tomo tomo5 --out tomo5_obs --maxdist 250
"""
import argparse, os, sys, time

import numpy as np
import pandas as pd
from scipy.sparse import load_npz, save_npz, csr_matrix, vstack

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from raytrace1d import load_layers, Model1D, Tracer
from tomo_rays import Grid

# folder naming is authoritative: OS = three-component, OH = hydrophone
OS_SITES = ["OS02", "OS06", "OS07", "OS08", "OS09"]     # OS12 excluded, see OBS_QC.md
OH_SITES = ["OH01", "OH03", "OH04", "OH05", "OH10", "OH11", "OH13", "OH14"]
CODE_OF = {"OS02": "OB2", "OS06": "OB6", "OS07": "OB7", "OS08": "OB8", "OS09": "OB9",
           "OH01": "OB1", "OH03": "OB3", "OH04": "OB4", "OH05": "OB5",
           "OH10": "O10", "OH11": "O11", "OH13": "O13", "OH14": "O14"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tomo", default="tomo5")
    ap.add_argument("--out", default="tomo5_obs")
    ap.add_argument("--info", default=os.path.join(HERE, "..", "config", "stations_info.csv"))
    ap.add_argument("--vel", default=os.path.join(HERE, "..", "nll", "nll_vel_full.in"))
    ap.add_argument("--maxdist", type=float, default=250.0)
    ap.add_argument("--seg", type=float, default=1.0)
    a = ap.parse_args()

    g = np.load(os.path.join(a.tomo, "coverage.npz"))
    grid = Grid(float(g["lat0"]), float(g["lat0"]) + int(g["nlat"]) * float(g["dlat"]),
                float(g["lon0"]), float(g["lon0"]) + int(g["nlon"]) * float(g["dlon"]),
                int(g["nz"]) * float(g["dz"]), float(g["dh_km"]), float(g["dz"]))
    # the saved grid is authoritative; re-deriving the cell counts from the
    # spacing can round differently and silently shift every cell index
    grid.nlat, grid.nlon, grid.nz = int(g["nlat"]), int(g["nlon"]), int(g["nz"])
    grid.dlat, grid.dlon, grid.dz = float(g["dlat"]), float(g["dlon"]), float(g["dz"])
    grid.n = grid.nlat * grid.nlon * grid.nz
    events = np.load(os.path.join(a.tomo, "events.npy"))

    info = pd.read_csv(a.info)
    coord = {r.sta: (r.lat, r.lon) for r in info.itertuples()}
    sites = [(s, "P") for s in OH_SITES] + [(s, "PS") for s in OS_SITES]
    missing = [s for s, _ in sites if CODE_OF[s] not in coord]
    if missing:
        print("no coordinates for", missing)
    print(f"{len(sites)} OBS sites, {len(events)} events, cutoff {a.maxdist:g} km")

    os.makedirs(a.out, exist_ok=True)
    total_hits = np.zeros(grid.n)
    for phase in ("P", "S"):
        G = load_npz(os.path.join(a.tomo, f"G_{phase}.npz")).tocsr()
        meta = np.load(os.path.join(a.tomo, f"meta_{phase}.npy"))
        tracer = Tracer(Model1D(load_layers(a.vel, phase), zmax=grid.zmax + 150.0))
        rows, cols, vals, newmeta = [], [], [], []
        r = 0
        t0 = time.time()
        for site, kind in sites:
            if phase == "S" and kind == "P":
                continue
            code = CODE_OF[site]
            if code not in coord:
                continue
            slat, slon = coord[code]
            for ei, (elat, elon, edep) in enumerate(events):
                dy = (slat - elat) * 111.19
                dx = (slon - elon) * 111.19 * np.cos(np.radians(elat))
                dist = np.hypot(dx, dy)
                if dist > a.maxdist:
                    continue
                got = tracer.ray_path(edep, max(dist, 0.05))
                if got is None:
                    continue
                x, z, tt = got
                s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(z)))])
                ns = max(int(s[-1] / a.seg), 2)
                su = np.linspace(0, s[-1], ns + 1)
                xu = np.interp(su, s, x); zu = np.interp(su, s, z)
                seg = np.diff(su)
                xm = 0.5 * (xu[:-1] + xu[1:]); zm = 0.5 * (zu[:-1] + zu[1:])
                f = xm / max(dist, 1e-6)
                idx, ok = grid.index(elat + f * (slat - elat),
                                     elon + f * (slon - elon), zm)
                idx, seg = idx[ok], seg[ok]
                if not len(idx):
                    continue
                order = np.argsort(idx)
                idx, seg = idx[order], seg[order]
                uniq, start = np.unique(idx, return_index=True)
                rows.append(np.full(len(uniq), r))
                cols.append(uniq)
                vals.append(np.add.reduceat(seg, start))
                newmeta.append((ei, r, 0.0, tt))
                r += 1
        print(f"  {phase}: {r} synthetic OBS rays in {time.time()-t0:.0f}s")
        if r:
            Gobs = csr_matrix((np.concatenate(vals),
                               (np.concatenate(rows), np.concatenate(cols))),
                              shape=(r, grid.n))
            Gall = vstack([G, Gobs]).tocsr()
            meta_all = np.vstack([meta, np.array(newmeta, dtype=float)])
        else:
            Gall, meta_all = G, meta
        save_npz(os.path.join(a.out, f"G_{phase}.npz"), Gall)
        np.save(os.path.join(a.out, f"meta_{phase}.npy"), meta_all)
        total_hits += np.asarray((Gall > 0).sum(axis=0)).ravel()

    np.savez(os.path.join(a.out, "coverage.npz"), hits=total_hits,
             lens=total_hits, nlat=g["nlat"], nlon=g["nlon"], nz=g["nz"],
             lat0=g["lat0"], lon0=g["lon0"], dlat=g["dlat"], dlon=g["dlon"],
             dz=g["dz"], dh_km=g["dh_km"])
    np.save(os.path.join(a.out, "events.npy"), events)
    base = np.load(os.path.join(a.tomo, "coverage.npz"))["hits"]
    print(f"cells with >=5 rays: {(base>=5).sum()} land-only -> "
          f"{(total_hits>=5).sum()} with OBS "
          f"(+{100*((total_hits>=5).sum()/(base>=5).sum()-1):.0f}%)")


if __name__ == "__main__":
    main()
