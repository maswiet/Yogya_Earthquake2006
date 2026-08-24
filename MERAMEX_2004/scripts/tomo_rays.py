#!/usr/bin/env python3
"""Ray tracing and coverage for the MERAMEX catalogue.

Reads the NonLinLoc per-event solutions, traces every weighted arrival through
the 1-D model, and writes:

  tomo/rays_<phase>.npz   sparse path-length matrix G (rays x cells) + ray metadata
  tomo/coverage.npz       hit count and total ray length per grid cell

The path-length matrix is what both the coverage maps and the checkerboard test
are built from, so they cannot disagree with each other.

Usage:
  tomo_rays.py --hypdir ../nll/loc --tag full --out ../tomo
"""
import argparse, glob, os, re, sys, time

import numpy as np
from scipy.sparse import csr_matrix, save_npz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from raytrace1d import load_layers, Model1D, Tracer

HERE = os.path.dirname(os.path.abspath(__file__))
GEO = re.compile(r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
QUAL = re.compile(r"RMS\s+([\d.eE+-]+)\s+Nphs\s+(\d+)\s+Gap\s+([\d.]+)")


class Grid:
    """Regular lat/lon/depth grid; cells are indexed C-order (lat, lon, depth)."""

    def __init__(self, lat0, lat1, lon0, lon1, zmax, dh_km=10.0, dz_km=10.0):
        self.lat0, self.lat1, self.lon0, self.lon1 = lat0, lat1, lon0, lon1
        self.latc = 0.5 * (lat0 + lat1)
        self.dlat = dh_km / 111.19
        self.dlon = dh_km / (111.19 * np.cos(np.radians(self.latc)))
        self.dz = dz_km
        self.nlat = int(np.ceil((lat1 - lat0) / self.dlat))
        self.nlon = int(np.ceil((lon1 - lon0) / self.dlon))
        self.nz = int(np.ceil(zmax / dz_km))
        self.n = self.nlat * self.nlon * self.nz
        self.dh_km, self.dz_km, self.zmax = dh_km, dz_km, zmax

    def index(self, lat, lon, z):
        i = np.floor((lat - self.lat0) / self.dlat).astype(int)
        j = np.floor((lon - self.lon0) / self.dlon).astype(int)
        k = np.floor(z / self.dz).astype(int)
        ok = (i >= 0) & (i < self.nlat) & (j >= 0) & (j < self.nlon) & (k >= 0) & (k < self.nz)
        return np.where(ok, (i * self.nlon + j) * self.nz + k, -1), ok

    def centres(self):
        i, j, k = np.unravel_index(np.arange(self.n), (self.nlat, self.nlon, self.nz))
        return (self.lat0 + (i + 0.5) * self.dlat,
                self.lon0 + (j + 0.5) * self.dlon,
                (k + 0.5) * self.dz)


def read_arrivals(hypdir, tag, stations, max_gap=360.0, max_rms=1.0):
    """(event index, station, phase, residual) for every weighted arrival."""
    files = sorted(f for f in glob.glob(f"{hypdir}/{tag}.*.grid0.loc.hyp") if "sum" not in f)
    events, arr = [], []
    for f in files:
        lat = lon = dep = None; rms = gap = None
        rows = []
        for line in open(f):
            if line.startswith("GEOGRAPHIC"):
                m = GEO.search(line)
                if m:
                    lat, lon, dep = (float(m.group(i)) for i in (1, 2, 3))
            elif line.startswith("QUALITY"):
                m = QUAL.search(line)
                if m:
                    rms, gap = float(m.group(1)), float(m.group(3))
            elif " > " in line and "GAU" in line:
                lf = line.split(" > ")[0].split()
                rf = line.split(" > ")[1].split()
                if float(rf[2]) <= 0:
                    continue
                rows.append((lf[0], lf[4], float(rf[1])))
        if lat is None or not rows:
            continue
        if gap is not None and gap > max_gap:
            continue
        if rms is not None and rms > max_rms:
            continue
        ei = len(events)
        events.append((lat, lon, max(dep, 0.0)))
        for sta, ph, res in rows:
            if sta in stations:
                arr.append((ei, sta, ph, res))
    return np.array(events), arr


def trace_all(events, arr, stations, grid, phase, vel_file, seg_km=1.0):
    layers = load_layers(vel_file, phase)
    tracer = Tracer(Model1D(layers, zmax=grid.zmax + 150.0))
    rows, cols, vals = [], [], []
    meta = []
    nfail = 0
    t0 = time.time()
    sel = [a for a in arr if a[2] == phase]
    for r, (ei, sta, ph, res) in enumerate(sel):
        elat, elon, edep = events[ei]
        slat, slon = stations[sta]
        dlat = (slat - elat) * 111.19
        dlon = (slon - elon) * 111.19 * np.cos(np.radians(elat))
        dist = np.hypot(dlat, dlon)
        got = tracer.ray_path(edep, max(dist, 0.05))
        if got is None:
            nfail += 1
            meta.append((ei, -1, res, np.nan))
            continue
        x, z, tt = got
        # resample the polyline at ~seg_km spacing
        s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(z)))])
        ns = max(int(s[-1] / seg_km), 2)
        su = np.linspace(0, s[-1], ns + 1)
        xu = np.interp(su, s, x); zu = np.interp(su, s, z)
        seglen = np.diff(su)
        xm = 0.5 * (xu[:-1] + xu[1:]); zm = 0.5 * (zu[:-1] + zu[1:])
        f = xm / max(dist, 1e-6)
        plat = elat + f * (slat - elat)
        plon = elon + f * (slon - elon)
        idx, ok = grid.index(plat, plon, zm)
        idx, seglen = idx[ok], seglen[ok]
        if not len(idx):
            nfail += 1
            meta.append((ei, -1, res, tt))
            continue
        order = np.argsort(idx)
        idx, seglen = idx[order], seglen[order]
        uniq, start = np.unique(idx, return_index=True)
        summed = np.add.reduceat(seglen, start)
        rows.append(np.full(len(uniq), len(meta)))
        cols.append(uniq)
        vals.append(summed)
        meta.append((ei, len(meta), res, tt))
        if (r + 1) % 5000 == 0:
            print(f"  {phase}: {r+1}/{len(sel)}  {time.time()-t0:.0f}s", flush=True)
    G = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(len(meta), grid.n))
    return G, np.array([(m[0], m[1], m[2], m[3]) for m in meta], dtype=float), nfail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hypdir", default=os.path.join(HERE, "..", "nll", "loc"))
    ap.add_argument("--tag", default="full")
    ap.add_argument("--stations", default=os.path.join(HERE, "..", "full",
                                                       "events_land_stations.csv"))
    ap.add_argument("--velP", default=os.path.join(HERE, "..", "nll", "nll_vel_full.in"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "tomo"))
    ap.add_argument("--lat0", type=float, default=-10.0)
    ap.add_argument("--lat1", type=float, default=-6.0)
    ap.add_argument("--lon0", type=float, default=108.8)
    ap.add_argument("--lon1", type=float, default=112.2)
    ap.add_argument("--zmax", type=float, default=200.0)
    ap.add_argument("--dh", type=float, default=10.0)
    ap.add_argument("--dz", type=float, default=10.0)
    ap.add_argument("--max-gap", type=float, default=360.0)
    ap.add_argument("--max-rms", type=float, default=1.0)
    a = ap.parse_args()

    import pandas as pd
    sdf = pd.read_csv(a.stations)
    stations = {r.id: (r.latitude, r.longitude) for r in sdf.itertuples()}
    grid = Grid(a.lat0, a.lat1, a.lon0, a.lon1, a.zmax, a.dh, a.dz)
    print(f"grid {grid.nlat} x {grid.nlon} x {grid.nz} = {grid.n} cells "
          f"({a.dh:g} km horizontal, {a.dz:g} km vertical)")

    events, arr = read_arrivals(a.hypdir, a.tag, stations, a.max_gap, a.max_rms)
    print(f"{len(events)} events, {len(arr)} weighted arrivals on {len(stations)} sites")

    os.makedirs(a.out, exist_ok=True)
    total_hits = np.zeros(grid.n); total_len = np.zeros(grid.n)
    for phase in ("P", "S"):
        G, meta, nfail = trace_all(events, arr, stations, grid, phase, a.velP)
        keep = meta[:, 1] >= 0
        print(f"{phase}: {int(keep.sum())} rays traced, {nfail} failed")
        save_npz(os.path.join(a.out, f"G_{phase}.npz"), G)
        np.save(os.path.join(a.out, f"meta_{phase}.npy"), meta)
        hits = np.asarray((G > 0).sum(axis=0)).ravel()
        lens = np.asarray(G.sum(axis=0)).ravel()
        np.savez(os.path.join(a.out, f"cov_{phase}.npz"), hits=hits, lens=lens)
        total_hits += hits; total_len += lens
    np.savez(os.path.join(a.out, "coverage.npz"),
             hits=total_hits, lens=total_len,
             nlat=grid.nlat, nlon=grid.nlon, nz=grid.nz,
             lat0=grid.lat0, lon0=grid.lon0, dlat=grid.dlat, dlon=grid.dlon,
             dz=grid.dz, dh_km=grid.dh_km)
    np.save(os.path.join(a.out, "events.npy"), events)
    print(f"cells with >=1 ray: {(total_hits>0).sum()} of {grid.n} "
          f"({100*(total_hits>0).mean():.1f}%)")
    print(f"cells with >=10 rays: {(total_hits>=10).sum()} "
          f"({100*(total_hits>=10).mean():.1f}%)")


if __name__ == "__main__":
    main()
