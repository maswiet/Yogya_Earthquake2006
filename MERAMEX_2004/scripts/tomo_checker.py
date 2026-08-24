#!/usr/bin/env python3
"""Checkerboard resolution test on the MERAMEX ray set.

Uses the same path-length matrices the coverage maps are built from, so the two
cannot disagree. The inversion carries the terms that actually compete with
velocity structure in a real local-earthquake tomography:

  * one origin-time parameter per event (the dominant trade-off with velocity),
  * separate P and S slowness fields,
  * damping and first-difference smoothing,
  * Gaussian noise at the pick uncertainties used in the location (0.1 s P,
    0.2 s S).

A noise-free test with fixed sources would look far better than anything the
data can actually deliver.

Usage:
  tomo_checker.py --cell 30 --amp 7
"""
import argparse, os, sys

import numpy as np
from scipy.sparse import load_npz, hstack, vstack, csr_matrix, eye
from scipy.sparse.linalg import lsqr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from raytrace1d import load_layers


def load_grid(tomo):
    z = np.load(os.path.join(tomo, "coverage.npz"))
    return {k: z[k] for k in z.files}


def cell_centres(g):
    nlat, nlon, nz = int(g["nlat"]), int(g["nlon"]), int(g["nz"])
    i, j, k = np.unravel_index(np.arange(nlat * nlon * nz), (nlat, nlon, nz))
    lat = g["lat0"] + (i + 0.5) * g["dlat"]
    lon = g["lon0"] + (j + 0.5) * g["dlon"]
    dep = (k + 0.5) * g["dz"]
    return lat, lon, dep, i, j, k


def checkerboard(g, cell_km, dz_cell_km, amp_pct):
    lat, lon, dep, i, j, k = cell_centres(g)
    x = (lon - g["lon0"]) * 111.19 * np.cos(np.radians(lat))
    y = (lat - g["lat0"]) * 111.19
    sx = np.floor(x / cell_km).astype(int)
    sy = np.floor(y / cell_km).astype(int)
    sz = np.floor(dep / dz_cell_km).astype(int)
    return (amp_pct / 100.0) * np.where((sx + sy + sz) % 2 == 0, 1.0, -1.0)


def smoothing_operator(g, mask_idx):
    nlat, nlon, nz = int(g["nlat"]), int(g["nlon"]), int(g["nz"])
    pos = -np.ones(nlat * nlon * nz, dtype=int)
    pos[mask_idx] = np.arange(len(mask_idx))
    i, j, k = np.unravel_index(mask_idx, (nlat, nlon, nz))
    rows, cols, vals = [], [], []
    r = 0
    for di, dj, dk in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        ii, jj, kk = i + di, j + dj, k + dk
        ok = (ii < nlat) & (jj < nlon) & (kk < nz)
        flat = np.where(ok, ((ii * nlon + jj) * nz + kk), 0)
        nb = np.where(ok, pos[flat], -1)
        good = ok & (nb >= 0)
        n = int(good.sum())
        if not n:
            continue
        rows.append(np.repeat(np.arange(r, r + n), 2))
        cols.append(np.column_stack([np.arange(len(mask_idx))[good], nb[good]]).ravel())
        vals.append(np.tile([1.0, -1.0], n))
        r += n
    return csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(r, len(mask_idx)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tomo", default=os.path.join(HERE, "..", "tomo"))
    ap.add_argument("--cell", type=float, default=30.0)
    ap.add_argument("--zcell", type=float, default=30.0)
    ap.add_argument("--amp", type=float, default=7.0)
    ap.add_argument("--minhit", type=int, default=5)
    ap.add_argument("--damp", type=float, default=1.0)
    ap.add_argument("--smooth", type=float, default=3.0)
    ap.add_argument("--noiseP", type=float, default=0.10)
    ap.add_argument("--noiseS", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--nray", type=int, default=0,
                    help="randomly keep this many rays in total; used to emulate "
                         "a smaller data set (Koulakov 2007 had ~13,000)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    g = load_grid(a.tomo)
    nz = int(g["nz"])
    ncell = int(g["nlat"]) * int(g["nlon"]) * nz
    lat, lon, dep, _, _, _ = cell_centres(g)

    GP = load_npz(os.path.join(a.tomo, "G_P.npz")).tocsr()
    GS = load_npz(os.path.join(a.tomo, "G_S.npz")).tocsr()
    mP = np.load(os.path.join(a.tomo, "meta_P.npy"))
    mS = np.load(os.path.join(a.tomo, "meta_S.npy"))
    GP, GS = GP[mP[:, 1] >= 0], GS[mS[:, 1] >= 0]
    evP = mP[mP[:, 1] >= 0, 0].astype(int)
    evS = mS[mS[:, 1] >= 0, 0].astype(int)

    if a.nray and a.nray < GP.shape[0] + GS.shape[0]:
        rs = np.random.default_rng(a.seed + 1)
        frac = a.nray / (GP.shape[0] + GS.shape[0])
        kp = rs.random(GP.shape[0]) < frac
        ks = rs.random(GS.shape[0]) < frac
        GP, evP = GP[kp], evP[kp]
        GS, evS = GS[ks], evS[ks]
        print(f"subsampled to {GP.shape[0]} P + {GS.shape[0]} S = "
              f"{GP.shape[0]+GS.shape[0]} rays")

    hits = np.asarray((GP > 0).sum(axis=0)).ravel() + np.asarray((GS > 0).sum(axis=0)).ravel()
    active = np.nonzero(hits >= a.minhit)[0]
    print(f"{len(active)} of {ncell} cells carry >= {a.minhit} rays "
          f"({100*len(active)/ncell:.1f}%)")

    vel = os.path.join(HERE, "..", "nll", "nll_vel_full.in")
    def u0(phase):
        lay = load_layers(vel, phase)
        tops = np.array([l[0] for l in lay]); v = np.array([l[1] for l in lay])
        return 1.0 / np.interp(dep, tops, v)
    u0P, u0S = u0("P"), u0("S")

    true_v = checkerboard(g, a.cell, a.zcell, a.amp)
    mtrueP = (-u0P * true_v)[active]
    mtrueS = (-u0S * true_v)[active]

    GPa, GSa = GP[:, active], GS[:, active]
    nev = int(max(evP.max(), evS.max())) + 1
    OP = csr_matrix((np.ones(len(evP)), (np.arange(len(evP)), evP)), shape=(len(evP), nev))
    OS = csr_matrix((np.ones(len(evS)), (np.arange(len(evS)), evS)), shape=(len(evS), nev))
    nact = len(active)
    A = vstack([
        hstack([GPa, csr_matrix((GPa.shape[0], nact)), OP]),
        hstack([csr_matrix((GSa.shape[0], nact)), GSa, OS]),
    ]).tocsr()
    print(f"design matrix {A.shape[0]} x {A.shape[1]} ({A.nnz/1e6:.1f}M nonzeros), "
          f"{nev} origin-time terms")

    rng = np.random.default_rng(a.seed)
    d = np.concatenate([GPa @ mtrueP, GSa @ mtrueS])
    d += np.concatenate([rng.normal(0, a.noiseP, GPa.shape[0]),
                         rng.normal(0, a.noiseS, GSa.shape[0])])

    L = smoothing_operator(g, active)
    reg = vstack([
        hstack([a.damp * eye(nact), csr_matrix((nact, nact + nev))]),
        hstack([csr_matrix((nact, nact)), a.damp * eye(nact), csr_matrix((nact, nev))]),
        hstack([a.smooth * L, csr_matrix((L.shape[0], nact + nev))]),
        hstack([csr_matrix((L.shape[0], nact)), a.smooth * L, csr_matrix((L.shape[0], nev))]),
    ]).tocsr()
    Afull = vstack([A, reg]).tocsr()
    dfull = np.concatenate([d, np.zeros(reg.shape[0])])

    print("running LSQR ...", flush=True)
    m = lsqr(Afull, dfull, atol=1e-8, btol=1e-8, iter_lim=400)[0]
    outP = np.full(ncell, np.nan); outS = np.full(ncell, np.nan)
    outP[active] = -m[:nact] / u0P[active]
    outS[active] = -m[nact:2 * nact] / u0S[active]

    path = a.out or os.path.join(a.tomo, f"checker_{int(a.cell)}km.npz")
    np.savez(path, true_v=true_v, recP=outP, recS=outS, active=active,
             lat=lat, lon=lon, dep=dep, cell=a.cell, zcell=a.zcell, amp=a.amp,
             nlat=g["nlat"], nlon=g["nlon"], nz=g["nz"], dz=g["dz"],
             lat0=g["lat0"], lon0=g["lon0"], dlat=g["dlat"], dlon=g["dlon"],
             hits=hits)

    print(f"\ncheckerboard {a.cell:g} km / {a.zcell:g} km vertical, +/-{a.amp:g}%")
    print("%10s %8s %8s %8s %8s %8s" % ("depth(km)", "ncell", "corrP", "ampP", "corrS", "ampS"))
    for k in range(nz):
        z0, z1 = k * g["dz"], (k + 1) * g["dz"]
        sel = active[(dep[active] >= z0) & (dep[active] < z1)]
        if len(sel) < 20:
            continue
        t = true_v[sel]
        row = [f"{int(z0)}-{int(z1)}".rjust(10), f"{len(sel):8d}"]
        for rec in (outP, outS):
            r = rec[sel]; good = np.isfinite(r)
            if good.sum() < 20:
                row += ["     n/a", "     n/a"]; continue
            c = np.corrcoef(t[good], r[good])[0, 1]
            amp = np.sum(t[good] * r[good]) / np.sum(t[good] ** 2)
            row += [f"{c:8.2f}", f"{100*amp:7.0f}%"]
        print(" ".join(row))


if __name__ == "__main__":
    main()
