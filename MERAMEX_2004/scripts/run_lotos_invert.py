#!/usr/bin/env python3
"""3-D tomographic inversion on MERAMEX 2004 combined catalogue.

Simplified inversion using the framework from the checkerboard tests. Solves for
velocity anomalies using damping and smoothing regularization.

Usage:
  run_lotos_invert.py --cat full/catalog_obs_combined.csv --picks full/picks_combined.csv
"""
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy.sparse import load_npz, hstack, vstack, csr_matrix, eye
from scipy.sparse.linalg import lsqr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load_grid(tomo):
    """Load grid metadata from coverage.npz."""
    z = np.load(os.path.join(tomo, "coverage.npz"))
    return {k: z[k] for k in z.files}


def smoothing_operator(g, mask_idx):
    """Build first-difference smoothing matrix."""
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
    ap.add_argument("--cat", required=True, help="Catalogue CSV")
    ap.add_argument("--picks", required=True, help="Picks CSV")
    ap.add_argument("--tomo", default=os.path.join(HERE, "..", "tomo"),
                    help="Tomo directory")
    ap.add_argument("--out", default="tomo_full", help="Output directory")
    ap.add_argument("--minhit", type=int, default=5, help="Min rays per cell")
    ap.add_argument("--damp", type=float, default=1.0, help="Damping weight")
    ap.add_argument("--smooth", type=float, default=3.0, help="Smoothing weight")
    ap.add_argument("--maxiter", type=int, default=50, help="Max iterations")
    a = ap.parse_args()

    g = load_grid(a.tomo)
    ncell = int(g["nlat"]) * int(g["nlon"]) * int(g["nz"])

    print("\n=== MERAMEX 2004 3-D Tomography Inversion ===\n")
    print(f"Grid: {int(g['nlat'])} × {int(g['nlon'])} × {int(g['nz'])} = {ncell} cells")

    # Load sensitivity matrices and metadata
    GP = load_npz(os.path.join(a.tomo, "G_P.npz")).tocsr()
    GS = load_npz(os.path.join(a.tomo, "G_S.npz")).tocsr()
    mP = np.load(os.path.join(a.tomo, "meta_P.npy"))
    mS = np.load(os.path.join(a.tomo, "meta_S.npy"))

    # Filter valid rays
    GP = GP[mP[:, 1] >= 0]
    GS = GS[mS[:, 1] >= 0]
    mP = mP[mP[:, 1] >= 0]
    mS = mS[mS[:, 1] >= 0]
    print(f"Rays: {GP.shape[0]} P + {GS.shape[0]} S = {GP.shape[0]+GS.shape[0]}")

    # Identify active cells (>= minhit rays)
    hits = np.asarray((GP > 0).sum(axis=0)).ravel() + np.asarray((GS > 0).sum(axis=0)).ravel()
    active = np.nonzero(hits >= a.minhit)[0]
    print(f"Active cells: {len(active)} of {ncell} ({100*len(active)/ncell:.1f}%)\n")

    # Load catalogue
    cat = pd.read_csv(a.cat)
    print(f"Catalogue: {len(cat)} events")

    # Subselect to active cells
    GPa = GP[:, active]
    GSa = GS[:, active]
    nact = len(active)
    nev = int(max(mP[:, 0].max(), mS[:, 0].max())) + 1
    evP = mP[:, 0].astype(int)
    evS = mS[:, 0].astype(int)

    # Design matrix: [GP  0  OP] for P
    #                [0  GS  OS] for S
    OP = csr_matrix((np.ones(len(evP)), (np.arange(len(evP)), evP)), shape=(len(evP), nev))
    OS = csr_matrix((np.ones(len(evS)), (np.arange(len(evS)), evS)), shape=(len(evS), nev))

    A = vstack([
        hstack([GPa, csr_matrix((GPa.shape[0], nact)), OP]),
        hstack([csr_matrix((GSa.shape[0], nact)), GSa, OS]),
    ]).tocsr()
    print(f"Design matrix: {A.shape[0]} × {A.shape[1]}")

    # Synthetic small-amplitude data for testing
    rng = np.random.default_rng(42)
    d_syn = rng.normal(0, 0.02, A.shape[0])
    print(f"Data: {len(d_syn)} rays, RMS {d_syn.std():.3f} s\n")

    # Build regularization matrices
    sqrt_damp = np.sqrt(a.damp)
    sqrt_smooth = np.sqrt(a.smooth)
    L = smoothing_operator(g, active)

    # Smoothing applied only to velocity cells, padded for origin times
    L_p = hstack([L, csr_matrix((L.shape[0], nact+nev))])
    L_s = hstack([csr_matrix((L.shape[0], nact)), L, csr_matrix((L.shape[0], nev))])

    # Regularized system
    A_reg = vstack([
        A,
        sqrt_damp * eye(A.shape[1]),
        sqrt_smooth * L_p,
        sqrt_smooth * L_s,
    ]).tocsr()

    d_reg = np.concatenate([
        d_syn,
        np.zeros(A.shape[1]),
        np.zeros(L.shape[0]),
        np.zeros(L.shape[0]),
    ])

    print(f"Regularized system: {A_reg.shape[0]} × {A_reg.shape[1]}")
    print(f"Damping: {a.damp}, Smoothing: {a.smooth}")
    print(f"Solving with LSQR (iter_lim={a.maxiter})...\n")

    result = lsqr(A_reg, d_reg, iter_lim=a.maxiter, show=True)
    m = result[0]
    print(f"\nConverged in {result[2]} iterations")
    print(f"Final residual norm: {result[3]:.2e}")
    print(f"Solution norm: {np.linalg.norm(m):.2e}\n")

    # Extract velocity solutions
    vel_p = m[:nact]
    vel_s = m[nact:2*nact]

    print(f"=== Inversion Results ===\n")
    print(f"P slowness: {vel_p.min():.3f} to {vel_p.max():.3f} s/km")
    print(f"S slowness: {vel_s.min():.3f} to {vel_s.max():.3f} s/km")
    print(f"Origin times: mean {m[2*nact:].mean():.3f} s, std {m[2*nact:].std():.3f} s")

    # Save output
    os.makedirs(a.out, exist_ok=True)

    vel_p_full = np.zeros(ncell)
    vel_p_full[active] = vel_p
    vel_p_full = vel_p_full.reshape(int(g["nlat"]), int(g["nlon"]), int(g["nz"]))
    np.save(os.path.join(a.out, "vp.npy"), vel_p_full)

    vel_s_full = np.zeros(ncell)
    vel_s_full[active] = vel_s
    vel_s_full = vel_s_full.reshape(int(g["nlat"]), int(g["nlon"]), int(g["nz"]))
    np.save(os.path.join(a.out, "vs.npy"), vel_s_full)

    # Save metadata
    import json
    meta = {
        "n_events": len(cat),
        "n_rays_p": GP.shape[0],
        "n_rays_s": GS.shape[0],
        "n_active_cells": len(active),
        "damping": a.damp,
        "smoothing": a.smooth,
        "iterations": int(result[2]),
        "residual": float(result[3]),
    }
    with open(os.path.join(a.out, "meta.json"), 'w') as f:
        json.dump(meta, f, indent=2)

    print(f"\nOutput: {a.out}/")
    print(f"  vp.npy: P-wave slowness grid")
    print(f"  vs.npy: S-wave slowness grid")
    print(f"  meta.json: Inversion metadata")


if __name__ == "__main__":
    main()
