#!/usr/bin/env python3
"""Prepare compact JSON payloads for the interactive 3-D velocity artifact.

Writes three flat-array JSON files (cells, events, coastline) sized for
inline embedding in a self-contained HTML/canvas visualization — see
build_velocity_3d.py.

Usage:
  gen_velocity_3d_data.py --out ../figures/velocity_3d_data
"""
import argparse, os, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

BOX_LAT = (-8.01, -7.87)
BOX_LON = (110.35, 110.53)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tomo", default=os.path.join(ROOT, "tomo"))
    ap.add_argument("--vp", default=os.path.join(ROOT, "tomo_full", "vp.npy"))
    ap.add_argument("--vs", default=os.path.join(ROOT, "tomo_full", "vs.npy"))
    ap.add_argument("--cat", default=os.path.join(ROOT, "full", "catalog_obs_combined.csv"))
    ap.add_argument("--coast", default=os.path.join(ROOT, "data", "coastline.xy"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "figures", "velocity_3d_data"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    z = np.load(os.path.join(a.tomo, "coverage.npz"))
    g = {k: z[k] for k in z.files}
    lat0, lon0, dlat, dlon, dz = g["lat0"], g["lon0"], g["dlat"], g["dlon"], g["dz"]
    nlat, nlon, nz = int(g["nlat"]), int(g["nlon"]), int(g["nz"])

    lat = lat0 + (np.arange(nlat) + 0.5) * dlat
    lon = lon0 + (np.arange(nlon) + 0.5) * dlon
    dep = (np.arange(nz) + 0.5) * dz

    vp = np.load(a.vp)
    vs = np.load(a.vs)
    u0p = 1.0 / np.interp(dep, [0, 200], [5.8, 8.0])
    u0s = 1.0 / np.interp(dep, [0, 200], [3.3, 4.6])

    LA, LO, DE = np.meshgrid(lat, lon, dep, indexing="ij")
    active = vp != 0
    print("active cells:", active.sum())

    vp_pct = np.zeros_like(vp)
    vs_pct = np.zeros_like(vs)
    for k in range(nz):
        m = vp[:, :, k] != 0
        vp_pct[:, :, k][m] = -100.0 * vp[:, :, k][m] / u0p[k]
        m2 = vs[:, :, k] != 0
        vs_pct[:, :, k][m2] = -100.0 * vs[:, :, k][m2] / u0s[k]

    cells = np.column_stack([
        LA[active].round(3), LO[active].round(3), DE[active].round(1),
        vp_pct[active].round(2), vs_pct[active].round(2),
    ])
    with open(os.path.join(a.out, "cells.json"), "w") as f:
        json.dump(cells.flatten().tolist(), f, separators=(",", ":"))
    print(f"cells: {cells.shape[0]} rows -> cells.json "
          f"({os.path.getsize(os.path.join(a.out,'cells.json'))/1024:.1f} KB)")

    cat = pd.read_csv(a.cat)
    cat = cat[cat.quality_pass | (cat.errh_km <= 15)]
    cat = cat[cat.latitude.between(lat0, lat0 + nlat * dlat) &
              cat.longitude.between(lon0, lon0 + nlon * dlon)]
    inbox = cat.latitude.between(*BOX_LAT) & cat.longitude.between(*BOX_LON)
    ev = np.column_stack([
        cat.latitude.round(3), cat.longitude.round(3), cat.depth.round(2), inbox.astype(int),
    ])
    with open(os.path.join(a.out, "events.json"), "w") as f:
        json.dump(ev.flatten().tolist(), f, separators=(",", ":"))
    print(f"events: {len(cat)} ({inbox.sum()} in Opak box) -> events.json "
          f"({os.path.getsize(os.path.join(a.out,'events.json'))/1024:.1f} KB)")

    lat_max, lon_max = lat0 + nlat * dlat, lon0 + nlon * dlon
    pad = 0.3
    segments, cur = [], []
    if os.path.exists(a.coast):
        with open(a.coast) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if len(cur) > 1:
                        segments.append(cur)
                    cur = []
                    continue
                parts = line.split()
                if len(parts) == 2:
                    lo, la = float(parts[0]), float(parts[1])
                    cur.append((la, lo))
        if len(cur) > 1:
            segments.append(cur)

    clipped = []
    for seg in segments:
        arr = np.array(seg)
        m = ((arr[:, 0] >= lat0 - pad) & (arr[:, 0] <= lat_max + pad) &
             (arr[:, 1] >= lon0 - pad) & (arr[:, 1] <= lon_max + pad))
        if m.sum() < 2:
            continue
        sub = arr[m][::2]
        if len(sub) >= 2:
            clipped.append(sub.round(3).tolist())
    with open(os.path.join(a.out, "coast.json"), "w") as f:
        json.dump(clipped, f, separators=(",", ":"))
    print(f"coastline: {len(clipped)} segments -> coast.json "
          f"({os.path.getsize(os.path.join(a.out,'coast.json'))/1024:.1f} KB)")

    print("\ngrid bounds (for the HTML build):")
    print(f"  LAT0={lat0}, LAT1={lat_max}, LON0={lon0}, LON1={lon_max}, DMAX={nz*dz}")


if __name__ == "__main__":
    main()
