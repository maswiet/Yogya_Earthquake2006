#!/usr/bin/env python3
"""Plot depth slices from a SIMULPS (simul2023) vlxyzltln.out velocity table.

Pivots on the native X/Y Cartesian grid (km from the origin), not on
lat/lon: simul2023's short-distance conversion scales east-west distance by
cos(latitude), so the same X column maps to a slightly different longitude
at every row -- pivoting on lat/lon directly fragments each column into as
many spurious columns as there are rows. Axis ticks are relabelled in
lat/lon using the corner values only, which is exact since X and Y are the
true rectangular grid.

Usage:
  plot_simulps_model.py --vlxyz ../simulps/inversion/vlxyzltln.out
      --out ../figures/simulps_velocity_slices.png
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLS = ["lat", "lon", "x", "y", "z", "vp", "vpvs", "vs",
        "vp_dws", "vp_dre", "vpvs_dws", "vpvs_dre"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vlxyz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-dws", type=float, default=5.0,
                    help="grey out nodes with fewer hits than this (unconstrained)")
    a = ap.parse_args()

    df = pd.read_csv(a.vlxyz, skiprows=2, sep=r"\s+", names=COLS)
    depths = sorted(df.z.unique())
    depths = [d for d in depths if d > 0]  # skip the surface/edge layer

    x_u = np.sort(df.x.unique())
    y_u = np.sort(df.y.unique())

    # lat/lon corner lookup for axis labelling only (exact at the grid corners)
    corner = df[(df.z == depths[0])]
    lon_at = {x: corner.loc[(corner.x == x) & (corner.y == y_u[0]), "lon"].iloc[0] for x in x_u}
    lat_at = {y: corner.loc[(corner.y == y) & (corner.x == x_u[0]), "lat"].iloc[0] for y in y_u}

    fig, axes = plt.subplots(2, len(depths), figsize=(3.1 * len(depths), 6.4))

    im = imD = None
    for i, z in enumerate(depths):
        sub = df[df.z == z]
        piv_vp = sub.pivot(index="y", columns="x", values="vp").reindex(index=y_u, columns=x_u)
        piv_dws = sub.pivot(index="y", columns="x", values="vp_dws").reindex(index=y_u, columns=x_u)
        vp0 = piv_vp.values[piv_dws.values <= a.min_dws]
        vp0 = np.nanmedian(vp0) if vp0.size else np.nanmedian(piv_vp.values)
        pct = 100.0 * (piv_vp.values - vp0) / vp0
        mask = piv_dws.values <= a.min_dws

        extent = [x_u[0], x_u[-1], y_u[-1], y_u[0]]

        axT = axes[0, i]
        im = axT.imshow(np.ma.masked_where(mask, pct), extent=extent,
                        aspect="equal", cmap="RdBu_r", vmin=-6, vmax=6,
                        interpolation="nearest", origin="upper")
        axT.set_facecolor("0.85")
        axT.set_title(f"{z:.0f} km", fontsize=11, fontweight="bold")
        axT.set_xticks(x_u); axT.set_xticklabels([f"{lon_at[x]:.1f}" for x in x_u], rotation=45, fontsize=7)
        axT.set_yticks(y_u); axT.set_yticklabels([f"{lat_at[y]:.1f}" for y in y_u], fontsize=7)
        if i == 0:
            axT.set_ylabel("Vp anomaly (%)\nLatitude (°N)")

        axB = axes[1, i]
        imD = axB.imshow(piv_dws.values, extent=extent, aspect="equal", cmap="viridis",
                         norm=matplotlib.colors.LogNorm(vmin=1, vmax=max(df.vp_dws.max(), 2)),
                         origin="upper")
        axB.set_xticks(x_u); axB.set_xticklabels([f"{lon_at[x]:.1f}" for x in x_u], rotation=45, fontsize=7)
        axB.set_yticks(y_u); axB.set_yticklabels([f"{lat_at[y]:.1f}" for y in y_u], fontsize=7)
        axB.set_xlabel("Longitude (°E)")
        if i == 0:
            axB.set_ylabel("Ray density (DWS)\nLatitude (°N)")

    fig.subplots_adjust(right=0.90, top=0.86, wspace=0.35, hspace=0.45)
    cbar_ax = fig.add_axes([0.92, 0.55, 0.015, 0.32])
    fig.colorbar(im, cax=cbar_ax, label="Vp anomaly (%)")
    cbar_ax2 = fig.add_axes([0.92, 0.12, 0.015, 0.32])
    fig.colorbar(imD, cax=cbar_ax2, label="DWS (ray hits)")

    fig.suptitle("SIMULPS (simul2023) Tomographic Inversion — MERAMEX 2004\n"
                 f"{len(x_u)}×{len(y_u)}×{len(depths)} grid (native X/Y, ~94 km spacing), "
                 "grey = fewer than 5 ray hits (unconstrained)",
                 fontsize=12, fontweight="bold")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=160, bbox_inches="tight")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
