#!/usr/bin/env python3
"""Plot depth slices from a SIMULPS (simul2023) vlxyzltln.out velocity table.

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
    depths = [d for d in depths if d > 0]  # skip the -5/edge layer if present

    nlat = df.lat.nunique()
    nlon = df.lon.nunique()
    lat_u = np.sort(df.lat.unique())
    lon_u = np.sort(df.lon.unique())

    fig, axes = plt.subplots(2, len(depths), figsize=(3.1 * len(depths), 6.4))

    for i, z in enumerate(depths):
        sub = df[df.z == z]
        piv_vp = sub.pivot(index="lat", columns="lon", values="vp").reindex(
            index=lat_u, columns=lon_u)
        piv_dws = sub.pivot(index="lat", columns="lon", values="vp_dws").reindex(
            index=lat_u, columns=lon_u)
        vp0 = piv_vp.values[piv_dws.values <= a.min_dws]
        vp0 = np.nanmedian(vp0) if vp0.size else np.nanmedian(piv_vp.values)
        pct = 100.0 * (piv_vp.values - vp0) / vp0
        mask = piv_dws.values <= a.min_dws

        axT = axes[0, i]
        im = axT.imshow(np.ma.masked_where(mask, pct), extent=[lon_u[0], lon_u[-1], lat_u[-1], lat_u[0]],
                        aspect="auto", cmap="RdBu_r", vmin=-6, vmax=6, interpolation="nearest")
        axT.set_facecolor("0.85")
        axT.set_title(f"{z:.0f} km", fontsize=11, fontweight="bold")
        if i == 0:
            axT.set_ylabel("Vp anomaly (%)\nLatitude")
        axT.set_xlabel("")

        axB = axes[1, i]
        imD = axB.imshow(piv_dws.values, extent=[lon_u[0], lon_u[-1], lat_u[-1], lat_u[0]],
                         aspect="auto", cmap="viridis", norm=matplotlib.colors.LogNorm(
                             vmin=1, vmax=max(df.vp_dws.max(), 2)))
        axB.set_xlabel("Longitude (°E)")
        if i == 0:
            axB.set_ylabel("Ray density (DWS)\nLatitude")

    cbar_ax = fig.add_axes([0.92, 0.55, 0.015, 0.35])
    fig.colorbar(im, cax=cbar_ax, label="Vp anomaly (%)")
    cbar_ax2 = fig.add_axes([0.92, 0.1, 0.015, 0.35])
    fig.colorbar(imD, cax=cbar_ax2, label="DWS (ray hits)")

    fig.suptitle("SIMULPS (simul2023) Tomographic Inversion — MERAMEX 2004\n"
                 f"{df.lat.nunique()}×{df.lon.nunique()}×{len(depths)} grid, "
                 "grey = fewer than 5 ray hits (unconstrained)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.90, 0.93])

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=160, bbox_inches="tight")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
