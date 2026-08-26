#!/usr/bin/env python3
"""Vertical cross-sections through the LOTOS velocity model.

Trench-perpendicular (N-S) profiles at fixed longitudes, showing P-wave
velocity anomaly with hypocentres from the combined catalogue projected onto
each swath. Companion to plot_slab_section.py (which shows raw hypocentres
only); this adds the tomographic structure the OBS network buys.

Usage:
  plot_velocity_sections.py --vp ../tomo_full/vp.npy --cat ../full/catalog_obs_combined.csv
      --profiles 109.6,110.5,111.4 --halfwidth 0.5 --out ../figures/velocity_sections.png
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))


def load_grid(tomo):
    z = np.load(os.path.join(tomo, "coverage.npz"))
    return {k: z[k] for k in z.files}


def cell_axes(g):
    """1-D coordinate arrays for cell centres along each grid axis."""
    lat = g["lat0"] + (np.arange(int(g["nlat"])) + 0.5) * g["dlat"]
    lon = g["lon0"] + (np.arange(int(g["nlon"])) + 0.5) * g["dlon"]
    dep = (np.arange(int(g["nz"])) + 0.5) * g["dz"]
    return lat, lon, dep


def to_pct(u, u0):
    """Slowness perturbation -> velocity percent anomaly, robust to u=0 (masked cell)."""
    out = np.full_like(u, np.nan)
    m = u0 > 0
    out[m] = -100.0 * u[m] / u0
    return out


def u0_profile(dep, v_top, v_bot, z_top, z_bot):
    """Piecewise-linear background slowness (crust->mantle) at each depth."""
    v = np.interp(dep, [z_top, z_bot], [v_top, v_bot])
    return 1.0 / v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tomo", default=os.path.join(ROOT, "tomo"))
    ap.add_argument("--vp", default=os.path.join(ROOT, "tomo_full", "vp.npy"))
    ap.add_argument("--cat", default=os.path.join(ROOT, "full", "catalog_obs_combined.csv"))
    ap.add_argument("--profiles", default="109.6,110.5,111.4")
    ap.add_argument("--halfwidth", type=float, default=0.5)
    ap.add_argument("--zmax", type=float, default=200.0)
    ap.add_argument("--out", default=os.path.join(ROOT, "figures", "velocity_sections.png"))
    a = ap.parse_args()

    g = load_grid(a.tomo)
    lat, lon, dep = cell_axes(g)
    vp = np.load(a.vp)  # (nlat, nlon, nz) slowness perturbation

    u0 = u0_profile(dep, 1/5.8, 1/8.0, 0, 200)  # ~crust to upper-mantle P velocity
    vp_pct = np.stack([to_pct(vp[:, :, k], u0[k]) for k in range(len(dep))], axis=-1)

    cat = pd.read_csv(a.cat)
    cat = cat[cat.quality_pass | (cat.errh_km <= 15)]  # keep reasonably located events

    lons = [float(v) for v in a.profiles.split(",")]
    fig, axes = plt.subplots(len(lons), 1, figsize=(9, 4.0 * len(lons)), sharex=True)
    if len(lons) == 1:
        axes = [axes]

    vmax = 4.0
    for i, (ax, lon0) in enumerate(zip(axes, lons)):
        # nearest grid column for this profile longitude
        j = int(np.argmin(np.abs(lon - lon0)))
        section = vp_pct[:, j, :].T  # (nz, nlat)

        im = ax.imshow(section, extent=[lat[0], lat[-1], dep[-1], dep[0]],
                        aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        interpolation="bilinear")

        band = cat[(cat.longitude - lon0).abs() <= a.halfwidth]
        ax.scatter(band.latitude, band.depth, s=10, facecolor="none",
                   edgecolor="k", lw=0.4, alpha=0.75, zorder=5)

        ax.set_ylim(a.zmax, 0)
        ax.set_xlim(lat[0], lat[-1])
        ax.set_ylabel("Depth (km)")
        ax.text(0.02, 0.92, f"{lon0:g}°E  ±{a.halfwidth:g}°   n={len(band)} events",
                transform=ax.transAxes, fontsize=10, fontweight="bold",
                va="top", bbox=dict(fc="white", alpha=0.75, ec="none", pad=2))
        ax.grid(alpha=0.25, lw=0.4)

    axes[-1].set_xlabel("Latitude (°N)")
    cbar = fig.colorbar(im, ax=axes, orientation="vertical", pad=0.015, aspect=35)
    cbar.set_label("P-wave velocity anomaly (%)")
    fig.suptitle("MERAMEX 2004 — Trench-Perpendicular Velocity Sections\n"
                 "Wadati–Benioff zone imaged with land + OBS rays", fontsize=13, fontweight="bold")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=170, bbox_inches="tight")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
