#!/usr/bin/env python3
"""Plot depth slices of P and S velocity anomalies from LOTOS inversion.

Generates a multi-panel figure showing velocity structure at key depths:
5, 10, 15, 20, 30, 40, 50 km. Each depth shows the P-wave anomaly with
S-wave comparison.

Usage:
  plot_velocity_slices.py --vp tomo_full/vp.npy --vs tomo_full/vs.npy
                          --out ../figures/velocity_slices.png
"""
import argparse, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import CenteredNorm

HERE = os.path.dirname(os.path.abspath(__file__))


def load_grid_meta():
    """Load grid metadata from tomo directory."""
    import sys
    sys.path.insert(0, HERE)

    # Hardcoded from the inversion setup
    return {
        'lat0': -10.0,
        'lon0': 108.8,
        'dlat': 111.19 / 900,  # degrees to km conversion
        'dlon': 111.19 * np.cos(np.radians(-8)) / 900,
        'dz': 10,
        'nlat': 45,
        'nlon': 38,
        'nz': 20,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vp", default="tomo_full/vp.npy", help="P-wave grid")
    ap.add_argument("--vs", default="tomo_full/vs.npy", help="S-wave grid")
    ap.add_argument("--out", default="../figures/velocity_slices.png",
                    help="Output figure")
    ap.add_argument("--dpi", type=int, default=150, help="Figure DPI")
    a = ap.parse_args()

    # Load velocity grids
    vp = np.load(a.vp)
    vs = np.load(a.vs)

    print(f"Loaded velocity grids: {vp.shape}")
    print(f"P slowness range: {vp[vp != 0].min():.5f} to {vp[vp != 0].max():.5f} s/km")
    print(f"S slowness range: {vs[vs != 0].min():.5f} to {vs[vs != 0].max():.5f} s/km")

    g = load_grid_meta()
    depths_km = np.arange(0, 20) * g['dz']  # 0, 10, 20, ..., 190 km
    lat = g['lat0'] + np.arange(g['nlat']) * 111.19 / 900
    lon = g['lon0'] + np.arange(g['nlon']) * 111.19 * np.cos(np.radians(-8)) / 900

    # Select depths to plot
    plot_depths = [5, 10, 15, 20, 30, 40, 50]
    depth_indices = [d // 10 for d in plot_depths if d < 200]

    nrows, ncols = len(depth_indices), 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4*nrows), dpi=a.dpi)
    if len(depth_indices) == 1:
        axes = axes.reshape(1, -1)

    # Color limits (percent velocity anomaly)
    vmin, vmax = -4, 4
    cmap = 'RdBu_r'  # Red = faster, Blue = slower
    norm = CenteredNorm(vcenter=0, halfrange=4)

    for i, depth_idx in enumerate(depth_indices):
        depth_km = plot_depths[i]

        # Extract slices
        vp_slice = vp[:, :, depth_idx]
        vs_slice = vs[:, :, depth_idx]

        # Convert slowness to velocity anomaly (%)
        # Reference slowness from 1-D model (approximate)
        u0p = 0.13  # 1/7.7 km/s
        u0s = 0.22  # 1/4.6 km/s
        vp_anom = (-vp_slice / u0p) * 100 if vp_slice.max() != 0 else vp_slice
        vs_anom = (-vs_slice / u0s) * 100 if vs_slice.max() != 0 else vs_slice

        # P-wave panel
        ax = axes[i, 0]
        im = ax.imshow(vp_anom, extent=[lon[0], lon[-1], lat[-1], lat[0]],
                        aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(f'P-wave anomaly, {depth_km} km', fontsize=12, fontweight='bold')
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.grid(True, alpha=0.3)

        # Add colorbar for first row
        if i == 0:
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02)
            cbar.set_label('Velocity anomaly (%)', fontsize=10)

        # S-wave panel
        ax = axes[i, 1]
        im = ax.imshow(vs_anom, extent=[lon[0], lon[-1], lat[-1], lat[0]],
                        aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(f'S-wave anomaly, {depth_km} km', fontsize=12, fontweight='bold')
        ax.set_xlabel('Longitude (°E)')
        ax.set_ylabel('Latitude (°N)')
        ax.grid(True, alpha=0.3)

        if i == 0:
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02)
            cbar.set_label('Velocity anomaly (%)', fontsize=10)

    plt.suptitle('MERAMEX 2004 Velocity Structure: Depth Slices', fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    plt.savefig(a.out, dpi=a.dpi, bbox_inches='tight')
    print(f"\nSaved: {a.out}")
    plt.close()

    # Also generate a single-panel high-resolution P-wave slice at 20 km
    fig, ax = plt.subplots(figsize=(12, 8), dpi=a.dpi)

    vp_slice_20km = vp[:, :, 2]  # depth index 2 = 20 km
    u0p = 0.13
    vp_anom_20km = (-vp_slice_20km / u0p) * 100

    im = ax.imshow(vp_anom_20km, extent=[lon[0], lon[-1], lat[-1], lat[0]],
                    aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title('P-wave Velocity Anomaly at 20 km Depth\nMERAMEX 2004 Full Campaign',
                  fontsize=14, fontweight='bold')
    ax.set_xlabel('Longitude (°E)', fontsize=12)
    ax.set_ylabel('Latitude (°N)', fontsize=12)
    ax.grid(True, alpha=0.3)
    cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label('Velocity anomaly (%)', fontsize=11)

    plt.tight_layout()
    fig_20km = a.out.replace('.png', '_20km.png')
    plt.savefig(fig_20km, dpi=a.dpi, bbox_inches='tight')
    print(f"Saved: {fig_20km}")
    plt.close()

    print("\n=== Velocity Slice Summary ===")
    print(f"Grid extent: {lat[0]:.2f}–{lat[-1]:.2f}°S, {lon[0]:.2f}–{lon[-1]:.2f}°E")
    print(f"Depth range: 0–190 km (10 km intervals)")
    print(f"Velocity anomaly scale: {vmin}% to {vmax}%")
    print(f"Red: faster velocities (denser/colder)")
    print(f"Blue: slower velocities (less dense/hotter)")


if __name__ == "__main__":
    main()
