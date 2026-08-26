#!/usr/bin/env python3
"""Compare this study's recovered slab structure to Koulakov et al. (2007).

Left panel: our 110.5°E P-velocity cross-section (land+OBS rays), the profile
closest to Koulakov's central-Java transect.

Right panel: a schematic redrawn from Koulakov et al. (2007)'s reported
parameters (NOT a reproduction of their figure) — slab dip, Wadati-Benioff
extent and forearc-mantle-wedge geometry, for a like-for-like read of dip
angle and wedge thickness against our own recovery.

Bottom: data-volume and resolution metrics side by side (from tomo/ASSESSMENT.md).

Usage:
  plot_koulakov_comparison.py --out ../figures/koulakov_comparison.png
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

# Koulakov et al. (2007), JGR 112, B05310 — reported parameters for the
# central Java transect (approximate, digitised from the text/figures'
# reported dip and depth ranges, not the gridded velocity model itself).
KOULAKOV = dict(
    dip_deg=50.0,           # slab dip beneath central Java
    slab_top_km=(0, 15),    # trench-proximal top depth
    wedge_thickness_km=80,  # forearc mantle wedge thickness
    max_depth_km=300,       # deepest resolved slab
    events=292,
    rays=13000,
)


def load_grid(tomo):
    z = np.load(os.path.join(tomo, "coverage.npz"))
    return {k: z[k] for k in z.files}


def cell_axes(g):
    lat = g["lat0"] + (np.arange(int(g["nlat"])) + 0.5) * g["dlat"]
    lon = g["lon0"] + (np.arange(int(g["nlon"])) + 0.5) * g["dlon"]
    dep = (np.arange(int(g["nz"])) + 0.5) * g["dz"]
    return lat, lon, dep


def to_pct(u, u0):
    out = np.full_like(u, np.nan)
    m = u0 > 0
    out[m] = -100.0 * u[m] / u0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tomo", default=os.path.join(ROOT, "tomo"))
    ap.add_argument("--vp", default=os.path.join(ROOT, "tomo_full", "vp.npy"))
    ap.add_argument("--cat", default=os.path.join(ROOT, "full", "catalog_obs_combined.csv"))
    ap.add_argument("--profile_lon", type=float, default=110.5)
    ap.add_argument("--halfwidth", type=float, default=0.5)
    ap.add_argument("--zmax", type=float, default=250.0)
    ap.add_argument("--out", default=os.path.join(ROOT, "figures", "koulakov_comparison.png"))
    a = ap.parse_args()

    g = load_grid(a.tomo)
    lat, lon, dep = cell_axes(g)
    vp = np.load(a.vp)
    u0 = 1.0 / np.interp(dep, [0, 200], [5.8, 8.0])
    vp_pct = np.stack([to_pct(vp[:, :, k], u0[k]) for k in range(len(dep))], axis=-1)

    cat = pd.read_csv(a.cat)
    cat = cat[cat.quality_pass | (cat.errh_km <= 15)]

    j = int(np.argmin(np.abs(lon - a.profile_lon)))
    section = vp_pct[:, j, :].T

    fig = plt.figure(figsize=(14, 7.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.1], hspace=0.35, wspace=0.18)
    axL = fig.add_subplot(gs[0, 0])
    axR = fig.add_subplot(gs[0, 1], sharey=axL)
    axT = fig.add_subplot(gs[1, :])

    # ---- left: our recovered section ----
    vmax = 4.0
    im = axL.imshow(section, extent=[lat[0], lat[-1], dep[-1], dep[0]],
                    aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    interpolation="bilinear")
    band = cat[(cat.longitude - a.profile_lon).abs() <= a.halfwidth]
    axL.scatter(band.latitude, band.depth, s=9, facecolor="none",
               edgecolor="k", lw=0.4, alpha=0.75, zorder=5)
    axL.set_ylim(a.zmax, 0)
    axL.set_xlabel("Latitude (°N)")
    axL.set_ylabel("Depth (km)")
    axL.set_title(f"This study — {a.profile_lon:g}°E section\n"
                  f"{len(band)} events, land + OBS rays ({KOULAKOV['events']*0+1005} total catalogue)",
                  fontsize=11, fontweight="bold")
    cbar = fig.colorbar(im, ax=axL, orientation="vertical", pad=0.02, shrink=0.9)
    cbar.set_label("Vp anomaly (%)")

    # ---- right: Koulakov (2007) schematic ----
    axR.set_facecolor("#f7f7f7")
    x0, x1 = lat[0], lat[-1]
    trench_x = x0 + 0.15 * (x1 - x0)
    xs = np.linspace(trench_x, x1, 200)
    slab_top = (xs - trench_x) * 111.19 * np.tan(np.radians(KOULAKOV["dip_deg"]))
    slab_bot = slab_top + 25  # ~25 km slab thickness, schematic
    axR.fill_between(xs, slab_top, slab_bot, color="#2166ac", alpha=0.55,
                     label=f"slab core (dip {KOULAKOV['dip_deg']:g}°, schematic)")
    wedge_bot = slab_top - KOULAKOV["wedge_thickness_km"] * 0.6
    axR.fill_between(xs, np.clip(wedge_bot, 0, None), slab_top, color="#b2182b", alpha=0.30,
                     label="mantle wedge (slow, schematic)")
    axR.axhline(KOULAKOV["max_depth_km"], color="0.4", ls=":", lw=1)
    axR.text(x1 - 0.05, KOULAKOV["max_depth_km"] - 6,
             f"resolved to ~{KOULAKOV['max_depth_km']:g} km (Koulakov 2007)",
             ha="right", fontsize=8, color="0.35")
    axR.set_ylim(a.zmax, 0)
    axR.set_xlim(x0, x1)
    axR.set_xlabel("Latitude (°N)")
    axR.set_title(f"Koulakov et al. (2007) — reported geometry\n"
                  f"schematic redrawn from dip/depth values, not their gridded model",
                  fontsize=11, fontweight="bold")
    axR.legend(loc="lower left", fontsize=8, framealpha=0.9)
    plt.setp(axR.get_yticklabels(), visible=False)

    # ---- bottom: metrics table ----
    axT.axis("off")
    rows = [
        ["", "Koulakov et al. 2007", "This study (land+OBS)"],
        ["Events", "292", "1,005"],
        ["P+S arrivals", "~13,000", "21,363"],
        ["Arrivals / event", "~45", "21.3"],
        ["Resolved cells (5 km grid)", "—", "39,869 (of ~267k)"],
        ["Pattern correlation @ 20 km, 15–20 km depth", "—", "0.78"],
        ["Pattern correlation @ 20 km, 55–60 km depth", "—", "0.56"],
        ["Reported slab dip", "~50°", "see left panel"],
        ["Max resolved depth", "~300 km", "~150–200 km (this grid)"],
    ]
    tbl = axT.table(cellText=rows, loc="center", cellLoc="left",
                    colWidths=[0.42, 0.29, 0.29])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("#dde6f0")
        cell.set_edgecolor("0.75")

    fig.suptitle("MERAMEX 2004 vs Koulakov et al. (2007): Slab Geometry & Data Volume",
                 fontsize=14, fontweight="bold", y=0.995)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=170, bbox_inches="tight")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
