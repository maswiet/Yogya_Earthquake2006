#!/usr/bin/env python3
"""Velocity structure beneath the future 2006 Opak-fault rupture zone.

Three panels:
  A. Map view at 15 km depth (upper-plate seismogenic depth for the 2006
     Mw 6.3), zoomed on the rupture box, with the 2004 catalogue and the
     2006 aftershock cloud (this repo's GrowClust relocation) overlaid.
  B. Vertical section through the rupture centroid (110.46°E), showing
     velocity structure and the striking absence of 2004 hypocentres in
     the box that failed two years later.
  C. Ray/station geometry at the rupture depth — this was the single
     best-instrumented patch of the whole network, so the empty box is
     a detection result, not a detection gap.

Usage:
  plot_opak_velocity.py --out ../figures/opak_velocity_detail.png
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
from plot_pilot import coast_segments

BOX = dict(lat=(-8.01, -7.87), lon=(110.35, 110.53))
CENTRE = (-7.940, 110.460)


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


def km_grid(sta, lat_lo, lat_hi, lon_lo, lon_hi, n=140, k=4):
    la = np.linspace(lat_lo, lat_hi, n)
    lo = np.linspace(lon_lo, lon_hi, n)
    LO, LA = np.meshgrid(lo, la)
    d = np.stack([np.hypot((LA - r.latitude) * 111.19,
                           (LO - r.longitude) * 111.19 * np.cos(np.radians(LA)))
                  for r in sta.itertuples()])
    return LO, LA, np.sort(d, axis=0)[k - 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tomo", default=os.path.join(ROOT, "tomo"))
    ap.add_argument("--vp", default=os.path.join(ROOT, "tomo_full", "vp.npy"))
    ap.add_argument("--cat", default=os.path.join(ROOT, "full", "catalog_obs_combined.csv"))
    ap.add_argument("--stations", default=os.path.join(ROOT, "full", "events_combined_stations.csv"))
    ap.add_argument("--aftershocks", default=os.path.join(ROOT, "..", "eqt", "full", "catalog_growclust.csv"))
    ap.add_argument("--coast", default=os.path.join(ROOT, "data", "coastline.xy"))
    ap.add_argument("--out", default=os.path.join(ROOT, "figures", "opak_velocity_detail.png"))
    a = ap.parse_args()

    g = load_grid(a.tomo)
    lat, lon, dep = cell_axes(g)
    vp = np.load(a.vp)
    u0 = 1.0 / np.interp(dep, [0, 200], [5.8, 8.0])
    vp_pct = np.stack([to_pct(vp[:, :, k], u0[k]) for k in range(len(dep))], axis=-1)

    cat = pd.read_csv(a.cat)
    cat = cat[cat.quality_pass | (cat.errh_km <= 15)]
    sta = pd.read_csv(a.stations).drop_duplicates("id") if os.path.exists(a.stations) else None
    aft = pd.read_csv(a.aftershocks) if os.path.exists(a.aftershocks) else pd.DataFrame(columns=["lat", "lon"])
    coast = coast_segments(a.coast) if os.path.exists(a.coast) else []

    lat_lo, lat_hi, lon_lo, lon_hi = -8.25, -7.60, 110.10, 110.85

    fig = plt.figure(figsize=(16, 5.8))
    axA = fig.add_subplot(1, 3, 1)
    axB = fig.add_subplot(1, 3, 2)
    axC = fig.add_subplot(1, 3, 3)

    # ---- A: map view at 15 km with velocity + seismicity ----
    zi = int(np.argmin(np.abs(dep - 15)))
    ii = (lat >= lat_lo - 0.3) & (lat <= lat_hi + 0.3)
    jj = (lon >= lon_lo - 0.3) & (lon <= lon_hi + 0.3)
    sub = vp_pct[np.ix_(ii, jj)][:, :, zi]
    vmax = 4.0
    im = axA.imshow(sub, extent=[lon[jj][0], lon[jj][-1], lat[ii][-1], lat[ii][0]],
                    aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                    interpolation="bilinear", zorder=1)
    for seg in coast:
        axA.plot(seg[:, 0], seg[:, 1], color="0.3", lw=0.7, zorder=2)
    if len(aft):
        axA.scatter(aft.lon, aft.lat, s=1.0, color="#333333", alpha=0.12, lw=0, zorder=3,
                   label=f"{len(aft):,} 2006 aftershocks")
    inwindow = cat[cat.latitude.between(lat_lo, lat_hi) & cat.longitude.between(lon_lo, lon_hi)]
    inbox = cat[cat.latitude.between(*BOX["lat"]) & cat.longitude.between(*BOX["lon"])]
    axA.scatter(inwindow.longitude, inwindow.latitude, s=45, marker="*", color="#7fc97f",
               edgecolor="k", lw=0.4, zorder=5, label=f"2004 events, wider window: {len(inwindow)}")
    axA.scatter(inbox.longitude, inbox.latitude, s=110, marker="*", color="gold",
               edgecolor="k", lw=0.7, zorder=6, label=f"2004 events IN rupture box: {len(inbox)}")
    axA.add_patch(plt.Rectangle((BOX["lon"][0], BOX["lat"][0]),
                                BOX["lon"][1] - BOX["lon"][0], BOX["lat"][1] - BOX["lat"][0],
                                fill=False, ec="k", lw=2, zorder=5))
    axA.set_xlim(lon_lo, lon_hi); axA.set_ylim(lat_lo, lat_hi)
    axA.set_xlabel("Longitude (°E)"); axA.set_ylabel("Latitude (°N)")
    axA.set_title("A. Vp anomaly at 15 km + seismicity", fontsize=11, fontweight="bold")
    axA.legend(loc="lower left", fontsize=7.5, framealpha=0.9)
    axA.set_aspect(1 / np.cos(np.radians(7.94)))
    plt.colorbar(im, ax=axA, label="Vp anomaly (%)", shrink=0.8)

    # ---- B: vertical section through rupture centroid ----
    j = int(np.argmin(np.abs(lon - CENTRE[1])))
    section = vp_pct[:, j, :].T
    im2 = axB.imshow(section, extent=[lat[0], lat[-1], dep[-1], dep[0]],
                     aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                     interpolation="bilinear")
    band = cat[(cat.longitude - CENTRE[1]).abs() <= 0.15]
    axB.scatter(band.latitude, band.depth, s=14, facecolor="none",
               edgecolor="k", lw=0.5, zorder=5)
    band_shallow = band[(band.latitude.between(*BOX["lat"])) & (band.depth <= 25)]
    axB.scatter(band_shallow.latitude, band_shallow.depth, s=60, marker="*",
               facecolor="gold", edgecolor="k", lw=0.7, zorder=6)
    axB.axvspan(BOX["lat"][0], BOX["lat"][1], color="k", alpha=0.08, zorder=1)
    axB.axvline(CENTRE[0], color="k", ls="--", lw=1)
    axB.set_ylim(60, 0)
    axB.set_xlim(lat_lo, lat_hi)
    axB.set_xlabel("Latitude (°N)"); axB.set_ylabel("Depth (km)")
    axB.set_title(f"B. Section at {CENTRE[1]:.2f}°E ±0.15°\n"
                  f"n={len(band)} nearby 2004 events, {len(inbox)} inside rupture box",
                  fontsize=11, fontweight="bold")
    plt.colorbar(im2, ax=axB, label="Vp anomaly (%)", shrink=0.8)

    # ---- C: detection capability at 2006 hypocentral depth ----
    if sta is not None:
        LO, LA, D4 = km_grid(sta, lat_lo, lat_hi, lon_lo, lon_hi)
        cf = axC.contourf(LO, LA, D4, levels=[0, 10, 15, 20, 30, 50], cmap="Greens_r", alpha=0.55)
        plt.colorbar(cf, ax=axC, label="Distance to 4th-nearest station (km)", shrink=0.8)
        axC.scatter(sta.longitude, sta.latitude, marker="^", s=70, facecolor="w",
                   edgecolor="k", lw=0.8, zorder=4, label="MERAMEX stations")
    axC.add_patch(plt.Rectangle((BOX["lon"][0], BOX["lat"][0]),
                                BOX["lon"][1] - BOX["lon"][0], BOX["lat"][1] - BOX["lat"][0],
                                fill=False, ec="#d62728", lw=2, zorder=5))
    axC.scatter(inbox.longitude, inbox.latitude, s=110, marker="*", color="gold",
               edgecolor="k", lw=0.7, zorder=6)
    axC.set_xlim(lon_lo, lon_hi); axC.set_ylim(lat_lo, lat_hi)
    axC.set_xlabel("Longitude (°E)"); axC.set_ylabel("Latitude (°N)")
    axC.set_title("C. Network geometry — best coverage\nin the whole 2004 deployment", fontsize=11, fontweight="bold")
    axC.legend(loc="lower left", fontsize=8, framealpha=0.9)
    axC.set_aspect(1 / np.cos(np.radians(7.94)))

    fig.suptitle("Near-Silent 2004 Seismicity at the Future Opak Fault Rupture\n"
                 "(2006 Mw 6.3 Yogyakarta earthquake source zone)",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=170, bbox_inches="tight")
    print("wrote", a.out)
    print(f"  2004 events inside rupture box (strict): {len(inbox)}")
    print(f"  2004 events in wider window: {len(inwindow)}")
    print(f"  2006 aftershocks in same box: {len(aft[aft.lat.between(*BOX['lat']) & aft.lon.between(*BOX['lon'])]) if len(aft) else 'N/A'}")


if __name__ == "__main__":
    main()
