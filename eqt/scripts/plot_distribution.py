#!/usr/bin/env python3
"""Aftershock distribution figure from the NLLoc well-constrained catalog:
map (colour=depth) + along/across-strike cross-sections + depth histogram.
Strike is estimated by PCA of the epicentres.
"""
import os, json, argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", default=f"{ROOT}/full/catalog_nll_good.csv")
    ap.add_argument("--out", default=f"{ROOT}/figures/aftershock_distribution.png")
    ap.add_argument("--title", default="Yogyakarta 2006 aftershocks — EQTransformer + NonLinLoc")
    a = ap.parse_args()
    e = pd.read_csv(a.catalog)
    lat0, lon0 = e.latitude.mean(), e.longitude.mean()
    kx = 111.0*np.cos(np.radians(lat0)); ky = 111.0
    x = (e.longitude-lon0)*kx; y = (e.latitude-lat0)*ky        # km E, N
    # PCA -> strike direction
    C = np.cov(np.vstack([x, y])); w, V = np.linalg.eigh(C)
    strike_vec = V[:, np.argmax(w)]                            # along-strike unit
    across_vec = np.array([-strike_vec[1], strike_vec[0]])
    az = (np.degrees(np.arctan2(strike_vec[0], strike_vec[1]))) % 180
    al = x*strike_vec[0] + y*strike_vec[1]                     # along-strike km
    ac = x*across_vec[0] + y*across_vec[1]                     # across-strike km

    periods = json.load(open(f"{ROOT}/config/stations_periods.json"))
    slon=[]; slat=[]
    for info in periods.values():
        for s in info.get("sites", []): slon.append(s["lon"]); slat.append(s["lat"])

    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.4, 1])
    # map
    axm = fig.add_subplot(gs[0, 0])
    sc = axm.scatter(e.longitude, e.latitude, s=5, c=e.depth, cmap="turbo",
                     vmin=0, vmax=20, alpha=0.5)
    axm.scatter(slon, slat, marker="^", s=120, c="k", edgecolor="w", zorder=6)
    # draw strike line through centroid
    t = np.linspace(al.min(), al.max(), 2)
    axm.plot(lon0 + t*strike_vec[0]/kx, lat0 + t*strike_vec[1]/ky, "k--", lw=1)
    axm.annotate("A", (lon0+al.min()*strike_vec[0]/kx, lat0+al.min()*strike_vec[1]/ky))
    axm.annotate("A'", (lon0+al.max()*strike_vec[0]/kx, lat0+al.max()*strike_vec[1]/ky))
    plt.colorbar(sc, ax=axm, label="depth (km)")
    axm.set_xlabel("Longitude"); axm.set_ylabel("Latitude"); axm.set_aspect("equal","box")
    axm.set_title(f"Epicentres (n={len(e)}); strike ≈ N{az:.0f}°E")
    # depth histogram
    axh = fig.add_subplot(gs[0, 1])
    axh.hist(e.depth, bins=np.arange(0, 26, 1), color="steelblue", edgecolor="k",
             orientation="horizontal")
    axh.set_ylim(25, 0); axh.set_ylabel("depth (km)"); axh.set_xlabel("events")
    axh.set_title(f"Depth (median {e.depth.median():.1f} km)")
    # along-strike section A-A'
    axa = fig.add_subplot(gs[1, 0])
    axa.scatter(al, e.depth, s=4, c=e.depth, cmap="turbo", vmin=0, vmax=20, alpha=0.4)
    axa.set_ylim(25, 0); axa.set_xlabel("along-strike distance A–A' (km)")
    axa.set_ylabel("depth (km)"); axa.set_title("Along-strike cross-section")
    # across-strike section
    axc = fig.add_subplot(gs[1, 1])
    axc.scatter(ac, e.depth, s=4, c=e.depth, cmap="turbo", vmin=0, vmax=20, alpha=0.4)
    axc.set_ylim(25, 0); axc.set_xlabel("across-strike distance (km)")
    axc.set_ylabel("depth (km)"); axc.set_title("Across-strike cross-section (dip)")
    fig.suptitle(a.title, fontsize=14, y=0.99)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    plt.tight_layout(rect=[0,0,1,0.98]); plt.savefig(a.out, dpi=120)
    print(f"n={len(e)}  strike N{az:.0f}E  depth median {e.depth.median():.1f} km")
    print("wrote", a.out)

if __name__ == "__main__":
    main()
