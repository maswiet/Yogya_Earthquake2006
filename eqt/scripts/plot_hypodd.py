#!/usr/bin/env python3
"""Compare single-event (hypoDD.loc) vs double-difference (hypoDD.reloc)
locations: maps + across-strike cross-sections showing the fault sharpening.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HD = f"{ROOT}/hypodd"
COLS = ["id","lat","lon","depth","x","y","z","ex","ey","ez","yr","mo","dy",
        "hr","mi","sc","mag","nccp","nccs","nctp","ncts","rcc","rct","cid"]

def load(f):
    d = pd.read_csv(f, sep=r"\s+", header=None, names=COLS, engine="python")
    return d

def strike_proj(lat, lon, lat0, lon0, vec):
    kx = 111.0*np.cos(np.radians(lat0)); ky = 111.0
    x = (lon-lon0)*kx; y = (lat-lat0)*ky
    al = x*vec[0]+y*vec[1]; ac = x*(-vec[1])+y*vec[0]
    return al, ac

def main():
    loc = load(f"{HD}/hypoDD.loc"); rel = load(f"{HD}/hypoDD.reloc")
    lat0, lon0 = rel.lat.mean(), rel.lon.mean()
    kx = 111.0*np.cos(np.radians(lat0))
    x = (rel.lon-lon0)*kx; y = (rel.lat-lat0)*111.0
    C = np.cov(np.vstack([x, y])); w, V = np.linalg.eigh(C)
    vec = V[:, np.argmax(w)]; az = np.degrees(np.arctan2(vec[0], vec[1])) % 180
    print(f"hypoDD: {len(rel)} relocated of {len(loc)} | strike N{az:.0f}E")

    fig, ax = plt.subplots(2, 2, figsize=(14, 12))
    for col, d, tag in [(0, loc, "single-event (NLLoc/VELEST)"), (1, rel, "double-difference (HypoDD)")]:
        sc = ax[0, col].scatter(d.lon, d.lat, s=4, c=d.depth, cmap="turbo", vmin=0, vmax=18, alpha=0.5)
        ax[0, col].set_title(f"Epicentres — {tag} (n={len(d)})")
        ax[0, col].set_xlabel("Lon"); ax[0, col].set_ylabel("Lat"); ax[0, col].set_aspect("equal", "box")
        ax[0, col].set_xlim(110.30, 110.60); ax[0, col].set_ylim(-8.05, -7.82)
        al, ac = strike_proj(d.lat, d.lon, lat0, lon0, vec)
        ax[1, col].scatter(ac, d.depth, s=4, c=d.depth, cmap="turbo", vmin=0, vmax=18, alpha=0.5)
        ax[1, col].set_ylim(20, 0); ax[1, col].set_xlim(-15, 15)
        ax[1, col].set_xlabel("across-strike distance (km)"); ax[1, col].set_ylabel("depth (km)")
        ax[1, col].set_title(f"Across-strike section — {tag}")
    plt.colorbar(sc, ax=ax[:, 1], label="depth (km)", shrink=0.6)
    plt.savefig(f"{ROOT}/figures/hypodd_relocation.png", dpi=120, bbox_inches="tight")
    # spread stats (scatter about local trend as a sharpness proxy)
    print("wrote figures/hypodd_relocation.png")

if __name__ == "__main__":
    main()
