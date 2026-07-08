#!/usr/bin/env python3
"""Parse GrowClust output, save catalog, and plot single-event vs GrowClust
(maps + across-strike cross-sections) for the full catalog."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CAT = f"{ROOT}/growclust/OUT/out.growclust_cat"
cols = ["yr","mo","dy","hr","mi","sc","evid","lat","lon","dep","mag","qID","cID",
        "nbranch","npair","ndifP","ndifS","rmsP","rmsS","eh","ez","et",
        "latC","lonC","depC"]
d = pd.read_csv(CAT, sep=r"\s+", header=None, names=cols, engine="python")
# save relocated catalog
d[["evid","lat","lon","dep","mag","nbranch","cID","rmsP","rmsS"]].to_csv(
    f"{ROOT}/full/catalog_growclust.csv", index=False)
reloc = d[d.nbranch >= 2]     # events actually relocated (in multi-event clusters)
print(f"GrowClust: {len(d)} events, {len(reloc)} in multi-event clusters "
      f"({d.cID.nunique()} clusters); largest cluster {d.nbranch.max()} events")
print(f"depth median: catalog {d.depC.median():.1f} km -> GrowClust {d.dep.median():.1f} km")

# strike from relocated
lat0, lon0 = reloc.lat.mean(), reloc.lon.mean(); kx = 111*np.cos(np.radians(lat0))
X = (reloc.lon-lon0)*kx; Y = (reloc.lat-lat0)*111
w, V = np.linalg.eigh(np.cov(np.vstack([X, Y]))); vec = V[:, np.argmax(w)]
az = np.degrees(np.arctan2(vec[0], vec[1])) % 180
def across(la, lo): return ((lo-lon0)*kx)*(-vec[1]) + ((la-lat0)*111)*vec[0]

fig, ax = plt.subplots(2, 2, figsize=(14, 12))
for col, (la, lo, dp, tag) in enumerate([
        (reloc.latC, reloc.lonC, reloc.depC, "single-event (NLLoc/VELEST)"),
        (reloc.lat,  reloc.lon,  reloc.dep,  "GrowClust")]):
    sc = ax[0, col].scatter(lo, la, s=3, c=dp, cmap="turbo", vmin=0, vmax=18, alpha=0.45)
    ax[0, col].set_title(f"Epicentres — {tag} (n={len(reloc)})")
    ax[0, col].set_xlabel("Lon"); ax[0, col].set_ylabel("Lat"); ax[0, col].set_aspect("equal", "box")
    ax[0, col].set_xlim(110.30, 110.60); ax[0, col].set_ylim(-8.05, -7.82)
    ax[1, col].scatter(across(la, lo), dp, s=3, c=dp, cmap="turbo", vmin=0, vmax=18, alpha=0.45)
    ax[1, col].set_ylim(20, 0); ax[1, col].set_xlim(-15, 15)
    ax[1, col].set_xlabel("across-strike distance (km)"); ax[1, col].set_ylabel("depth (km)")
    ax[1, col].set_title(f"Across-strike section — {tag}")
plt.colorbar(sc, ax=ax[:, 1], label="depth (km)", shrink=0.6)
fig.suptitle(f"GrowClust relocation of the full catalog (strike ≈ N{az:.0f}°E)", fontsize=14)
plt.savefig(f"{ROOT}/figures/growclust_relocation.png", dpi=120, bbox_inches="tight")
print("wrote figures/growclust_relocation.png")
