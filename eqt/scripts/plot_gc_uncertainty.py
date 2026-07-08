#!/usr/bin/env python3
"""Parse GrowClust bootstrap output: relative-location uncertainties (eh/ez/et),
update catalog, and plot distributions + map coloured by horizontal error."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
cols = ["yr","mo","dy","hr","mi","sc","evid","lat","lon","dep","mag","qID","cID",
        "nbranch","npair","ndifP","ndifS","rmsP","rmsS","eh","ez","et",
        "latC","lonC","depC"]
d = pd.read_csv(f"{ROOT}/growclust/OUT/out.growclust_cat", sep=r"\s+", header=None,
                names=cols, engine="python")
d[["evid","lat","lon","dep","mag","nbranch","cID","rmsP","rmsS","eh","ez","et"]].to_csv(
    f"{ROOT}/full/catalog_growclust.csv", index=False)

# bootstrap errors are set for relocated events (nbranch>=2 & eh>=0)
b = d[(d.nbranch >= 2) & (d.eh >= 0) & (d.ez >= 0)].copy()
print(f"{len(d)} events, {len(b)} with bootstrap errors")
for k in ("eh","ez","et"):
    print(f"  {k}: median {b[k].median()*1000:.0f} m  (90th pct {b[k].quantile(.9)*1000:.0f} m)")

fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1, 3, 1)
ax1.hist(b.eh*1000, bins=np.arange(0, 400, 15), color="tab:blue", edgecolor="k")
ax1.axvline(b.eh.median()*1000, color="k", ls="--")
ax1.set_xlabel("horizontal error eh (m)"); ax1.set_ylabel("events")
ax1.set_title(f"Relative horizontal uncertainty (median {b.eh.median()*1000:.0f} m)")
ax2 = fig.add_subplot(1, 3, 2)
ax2.hist(b.ez*1000, bins=np.arange(0, 600, 20), color="tab:red", edgecolor="k")
ax2.axvline(b.ez.median()*1000, color="k", ls="--")
ax2.set_xlabel("vertical error ez (m)"); ax2.set_ylabel("events")
ax2.set_title(f"Relative vertical uncertainty (median {b.ez.median()*1000:.0f} m)")
ax3 = fig.add_subplot(1, 3, 3)
sc = ax3.scatter(b.lon, b.lat, s=6, c=b.eh*1000, cmap="viridis_r", vmin=0, vmax=200)
plt.colorbar(sc, ax=ax3, label="eh (m)"); ax3.set_aspect("equal", "box")
ax3.set_xlabel("Lon"); ax3.set_ylabel("Lat"); ax3.set_title("Epicentres coloured by horizontal error")
fig.suptitle(f"GrowClust bootstrap (nboot=100) relative-location uncertainties (n={len(b)})", fontsize=13)
plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig(f"{ROOT}/figures/growclust_uncertainty.png", dpi=120)
print("wrote figures/growclust_uncertainty.png")
