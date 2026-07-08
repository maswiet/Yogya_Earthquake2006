#!/usr/bin/env python3
"""Temporal aftershock analysis: map coloured by time + time series
(along-strike migration, depth vs time, daily rate & cumulative)."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAIN = pd.Timestamp("2006-05-27", tz="UTC")

e = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
e = e[e.gap < 180].copy()
e["t"] = pd.to_datetime(e["time"], utc=True)
e = e.sort_values("t").reset_index(drop=True)
e["days"] = (e["t"] - MAIN).dt.total_seconds()/86400.0

# along-strike distance (PCA of epicentres)
lat0, lon0 = e.latitude.mean(), e.longitude.mean()
kx = 111.0*np.cos(np.radians(lat0))
x = (e.longitude-lon0)*kx; y = (e.latitude-lat0)*111.0
w, V = np.linalg.eigh(np.cov(np.vstack([x, y]))); vec = V[:, np.argmax(w)]
az = np.degrees(np.arctan2(vec[0], vec[1])) % 180
e["al"] = x*vec[0] + y*vec[1]

fig = plt.figure(figsize=(16, 9))
gs = fig.add_gridspec(3, 2, width_ratios=[1.25, 1])
# --- temporal map (left, spans 3 rows) ---
axm = fig.add_subplot(gs[:, 0])
sc = axm.scatter(e.longitude, e.latitude, s=6, c=e.days, cmap="viridis", alpha=0.55)
axm.set_aspect("equal", "box"); axm.set_xlabel("Longitude"); axm.set_ylabel("Latitude")
axm.set_title(f"Aftershock epicentres coloured by time (strike ≈ N{az:.0f}°E)")
cb = plt.colorbar(sc, ax=axm, label="days since mainshock (27 May 2006)", shrink=0.8)
# --- along-strike migration ---
ax1 = fig.add_subplot(gs[0, 1])
ax1.scatter(e.t, e.al, s=4, c=e.days, cmap="viridis", alpha=0.4)
ax1.set_ylabel("along-strike (km)\nSW  →  NE"); ax1.set_title("Along-strike position vs time (migration)")
# --- depth vs time ---
ax2 = fig.add_subplot(gs[1, 1], sharex=ax1)
ax2.scatter(e.t, e.depth, s=4, c=e.depth, cmap="turbo_r", vmin=0, vmax=18, alpha=0.4)
ax2.set_ylim(20, 0); ax2.set_ylabel("depth (km)"); ax2.set_title("Depth vs time")
# --- daily rate + cumulative ---
ax3 = fig.add_subplot(gs[2, 1], sharex=ax1)
daily = e.set_index("t").resample("D").size()
ax3.bar(daily.index, daily.values, width=1.0, color="tab:red", alpha=0.7)
ax3.set_ylabel("events/day", color="tab:red"); ax3.set_yscale("log")
ax3b = ax3.twinx()
ax3b.plot(e.t, np.arange(1, len(e)+1), color="navy", lw=2)
ax3b.set_ylabel("cumulative", color="navy")
ax3.set_title("Daily rate (log) & cumulative")
for ax in (ax1, ax2, ax3):
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
fig.suptitle(f"Yogyakarta 2006 temporal aftershock evolution (n={len(e)})", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.98])
out = f"{ROOT}/figures/temporal_aftershocks.png"
plt.savefig(out, dpi=120); print(f"n={len(e)}  strike N{az:.0f}E; wrote {out}")
