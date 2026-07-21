#!/usr/bin/env python3
"""Aftershock event-rate vs time (Omori decay) + cumulative count."""
import os
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
e = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
e["t"] = pd.to_datetime(e["time"], utc=True)
e = e[e.gap < 180]
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compute_magnitudes import mbs_mc
mc = float(mbs_mc(e.ML.values))          # b-stability completeness
daily = e.set_index("t").resample("D").size()
daily_mc = e[e.ML >= mc].set_index("t").resample("D").size().reindex(daily.index, fill_value=0)

# Omori: n(t) = K/(t+c)^p  (t in days since mainshock 2006-05-27)
mainshock = pd.Timestamp("2006-05-27", tz="UTC")
tdays = np.array([(d - mainshock).days + 0.5 for d in daily.index], dtype=float)
y = daily.values.astype(float)
m = tdays > 0
def omori(t, K, c, p): return K/np.power(t+c, p)
try:
    popt,_ = curve_fit(omori, tdays[m], y[m], p0=[400,0.5,1.0],
                       bounds=([1,0.01,0.3],[1e5,20,2.5]), maxfev=20000)
    K,c,p = popt
    fit_ok = True
except Exception as ex:
    print("Omori fit failed:", ex); fit_ok = False

fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
ax[0].bar(daily.index, daily.values, width=1.0, color="lightgray", label="all events")
ax[0].bar(daily_mc.index, daily_mc.values, width=1.0, color="tab:red", alpha=0.8, label=f"ML ≥ {mc:g}")
if fit_ok:
    ax[0].plot(daily.index, omori(tdays, *popt), "b-", lw=2,
               label=f"Omori: p={p:.2f}, c={c:.2f} d")
ax[0].set_ylabel("events / day"); ax[0].set_title("Daily aftershock rate (Omori decay)")
ax[0].legend(); plt.setp(ax[0].get_xticklabels(), rotation=45, ha="right", fontsize=8)
ax[0].set_yscale("log"); ax[0].set_ylim(1, None)

axc = ax[1]
axc.plot(e.sort_values("t")["t"], np.arange(1, len(e)+1), color="navy", lw=2)
axc.set_ylabel("cumulative number of events"); axc.set_title(f"Cumulative aftershocks (n={len(e)})")
plt.setp(axc.get_xticklabels(), rotation=45, ha="right", fontsize=8)
plt.tight_layout(); plt.savefig(f"{ROOT}/figures/aftershock_rate_decay.png", dpi=120)
print(f"n={len(e)}  peak/day={daily.max()} on {daily.idxmax().date()}")
if fit_ok: print(f"Omori p={p:.2f}, c={c:.2f} d, K={K:.0f}")
print("wrote figures/aftershock_rate_decay.png")
