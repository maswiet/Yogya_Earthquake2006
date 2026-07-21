#!/usr/bin/env python3
"""Detection limit of the XN temporary array: is the small-magnitude tail real?

Four panels establishing, independently of the Gutenberg-Richter statistics,
what the physical ML detection floor of a *surface* 1 Hz network at 16-20 km
hypocentral distance is, and hence which part of the catalogue is usable:

  A  Per-reading ML vs hypocentral distance, with instrument noise-floor curves.
  B  Nearest-station distance vs event ML. A genuinely magnitude-limited
     catalogue detects its smallest events only where a station happens to be
     close, so this curve should decrease towards small ML. It does (12.7 km
     at ML 1 down to 9.0 km at ML -1), which is the expected behaviour.
     NOTE: before the simulate()-taper fix in build_amplitudes.py this curve
     was flat, because events near a day boundary had their amplitudes
     suppressed by up to ~200x and were scattered into the small-ML bins
     irrespective of geometry. The flatness was the artifact, not a result.
  C  Gutenberg-Richter, with Mc compared against the predicted noise floor.
     These do NOT coincide: the single-station ambient floor at the median
     nearest-station distance is ML -0.11, while network completeness (by
     b-value stability) is ML +0.50. The 0.6-unit gap is the cost of the
     network requirement -- an event must clear the noise at enough stations
     to be associated and located (>=8 phases, >=3 amplitudes), not at one.
     A single-station floor is therefore a lower bound on Mc, never an
     estimate of it.
  D  Why geothermal MEQ networks routinely reach ML -2 and this one does not:
     the same ML floor as a function of recording distance.

Noise floor: n counts RMS -> ground velocity n/S m/s (S = L4-3D sensitivity),
displacement v/(2*pi*f), Wood-Anderson magnification 2080, then the standard
Hutton & Boore (1987) distance terms.
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

SENS_L4 = 1.7e8    # counts/(m/s), L4-3D geophone + EDL digitiser
WA_GAIN = 2080.0   # Wood-Anderson static magnification (IASPEI)
F_CHAR  = 8.0      # Hz, characteristic frequency of the small-event S phase

def _mc(mags):
    """Completeness by b-value stability, as compute_magnitudes.py defines it."""
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from compute_magnitudes import mbs_mc
    return float(mbs_mc(mags))

# Noise levels to draw, in digitiser counts RMS.
NOISE_LEVELS = [(1, "1 count - digitiser floor", "0.35"),
                (5, "5 counts", "0.55"),
                (50, "50 counts - ambient (tropical, populated)", "0.15")]


def ml_from_amp(A_mm, R_km):
    """Hutton & Boore (1987) local magnitude."""
    return np.log10(A_mm) + 1.110*np.log10(R_km/100.) + 0.00189*(R_km-100.) + 3.0


def wa_amp_of_noise(counts, f=F_CHAR):
    """Wood-Anderson displacement (mm) produced by `counts` RMS of noise."""
    v = counts/SENS_L4                 # m/s
    d = v/(2*np.pi*f)                  # m
    return d*WA_GAIN*1000.             # mm


def ml_floor(R_km, counts, f=F_CHAR):
    return ml_from_amp(wa_amp_of_noise(counts, f), R_km)


def b_value(mags, mc, dM=0.1):
    """Aki-Utsu maximum likelihood, identical to compute_magnitudes.py."""
    m = mags[mags >= mc-1e-9]
    b = np.log10(np.e)/(m.mean()-(mc-dM/2))
    return b, b/np.sqrt(len(m)), len(m)


amp = pd.read_csv(f"{ROOT}/full/amplitudes.csv")
cat = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
MC = _mc(cat.ML.values)
amp = amp[amp.amp_mm > 0].copy()
amp["ML_pick"] = ml_from_amp(amp.amp_mm.values, amp.hypo_km.values)

fig, ax = plt.subplots(2, 2, figsize=(14, 10))

# ---- A: reading-level detections against the noise floor ---------------------
a = ax[0, 0]
hb = a.hexbin(amp.hypo_km, amp.ML_pick, gridsize=55, bins="log",
              cmap="Blues", mincnt=1, extent=(0, 60, -4, 4))
fig.colorbar(hb, ax=a, label="readings per cell (log)")
rr = np.linspace(1, 60, 200)
for counts, label, col in NOISE_LEVELS:
    a.plot(rr, ml_floor(rr, counts), color=col, lw=2,
           ls="--" if counts != 50 else "-", label=label)
a.axhline(MC, color="tab:red", lw=1.5, ls=":", label=f"Mc = {MC:g} (Gutenberg-Richter)")
a.set_xlim(0, 60); a.set_ylim(-4, 4)
a.set_xlabel("hypocentral distance (km)"); a.set_ylabel("per-reading ML")
a.set_title("A  Amplitude readings vs instrumental noise floor")
a.legend(fontsize=7.5, loc="lower right")

# ---- B: nearest-station distance vs event ML (the diagnostic) ---------------
b = ax[0, 1]
near = amp.groupby("evid").hypo_km.min().rename("r_min")
ev = cat.join(near, on="evid") if "evid" in cat.columns else cat.copy()
ev = ev.dropna(subset=["r_min"])
edges = np.arange(-2.0, 2.01, 0.25)
mid, med, q1, q3, n = [], [], [], [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    s = ev.r_min[(ev.ML >= lo) & (ev.ML < hi)]
    if len(s) >= 20:
        mid.append((lo+hi)/2); med.append(s.median())
        q1.append(s.quantile(0.25)); q3.append(s.quantile(0.75)); n.append(len(s))
mid, med = np.array(mid), np.array(med)
b.fill_between(mid, q1, q3, color="tab:blue", alpha=0.2, label="25-75th percentile")
b.plot(mid, med, "o-", color="tab:blue", lw=2, label="median nearest station")
trend = np.polyfit(mid, med, 1)
b.plot(mid, np.polyval(trend, mid), "--", color="0.4", lw=1.2,
       label=f"trend {trend[0]:+.1f} km per ML unit")
b.axvline(MC, color="tab:red", lw=1.5, ls=":", label=f"Mc = {MC:g}")
b.set_xlabel("event ML"); b.set_ylabel("distance to nearest recording station (km)")
b.set_title("B  Smaller events are detected closer to a station\n"
            "(expected for a magnitude-limited catalogue)")
b.set_ylim(0, max(q3)*1.15); b.legend(fontsize=8, loc="lower right")
for x, y, c in zip(mid[::3], med[::3], n[::3]):
    b.annotate(f"n={c}", (x, y), textcoords="offset points", xytext=(0, 7),
               ha="center", fontsize=6.5, color="0.35")

# ---- C: Gutenberg-Richter with the predicted floor --------------------------
c = ax[1, 0]
m = cat.ML.values
edges = np.arange(np.floor(m.min()*10)/10, m.max()+0.1, 0.1)
centres = edges[:-1]+0.05
cum = np.array([(m >= e).sum() for e in edges[:-1]])
c.semilogy(centres, cum, "o", ms=3.5, color="0.25", label="cumulative N(>=ML)")
bval, berr, nfit = b_value(m, MC)
fitm = (centres >= MC) & (cum > 0)
n_mc = (m >= MC).sum()
c.semilogy(centres[fitm], n_mc*10**(-bval*(centres[fitm]-MC)), "-",
           color="tab:red", lw=2,
           label=f"b = {bval:.2f} $\\pm$ {berr:.2f} (Aki-Utsu, ML $\\geq$ {MC:g})")
c.axvline(MC, color="tab:red", ls=":", lw=1.5, label=f"network Mc = {MC:+.2f}")
r_typ = float(near.median())
c.axvspan(ml_floor(r_typ, 5), ml_floor(r_typ, 50), color="tab:orange", alpha=0.18,
          label=f"1-station noise floor @ {r_typ:.0f} km")
c.annotate("", xy=(MC, 3e3), xytext=(ml_floor(r_typ, 50), 3e3),
           arrowprops=dict(arrowstyle="<->", color="0.3", lw=1.4))
c.text((MC+ml_floor(r_typ, 50))/2, 4.2e3,
       f"{MC-ml_floor(r_typ, 50):.2f} ML\nnetwork requirement",
       ha="center", fontsize=7.5, color="0.25")
c.set_xlabel("ML"); c.set_ylabel("cumulative number")
c.set_title("C  Network Mc sits 0.6 units ABOVE the 1-station noise floor")
c.legend(fontsize=8); c.set_xlim(-2.5, 4)

# ---- D: why geothermal networks reach ML -2 ---------------------------------
d = ax[1, 1]
rr = np.logspace(np.log10(0.5), np.log10(60), 200)
d.plot(rr, ml_floor(rr, 50), color="tab:blue", lw=2.5,
       label="surface network, ambient noise (this study)")
d.plot(rr, ml_floor(rr, 5), color="tab:green", lw=2.5,
       label="borehole siting (-20 dB ambient)")
r_lo, r_hi = float(near.quantile(0.25)), float(near.quantile(0.75))
d.axvspan(1, 3, color="tab:green", alpha=0.12)
d.axvspan(r_lo, r_hi, color="tab:blue", alpha=0.12)
d.annotate("geothermal MEQ\n(borehole, 1-3 km)", (1.8, ml_floor(1.8, 5)),
           textcoords="offset points", xytext=(6, -34), fontsize=8,
           color="tab:green", ha="left",
           arrowprops=dict(arrowstyle="->", color="tab:green", lw=1.2))
d.annotate(f"Yogya XN array\n(surface, {r_lo:.0f}-{r_hi:.0f} km)",
           (r_typ, ml_floor(r_typ, 50)),
           textcoords="offset points", xytext=(-14, 26), fontsize=8,
           color="tab:blue", ha="right",
           arrowprops=dict(arrowstyle="->", color="tab:blue", lw=1.2))
gap = ml_floor(r_typ, 50) - ml_floor(1.8, 5)
d.set_xscale("log")
d.set_xlabel("recording distance (km)"); d.set_ylabel("ML detection floor")
d.set_title(f"D  Distance + siting explain the {gap:.1f} magnitude-unit offset\n"
            "between geothermal MEQ practice and this deployment")
d.legend(fontsize=8, loc="upper left"); d.grid(alpha=0.3, which="both")

plt.tight_layout()
out = f"{ROOT}/figures/detection_noise_floor.png"
plt.savefig(out, dpi=140)

# ---- numbers for the manuscript --------------------------------------------
print(f"WA amp of 1/5/50 counts @ {F_CHAR:g} Hz (mm): "
      + ", ".join(f"{wa_amp_of_noise(c):.2e}" for c, _, _ in NOISE_LEVELS))
print(f"median nearest-station distance: {r_typ:.1f} km")
for counts, label, _ in NOISE_LEVELS:
    print(f"  ML floor @ {r_typ:.0f} km, {label:<44s} = {ml_floor(r_typ, counts):+.2f}")
print(f"empirical Mc = {MC:+.2f}   b = {bval:.2f} +/- {berr:.2f} (n={nfit})")
print(f"nearest-station distance median across ML bins: "
      f"{med.min():.1f}-{med.max():.1f} km (spread {med.max()-med.min():.1f} km)")
for thr in (-1.0, -0.5, 0.0):
    print(f"  events with ML >= {thr:+.1f}: {(cat.ML >= thr).sum():6d}")
print(f"geothermal-vs-here offset: {gap:.2f} magnitude units")
print(f"wrote {os.path.relpath(out, ROOT)}")
