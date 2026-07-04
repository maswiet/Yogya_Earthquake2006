#!/usr/bin/env python3
"""Extract per-pick (distance, travel-time, phase) from NLLoc phase lines and
do the direct travel-time analysis:
  - Wadati diagram (ts-tp vs tp)  -> Vp/Vs
  - travel-time vs hypocentral distance (P & S) -> apparent velocities
Writes full/arrivals.csv and full/ttime_analysis.png
"""
import os, re, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HYPS = sorted(glob.glob(f"{ROOT}/nll/loc/yogya.2*.grid0.loc.hyp"))

rows = []
for evid, path in enumerate(HYPS):
    depth = None; gap = 9999.0; rms = 9.0
    for line in open(path):
        if line.startswith("GEOGRAPHIC"):
            m = re.search(r"Depth\s+(-?[\d.]+)", line); depth = float(m.group(1)) if m else None
        elif line.startswith("QUALITY"):
            mg = re.search(r"RMS\s+([\d.eE+-]+).*Gap\s+([\d.]+)", line)
            if mg: rms = float(mg.group(1)); gap = float(mg.group(2))
        elif " > " in line and "GAU" in line:
            left, right = line.split(" > ")
            lf = left.split(); rf = right.split()
            sta, phase = lf[0], lf[4]
            ttpred, res = float(rf[0]), float(rf[1])
            # after '>': TTpred Res Weight StaLoc(X Y Z) SDist SAzim ...
            sdist = float(rf[6])                 # epicentral distance (km)
            if depth is None:
                continue
            obs_tt = ttpred + res
            hypo = float(np.sqrt(sdist**2 + depth**2))
            rows.append((evid, sta, phase, obs_tt, sdist, hypo, depth, gap, rms))

df = pd.DataFrame(rows, columns=["evid","sta","phase","tt","epi","hypo","depth","gap","rms"])
df.to_csv(f"{ROOT}/full/arrivals.csv", index=False)
print(f"{len(df)} phase observations from {df.evid.nunique()} events")

# quality subset
q = df[(df.gap < 180) & (df.rms < 0.4) & (df.tt > 0)]
P = q[q.phase == "P"]; S = q[q.phase == "S"]
print(f"quality picks: P={len(P)} S={len(S)}")

# ---- Wadati: pair P & S per (evid, sta) ----
pv = P[["evid","sta","tt"]].rename(columns={"tt":"tp"})
sv = S[["evid","sta","tt"]].rename(columns={"tt":"ts"})
w = pv.merge(sv, on=["evid","sta"])
w["tsp"] = w["ts"] - w["tp"]
w = w[(w.tsp > 0) & (w.tp > 0) & (w.tsp < 30)]
# robust linear fit tsp = (VpVs-1)*tp  (through slope)
A = np.polyfit(w["tp"], w["tsp"], 1)
vpvs = A[0] + 1.0
print(f"Wadati: n={len(w)} pairs, slope={A[0]:.3f} -> Vp/Vs = {vpvs:.3f}")

# ---- tt vs hypocentral distance: apparent velocity via robust linear fit ----
def app_vel(d):
    # fit tt = a*dist + b ; apparent velocity = 1/a
    a, b = np.polyfit(d["hypo"], d["tt"], 1)
    return 1.0/a, b, a
vP, bP, aP = app_vel(P); vS, bS, aS = app_vel(S)
print(f"apparent Vp = {vP:.2f} km/s (intercept {bP:.2f}s) ; Vs = {vS:.2f} km/s ; Vp/Vs(app)={vP/vS:.2f}")

# ---- figure ----
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
# Wadati
ws = w.sample(min(4000, len(w)), random_state=1)
ax[0].scatter(ws.tp, ws.tsp, s=3, alpha=0.25, color="navy")
x = np.linspace(0, w.tp.quantile(.99), 50)
ax[0].plot(x, A[0]*x + A[1], "r-", lw=2, label=f"Vp/Vs = {vpvs:.2f}")
ax[0].set_xlabel("tP  (s)"); ax[0].set_ylabel("tS - tP  (s)")
ax[0].set_title(f"Wadati diagram (n={len(w)})"); ax[0].legend()
# tt-distance P
for a_, d_, v_, c_, nm in [(ax[1], P, vP, "tab:blue", "P"), (ax[2], S, vS, "tab:red", "S")]:
    ds = d_.sample(min(6000, len(d_)), random_state=1)
    a_.scatter(ds.hypo, ds.tt, s=3, alpha=0.2, color=c_)
    xx = np.linspace(0, d_.hypo.quantile(.99), 50)
    slope = 1.0/v_
    b = (d_.tt - slope*d_.hypo).median()
    a_.plot(xx, slope*xx + b, "k-", lw=2, label=f"apparent V{nm} = {v_:.2f} km/s")
    a_.set_xlabel("hypocentral distance (km)"); a_.set_ylabel(f"t{nm} (s)")
    a_.set_title(f"{nm} travel-time vs distance"); a_.legend()
plt.tight_layout(); plt.savefig(f"{ROOT}/full/ttime_analysis.png", dpi=110)
print("wrote full/arrivals.csv and full/ttime_analysis.png")
