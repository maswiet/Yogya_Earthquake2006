"""Validate EQTransformer picks + magnitudes against the Anggraini (2013) manual
catalogue -- the pick-level ground truth a reviewer will demand for an ML-picker
paper.

Inputs (from the group, June 3-7 2006, same GFZ deployment):
  phase_300.dat  event headers (origin + hypo + ML) followed by per-station
                 travel times: "STA  tt  weight  phase"
  station.dat    manual station codes + lat/lon
Our side:
  nll/loc_v/*.hyp   NonLinLoc picks (absolute arrival times per TF station/phase)
  full/catalog_magnitude.csv   our origins + ML

Method:
  1. Map manual station codes -> our TF codes by June-3-7 coordinates.
  2. Match events by origin time (<=5 s).
  3. For each matched event + common (station, phase): compare travel times
     referenced to each catalogue's OWN origin (removes any origin/clock offset)
     -> pick residual = our_tt - manual_tt.
  4. Event-matched ML comparison (resolves the ~1-unit scale offset per event).

Outputs: figures/pick_validation.png + printed statistics.
"""
import os, re, glob, json
import numpy as np, pandas as pd
from obspy import UTCDateTime
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PHASE = os.environ.get("ANGGRAINI_PHASE", os.path.expanduser("~/Downloads/phase_300.dat"))
MSTA = os.environ.get("ANGGRAINI_STA", os.path.expanduser("~/Downloads/station.dat"))
JDAY_MID = 156          # representative julian day for the 3-7 Jun window
MATCH_S = 5.0           # event origin-time match tolerance


def station_map():
    """manual code -> TF code, by active June-3-7 position (handles movers)."""
    man = {}
    for line in open(MSTA):
        f = line.split()
        if len(f) >= 3:
            man[f[0]] = (float(f[1]), float(f[2]))
    raw = json.load(open(f"{ROOT}/config/stations_raw.json"))
    per = json.load(open(f"{ROOT}/config/stations_periods.json"))
    tfpos = {}
    for tf, v in raw.items():
        pos = (v["latitude"], v["longitude"])
        if tf in per:                       # use the site active in the window
            for s in per[tf]["sites"]:
                if s["day_start"] <= JDAY_MID <= s["day_end"]:
                    pos = (s["lat"], s["lon"]); break
        tfpos[tf] = pos
    def km(a, b):
        return np.hypot((a[0]-b[0])*111.2, (a[1]-b[1])*111.2*np.cos(np.radians(-7.9)))
    mp = {}
    for mc, p in man.items():
        best = min(tfpos, key=lambda t: km(p, tfpos[t]))
        if km(p, tfpos[best]) < 0.5:
            mp[mc] = best
    return mp


def read_manual():
    evs = []
    cur = None
    for line in open(PHASE):
        if line.startswith("#"):
            f = line[1:].split()
            try:
                ot = UTCDateTime(int(f[0]), int(f[1]), int(f[2]), int(f[3]),
                                 int(f[4]), float(f[5]))
            except Exception:
                cur = None; continue
            cur = dict(ot=ot, lat=float(f[6]), lon=float(f[7]), dep=float(f[8]),
                       ml=float(f[9]), ph=[])
            if -1 <= cur["ml"] <= 6:        # drop the ML 10.89 / bad rows
                evs.append(cur)
            else:
                cur = None
        elif line.strip() and cur is not None:
            f = line.split()
            cur["ph"].append((f[0], float(f[1]), f[3]))   # sta, tt, phase
    return evs


def read_ours():
    """Our events: origin + {(TF,phase): arrival} from NLLoc hyp files."""
    geo = re.compile(r"GEOGRAPHIC\s+OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)")
    evs = []
    for path in sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp")):
        ot = None; ph = {}
        for line in open(path):
            if line.startswith("GEOGRAPHIC"):
                m = geo.search(line)
                if m:
                    y, mo, d, h, mi = map(int, m.groups()[:5])
                    ot = UTCDateTime(y, mo, d, h, mi, 0) + float(m.group(6))
            elif " > " in line and "GAU" in line:
                lf = line.split(" > ")[0].split()
                sta, phase = lf[0], lf[4]
                if phase in ("P", "S"):
                    t = UTCDateTime(f"{lf[6][:4]}-{lf[6][4:6]}-{lf[6][6:8]}T"
                                    f"{lf[7][:2]}:{lf[7][2:4]}:00") + float(lf[8])
                    ph[(sta, phase)] = t
        if ot is not None and "2006-06-03" <= str(ot)[:10] <= "2006-06-07":
            evs.append((ot, ph))
    return evs


def main():
    mp = station_map()
    man = read_manual()
    ours = read_ours()
    print(f"station map ({len(mp)}): " + ", ".join(f"{k}->{v}" for k, v in sorted(mp.items())))
    print(f"manual events {len(man)}   our events (3-7 Jun) {len(ours)}")

    ot_arr = np.array([e[0].timestamp for e in ours])
    order = np.argsort(ot_arr); ot_arr = ot_arr[order]
    ours_s = [ours[i] for i in order]

    resid = []          # per matched pick: phase, station, raw dt, event index
    mlpair = []
    ourml = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    ourml["ts"] = (pd.to_datetime(ourml.time) - pd.Timestamp("1970-01-01")).dt.total_seconds()
    matched = 0
    for ei, e in enumerate(man):
        ts = e["ot"].timestamp
        i = np.searchsorted(ot_arr, ts)
        cand = [j for j in (i-1, i) if 0 <= j < len(ot_arr) and abs(ot_arr[j]-ts) <= MATCH_S]
        if not cand:
            continue
        j = min(cand, key=lambda k: abs(ot_arr[k]-ts))
        oot, oph = ours_s[j]
        matched += 1
        for sta, tt, phase in e["ph"]:
            tf = mp.get(sta)
            if tf and (tf, phase) in oph:
                resid.append((phase, sta, (oph[(tf, phase)] - oot) - tt, ei))
        k = np.argmin(np.abs(ourml.ts.values - oot.timestamp))
        if abs(ourml.ts.values[k]-oot.timestamp) < 1:
            mlpair.append((e["ml"], ourml.ML.values[k]))

    resid = pd.DataFrame(resid, columns=["phase", "sta", "dt", "ei"])
    # demean per event: removes the origin-time/location offset, leaving pick scatter
    resid["dt_dm"] = resid.dt - resid.groupby("ei").dt.transform("median")
    mlp = pd.DataFrame(mlpair, columns=["ml_man", "ml_us"])
    rP = resid.dt_dm[resid.phase == "P"]; rS = resid.dt_dm[resid.phase == "S"]
    print(f"\nmatched events: {matched}/{len(man)} ({100*matched/len(man):.0f}%)")
    print("pick residuals, per-event demeaned (pure pick scatter):")
    for lab, r in [("P", rP), ("S", rS)]:
        print(f"  {lab}: n={len(r)}  median {r.median():+.2f}s  MAD {(r-r.median()).abs().median():.2f}s"
              f"  |dt|<0.3s: {100*(r.abs()<0.3).mean():.0f}%")
    print("per-station median demeaned residual (flags any mis-mapped station):")
    st = resid.groupby("sta").dt_dm.agg(["median", "size"]).sort_values("median")
    print("  " + "  ".join(f"{s}:{r['median']:+.2f}({int(r['size'])})" for s, r in st.iterrows()))
    off = mlp.ml_man - mlp.ml_us
    print(f"\nevent-matched ML: n={len(mlp)}  offset (manual-ours) median {off.median():+.2f}"
          f"  mean {off.mean():+.2f}  std {off.std():.2f}  r={np.corrcoef(mlp.ml_man, mlp.ml_us)[0,1]:.2f}")

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for r, lab, col in [(rP, "P", "tab:blue"), (rS, "S", "tab:red")]:
        ax[0].hist(r.clip(-1, 1), bins=41, alpha=0.6, color=col,
                   label=f"{lab} (n={len(r)}, MAD {(r-r.median()).abs().median():.2f}s)")
    ax[0].axvline(0, color="0.3", lw=1)
    ax[0].set_xlabel("pick residual, per-event demeaned  (s)"); ax[0].set_ylabel("count")
    ax[0].set_title("A  Pick-time scatter vs manual"); ax[0].legend(fontsize=8)

    ax[1].scatter(mlp.ml_man, mlp.ml_us, s=12, alpha=0.5, color="tab:purple")
    lim = [-1, 4]
    ax[1].plot(lim, lim, "k--", lw=1, label="1:1")
    ax[1].plot(lim, [x-off.median() for x in lim], "r--", lw=1,
               label=f"offset {off.median():+.2f}")
    ax[1].set_xlim(lim); ax[1].set_ylim(-2, 3)
    ax[1].set_xlabel("Anggraini ML (manual)"); ax[1].set_ylabel("our ML")
    ax[1].set_title(f"B  Event-matched ML (n={len(mlp)})"); ax[1].legend(fontsize=8)

    ax[2].hist(off, bins=31, color="tab:green", alpha=0.7)
    ax[2].axvline(off.median(), color="tab:red", lw=1.5,
                  label=f"median {off.median():+.2f}")
    ax[2].set_xlabel("ML offset  manual − ours"); ax[2].set_ylabel("count")
    ax[2].set_title("C  ML offset distribution"); ax[2].legend(fontsize=8)

    fig.suptitle("EQTransformer vs Anggraini (2013) manual picks, 3-7 June 2006", fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = f"{ROOT}/figures/pick_validation.png"
    plt.savefig(out, dpi=140)
    print(f"wrote {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
