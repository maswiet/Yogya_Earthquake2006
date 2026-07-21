#!/usr/bin/env python3
"""Derive station magnitude corrections and compare them with the VELOCITY
(travel-time) station corrections from VELEST.

Per-station ML for the same event scatters by up to +/-0.5 magnitude units
(TF16 reads 0.50 low, TF17 0.26 high), which is site response, not measurement
noise: it is systematic per station and stable across the deployment. Solving
for it both tightens the catalogue magnitudes and provides an independent
probe of the same sediment/limestone contrast that the VELEST travel-time
corrections map across the Opak fault.

Model:  ML_ij = M_j + c_i + e_ij     (station i, event j)
solved by alternating medians (robust to the outliers a picked catalogue
always contains), with the datum fixed by forcing the amplitude-weighted mean
correction to zero so the overall magnitude scale is unchanged.

Note the sign relationship with VELEST is NOT assumed. Thick sediment slows P
arrivals (positive travel-time correction) and would normally also amplify
ground motion (negative ML correction, since the station reads high), so a
negative correlation is the naive expectation -- but TF16 is slow (+0.40 s)
AND reads 0.50 low, which that story does not explain. Panel C reports the
correlation and leaves the interpretation to the data.

Writes config/station_ml_corrections.json and figures/station_ml_corrections.png.
"""
import os, re, json, argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VELEST_OUT = f"{ROOT}/velest/yogya.OUT"


def ml_pick(A_mm, R_km):
    return np.log10(A_mm) + 1.110*np.log10(R_km/100.) + 0.00189*(R_km-100.) + 3.0


def solve(a, n_iter=30, min_obs=30, tol=1e-4):
    """Alternating-median solution for event magnitudes + station corrections."""
    counts = a.groupby("sta").size()
    keep = counts[counts >= min_obs].index
    a = a[a.sta.isin(keep)].copy()
    corr = pd.Series(0.0, index=sorted(a.sta.unique()))
    hist = []
    for it in range(n_iter):
        a["ml_c"] = a.ml - a.sta.map(corr)
        ev = a.groupby("evid").ml_c.median().rename("M")
        a["M"] = a.evid.map(ev)
        new = a.groupby("sta").apply(lambda g: (g.ml - g.M).median(),
                                     include_groups=False)
        new = new - np.average(new, weights=counts[new.index])   # fix the datum
        shift = float(np.max(np.abs(new - corr)))
        corr = new
        scatter = float((a.ml - a.M - a.sta.map(corr)).abs().median())
        hist.append(scatter)
        if shift < tol:
            break
    a["ml_c"] = a.ml - a.sta.map(corr)
    ev = a.groupby("evid").ml_c.median().rename("ML_corr")
    return corr, ev, a, hist, it+1


def velest_corrections():
    """P-wave station corrections from the VELEST output, if available.

    Same block and triplet layout that interpret_stacorr.py reads.
    """
    if not os.path.exists(VELEST_OUT):
        return {}
    lines = open(VELEST_OUT, errors="ignore").read().splitlines()
    idx = [i for i, l in enumerate(lines) if "Adjusted station corrections:" in l]
    if not idx:
        return {}
    out = {}
    for l in lines[max(idx)+2: max(idx)+8]:
        t = l.split(); k = 0
        while k+1 < len(t):
            try:
                out[t[k]] = float(t[k+1]); k += 3
            except (ValueError, IndexError):
                break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--amp", default=f"{ROOT}/full/amplitudes.csv")
    ap.add_argument("--out", default=f"{ROOT}/config/station_ml_corrections.json")
    ap.add_argument("--fig", default=f"{ROOT}/figures/station_ml_corrections.png")
    ap.add_argument("--min_obs", type=int, default=30)
    a_ = ap.parse_args()

    a = pd.read_csv(a_.amp)
    a = a[a.amp_mm > 0].copy()
    a["ml"] = ml_pick(a.amp_mm.values, a.hypo_km.values)

    # scatter before correction, for the before/after statement
    ev0 = a.groupby("evid").ml.median().rename("M0")
    a["M0"] = a.evid.map(ev0)
    before = float((a.ml - a.M0).abs().median())
    spread_before = a.groupby("sta").apply(lambda g: (g.ml-g.M0).median(),
                                           include_groups=False)

    corr, ev, a, hist, nit = solve(a, min_obs=a_.min_obs)
    after = hist[-1]

    json.dump({k: round(float(v), 3) for k, v in corr.items()},
              open(a_.out, "w"), indent=1, sort_keys=True)

    vel = velest_corrections()
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))

    s = corr.sort_values()
    ax[0].barh(range(len(s)), s.values,
               color=["tab:red" if v > 0 else "tab:blue" for v in s.values])
    ax[0].set_yticks(range(len(s))); ax[0].set_yticklabels(s.index, fontsize=8)
    ax[0].axvline(0, color="0.3", lw=1)
    ax[0].set_xlabel("station ML correction (magnitude units)")
    ax[0].set_title("A  Station magnitude corrections\n"
                    "(positive = station reads high, subtracted from ML)")

    ax[1].plot(range(1, len(hist)+1), hist, "o-", color="tab:purple")
    ax[1].axhline(before, color="0.4", ls="--",
                  label=f"before = {before:.3f}")
    ax[1].set_xlabel("iteration"); ax[1].set_ylabel("median |residual| (ML units)")
    ax[1].set_title(f"B  Scatter reduced {before:.3f} -> {after:.3f}\n"
                    f"({100*(1-after/before):.0f}% tighter, {nit} iterations)")
    ax[1].legend(fontsize=8)

    common = [s for s in corr.index if s in vel]
    if len(common) >= 4:
        x = np.array([vel[s] for s in common]); y = corr[common].values
        ax[2].scatter(x, y, s=45, color="tab:green", zorder=3)
        for s_, xi, yi in zip(common, x, y):
            ax[2].annotate(s_, (xi, yi), fontsize=7,
                           textcoords="offset points", xytext=(4, 4))
        r = np.corrcoef(x, y)[0, 1]
        k = np.polyfit(x, y, 1)
        xr = np.linspace(x.min(), x.max(), 10)
        ax[2].plot(xr, np.polyval(k, xr), "--", color="0.4",
                   label=f"r = {r:+.2f}")
        ax[2].axhline(0, color="0.8", lw=0.8); ax[2].axvline(0, color="0.8", lw=0.8)
        ax[2].set_xlabel("VELEST P travel-time correction (s)")
        ax[2].set_ylabel("station ML correction")
        ax[2].set_title("C  Amplitude vs travel-time site response\n"
                        "(do the two probes see the same sites?)")
        ax[2].legend(fontsize=8)
    else:
        ax[2].text(0.5, 0.5, "VELEST station corrections\nnot parsed",
                   ha="center", va="center", transform=ax[2].transAxes,
                   fontsize=10, color="0.5")
        ax[2].axis("off")

    plt.tight_layout(); plt.savefig(a_.fig, dpi=140)

    print(f"stations solved      : {len(corr)}  (>= {a_.min_obs} readings)")
    print(f"scatter before/after : {before:.3f} -> {after:.3f} ML "
          f"({100*(1-after/before):.0f}% tighter)")
    print(f"correction range     : {corr.min():+.2f} .. {corr.max():+.2f}")
    print("\nper station (raw median residual -> solved correction):")
    for st in corr.sort_values().index:
        print(f"  {st:8s} {spread_before.get(st, np.nan):+6.2f} -> {corr[st]:+6.2f}")
    if len(common) >= 4:
        print(f"\ncorrelation with VELEST P corrections: r = "
              f"{np.corrcoef([vel[s] for s in common], corr[common].values)[0,1]:+.2f} "
              f"(n={len(common)})")
    print(f"\nwrote {os.path.relpath(a_.out, ROOT)} and {os.path.relpath(a_.fig, ROOT)}")


if __name__ == "__main__":
    main()
