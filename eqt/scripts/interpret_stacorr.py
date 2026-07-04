#!/usr/bin/env python3
"""Interpret VELEST station corrections in terms of the Opak-fault geology:
western stations sit on Merapi sediment (slow, +delay), eastern on limestone
(fast, -delay). Plots corrections vs longitude and quantifies the W-E contrast.
"""
import os, re, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = f"{ROOT}/velest/yogya.OUT"
OPAK_LON = 110.43          # approx Opak-fault longitude in the study area

# our per-period code -> official XN name
TF2XN = {"TF12":"WON","TF13":"PEL","TF14":"RAT","TF16":"WAN","TF18":"BUM",
         "TF19":"PAL","TF15a":"NGL","TF10a":"PRA","TF10b":"IMO","TF17":"BOG",
         "TF11b":"KRI","TF07b":"KEM","TF09b":"KARA"}

def station_lons():
    per = json.load(open(f"{ROOT}/config/stations_periods.json"))
    import string
    lon = {}
    for code, info in per.items():
        sites = info.get("sites", [])
        if len(sites) <= 1:
            if sites: lon[code] = sites[0]["lon"]
        else:
            for suf, s in zip(string.ascii_lowercase, sites):
                lon[f"{code}{suf}"] = s["lon"]
    return lon

def main():
    lines = open(OUT).read().splitlines()
    i = max(i for i,l in enumerate(lines) if "Adjusted station corrections:" in l)
    cor = {}
    for l in lines[i+2:i+8]:
        t = l.split(); k = 0
        while k+1 < len(t):
            try: cor[t[k]] = float(t[k+1]); k += 3
            except (ValueError, IndexError): break
    lon = station_lons()
    rows = []
    for code, c in cor.items():
        if code not in lon: continue
        rows.append({"code": code, "xn": TF2XN.get(code, code), "lon": lon[code],
                     "stcor": c, "side": "West (sediment)" if lon[code] < OPAK_LON else "East (limestone)"})
    d = pd.DataFrame(rows).sort_values("lon")
    w = d[d.side.str.startswith("West")]; e = d[d.side.str.startswith("East")]
    print("Station P-corrections vs Opak fault:")
    for _, r in d.iterrows():
        print(f"  {r.xn:5} ({r.code:6}) lon {r.lon:.3f}  stcor {r.stcor:+.2f} s  [{r.side}]")
    print(f"\nWest (sediment) mean: {w.stcor.mean():+.2f} s ; East (limestone) mean: {e.stcor.mean():+.2f} s")
    print(f"W-E contrast: {w.stcor.mean()-e.stcor.mean():+.2f} s")

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.axvline(OPAK_LON, color="k", ls="--", lw=1.5, label="Opak fault (~110.43°E)")
    ax.scatter(w.lon, w.stcor, s=140, c="tab:blue", edgecolor="k", zorder=5, label="West: Merapi sediment")
    ax.scatter(e.lon, e.stcor, s=140, c="tab:red", edgecolor="k", zorder=5, label="East: limestone")
    for _, r in d.iterrows():
        ax.annotate(r.xn, (r.lon, r.stcor), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=9, fontweight="bold")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("P station correction (s)")
    ax.set_title("VELEST station corrections vs Opak fault\n"
                 "West (sediment) = slow / +delay ; East (limestone) = fast / −delay")
    ax.legend(); ax.grid(alpha=0.3)
    out = f"{ROOT}/figures/station_corrections_opak.png"
    plt.tight_layout(); plt.savefig(out, dpi=120)
    print("wrote", out)

if __name__ == "__main__":
    main()
