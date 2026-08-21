#!/usr/bin/env python3
"""Write the NonLinLoc control files + phase file for a MERAMEX pick set.

Uses GRID2D travel-time grids (one 2-D grid per station per phase) because the
velocity model is 1-D and the array has tens of stations -- far cheaper than the
3-D grids used for the 12-station 2006 network.

Station codes are per-deployment-period (<CODE>a/b/...), matching associate.py.

Usage:
  gen_nll.py --events ../pilot/events.csv --assignments ../pilot/events_assignments.csv \
      --outdir ../nll --tag pilot --min-phases 6
"""
import argparse, json, os, string

import pandas as pd
from obspy import UTCDateTime

# depth_km, Vp, Vs.  Crust = VELEST minimum-1D from the 2006 Yogya aftershocks;
# below 40 km the ak135 mantle is appended so slab events are not squeezed into
# a crustal half-space (ak135 mantle Vp/Vs ~1.80, not the crustal 1.73).
MODEL = [(-3.0, 2.50, 1.44), (0.0, 2.90, 1.67), (0.7, 4.30, 2.47),
         (2.0, 4.65, 2.56), (4.0, 5.49, 3.17), (7.0, 5.49, 3.34),
         (10.0, 6.30, 3.64), (13.0, 6.39, 3.73), (16.0, 6.55, 3.78),
         (22.0, 6.80, 3.92), (30.0, 7.20, 4.15),
         (40.0, 8.04, 4.48), (80.0, 8.05, 4.49), (120.0, 8.18, 4.51),
         (165.0, 8.30, 4.52), (210.0, 8.48, 4.61)]
ERR = {"P": 0.10, "S": 0.20}


def gtsrce_lines(periods, keep):
    out, coords = [], {}
    for code, info in periods.items():
        sites = info["sites"]
        sufs = [""] if len(sites) == 1 else list(string.ascii_lowercase)
        for suf, s in zip(sufs, sites):
            pid = f"{code}{suf}"
            if keep is not None and pid not in keep:
                continue
            coords[pid] = s
            out.append(f"GTSRCE {pid:<6} LATLON {s['lat']:.5f} {s['lon']:.5f} "
                       f"0.0 {s['elev_m']/1000.0:.3f}")
    return out, coords


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--events", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--outdir", default=os.path.join(here, "..", "nll"))
    ap.add_argument("--tag", default="pilot")
    ap.add_argument("--periods", default=os.path.join(here, "..", "config",
                                                      "stations_periods.json"))
    ap.add_argument("--min-phases", type=int, default=6)
    ap.add_argument("--nxy", type=int, default=401,
                    help="location-grid horizontal nodes at 1 km (default +-200 km)")
    ap.add_argument("--nz", type=int, default=254,
                    help="grid nodes in depth at 1 km (default -3 to 250 km)")
    ap.add_argument("--ndist", type=int, default=601,
                    help="2-D travel-time grid nodes in epicentral distance (km)")
    a = ap.parse_args()

    out = os.path.abspath(a.outdir)
    for sub in ("model", "time", "obs", "loc"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    periods = json.load(open(a.periods))
    asg = pd.read_csv(a.assignments)
    ev = pd.read_csv(a.events)
    keep_ev = set(ev["idx"])
    asg = asg[asg["event_idx"].isin(keep_ev)]
    used = set(asg["station"].unique())

    lines, coords = gtsrce_lines(periods, used)
    lat0 = round(sum(s["lat"] for s in coords.values()) / len(coords), 3)
    lon0 = round(sum(s["lon"] for s in coords.values()) / len(coords), 3)
    half = (a.nxy - 1) // 2
    trans = f"TRANS SIMPLE {lat0} {lon0} 0.0"
    # NLL requires xOrig = yOrig = 0 for a 2-D (2 x nDist x nZ) grid; the y axis
    # is epicentral distance, so it must reach the farthest station-event pair.
    vggrid_2d = f"VGGRID 2 {a.ndist} {a.nz} 0.0 0.0 -3.0 1.0 1.0 1.0 SLOW_LEN"
    locgrid = f"LOCGRID {a.nxy} {a.nxy} {a.nz} -{half}.0 -{half}.0 -3.0 1.0 1.0 1.0 PROB_DENSITY SAVE"

    def layers(phase):
        return "\n".join(
            f"LAYER {d:6.1f}  {vp:.2f} 0.00  {vs:.2f} 0.00  2.70 0.00"
            for d, vp, vs in MODEL)

    tag = a.tag
    with open(f"{out}/nll_vel_{tag}.in", "w") as fh:
        fh.write(f"CONTROL 1 54321\n{trans}\nVGOUT ./model/{tag}\nVGTYPE P\n"
                 f"{vggrid_2d}\n# VELEST minimum-1D (Yogya 2006) extended to the mantle\n"
                 f"{layers('P')}\n")
    with open(f"{out}/nll_vel_{tag}_S.in", "w") as fh:
        fh.write(f"CONTROL 1 54321\n{trans}\nVGOUT ./model/{tag}\nVGTYPE S\n"
                 f"{vggrid_2d}\n{layers('S')}\n")
    for ph in ("P", "S"):
        with open(f"{out}/nll_time_{tag}_{ph}.in", "w") as fh:
            fh.write(f"CONTROL 1 54321\n{trans}\n"
                     f"GTFILES ./model/{tag} ./time/{tag} {ph}\n"
                     f"GTMODE GRID2D ANGLES_NO\n" + "\n".join(lines) +
                     "\n\nGT_PLFD 1.0e-3 0\n")
    with open(f"{out}/nll_loc_{tag}.in", "w") as fh:
        fh.write(f"""CONTROL 1 54321
{trans}
LOCSIG "MERAMEX 2004 - EQTransformer + PyOcto, VELEST 1-D"
LOCFILES ./obs/{tag}.obs NLLOC_OBS ./time/{tag} ./loc/{tag}
LOCHYPOUT SAVE_NLLOC_ALL NLL_FORMAT_VER_2
LOCSEARCH OCT 16 16 12 0.01 50000 5000 0 1
{locgrid}
LOCMETH EDT_OT_WT 9999.0 4 -1 -1 -1 -1 -1 1
LOCGAU 0.2 0.0
LOCGAU2 0.02 0.05 2.0
LOCPHASEID P P
LOCPHASEID S S
LOCQUAL2ERR 0.1 0.2 0.5 1.0 2.0
LOCANGLES ANGLES_NO 5
""")

    nev = 0
    with open(f"{out}/obs/{tag}.obs", "w") as fh:
        for eidx, g in asg.groupby("event_idx"):
            wrote = 0
            for _, r in g.iterrows():
                sta = r["station"]
                if sta not in coords:
                    continue
                ph = str(r["phase"]).upper()
                t = UTCDateTime(float(r["time"]))
                fh.write(f"{sta:<6} ?    ?    ? {ph:<6} ? "
                         f"{t.year:04d}{t.month:02d}{t.day:02d} {t.hour:02d}{t.minute:02d} "
                         f"{t.second + t.microsecond*1e-6:07.4f} GAU "
                         f"{ERR.get(ph, 0.15):9.2e} -1.00e+00 -1.00e+00 -1.00e+00\n")
                wrote += 1
            if wrote >= a.min_phases:
                fh.write("\n")
                nev += 1
    print(f"grid origin {lat0} {lon0} | {len(lines)} station-sites | "
          f"{nev} events -> {out}/obs/{tag}.obs")


if __name__ == "__main__":
    main()
