#!/usr/bin/env python3
"""Associate MERAMEX EQTransformer picks into events with PyOcto (env `assoc`).

Geometry is time-dependent: 21 MERAMEX codes were re-occupied at a new site
mid-campaign, so each deployment period becomes its own PyOcto station
(<CODE>a, <CODE>b, ...) and every pick is relabelled to the site that was
actually recording on that julian day.

Velocity model defaults to the VELEST minimum-1D model derived from the 2006
Yogyakarta aftershocks (eqt/nll/nll_vel_v.in), extended to mantle depths so
Wadati-Benioff events are not forced into the crust.

Usage:
  associate.py --picks ../pilot/picks.csv --out ../pilot/events.csv \
      --tmpdir ../pilot/octo_tmp --n_picks 8 --n_p_and_s 3 --zmax 200
"""
import argparse, json, os, string, tempfile

import numpy as np
import pandas as pd
import pyocto
from obspy import UTCDateTime

# depth_km, Vp, Vs -- VELEST minimum-1D (Yogya 2006) crust + ak135 mantle
MODEL = [(-2.0, 2.50, 1.44), (0.0, 2.90, 1.67), (0.7, 4.30, 2.47),
         (2.0, 4.65, 2.56), (4.0, 5.49, 3.17), (7.0, 5.49, 3.34),
         (10.0, 6.30, 3.64), (13.0, 6.39, 3.73), (16.0, 6.55, 3.78),
         (22.0, 6.80, 3.92), (30.0, 7.20, 4.15),
         (40.0, 8.04, 4.48), (80.0, 8.05, 4.49), (120.0, 8.18, 4.51),
         (165.0, 8.30, 4.52), (210.0, 8.48, 4.61)]


def build_velocity_model(tmpdir, zmax, xmax, tolerance, cutoff):
    layers = [(max(d, 0.0), vp, vs) for d, vp, vs in MODEL if d <= zmax]
    df = pd.DataFrame({"depth": [d for d, _, _ in layers],
                       "vp": [p for _, p, _ in layers],
                       "vs": [s for _, _, s in layers]}).drop_duplicates("depth")
    path = os.path.join(tmpdir, "meramex_vmodel")
    pyocto.VelocityModel1D.create_model(df, 1.0, xmax, zmax, path)
    return pyocto.VelocityModel1D(path, tolerance=tolerance,
                                  association_cutoff_distance=cutoff)


def build_station_table(periods, only=None):
    """(stations_df, mapper) -- mapper(code, julday) -> period-specific id."""
    rows, ranges = [], {}
    for code, info in periods.items():
        if only is not None and code not in only:
            continue
        sites = info["sites"]
        suffixes = [""] if len(sites) == 1 else list(string.ascii_lowercase)
        for suf, s in zip(suffixes, sites):
            pid = f"{code}{suf}"
            rows.append({"id": pid, "latitude": s["lat"], "longitude": s["lon"],
                         "elevation": float(s["elev_m"])})
            ranges.setdefault(code, []).append((s["day_start"], s["day_end"], pid))

    def mapper(code, jday):
        cand = ranges.get(code)
        if not cand:
            return None
        for d0, d1, pid in cand:
            if d0 <= jday <= d1:
                return pid
        return min(cand, key=lambda c: min(abs(jday - c[0]), abs(jday - c[1])))[2]

    return pd.DataFrame(rows), mapper


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--picks", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--periods", default=os.path.join(here, "..", "config",
                                                      "stations_periods.json"))
    ap.add_argument("--tmpdir", default=None)
    ap.add_argument("--n_picks", type=int, default=8)
    ap.add_argument("--n_p_and_s", type=int, default=3)
    ap.add_argument("--pmin", type=float, default=0.0)
    ap.add_argument("--zmax", type=float, default=200.0)
    ap.add_argument("--tolerance", type=float, default=2.0)
    ap.add_argument("--cutoff", type=float, default=200.0)
    ap.add_argument("--pad", type=float, default=0.6, help="degrees around array")
    ap.add_argument("--pad-south", type=float, default=2.5,
                    help="extra degrees south of the array: the Java trench lies "
                         "~2.5 deg south of the southernmost station, and events "
                         "pile up on the association-box edge if it is too tight")
    ap.add_argument("--area", default=None,
                    help="explicit lat0,lat1,lon0,lon1 override")
    a = ap.parse_args()

    tmpdir = os.path.abspath(a.tmpdir) if a.tmpdir else tempfile.mkdtemp()
    os.makedirs(tmpdir, exist_ok=True)

    p = pd.read_csv(a.picks)
    if a.pmin > 0:
        p = p[p["probability"] >= a.pmin]
    codes = set(p["station"].unique())
    periods = json.load(open(a.periods))
    stations, mapper = build_station_table(periods, only=codes)

    t = p["peak_time"].apply(UTCDateTime)
    pid = [mapper(c, tt.julday) for c, tt in zip(p["station"], t)]
    picks = pd.DataFrame({"station": pid, "phase": p["phase"].str.upper(),
                          "time": [tt.timestamp for tt in t]}).dropna()
    picks = picks[picks["station"].isin(set(stations["id"]))]
    print(f"{len(picks)} picks on {picks['station'].nunique()}/{len(stations)} "
          f"station-sites ({len(codes)} codes)")

    lat, lon = stations["latitude"], stations["longitude"]
    if a.area:
        lat0, lat1, lon0, lon1 = (float(x) for x in a.area.split(","))
    else:
        lat0, lat1 = lat.min() - a.pad - a.pad_south, lat.max() + a.pad
        lon0, lon1 = lon.min() - a.pad, lon.max() + a.pad
    print(f"association area: lat {lat0:.2f}..{lat1:.2f}  lon {lon0:.2f}..{lon1:.2f}  "
          f"depth 0..{a.zmax:.0f} km")
    xmax = 1.5 * 111.19 * float(max(lat1 - lat0, lon1 - lon0))
    vmodel = build_velocity_model(tmpdir, a.zmax, xmax, a.tolerance, a.cutoff)
    assoc = pyocto.OctoAssociator.from_area(
        lat=(lat0, lat1), lon=(lon0, lon1),
        zlim=(0, a.zmax), time_before=300, velocity_model=vmodel,
        n_picks=a.n_picks, n_p_and_s_picks=a.n_p_and_s)
    assoc.transform_stations(stations)
    events, assignments = assoc.associate(picks, stations)
    print("associated events:", len(events))

    if len(events):
        assoc.transform_events(events)
        events["time_utc"] = events["time"].apply(lambda x: str(UTCDateTime(x)))
        npick = assignments.groupby("event_idx").size().rename("n_picks")
        nsta = assignments.groupby("event_idx")["station"].nunique().rename("n_stations")
        nP = assignments[assignments.phase == "P"].groupby("event_idx").size().rename("n_P")
        nS = assignments[assignments.phase == "S"].groupby("event_idx").size().rename("n_S")
        for s in (npick, nsta, nP, nS):
            events = events.merge(s, left_on="idx", right_index=True, how="left")
        events[["n_P", "n_S"]] = events[["n_P", "n_S"]].fillna(0).astype(int)

    out = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cols = [c for c in ["idx", "time_utc", "latitude", "longitude", "depth",
                        "n_picks", "n_stations", "n_P", "n_S"] if c in events.columns]
    (events[cols] if len(events) else events).to_csv(out, index=False)
    assignments.to_csv(out.replace(".csv", "_assignments.csv"), index=False)
    stations.to_csv(out.replace(".csv", "_stations.csv"), index=False)
    print("wrote", out, "| events:", len(events))


if __name__ == "__main__":
    main()
