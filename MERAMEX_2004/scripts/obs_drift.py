#!/usr/bin/env python3
"""Measure and correct OBS clock drift from residuals vs catalogue events.

For each OBS pick, find the nearest catalogue event in time, predict the arrival
time at that OBS station using the catalogue location, and compute the residual.
Fit a linear drift per station and apply the correction.

Usage:
  obs_drift.py --catalog ../full/catalog_nll.csv --picks ../full/picks_obs.csv \
      --stations ../full/events_land_stations.csv \
      --out ../full/picks_obs_corrected.csv --drift ../full/obs_drift.txt
"""
import argparse, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from raytrace1d import load_layers, Model1D, Tracer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--picks", required=True)
    ap.add_argument("--stations", required=True)
    ap.add_argument("--vel", default=os.path.join(os.path.dirname(__file__), "..", "nll", "nll_vel_full.in"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--drift", required=True)
    ap.add_argument("--maxdist-to-event", type=float, default=5.0, help="max seconds from pick to event")
    ap.add_argument("--outlier-sigma", type=float, default=2.0)
    a = ap.parse_args()

    cat = pd.read_csv(a.catalog)
    picks = pd.read_csv(a.picks)
    stations_df = pd.read_csv(a.stations).drop_duplicates("id")
    stations = {r.id: (r.latitude, r.longitude) for r in stations_df.itertuples()}

    # filter to OBS picks only
    picks_obs = picks[picks.station.str.startswith(("OS", "OH"))].copy()
    print(f"processing {len(picks_obs)} OBS picks from {picks_obs.station.nunique()} stations")

    # set up ray tracer for P waves
    tracer = Tracer(Model1D(load_layers(a.vel, "P")))

    # convert times to numpy datetime64 for fast comparison
    cat["t_ns"] = pd.to_datetime(cat.time_utc, utc=True).values
    picks_obs["t_ns"] = pd.to_datetime(picks_obs.peak_time, utc=True).values
    cat_t_array = cat.t_ns.values
    cat_t_sec = (cat.t_ns - cat.t_ns[0]).dt.total_seconds().values

    residuals_by_sta = {}
    for sta_code in picks_obs.station.unique():
        picks_sta = picks_obs[picks_obs.station == sta_code]
        if sta_code not in stations:
            print(f"  {sta_code}: not in stations file, skipping")
            continue

        slat, slon = stations[sta_code]
        residuals = []; event_times = []

        for _, p in picks_sta.iterrows():
            # find nearest event
            dt_ns = np.abs((cat_t_array - p.t_ns) / np.timedelta64(1, "s"))
            if dt_ns.min() > a.maxdist_to_event:
                continue
            i = dt_ns.argmin()
            ev = cat.iloc[i]

            # ray parameters
            dlat = (slat - ev.latitude) * 111.19
            dlon = (slon - ev.longitude) * 111.19 * np.cos(np.radians(ev.latitude))
            dist = np.hypot(dlat, dlon)

            # predict P arrival
            got = tracer.trace(ev.depth, dist)
            if got is None:
                continue
            tt_pred, _, _ = got

            # residual: pick_time - (event_time + predicted_travel_time)
            pick_sec = (p.t_ns - cat.t_ns[0]) / np.timedelta64(1, "s")
            event_sec = cat_t_sec[i]
            res = pick_sec - event_sec - tt_pred
            residuals.append(res); event_times.append(event_sec)

        if len(residuals) < 5:
            print(f"  {sta_code}: only {len(residuals)} pairs, skipping")
            continue

        residuals = np.array(residuals); event_times = np.array(event_times)

        # linear fit
        A = np.column_stack([event_times, np.ones(len(event_times))])
        m, b = np.linalg.lstsq(A, residuals, rcond=None)[0]
        pred = A @ np.array([m, b])
        rmse = np.sqrt(np.mean((residuals - pred) ** 2))

        # outlier screening: refit without points >2σ away
        good = np.abs(residuals - pred) < a.outlier_sigma * rmse
        if good.sum() >= 5:
            A_good = A[good]; res_good = residuals[good]
            m, b = np.linalg.lstsq(A_good, res_good, rcond=None)[0]
            pred_good = A_good @ np.array([m, b])
            rmse = np.sqrt(np.mean((res_good - pred_good) ** 2))
            n_use = good.sum()
        else:
            n_use = len(residuals)

        residuals_by_sta[sta_code] = (m, b, rmse)
        # drift in seconds per day (event_times are in seconds from t0)
        t0 = cat.t_ns.iloc[0]
        t1 = cat.t_ns.iloc[-1]
        days_span = (t1 - t0) / np.timedelta64(1, "D")
        drift_per_day = m * days_span / (event_times[-1] - event_times[0]) if event_times[-1] > event_times[0] else 0
        print(f"  {sta_code}: {n_use}/{len(residuals)} pairs | "
              f"drift {m:.2e} s/s ({drift_per_day:.3f} s/day) offset {b:.3f} s ± {rmse:.3f} s")

    # write drift summary
    with open(a.drift, "w") as f:
        f.write("station  slope_s_per_sec  intercept_s  rmse_s\n")
        for sta in sorted(residuals_by_sta):
            m, b, rmse = residuals_by_sta[sta]
            f.write(f"{sta:8s} {m:15.6e} {b:12.6f} {rmse:8.3f}\n")
    print(f"\nwrote {a.drift}")

    # apply corrections: dt_corrected = dt - (m * t_sec + b)
    picks_out = picks.copy()
    # add numeric time for all picks
    t0_ns = pd.to_datetime(cat.time_utc.iloc[0], utc=True)
    picks_out["t_sec"] = (pd.to_datetime(picks_out.peak_time, utc=True) - t0_ns).dt.total_seconds()

    for sta, (m, b, _) in residuals_by_sta.items():
        mask = picks_out.station == sta
        picks_out.loc[mask, "t_sec"] -= (m * picks_out.loc[mask, "t_sec"] + b)

    # convert back to ISO format
    picks_out["peak_time"] = (t0_ns + pd.to_timedelta(picks_out.t_sec, unit="s")).dt.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    picks_out = picks_out.drop(columns=["t_sec"])
    picks_out.to_csv(a.out, index=False)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
