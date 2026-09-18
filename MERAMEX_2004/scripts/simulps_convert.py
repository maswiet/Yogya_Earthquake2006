#!/usr/bin/env python3
"""Convert MERAMEX 2004 catalogue/picks/stations into SIMULPS (simul2023) input files.

Coordinate convention (verified against simul2023.f subroutines input2/input4 and
the real 'spi' example from the Zenodo package): origin at the SE (lower-right)
corner of the study box, X increases WEST, Y increases NORTH, longitude is
stored as NEGATIVE for East (classic Thurber/SIMULPS internal convention), and
station/event records use explicit N/S, E/W letters with unsigned degree/minute
magnitudes, which the program negates internally to match.

Writes, into --out:
  STNS  - station file (file02), format(2x,a4,i2,a1,f5.2,i4,a1,f5.2,i5,2f5.2,2i3)
  EQKS  - event/phase file (file04), kttfor=2 "phs" one-pick-per-line format
  MOD   - velocity model (file03), free-format grid + Vp then Vp/Vs blocks
  CNTL  - control file (file01), adapted from the Zenodo 'spi' example

Usage:
  simulps_convert.py --out ../simulps/run1 --nx 7 --ny 7 --nz 8
"""
import argparse, json, os, string
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))


def build_station_table(periods, only=None):
    """Copied verbatim from associate.py's build_station_table (avoids importing
    that module's pyocto/obspy dependency). Same a/b/c deployment-period
    disambiguation used when picks were associated into events, so that
    codes like 'AI4a' in events_combined_assignments.csv resolve to the
    right coordinates."""
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

# study-box origin (SE corner): southernmost lat, easternmost lon
LAT0 = -10.0
LON0 = 112.2  # deg E magnitude; stored internally as negative


def dm(val):
    """decimal degrees -> (unsigned integer degrees, unsigned minutes)."""
    a = abs(val)
    d = int(a)
    m = (a - d) * 60.0
    return d, m


def load_velest_layers(path):
    layers = []  # (top_km, vp, vs)
    with open(path) as f:
        for line in f:
            if line.startswith("LAYER"):
                p = line.split()
                top, vp, vs = float(p[1]), float(p[2]), float(p[4])
                layers.append((top, vp, vs))
    return layers


def write_stations(out_dir, periods_json, obs_csv, used_stations):
    periods = json.load(open(periods_json))
    land, _mapper = build_station_table(periods)
    obs = pd.read_csv(obs_csv)
    rows = []
    for r in land.itertuples():
        if r.id in used_stations:
            rows.append((r.id, r.latitude, r.longitude, float(r.elevation)))
    for r in obs.itertuples():
        if r.id in used_stations:
            rows.append((r.id, r.latitude, r.longitude, float(r.elev_m)))
    rows.sort(key=lambda x: x[0])

    missing = used_stations - {r[0] for r in rows}
    if missing:
        print(f"WARNING: {len(missing)} station codes used in assignments have no "
              f"coordinates: {sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}")

    lines = []
    ltdo, oltm = dm(LAT0)
    lndo, olnm = dm(LON0)
    # nzco=0/3 (Transverse Mercator) hit a real bug in simul2023's utm_geo:
    # the southern-hemisphere false-northing check is "if(north.eq.1)"
    # instead of "if(inorth.eq.1)", so the local variable `north` is never
    # set to 0 or 10,000,000 for any southern-hemisphere point (all of
    # Java) and utm_geo adds an uninitialised garbage value to y. nzco=5
    # (Thurber's original short-distance conversion) never calls
    # tmllxy/utm_geo at all -- correct and simpler for a ~400 km study box.
    lines.append(f"{-ltdo:4d}{oltm:6.2f} {-lndo:4d}{olnm:6.2f}{0.0:8.2f}{5:3d}{0.0:7.1f}"
                 "  origin,rota,nzco,cmerid")
    lines.append(f"{len(rows):5d}                      nsta")
    for sta, lat, lon, elev in rows:
        latd, latm = dm(lat)
        lond, lonm = dm(lon)
        name4 = sta[:4].ljust(4)
        line = (f"  {name4}{latd:2d}S{latm:5.2f}{lond:4d}E{lonm:5.2f}"
                f"{int(round(elev)):5d}{0.0:5.2f}{0.0:5.2f}{0:3d}{0:3d}")
        lines.append(line)

    with open(os.path.join(out_dir, "STNS"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"STNS: {len(rows)} stations")
    return {sta: (lat, lon) for sta, lat, lon, elev in rows}


def to_ddmm(deg):
    """signed decimal degrees -> (unsigned int degrees, unsigned minutes, hemisphere sign)"""
    d, m = dm(deg)
    return d, m


def write_events(out_dir, cat_csv, assign_csv, station_names, evid_prefix="ev",
                  max_errh_km=15.0, max_gap=300.0, max_rms=0.6, max_depth_km=180.0):
    cat = pd.read_csv(cat_csv)
    asg = pd.read_csv(assign_csv)
    epoch = pd.Timestamp("1970-01-01", tz="UTC")

    good = (cat.quality_pass | (cat.errh_km <= max_errh_km)) & \
           (cat.gap <= max_gap) & (cat.rms <= max_rms) & (cat.depth <= max_depth_km)
    n_all = len(cat)
    cat = cat[good]
    print(f"quality filter: kept {len(cat)} of {n_all} events "
          f"(errh<={max_errh_km}km or quality_pass, gap<={max_gap}, rms<={max_rms}, "
          f"depth<={max_depth_km}km)")

    cat_t = pd.to_datetime(cat.time_utc, utc=True)

    n_events = 0
    n_picks = 0
    lines = []

    for idx, ev in cat.iterrows():
        picks = asg[asg.event_idx == idx]
        picks = picks[picks.station.isin(station_names)]
        if len(picks) < 4:
            continue

        ot = cat_t.loc[idx]
        ot_sec = (ot - epoch).total_seconds()
        p_times = {p.station: p.time - ot_sec for p in picks.itertuples() if p.phase == "P"}

        pick_lines = []
        for p in picks.itertuples():
            phase = p.phase
            if phase == "P":
                val = p.time - ot_sec
                if val <= 0 or val > 120:
                    continue
            elif phase == "S":
                tp = p_times.get(p.station)
                if tp is None:
                    continue
                val = (p.time - ot_sec) - tp  # SIMULPS wants S as S-P differential
                if val <= 0 or val > 60:
                    continue
            else:
                continue
            sta4 = p.station[:4].ljust(4)
            rmk = f" {phase} 0"
            pick_lines.append(f"{sta4} {rmk} {val:9.4f}".ljust(110))

        if len(pick_lines) < 4:
            continue

        ltd, ltm = dm(ev.latitude)
        lnd, lnm = dm(ev.longitude)
        iyrmo = f"{ot.year % 100:02d}{ot.month:02d}"
        iday = f"{ot.day:02d}"
        ihr = f"{ot.hour:02d}"
        mino = ot.minute
        seco = ot.second + ot.microsecond / 1e6
        dep = max(ev.depth, 0.1)
        mag = 0.0
        evid = f"{evid_prefix}{idx:06d}"[:11].rjust(11)

        header = (f"{iyrmo}{iday} {ihr}{mino:02d} {seco:5.2f}"
                  f"{ltd:3d}S{ltm:5.2f} {lnd:3d}E{lnm:5.2f}"
                  f"{dep:7.2f}{mag:7.2f}")
        header = header.ljust(67) + evid
        header = header.ljust(110)

        lines.append(header)
        lines.extend(pick_lines)
        lines.append(" " * 110)  # blank line terminates event's pick list
        n_events += 1
        n_picks += len(pick_lines)

    with open(os.path.join(out_dir, "EQKS"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"EQKS: {n_events} events, {n_picks} picks written")
    return n_events


def write_velocity(out_dir, vel_file, nx, ny, nz, xn, yn, zn):
    layers = load_velest_layers(vel_file)
    tops = np.array([l[0] for l in layers])
    vps = np.array([l[1] for l in layers])
    vss = np.array([l[2] for l in layers])

    def vp_at(z):
        return float(np.interp(z, tops, vps))

    def vpvs_at(z):
        vp = np.interp(z, tops, vps)
        vs = np.interp(z, tops, vss)
        return float(vp / vs)

    lines = []
    lines.append(f"c MERAMEX 2004 initial model from VELEST 1-D (nll/nll_vel_obs.in)")
    lines.append(f" 1.0{nx:3d}{ny:3d}{nz:3d}")
    lines.append(" ".join(f"{v:6.1f}" for v in xn))
    lines.append(" ".join(f"{v:6.1f}" for v in yn))
    lines.append(" ".join(f"{v:6.1f}" for v in zn))
    lines.append("  0  0  0")  # no fixed nodes
    lines.append("  0  0  0")  # no linked nodes

    for k in range(nz):
        vp = vp_at(zn[k])
        row = " ".join(f"{vp:6.2f}" for _ in range(nx))
        for _j in range(ny):
            lines.append(row)
    for k in range(nz):
        vpvs = vpvs_at(zn[k])
        row = " ".join(f"{vpvs:6.2f}" for _ in range(nx))
        for _j in range(ny):
            lines.append(row)

    with open(os.path.join(out_dir, "MOD"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"MOD: {nx}x{ny}x{nz} grid, Vp {vp_at(zn[0]):.2f}-{vp_at(zn[-1]):.2f} km/s")


def write_control(out_dir, neqs, nitmax=4, ndip=9, iskip=3, nitloc=6, dxmax=3.0, eigtol=0.01):
    lines = [
        f"{neqs:5d}{0:4d}{0:4d}{1.0:6.1f}{4:3d}{3:3d}{0:3d}"
        "   neqs,nsht,nbls,wtsht,kout,kout2,kout3",
        f"{0:4d}{1.0:6.1f}{4:4d}{2.0:6.1f}{6.0:6.1f}{10.0:6.1f}"
        "   ntel,wttel,ntobmin,reste1,reste2,reste3",
        f"{0:3d}{50.0:6.1f}{1.5:5.1f}{1:3d}{2.0:5.1f}{0.80:6.2f}"
        "   iuserdt,edisrdt,facrayr,mxiqrdt,wtrdtsht,res4",
        f"{0:3d}{0:5d}{20.0:5.1f}{1.5:5.1f}{1:3d}{8:4d}{1.0:5.1f}{0.40:5.2f}"
        f"{75.0:6.1f}{25:4d}"
        "   iuseclu,nclu,edisedt,facraye,mxiqedt,nstaepair,wtepdt,wrmscemx,gminedt0,minedt",
        f"{0:4d}{5.0:5.1f}{40.0:6.1f}{200.0:6.1f}{8.0:6.1f}{0:4d}"
        "   ngsw,dhsw,boxmin,zswmax,wtgsw,itopo",
        f"{nitloc:3d}{0.62:5.2f}{eigtol:6.3f}{0.03:5.2f}{-0.5:6.2f}{dxmax:5.1f}{0.05:5.2f}{0.0:5.2f}"
        "   nitloc,wtsp,eigtol,rmscut,zmin,dxmax,rderr,ercof",
        f"{5:3d}{0.30:5.2f}{0.04:5.2f}{9.5:5.2f}{1.5:5.2f}"
        "   nhitct,dvpmx,dvpvsmx,corrvpmax,corrvpmin",
        f"{0:3d}{500.0:7.1f}{300.0:7.1f}{25.0:6.1f}{9500.0:7.1f}{2.0:5.2f}{1:3d}"
        "   idmp,vdamp(1..4),stepl,ieft",
        f"{0:3d}{0:3d}{nitmax:3d}{2:3d}{0.001:7.4f}{0:3d}{0.08:5.2f}{0:3d}{0:3d}"
        "   ires,i3d,nitmax,kspr,snrmct,ihomo,rmstop,ifixl,ifixz",
        f"{100.0:7.2f}{300.0:7.2f}{0.60:5.2f}{2.0:5.2f}{5.0:5.2f}"
        "   delt1,delt2,res1,res2,res3",
        f"{ndip:3d}{iskip:3d}{1.5:5.2f}{2.0:5.2f}"
        "   ndip,iskip,scale1,scale2",
        f"{1.4:5.2f}{0.004:7.4f}{8:3d}{15:3d}"
        "   xfac,tlim,nitpb(1),nitpb(2)",
        f"{2:3d}{1:3d}{1:3d}{0:3d}{0:3d}"
        "   kttfor,iusep,iuses,invdel,iuse2t",
        f"{0:3d}{75.0:7.1f}{0.07:6.2f}{1800.0:7.1f}{12.0:6.1f}"
        "   iuseq,dqmax,qdamp,qrmax,qmin",
    ]
    with open(os.path.join(out_dir, "CNTL"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"CNTL: neqs={neqs}, nitmax={nitmax}, ndip={ndip}, iskip={iskip}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cat", default=os.path.join(ROOT, "full", "catalog_obs_combined.csv"))
    ap.add_argument("--assign", default=os.path.join(ROOT, "full", "events_combined_assignments.csv"))
    ap.add_argument("--periods", default=os.path.join(ROOT, "config", "stations_periods.json"))
    ap.add_argument("--obs-sta", default=os.path.join(ROOT, "full", "obs_stations_elev.csv"))
    ap.add_argument("--vel", default=os.path.join(ROOT, "nll", "nll_vel_obs.in"))
    ap.add_argument("--picks", default=os.path.join(ROOT, "full", "picks_combined.csv"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "simulps", "run1"))
    ap.add_argument("--nx", type=int, default=7)
    ap.add_argument("--ny", type=int, default=7)
    ap.add_argument("--nz", type=int, default=8)
    ap.add_argument("--quasi1d", action="store_true",
                    help="collapse to nx=ny=3 (1-D starting model check, per docs Hint 1)")
    ap.add_argument("--nitmax", type=int, default=4)
    ap.add_argument("--nitloc", type=int, default=6)
    ap.add_argument("--dxmax", type=float, default=3.0)
    ap.add_argument("--eigtol", type=float, default=0.01)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    asg_all = pd.read_csv(a.assign)
    used_stations = set(asg_all.station.unique())

    coords = write_stations(a.out, a.periods, a.obs_sta, used_stations)
    neqs = write_events(a.out, a.cat, a.assign, set(coords.keys()))

    if a.quasi1d:
        nx, ny, nz = 3, 3, a.nz
        xn = [-999.0, 0.0, 999.0]
        yn = [-999.0, 0.0, 999.0]
        zn = list(np.linspace(-5, 200, nz))
        ndip, iskip = 9, 4
    else:
        nx, ny, nz = a.nx, a.ny, a.nz
        xn = np.concatenate([[-60.0], np.linspace(0, 375, nx - 2), [450.0]])
        yn = np.concatenate([[-60.0], np.linspace(0, 445, ny - 2), [505.0]])
        # outer z edge must be well past 2x the deepest real event (some
        # MERAMEX events reach ~240 km): "no part of a ray must reach over
        # halfway to an outer node" (simul2023.f_notes.txt, Hint 3).
        zn = np.concatenate([[-5.0], np.array([0, 10, 25, 50, 100, 150, 200])[:nz-2], [600.0]])
        ndip, iskip = 9, 3

    write_velocity(a.out, a.vel, nx, ny, nz, xn, yn, zn)
    write_control(a.out, neqs, nitmax=a.nitmax, ndip=ndip, iskip=iskip,
                  nitloc=a.nitloc, dxmax=a.dxmax, eigtol=a.eigtol)

    print(f"\nWrote SIMULPS inputs to {a.out}")


if __name__ == "__main__":
    main()
