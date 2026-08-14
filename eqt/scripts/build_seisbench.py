"""Build a SeisBench-format waveform archive from the located catalogue -- the
ML-ready data product for the Seismica Data Report (train/benchmark ML pickers).

One trace per (event, station): a 3-component (Z,N,E) 60 s window at 100 Hz
around the P pick, labelled with P and S arrival samples, plus event metadata
(origin, hypocentre, ML with the +0.41 local-scale note, gap, quality flag) and
station metadata. Raw counts (L4-3D response NOT removed) -- the STEAD/EQT
convention; users normalise/restitute as needed.

Efficient pass: group picks by (station folder, julian day), read each day once,
cut + resample only the short per-event windows.

Outputs (default eqt/seisbench/):  metadata.csv  +  waveforms.hdf5
Test:  build_seisbench.py --limit 3    Full:  build_seisbench.py
"""
import os, re, glob, json, hashlib, argparse
import numpy as np
import obspy
from obspy.core import UTCDateTime
import seisbench.data as sbd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_CANDS = ["/Volumes/Untitled/DATA-GFZ-Gempa-JOgja-tahap-2",
          "/Volumes/Untitled 1/DATA-GFZ-Gempa-JOgja-tahap-2"]
BASE = os.environ.get("EDL_BASE") or next((p for p in _CANDS if os.path.isdir(p)), _CANDS[0])
COMP = {"pri0": "Z", "pri1": "N", "pri2": "E"}
year_re = re.compile(r"e\d{4}(\d{2})\d+\.pri1$")
HYP = sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))

SR = 100.0           # target sampling rate (Hz)
PRE, LEN = 15.0, 60.0   # window: P-PRE .. P-PRE+LEN ; P at sample PRE*SR
NPTS = int(LEN*SR)
ML_TIE = 0.41        # add to our ML to tie to the local (Anggraini) scale


def parse_events():
    """evid -> dict(origin, lat, lon, dep, ml, gap, nsta, picks{sta:{P,S,dist}})."""
    import pandas as pd
    mag = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv").set_index("evid")
    q = pd.read_csv(f"{ROOT}/full/catalog_quality.csv").set_index("evid")
    geo = re.compile(r"GEOGRAPHIC\s+OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                     r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
    qual = re.compile(r"Gap\s+([\d.]+)")
    evs = {}
    for evid, path in enumerate(HYP):
        ot = lat = lon = dep = None; gap = np.nan; picks = {}
        for line in open(path):
            if line.startswith("GEOGRAPHIC"):
                m = geo.search(line)
                if m:
                    y, mo, d, h, mi = map(int, m.groups()[:5])
                    ot = UTCDateTime(y, mo, d, h, mi, 0) + float(m.group(6))
                    lat, lon, dep = float(m.group(7)), float(m.group(8)), float(m.group(9))
            elif line.startswith("QUALITY"):
                mm = qual.search(line); gap = float(mm.group(1)) if mm else np.nan
            elif " > " in line and "GAU" in line:
                lf = line.split(" > ")[0].split(); rf = line.split(" > ")[1].split()
                sta, ph = lf[0], lf[4]
                if ph not in ("P", "S"):
                    continue
                t = (UTCDateTime(f"{lf[6][:4]}-{lf[6][4:6]}-{lf[6][6:8]}T"
                                 f"{lf[7][:2]}:{lf[7][2:4]}:00") + float(lf[8]))
                picks.setdefault(sta, {})[ph] = t
                picks[sta]["dist"] = float(rf[6])
        if ot is None or not picks:
            continue
        evs[evid] = dict(ot=ot, lat=lat, lon=lon, dep=dep, gap=gap,
                         ml=float(mag.ML.get(evid, np.nan)),
                         nsta=int(mag.n_sta.get(evid, len(picks))),
                         passed=bool(q["pass"].get(evid, False)), picks=picks)
    return evs


def station_meta():
    st = {}
    for line in open(f"{ROOT}/nll/stations_gtsrce.txt"):
        f = line.split()
        if len(f) >= 5 and f[0] == "GTSRCE":
            st[f[1]] = (float(f[3]), float(f[4]))
    return st


def read_day(folder, jday):
    """3-component day stream (native rate), channels renamed HHZ/HHN/HHE."""
    d = os.path.join(BASE, folder, f"{jday:03d}")
    if not os.path.isdir(d):
        return None
    st = obspy.Stream()
    for ext, ch in COMP.items():
        s = obspy.Stream()
        for p1 in sorted(glob.glob(os.path.join(d, "*.pri1"))):
            m = year_re.search(os.path.basename(p1))
            if not (m and m.group(1) == "06"):
                continue
            f = p1.replace(".pri1", "." + ext)
            if os.path.exists(f):
                try:
                    s += obspy.read(f, format="MSEED")
                except Exception:
                    pass
        for tr in s:
            tr.stats.channel = "HH"+ch
        if len(s):
            s.merge(method=1, fill_value=0); st += s
    return st if len(st) else None


def split_of(evid):
    h = int(hashlib.md5(str(evid).encode()).hexdigest(), 16) % 100
    return "train" if h < 80 else ("dev" if h < 90 else "test")


def folder_for(sta):
    return f"tf30{sta[2:4]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{ROOT}/seisbench")
    ap.add_argument("--limit", type=int, default=0, help="process only N station-days (test)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    evs = parse_events()
    smeta = station_meta()
    # group (folder, jday) -> list of (evid, sta)
    groups = {}
    for evid, e in evs.items():
        for sta in e["picks"]:
            groups.setdefault((folder_for(sta), e["ot"].julday), []).append((evid, sta))
    keys = sorted(groups)
    if a.limit:
        keys = keys[:a.limit]
    print(f"events {len(evs)}  station-days {len(keys)}  BASE {BASE}")

    ntr = 0
    with sbd.WaveformDataWriter(f"{a.out}/metadata.csv", f"{a.out}/waveforms.hdf5") as w:
        w.data_format = {"dimension_order": "CW", "component_order": "ZNE",
                         "measurement": "velocity", "unit": "counts",
                         "instrument_response": "not restituted",
                         "sampling_rate": SR}
        for gi, (folder, jday) in enumerate(keys, 1):
            day = read_day(folder, jday)
            if day is None:
                continue
            for evid, sta in groups[(folder, jday)]:
                e = evs[evid]; pk = e["picks"][sta]
                if "P" not in pk:
                    continue
                p = pk["P"]; s = pk.get("S")
                t0 = p - PRE
                seg = day.slice(t0, t0+LEN+2).copy()
                seg = obspy.Stream([tr for tr in seg if tr.stats.npts > 10])
                if not len(seg):
                    continue
                seg.detrend("demean")
                try:
                    seg.resample(SR)
                except Exception:
                    continue
                arr = np.zeros((3, NPTS), np.float32)
                for ci, ch in enumerate("ZNE"):
                    tr = seg.select(channel="HH"+ch)
                    if len(tr):
                        d = tr[0].slice(t0, t0+LEN).data[:NPTS]
                        arr[ci, :len(d)] = d
                if not np.any(arr):
                    continue
                net, code = "YK", sta
                slat, slon = smeta.get(sta, (np.nan, np.nan))
                meta = {
                    "trace_name": f"{evid:05d}_{sta}",
                    "source_id": evid,
                    "source_origin_time": str(e["ot"]),
                    "source_latitude_deg": e["lat"], "source_longitude_deg": e["lon"],
                    "source_depth_km": e["dep"],
                    "source_magnitude": e["ml"], "source_magnitude_type": "ML",
                    "source_magnitude_local_tie": round(e["ml"]+ML_TIE, 2),
                    "source_gap_deg": e["gap"], "source_num_stations": e["nsta"],
                    "source_quality_passed": e["passed"],
                    "station_network_code": net, "station_code": code,
                    "station_latitude_deg": slat, "station_longitude_deg": slon,
                    "trace_channel": "HH", "trace_component_order": "ZNE",
                    "trace_sampling_rate_hz": SR,
                    "trace_start_time": str(t0),
                    "trace_p_arrival_sample": int(round(PRE*SR)),
                    "trace_p_status": "automatic",
                    "trace_s_arrival_sample": (int(round((PRE+(s-p))*SR))
                                               if s is not None and (s-p) < LEN-PRE else None),
                    "trace_s_status": "automatic" if s is not None else None,
                    "path_ep_distance_km": pk.get("dist", np.nan),
                    "split": split_of(evid),
                }
                w.add_trace(meta, arr)
                ntr += 1
            if gi % 50 == 0:
                print(f"  [{gi}/{len(keys)}] traces so far: {ntr}", flush=True)
    print(f"done: {ntr} traces -> {os.path.relpath(a.out, ROOT)}/")


if __name__ == "__main__":
    main()
