#!/usr/bin/env python3
"""Shared I/O for the MERAMEX 2004 archive.

The archive is already miniSEED, but three recorder families label it
differently and only the *directory* name carries the real station code:

  EDL  <STA>/<DDD>/E<serial><yymmddhhmmss>.PRI0/1/2   100 Hz, STEIM1, little-endian
                                                      header station = recorder serial
                                                      channels p0/p1/p2 -> Z/N/E
  SAM  <STA>/<DDD>/<STA>.IN.<CHA>.2004.<DDD>.<hhmmss>  100 Hz, STEIM1, big-endian
                                                      CHA = SPZ/SPN/SPE (Trillium T3/3T)
                                                         or BHZ/BHN/BHE (Guralp 3ESP)
  OBS  <STA>/<DDD>/<STA>.IN.P0.D.2004.<DDD>.<hhmm>      50 Hz, STEIM2, big-endian

`build_stream` returns a 3-component ObsPy Stream with SEED-style channel codes
and the folder-derived station code, ready to hand to SeisBench.
"""
import csv, os

import obspy

NETWORK = "XM"          # placeholder network code for the campaign
EDL_MAP = {"PRI0": "Z", "PRI1": "N", "PRI2": "E"}
OBS_MAP = {"P0": "Z", "P1": "N", "P2": "E", "P3": "H"}

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "..", "config")


def load_index(path=None):
    """Station-day work list produced by index_data.py."""
    path = path or os.path.join(CONFIG, "station_days.csv")
    rows = []
    for r in csv.DictReader(open(path)):
        r["day"] = int(r["day"])
        r["nfiles"] = int(r["nfiles"])
        r["nbytes"] = int(r["nbytes"])
        r["path"] = os.path.join(r["root"], r["doydir"], r["station"], "%03d" % r["day"])
        rows.append(r)
    return rows


def discover(path):
    """[(filename, component)] for one station-day directory.

    Handles both naming schemes without assuming a fixed channel code: the SAM
    recorders wrote SPZ/SPN/SPE at the Kiel Trillium sites and BHZ/BHN/BHE at
    the GFZ Guralp sites, and the OBS wrote P0/P1/P2/P3.
    """
    try:
        names = os.listdir(path)
    except OSError:
        return []
    found = []
    for n in names:
        ext = n.rsplit(".", 1)[-1].upper()
        if ext in EDL_MAP:                       # EDL raw: *.PRI0/1/2
            found.append((n, EDL_MAP[ext]))
            continue
        parts = n.split(".")
        if len(parts) >= 5:                      # <STA>.<NET>.<CHA>....
            cha = parts[2].upper()
            comp = OBS_MAP.get(cha) or (cha[-1] if cha[-1:] in "ZNEH" else None)
            if comp:
                found.append((n, comp))
    return sorted(found)


def _band_code(sampling_rate, kind):
    """SEED band code: E for 100 Hz short-period, H for 100 Hz broadband,
    S/B for the 50 Hz OBS records."""
    if sampling_rate >= 80:
        return "H" if kind in ("SAM", "OBS") else "E"
    return "B" if kind in ("SAM", "OBS") else "S"


def build_stream(path, station, kind, merge_fill=0):
    """Read one station-day directory into a relabelled, merged Stream.

    Returns None when the directory holds no readable waveform data.
    """
    files = discover(path)
    if not files:
        return None

    st = obspy.Stream()
    for name, comp in files:
        try:
            seg = obspy.read(os.path.join(path, name), format="MSEED")
        except Exception:
            continue          # a truncated 30-min segment must not kill the day
        for tr in seg:
            tr.stats.network = NETWORK
            tr.stats.station = station
            tr.stats.location = ""
            tr.stats.channel = f"{_band_code(tr.stats.sampling_rate, kind)}H{comp}"
        st += seg
    if not len(st):
        return None

    st = st.select(component="[ZNE]")
    if not len(st):
        return None
    st.merge(method=1, fill_value=merge_fill)
    st.sort()
    return st
