#!/usr/bin/env python3
"""Read NonLinLoc summary .hyp files, hypocentre *and* phase records.

`parse_nll.py` only needs the hypocentres; the relocation codes also need the
per-station travel times, so that parsing lives here and both use it.
"""
import math, re

from obspy import UTCDateTime

GEO = re.compile(r"^GEOGRAPHIC\s+OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+"
                 r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
QUAL = re.compile(r"RMS\s+([\d.eE+-]+)\s+Nphs\s+(\d+)\s+Gap\s+([\d.]+)\s+Dist\s+([\d.]+)")
STAT = re.compile(r"CovXX\s+(-?[\d.eE+-]+).*?YY\s+(-?[\d.eE+-]+)\s+YZ\s+-?[\d.eE+-]+\s+"
                  r"ZZ\s+(-?[\d.eE+-]+)")


def read_hyp(path, require_phases=True):
    """-> list of dicts with ot/lat/lon/dep/rms/nphs/gap/dist/eh/ez and
    `phs` = [(station, phase, travel_time_seconds), ...].

    NLLoc's `*.sum.*.hyp` holds hypocentres only; the PHASE blocks live in the
    per-event files, so relocation input must be built from those (see
    `read_events`).
    """
    events, cur = [], None
    for line in open(path):
        if line.startswith("GEOGRAPHIC"):
            m = GEO.search(line)
            if not m:
                continue
            y, mo, d, h, mi = map(int, m.groups()[:5])
            cur = {"ot": UTCDateTime(y, mo, d, h, mi, 0) + float(m.group(6)),
                   "lat": float(m.group(7)), "lon": float(m.group(8)),
                   "dep": float(m.group(9)), "phs": []}
        elif cur is None:
            continue
        elif line.startswith("QUALITY"):
            m = QUAL.search(line)
            if m:
                cur.update(rms=float(m.group(1)), nphs=int(m.group(2)),
                           gap=float(m.group(3)), dist=float(m.group(4)))
        elif line.startswith("STATISTICS"):
            m = STAT.search(line)
            if m:
                cxx, cyy, czz = map(float, m.groups())
                cur["eh"] = round(math.sqrt(max(cxx, 0) + max(cyy, 0)), 3)
                cur["ez"] = round(math.sqrt(max(czz, 0)), 3)
        elif " > " in line and "GAU" in line:
            left, right = line.split(" > ")[0].split(), line.split(" > ")[1].split()
            try:
                tt = float(right[0]) + float(right[1])      # TTpred + residual
            except (IndexError, ValueError):
                continue
            if tt > 0:
                cur["phs"].append((left[0], left[4], tt))
        elif line.startswith("END_NLLOC"):
            if cur and (cur["phs"] or not require_phases):
                events.append(cur)
            cur = None
    if cur and (cur["phs"] or not require_phases):
        events.append(cur)
    return events


def read_events(locdir, tag, require_phases=True):
    """All per-event .hyp files for one run, sorted by origin time."""
    import glob, os
    pattern = os.path.join(locdir, f"{tag}.2*.grid0.loc.hyp")
    events = []
    for path in sorted(glob.glob(pattern)):
        events.extend(read_hyp(path, require_phases=require_phases))
    events.sort(key=lambda e: e["ot"])
    return events
