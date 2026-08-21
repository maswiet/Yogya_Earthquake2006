#!/usr/bin/env python3
"""Turn config/stations_info.csv (parsed INFO.DAT) into stations_periods.json.

21 MERAMEX sites were re-occupied at a different location during the campaign,
so a station code alone does not define a geometry. Each deployment period gets
its own entry with the julian-day window it is valid for; the associator splits
multi-period codes into <CODE>a, <CODE>b, ... exactly like the 2006 pipeline.
"""
import csv, json, os, sys
from obspy import UTCDateTime

here = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(here, "..", "config", "stations_info.csv")
dst = os.path.join(here, "..", "config", "stations_periods.json")

periods = {}
for r in csv.DictReader(open(src)):
    t0, t1 = UTCDateTime(r["start"]), UTCDateTime(r["end"])
    site = dict(lat=float(r["lat"]), lon=float(r["lon"]), elev_m=int(r["elev_m"]),
                kind=r["kind"], sensor=r["sensor"], sensor_sn=r["sensor_sn"],
                recorder=r["recorder"], start=str(t0), end=str(t1),
                day_start=t0.julday, day_end=t1.julday)
    periods.setdefault(r["sta"], {"sites": []})["sites"].append(site)

for code, info in periods.items():
    info["sites"].sort(key=lambda s: s["day_start"])
    # close open-ended gaps: a period runs until the next one starts
    for a, b in zip(info["sites"], info["sites"][1:]):
        a["day_end"] = min(a["day_end"], b["day_start"])

json.dump(periods, open(dst, "w"), indent=1)
multi = {k: len(v["sites"]) for k, v in periods.items() if len(v["sites"]) > 1}
print(f"wrote {dst}: {len(periods)} station codes, "
      f"{sum(len(v['sites']) for v in periods.values())} periods")
print(f"multi-period sites ({len(multi)}):", multi)
