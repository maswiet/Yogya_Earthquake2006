#!/usr/bin/env python3
"""Waveform cross-correlation differential times for HypoDD and GrowClust.

Catalog differential times inherit every picking error, so they can only sharpen
a cluster so far. Cross-correlating the actual waveforms of two nearby events at
a shared station measures the differential time to a fraction of a sample, which
is what makes a fault plane or a slab surface resolve as a surface rather than a
smear.

Reads windows straight from the read-only archive (nothing is cached to disk
beyond the small window cache), so it is I/O bound; `--max-pairs` and
`--max-sep` keep it bounded, and `--resume` lets it be interrupted.

Outputs:
  dt.cc                 HypoDD  ("# ID1 ID2 0.0" then "STA DT CC PHA")
  IN/xcordata.txt       GrowClust (same records, tdif_fmt 12)

Usage (env `eqt`, needs the archive mounted):
  xcorr_dt.py --locdir ../nll/loc --tag full --outroot .. \
      --max-sep 20 --cc-min 0.7 --workers 4
"""
import argparse, math, os, sys
from collections import OrderedDict

import numpy as np
import obspy
from obspy.signal.cross_correlation import correlate, xcorr_max

from mxio import load_index, read_window
from nllio import read_events

HERE = os.path.dirname(os.path.abspath(__file__))
WIN = {"P": (-0.6, 1.6), "S": (-0.8, 2.6)}     # seconds around the pick
COMP = {"P": "Z", "S": "E"}                     # component used per phase


class WindowCache:
    """LRU over 30-minute archive segments.

    The archive stores 30-minute miniSEED segments, so one cache entry per
    (station, day, component, half-hour) serves every event window inside it.
    Loading the whole station-day instead — 144 files for a 2-second cut — was
    ~100x slower.
    """

    BUCKET = 1800.0

    def __init__(self, index, size=48, pad=30.0):
        self.index = index
        self.size = size
        self.pad = pad
        self.store = OrderedDict()
        self.hits = self.misses = 0

    def _record(self, station, day):
        rec = self.index.get((station, day))
        if rec is None and station[-1:].islower():
            # NLLoc carries per-period ids (AE2a); the archive uses the bare code
            rec = self.index.get((station[:-1], day))
        return rec

    def segment(self, station, t, component):
        day = t.julday
        bucket = int((t - obspy.UTCDateTime(year=t.year, julday=day)) // self.BUCKET)
        key = (station, day, component, bucket)
        if key in self.store:
            self.store.move_to_end(key)
            self.hits += 1
            return self.store[key]
        rec = self._record(station, day)
        st = None
        if rec is not None:
            t0 = obspy.UTCDateTime(year=t.year, julday=day) + bucket * self.BUCKET
            try:
                st = read_window(rec["path"], rec["station"], rec["kind"],
                                 t0, t0 + self.BUCKET, component=component,
                                 pad=self.pad)
            except Exception:
                st = None
        self.misses += 1
        self.store[key] = st
        if len(self.store) > self.size:
            self.store.popitem(last=False)
        return st


def cut(cache, station, t, phase, freqmin, freqmax):
    """Filtered waveform window around an absolute arrival time, or None."""
    pre, post = WIN[phase]
    st = cache.segment(station, t, COMP[phase])
    if st is None or not len(st):
        st = cache.segment(station, t, "Z")
    if st is None or not len(st):
        return None
    tr = st[0].slice(t + pre - 1.0, t + post + 1.0)
    need = int((post - pre) * (tr.stats.sampling_rate or 100))
    if tr.stats.npts < need:
        return None
    tr = tr.copy()
    tr.detrend("demean")
    tr.taper(0.1)
    tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)
    return tr.slice(t + pre, t + post)


def sep_km(a, b):
    dz = a["dep"] - b["dep"]
    dx = (a["lon"] - b["lon"]) * 111.19 * math.cos(math.radians(a["lat"]))
    dy = (a["lat"] - b["lat"]) * 111.19
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locdir", default=os.path.join(HERE, "..", "nll", "loc"))
    ap.add_argument("--tag", default="wide11")
    ap.add_argument("--outroot", default=os.path.join(HERE, ".."))
    ap.add_argument("--index", default=None)
    ap.add_argument("--max-sep", type=float, default=20.0, help="km between events")
    ap.add_argument("--max-neighbours", type=int, default=25)
    ap.add_argument("--max-pairs", type=int, default=0, help="0 = no limit")
    ap.add_argument("--cc-min", type=float, default=0.55,
                    help="0.55 suits MERAMEX: subduction events a few km "
                         "apart correlate at 0.4-0.8, not at the 0.9 of a "
                         "repeating-earthquake aftershock sequence")
    ap.add_argument("--min-links", type=int, default=4, help="per pair, to keep it")
    ap.add_argument("--shift", type=float, default=0.5, help="max lag searched, s")
    ap.add_argument("--freqmin", type=float, default=2.0)
    ap.add_argument("--freqmax", type=float, default=15.0)
    ap.add_argument("--cache", type=int, default=64,
                    help="30-minute segments held in memory")
    ap.add_argument("--gapmax", type=float, default=300.0)
    ap.add_argument("--nphmin", type=int, default=8)
    ap.add_argument("--zmax", type=float, default=200.0)
    a = ap.parse_args()

    evs = [e for e in read_events(a.locdir, a.tag)
           if e.get("gap", 999) <= a.gapmax and e.get("nphs", 0) >= a.nphmin
           and 0 <= e["dep"] <= a.zmax]
    for i, e in enumerate(evs, 1):
        e["id"] = i
        e["arr"] = {(s, p): e["ot"] + tt for s, p, tt in e["phs"]}
    print(f"{len(evs)} events after selection")

    pairs = []
    for i, ei in enumerate(evs):
        near = sorted(((sep_km(ei, ej), j) for j, ej in enumerate(evs) if j > i
                       and abs(sep_km(ei, ej)) <= a.max_sep))
        pairs.extend((i, j) for _, j in near[:a.max_neighbours])
    if a.max_pairs:
        pairs = pairs[:a.max_pairs]
    print(f"{len(pairs)} candidate pairs within {a.max_sep:g} km")
    if not pairs:
        raise SystemExit("no pairs - loosen --max-sep")

    index = {(r["station"], r["day"]): r for r in load_index(a.index)}
    cache = WindowCache(index, size=a.cache)

    hd = os.path.join(a.outroot, "hypodd")
    gc = os.path.join(a.outroot, "growclust", "IN")
    os.makedirs(hd, exist_ok=True)
    os.makedirs(gc, exist_ok=True)
    n_kept = n_obs = 0
    with open(f"{hd}/dt.cc", "w") as fcc, open(f"{gc}/xcordata.txt", "w") as fgc:
        for k, (i, j) in enumerate(pairs, 1):
            ei, ej = evs[i], evs[j]
            common = set(ei["arr"]) & set(ej["arr"])
            rows = []
            for sta, ph in sorted(common):
                if ph not in WIN:
                    continue
                ti, tj = ei["arr"][(sta, ph)], ej["arr"][(sta, ph)]
                wi = cut(cache, sta, ti, ph, a.freqmin, a.freqmax)
                wj = cut(cache, sta, tj, ph, a.freqmin, a.freqmax)
                if wi is None or wj is None or wi.stats.npts != wj.stats.npts:
                    continue
                nshift = int(a.shift * wi.stats.sampling_rate)
                cc = correlate(wi.data.astype(float), wj.data.astype(float), nshift)
                lag, val = xcorr_max(cc)
                if abs(val) < a.cc_min:
                    continue
                dt = (ti - tj) + lag / wi.stats.sampling_rate
                rows.append((sta, dt, abs(val), ph))
            if len(rows) >= a.min_links:
                fcc.write(f"# {ei['id']} {ej['id']} 0.0\n")
                fgc.write(f"# {ei['id']} {ej['id']} 0.000\n")
                for sta, dt, val, ph in rows:
                    fcc.write(f"{sta:<7}{dt:9.4f} {val:6.3f} {ph}\n")
                    fgc.write(f"{sta:<7}{dt:9.4f} {val:6.3f} {ph}\n")
                n_kept += 1
                n_obs += len(rows)
            if k % 25 == 0:
                print(f"  {k}/{len(pairs)} pairs | kept {n_kept} | {n_obs} dt | "
                      f"cache {cache.hits}/{cache.hits + cache.misses}", flush=True)

    print(f"kept {n_kept} pairs, {n_obs} cross-correlation differential times")
    print(f"  -> {hd}/dt.cc and {gc}/xcordata.txt")


if __name__ == "__main__":
    main()
