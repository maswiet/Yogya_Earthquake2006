#!/usr/bin/env python3
"""Streaming EQTransformer picker for the MERAMEX 2004 archive (env `eqt`).

Reads one station-day at a time straight off the read-only archive drives,
relabels it, runs EQTransformer, appends the picks to a CSV and records the
station-day in a progress file. No waveform data is ever written to disk, which
matters: the full archive is ~500 GB and the internal disk has ~80 GB free.

Usage:
  run_picks.py --out ../pilot/picks.csv --progress ../pilot/done.txt \
      --stations AI1,AI3,... --dmin 155 --dmax 165 --device mps
"""
import argparse, csv, os, sys, time

import obspy
import seisbench.models as sbm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mxio import build_stream, load_index

FIELDS = ["station", "kind", "day", "phase", "peak_time", "probability", "channel"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--progress", required=True)
    ap.add_argument("--index", default=None)
    ap.add_argument("--model", default="original")
    ap.add_argument("--device", default="mps")
    ap.add_argument("--pthr", type=float, default=0.3)
    ap.add_argument("--sthr", type=float, default=0.3)
    ap.add_argument("--dthr", type=float, default=0.3)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--dmin", type=int, default=127)
    ap.add_argument("--dmax", type=int, default=282)
    ap.add_argument("--stations", default=None, help="comma list; default = all land")
    ap.add_argument("--kinds", default="EDL,SAM")
    ap.add_argument("--min-bytes", type=int, default=5_000_000,
                    help="skip station-days smaller than this (transit / dead days)")
    a = ap.parse_args()

    out, prog = os.path.abspath(a.out), os.path.abspath(a.progress)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    done = set(open(prog).read().split()) if os.path.exists(prog) else set()

    kinds = set(a.kinds.split(","))
    want = set(a.stations.split(",")) if a.stations else None
    work = [r for r in load_index(a.index)
            if r["kind"] in kinds
            and a.dmin <= r["day"] <= a.dmax
            and (want is None or r["station"] in want)
            and r["nbytes"] >= a.min_bytes]
    work.sort(key=lambda r: (r["day"], r["station"]))
    todo = [r for r in work if f"{r['station']}/{r['day']}" not in done]
    print(f"{len(work)} station-days selected, {len(todo)} still to do "
          f"({sum(r['nbytes'] for r in todo)/1e9:.1f} GB to read)", flush=True)
    if not todo:
        return

    model = sbm.EQTransformer.from_pretrained(a.model)
    model.to(a.device)
    model.eval()

    new = not os.path.exists(out)
    fh = open(out, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
    pf = open(prog, "a")

    t_start = time.time()
    npick_tot = 0
    for i, r in enumerate(todo, 1):
        t0 = time.time()
        key = f"{r['station']}/{r['day']}"
        try:
            st = build_stream(r["path"], r["station"], r["kind"])
        except Exception as e:
            print(f"  {key}: READ FAILED {e}", flush=True)
            st = None
        n = 0
        if st is not None and len(st) >= 3:
            try:
                picks = model.classify(st, batch_size=a.batch,
                                       P_threshold=a.pthr, S_threshold=a.sthr,
                                       detection_threshold=a.dthr).picks
                ch = ",".join(sorted({t.stats.channel for t in st}))
                for p in picks:
                    w.writerow(dict(station=r["station"], kind=r["kind"], day=r["day"],
                                    phase=p.phase, peak_time=str(p.peak_time),
                                    probability=round(float(p.peak_value), 4),
                                    channel=ch))
                n = len(picks)
            except Exception as e:
                print(f"  {key}: PICK FAILED {e}", flush=True)
        fh.flush()
        pf.write(key + "\n")
        pf.flush()
        npick_tot += n
        el = time.time() - t_start
        eta = el / i * (len(todo) - i)
        print(f"[{i}/{len(todo)}] {key:12s} ntr={0 if st is None else len(st)} "
              f"picks={n:5d}  {time.time()-t0:5.1f}s  total={npick_tot}  "
              f"ETA {eta/60:.0f} min", flush=True)

    fh.close(); pf.close()
    print(f"DONE  {npick_tot} picks -> {out}  ({(time.time()-t_start)/60:.1f} min)")


if __name__ == "__main__":
    main()
