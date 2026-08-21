#!/usr/bin/env python3
"""Build the MERAMEX station-day index from the two read-only archive drives.

Walks <ROOT>/<DOY-folder>/<STA>/<DDD>/ and records, for every station-day, the
directory, the file counts per extension, and the total byte size. The result
(`config/station_days.csv`) is the work list for the picking driver.

Usage:  index_data.py [--out ../config/station_days.csv]
"""
import argparse, csv, os, sys, time

ROOTS = [
    "/Volumes/Untitled/MERAMEX DATA",
    "/Volumes/Untitled 1/MERAMEX_LANJUTAN",
]
DATA_EXT = ("PRI0", "PRI1", "PRI2")          # EDL short-period / T40
SOH_EXT = ("GPS", "GST", "MSG", "PLL")


def kind_of(exts):
    if any(e.startswith("PRI") for e in exts):
        return "EDL"
    if any(e.startswith("SP") for e in exts):
        return "SAM"
    return "OTHER"


def scan(roots):
    for root in roots:
        if not os.path.isdir(root):
            print(f"WARNING: {root} not mounted", file=sys.stderr)
            continue
        for doydir in sorted(os.listdir(root)):
            p1 = os.path.join(root, doydir)
            if doydir.startswith(("@", "$", ".")) or not os.path.isdir(p1):
                continue
            for sta in sorted(os.listdir(p1)):
                p2 = os.path.join(p1, sta)
                if not os.path.isdir(p2):
                    continue
                for day in sorted(os.listdir(p2)):
                    p3 = os.path.join(p2, day)
                    if not (day.isdigit() and os.path.isdir(p3)):
                        continue
                    exts, nbytes = {}, 0
                    try:
                        with os.scandir(p3) as it:
                            for e in it:
                                if not e.is_file():
                                    continue
                                ext = e.name.rsplit(".", 1)[-1].upper()
                                # SAM files end .2004.159.000000 -> use channel field
                                if ext.isdigit():
                                    parts = e.name.split(".")
                                    ext = parts[2].upper() if len(parts) > 3 else "UNK"
                                exts[ext] = exts.get(ext, 0) + 1
                                try:
                                    nbytes += e.stat().st_size
                                except OSError:
                                    pass
                    except OSError as err:
                        yield dict(root=root, doydir=doydir, station=sta, day=int(day),
                                   kind="ERR", nseg=0, nbytes=0, exts=str(err))
                        continue
                    nseg = sum(v for k, v in exts.items() if k not in SOH_EXT)
                    yield dict(root=root, doydir=doydir, station=sta, day=int(day),
                               kind=kind_of(exts), nseg=nseg, nbytes=nbytes,
                               exts=";".join(f"{k}={v}" for k, v in sorted(exts.items())))


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--out", default=os.path.join(here, "..", "config", "station_days.csv"))
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    t0 = time.time()
    n = 0
    with open(out, "w", newline="") as fh:
        w = None
        for row in scan(ROOTS):
            if w is None:
                w = csv.DictWriter(fh, fieldnames=list(row.keys()))
                w.writeheader()
            w.writerow(row)
            n += 1
            if n % 500 == 0:
                print(f"{n} station-days  {time.time()-t0:.0f}s", flush=True)
    print(f"wrote {n} station-days to {out} in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
