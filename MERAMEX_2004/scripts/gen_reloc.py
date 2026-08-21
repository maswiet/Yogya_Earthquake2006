#!/usr/bin/env python3
"""Build HypoDD and GrowClust inputs from a located MERAMEX run.

Relative relocation is what turns a cloud of absolute hypocentres into visible
structure: common path effects cancel between nearby event pairs, so a fault
plane or a slab surface sharpens from tens of kilometres of scatter down to the
scale of the differential-time precision.

Two caveats specific to this dataset, both configurable here:

* MERAMEX seismicity is spread over ~500 km and 0-250 km depth, far sparser than
  an aftershock sequence, so the pair-linking limits have to be looser than the
  2006 Yogyakarta settings or almost nothing links.
* Double-difference codes assume rays that are near-parallel for a pair. That
  holds for shallow crustal clusters and reasonably for the intermediate-depth
  slab, but degrades for the deepest events; `--zmax` keeps the relocation to
  the depth range where the assumption is defensible.

Usage:
  gen_reloc.py --locdir ../nll/loc --tag full --outroot .. \
      --gapmax 250 --nphmin 8 --rmsmax 0.6 --zmax 200
"""
import argparse, json, os, string

from nllio import read_events

HERE = os.path.dirname(os.path.abspath(__file__))

# Same 1-D model as gen_nll.py, clipped to z >= 0 (both codes want positive
# depths): VELEST minimum-1D crust from the 2006 aftershocks + ak135 mantle.
MODEL = [(0.0, 2.90), (0.7, 4.30), (2.0, 4.65), (4.0, 5.49), (7.0, 5.49),
         (10.0, 6.30), (13.0, 6.39), (16.0, 6.55), (22.0, 6.80), (30.0, 7.20),
         (40.0, 8.04), (80.0, 8.05), (120.0, 8.18), (165.0, 8.30), (210.0, 8.48)]
VPVS = 1.75


def station_table(periods_path, keep=None):
    """{period_id: (lat, lon)} using the same <CODE><suffix> ids as gen_nll.py."""
    periods = json.load(open(periods_path))
    out = {}
    for code, info in periods.items():
        sites = info["sites"]
        sufs = [""] if len(sites) == 1 else list(string.ascii_lowercase)
        for suf, s in zip(sufs, sites):
            pid = f"{code}{suf}"
            if keep is None or pid in keep:
                out[pid] = (s["lat"], s["lon"])
    return out


def write_hypodd(hd, events, stations, a):
    os.makedirs(hd, exist_ok=True)
    with open(f"{hd}/hypoDD.pha", "w") as f:
        for eid, e in enumerate(events, 1):
            ot = e["ot"]
            f.write(f"# {ot.year} {ot.month} {ot.day} {ot.hour} {ot.minute} "
                    f"{ot.second + ot.microsecond * 1e-6:6.3f} "
                    f"{e['lat']:.4f} {e['lon']:.4f} {e['dep']:.2f} "
                    f"0.0 {e.get('eh', 0):.2f} {e.get('ez', 0):.2f} "
                    f"{e.get('rms', 0):.2f} {eid}\n")
            for sta, ph, tt in e["phs"]:
                if sta in stations:
                    f.write(f"{sta:<7}{tt:8.3f} {1.0 if ph == 'P' else 0.5:.2f} {ph}\n")
    with open(f"{hd}/station.dat", "w") as f:
        for pid, (la, lo) in sorted(stations.items()):
            f.write(f"{pid:<7}{la:9.4f} {lo:10.4f}\n")

    open(f"{hd}/ph2dt.inp", "w").write(
        f"* ph2dt.inp - MERAMEX 2004\nstation.dat\nhypoDD.pha\n"
        f"*MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS\n"
        f"   0   {a.maxdist:6.0f} {a.maxsep:6.0f} {a.maxngh:6d} "
        f"{a.minlnk:6d} {a.minobs:6d} {a.maxobs:6d}\n")

    # IDAT 2 = catalog differential times only; 3 = also waveform dt.cc.
    # Asking for 3 without a dt.cc file makes hypoDD abort immediately.
    idat = 3 if os.path.exists(f"{hd}/dt.cc") else 2
    mod_top = " ".join(f"{d:.1f}" for d, _ in MODEL)
    mod_v = " ".join(f"{v:.2f}" for _, v in MODEL)
    open(f"{hd}/hypoDD.inp", "w").write(f"""* hypoDD.inp - MERAMEX 2004
*--- INPUT FILE SELECTION
dt.cc
dt.ct
event.sel
station.dat
*--- OUTPUT FILE SELECTION
hypoDD.loc
hypoDD.reloc
hypoDD.sta
hypoDD.res
hypoDD.src
*--- DATA TYPE SELECTION: IDAT IPHA DIST
    {idat}     3   {a.maxdist:.0f}
*--- EVENT CLUSTERING: OBSCC OBSCT
    0     {a.minobs}
*--- SOLUTION CONTROL: ISTART ISOLV NSET
    2     2     4
*--- DATA WEIGHTING AND REWEIGHTING
* NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP
    5    -9    -9   -9   -9   1.0   0.5   -9   -9   80
    5    -9    -9   -9   -9   1.0   0.5   6    {a.maxsep:.0f}  80
    5    -9    -9   -9   -9   1.0   0.8   6    {a.maxsep / 2:.0f}  60
    5    -9    -9   -9   -9   1.0   0.8   4    {a.maxsep / 4:.0f}  60
*--- 1D MODEL: NLAY RATIO TOP VEL
   {len(MODEL)}   {VPVS}
{mod_top}
{mod_v}
*--- CLUSTER/EVENT SELECTION: CID  ID
    0
""")


def write_growclust(gc, events, stations, a):
    for d in ("IN", "OUT", "TT"):
        os.makedirs(f"{gc}/{d}", exist_ok=True)
    with open(f"{gc}/IN/evlist.txt", "w") as f:
        for eid, e in enumerate(events, 1):
            ot = e["ot"]
            f.write(f"{ot.year:4d} {ot.month:2d} {ot.day:2d} {ot.hour:2d} "
                    f"{ot.minute:2d} {ot.second + ot.microsecond * 1e-6:6.3f} "
                    f"{e['lat']:9.5f} {e['lon']:10.5f} {e['dep']:8.3f} "
                    f"0.000 {e.get('eh', 0):6.3f} {e.get('ez', 0):6.3f} "
                    f"{e.get('rms', 0):6.3f} {eid:9d}\n")
    with open(f"{gc}/IN/stlist.txt", "w") as f:
        for pid, (la, lo) in sorted(stations.items()):
            f.write(f"{pid:<6}{la:9.4f} {lo:10.4f}\n")
    with open(f"{gc}/IN/vzmodel.txt", "w") as f:      # piecewise constant
        for i in range(len(MODEL) - 1):
            d0, v0 = MODEL[i]
            d1, _ = MODEL[i + 1]
            f.write(f"{d0:7.2f} {v0:5.2f} 0.00\n{d1:7.2f} {v0:5.2f} 0.00\n")
        d, v = MODEL[-1]
        f.write(f"{d:7.2f} {v:5.2f} 0.00\n{a.zmax + 60:7.2f} {v:5.2f} 0.00\n")

    open(f"{gc}/growclust.inp", "w").write(f"""* GrowClust control - MERAMEX 2004
* evlist_fmt
1
* fin_evlist
IN/evlist.txt
* stlist_fmt
1
* fin_stlist
IN/stlist.txt
* xcordat_fmt  tdif_fmt
1  12
* fin_xcordat
IN/xcordata.txt
* fin_vzmdl
IN/vzmodel.txt
* fout_vzfine
TT/vzfine.txt
* fout_pTT
TT/tt.pg
* fout_sTT
TT/tt.sg
* vpvs_factor  rayparam_min
  {VPVS}       0.0
* tt_dep0 tt_dep1 tt_ddep
  0.  {a.zmax + 60:.0f}.  2.
* tt_del0 tt_del1 tt_ddel
  0.  {a.maxdist:.0f}.  4.
* rmin delmax rmsmax
  0.0  {a.maxdist:.0f}  {a.rmsmax_gc}
* rpsavgmin rmincut ngoodmin iponly
  0  0  0  0
* nboot nbranch_min
  {a.nboot}  {a.nbranch_min}
* fout_cat
OUT/out.growclust_cat
* fout_clust
OUT/out.growclust_clust
* fout_log
OUT/out.growclust_log
* fout_boot
OUT/out.growclust_boot
""")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locdir", default=os.path.join(HERE, "..", "nll", "loc"))
    ap.add_argument("--tag", default="wide11")
    ap.add_argument("--outroot", default=os.path.join(HERE, ".."))
    ap.add_argument("--periods", default=os.path.join(HERE, "..", "config",
                                                      "stations_periods.json"))
    # event selection
    ap.add_argument("--gapmax", type=float, default=300.0)
    ap.add_argument("--nphmin", type=int, default=8)
    ap.add_argument("--rmsmax", type=float, default=0.6)
    ap.add_argument("--zmax", type=float, default=200.0)
    ap.add_argument("--ehmax", type=float, default=25.0)
    # pair linking
    ap.add_argument("--maxdist", type=float, default=400.0, help="km, event-station")
    ap.add_argument("--maxsep", type=float, default=20.0, help="km, event-event")
    ap.add_argument("--maxngh", type=int, default=30)
    ap.add_argument("--minlnk", type=int, default=6)
    ap.add_argument("--minobs", type=int, default=6)
    ap.add_argument("--maxobs", type=int, default=60)
    # growclust
    ap.add_argument("--rmsmax-gc", type=float, default=0.50)
    ap.add_argument("--nboot", type=int, default=100)
    ap.add_argument("--nbranch-min", type=int, default=2)
    a = ap.parse_args()

    evs = read_events(a.locdir, a.tag)
    sel = [e for e in evs
           if e.get("gap", 999) <= a.gapmax
           and e.get("nphs", 0) >= a.nphmin
           and e.get("rms", 9) <= a.rmsmax
           and e.get("eh", 999) <= a.ehmax
           and 0 <= e["dep"] <= a.zmax]
    if not sel:
        raise SystemExit("no events passed the selection - loosen the cuts")

    used = {s for e in sel for s, _, _ in e["phs"]}
    stations = station_table(a.periods, keep=used)

    hd = os.path.join(a.outroot, "hypodd")
    gc = os.path.join(a.outroot, "growclust")
    write_hypodd(hd, sel, stations, a)
    write_growclust(gc, sel, stations, a)

    nph = sum(1 for e in sel for s, _, _ in e["phs"] if s in stations)
    print(f"{len(evs)} located -> {len(sel)} selected "
          f"(gap<={a.gapmax:.0f}, nphs>={a.nphmin}, rms<={a.rmsmax}, "
          f"errH<={a.ehmax:.0f} km, z<={a.zmax:.0f} km)")
    print(f"  {nph} phases on {len(stations)} station-sites")
    print(f"  wrote {hd}/ and {gc}/")


if __name__ == "__main__":
    main()
