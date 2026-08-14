"""Proof-of-concept event browser for the catalogue data product.

Generates, per event: a multi-station Wood-Anderson waveform preview PNG, and a
single index.html with a sortable table (time, location, depth, ML, gap, quality)
where each row links to its preview. This is the "clickable event table" for the
Seismica Data Report's GitHub/Pages appendix; the PNGs + waveform data archive
themselves go to a DOI-citable Zenodo deposit (GitHub is not DOI-archival).

Usage:
  build_event_browser.py --sample 12         # POC: 12 events spread over ML
  build_event_browser.py --all               # full run (hours; reads raw volume)

Runs in the eqt env (obspy). Needs the raw EDL volume mounted for the PNGs.
"""
import os, re, glob, math, argparse
import numpy as np, pandas as pd
import obspy
from obspy.core import UTCDateTime
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_CANDS = ["/Volumes/Untitled/DATA-GFZ-Gempa-JOgja-tahap-2",
          "/Volumes/Untitled 1/DATA-GFZ-Gempa-JOgja-tahap-2"]
BASE = os.environ.get("EDL_BASE") or next((p for p in _CANDS if os.path.isdir(p)), _CANDS[0])
COMP = {"pri1": "N", "pri2": "E", "pri0": "Z"}
year_re = re.compile(r"e\d{4}(\d{2})\d+\.pri1$")
OUTDIR = f"{ROOT}/figures/events"

w0 = 2*math.pi*1.0; h = 0.707
PAZ_L4 = {"poles": [complex(-h*w0, w0*math.sqrt(1-h*h)), complex(-h*w0, -w0*math.sqrt(1-h*h))],
          "zeros": [0j, 0j], "gain": 1.0, "sensitivity": 1.7e8}
PAZ_WA = {"poles": [complex(-6.2832, -4.7124), complex(-6.2832, 4.7124)],
          "zeros": [0j, 0j], "gain": 1.0, "sensitivity": 2080.0}
HYP = sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))


def picks_of(evid):
    rows = []
    for line in open(HYP[evid]):
        if " > " in line and "GAU" in line:
            lf = line.split(" > ")[0].split(); rf = line.split(" > ")[1].split()
            if lf[4] in ("P", "S"):
                t = (UTCDateTime(f"{lf[6][:4]}-{lf[6][4:6]}-{lf[6][6:8]}T"
                                 f"{lf[7][:2]}:{lf[7][2:4]}:00") + float(lf[8]))
                rows.append((lf[0], lf[4], t, float(rf[6])))
    return rows


def read_wa(sta, t0, t1):
    folder = f"tf30{sta[2:4]}"
    st = obspy.Stream()
    for jd in sorted({t0.julday, t1.julday}):
        d = os.path.join(BASE, folder, f"{jd:03d}")
        if not os.path.isdir(d):
            continue
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
            st += s
    if not len(st):
        return None
    st.merge(method=1, fill_value=None); st.trim(t0-5, t1+5)
    for tr in st:
        if hasattr(tr.data, "filled"):
            tr.data = tr.data.filled(0.0)
    if not len(st):
        return None
    st.detrend("demean")
    try:
        st.simulate(paz_remove=PAZ_L4, paz_simulate=PAZ_WA, water_level=10)
        st.filter("bandpass", freqmin=1.0, freqmax=20.0, corners=4, zerophase=True)
    except Exception:
        return None
    return st


def make_preview(evid, ev, nmax=6):
    pk = picks_of(evid)
    sp = sorted([p for p in pk if p[1] == "S"], key=lambda p: p[3])[:nmax]
    if not sp:
        return False
    fig, axes = plt.subplots(len(sp), 1, figsize=(8, 1.3*len(sp)+0.5), sharex=True, squeeze=False)
    ot = pd.to_datetime(ev.time)
    ok = False
    for i, (sta, _, stime, dist) in enumerate(sp):
        pt = [p for p in pk if p[0] == sta and p[1] == "P"]
        ptime = pt[0][2] if pt else stime-2
        st = read_wa(sta, ptime-8, stime+18)
        ax = axes[i][0]
        if st:
            tr = max(st, key=lambda x: x.stats.npts)
            ax.plot(tr.times()+(tr.stats.starttime-stime), tr.data, lw=0.5, color="0.15")
            ax.axvline(float(ptime-stime), color="tab:blue", lw=1)
            ax.axvline(0, color="tab:red", lw=1)
            ok = True
        ax.text(0.01, 0.9, f"{sta} {dist:.0f} km", transform=ax.transAxes, fontsize=8,
                va="top", weight="bold")
        ax.set_yticks([]); ax.tick_params(labelsize=7)
    axes[-1][0].set_xlabel("time relative to S pick (s)", fontsize=8)
    fig.suptitle(f"evid {evid} — {ev.time[:19]} UTC — ML {ev.ML:+.2f} — "
                 f"{ev.latitude:.3f},{ev.longitude:.3f} — {ev.depth:.1f} km — gap {ev.gap:.0f}°",
                 fontsize=9)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs(OUTDIR, exist_ok=True)
    plt.savefig(f"{OUTDIR}/ev{evid:05d}.png", dpi=90); plt.close(fig)
    return ok


def write_index(cat):
    rows = []
    for _, e in cat.iterrows():
        img = f"ev{int(e.evid):05d}.png"
        exists = os.path.exists(f"{OUTDIR}/{img}")
        link = f'<a href="{img}">{int(e.evid)}</a>' if exists else str(int(e.evid))
        rows.append(f"<tr><td>{link}</td><td>{e.time[:19]}</td><td>{e.ML:+.2f}</td>"
                    f"<td>{e.latitude:.3f}</td><td>{e.longitude:.3f}</td>"
                    f"<td>{e.depth:.1f}</td><td>{int(e.n_sta)}</td><td>{e.gap:.0f}</td>"
                    f"<td>{'✓' if e.get('pass', True) else ''}</td></tr>")
    html = f"""<!doctype html><meta charset=utf-8>
<title>Yogyakarta 2006 aftershock catalogue — event browser</title>
<style>body{{font:14px system-ui;margin:1.5em}}table{{border-collapse:collapse}}
th,td{{padding:3px 9px;border-bottom:1px solid #ddd;text-align:right}}
th{{cursor:pointer;background:#f4f4f4}}td:nth-child(2){{text-align:left}}
a{{color:#06c;text-decoration:none}}</style>
<h2>Yogyakarta 2006 aftershock catalogue ({len(cat)} events shown)</h2>
<p>Click an event id to view its multi-station Wood-Anderson waveforms.
Full data + waveform archive: Zenodo DOI (placeholder).</p>
<table id=t><thead><tr>
<th>evid</th><th>time (UTC)</th><th>ML</th><th>lat</th><th>lon</th>
<th>depth</th><th>nsta</th><th>gap</th><th>QC</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table>
<script>
document.querySelectorAll('th').forEach((h,i)=>h.onclick=()=>{{
 const t=document.querySelector('#t tbody');
 [...t.rows].sort((a,b)=>{{const x=a.cells[i].innerText,y=b.cells[i].innerText;
 const nx=parseFloat(x),ny=parseFloat(y);
 return isNaN(nx)?x.localeCompare(y):nx-ny;}}).forEach(r=>t.appendChild(r));}});
</script>"""
    open(f"{OUTDIR}/index.html", "w").write(html)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=12)
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    cat = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    q = pd.read_csv(f"{ROOT}/full/catalog_quality.csv")
    cat = cat.merge(q[["evid", "pass"]], on="evid", how="left")
    if a.all:
        sel = cat
    else:
        sel = cat[cat['pass']].sort_values("ML").iloc[::max(1, len(cat)//a.sample)][:a.sample]
    print(f"generating previews for {len(sel)} events -> {OUTDIR}")
    made = 0
    for _, e in sel.iterrows():
        try:
            if make_preview(int(e.evid), e):
                made += 1
        except Exception as ex:
            print(f"  ev{int(e.evid)}: {ex}")
    write_index(cat if a.all else sel)
    print(f"previews made: {made}/{len(sel)}; wrote {os.path.relpath(OUTDIR, ROOT)}/index.html")


if __name__ == "__main__":
    main()
