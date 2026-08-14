"""Generate the per-event waveform gallery + browser DIRECTLY from the SeisBench
archive (fast, drive-independent, guaranteed consistent with the deposited data).

For each event: one PNG showing its stations' vertical-component windows with the
labelled P/S samples; plus a sortable index.html linking every event to its PNG.
Far faster than re-reading the raw volume, and the previews match exactly what is
archived in waveforms.hdf5.

Usage: build_gallery.py --out ~/Yogya2006_Zenodo/event_browser [--limit N]
"""
import os, argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import seisbench.data as sbd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", default=f"{ROOT}/seisbench")
    ap.add_argument("--out", default=os.path.expanduser("~/Yogya2006_Zenodo/event_browser"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--nmax", type=int, default=6, help="max stations per preview")
    a = ap.parse_args()
    pdir = f"{a.out}/previews"; os.makedirs(pdir, exist_ok=True)

    d = sbd.WaveformDataset(a.archive)
    m = d.metadata.reset_index(drop=True)
    m["_row"] = np.arange(len(m))
    cat = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
    q = pd.read_csv(f"{ROOT}/full/catalog_quality.csv")
    cat = cat.merge(q[["evid", "pass"]], on="evid", how="left")

    by_ev = m.groupby("source_id")
    evids = sorted(by_ev.groups)
    if a.limit:
        evids = evids[:a.limit]
    print(f"generating {len(evids)} event previews -> {pdir}")

    made = 0
    for k, evid in enumerate(evids, 1):
        rows = by_ev.get_group(evid)
        rows = rows.reindex(rows.path_ep_distance_km.sort_values().index)[:a.nmax]
        ev = cat[cat.evid == evid]
        if not len(ev):
            continue
        ev = ev.iloc[0]
        n = len(rows)
        fig, axes = plt.subplots(n, 1, figsize=(8, 1.15*n+0.5), sharex=True, squeeze=False)
        for i, (_, r) in enumerate(rows.iterrows()):
            w = d.get_waveforms(int(r._row))
            t = np.arange(w.shape[1])/r.trace_sampling_rate_hz
            ax = axes[i][0]
            ax.plot(t, w[2], lw=0.5, color="0.15")           # E component
            ax.axvline(r.trace_p_arrival_sample/r.trace_sampling_rate_hz,
                       color="tab:blue", lw=1)
            if pd.notna(r.trace_s_arrival_sample):
                ax.axvline(r.trace_s_arrival_sample/r.trace_sampling_rate_hz,
                           color="tab:red", lw=1)
            ax.text(0.01, 0.9, f"{r.station_code} {r.path_ep_distance_km:.0f} km",
                    transform=ax.transAxes, fontsize=8, va="top", weight="bold")
            ax.set_yticks([]); ax.tick_params(labelsize=7)
        axes[-1][0].set_xlabel("time in window (s)  —  P blue, S red", fontsize=8)
        fig.suptitle(f"evid {evid} — {ev.time[:19]} UTC — ML {ev.ML:+.2f} "
                     f"(tied {ev.ML+0.41:+.2f}) — {ev.latitude:.3f},{ev.longitude:.3f} — "
                     f"{ev.depth:.1f} km — gap {ev.gap:.0f}°", fontsize=8.5)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(f"{pdir}/ev{evid:05d}.png", dpi=85); plt.close(fig)
        made += 1
        if k % 500 == 0:
            print(f"  [{k}/{len(evids)}] made {made}", flush=True)

    # index.html
    rows = []
    for _, e in cat[cat.evid.isin(evids)].iterrows():
        img = f"previews/ev{int(e.evid):05d}.png"
        rows.append(f"<tr><td><a href='{img}'>{int(e.evid)}</a></td><td>{e.time[:19]}</td>"
                    f"<td>{e.ML:+.2f}</td><td>{e.ML+0.41:+.2f}</td><td>{e.latitude:.3f}</td>"
                    f"<td>{e.longitude:.3f}</td><td>{e.depth:.1f}</td><td>{int(e.n_sta)}</td>"
                    f"<td>{e.gap:.0f}</td><td>{'✓' if e['pass'] else ''}</td></tr>")
    html = f"""<!doctype html><meta charset=utf-8>
<title>Yogyakarta 2006 aftershock catalogue — event browser</title>
<style>body{{font:14px system-ui;margin:1.5em}}table{{border-collapse:collapse}}
th,td{{padding:3px 9px;border-bottom:1px solid #ddd;text-align:right}}
th{{cursor:pointer;background:#f4f4f4;position:sticky;top:0}}
td:nth-child(2){{text-align:left}}a{{color:#06c;text-decoration:none}}</style>
<h2>Yogyakarta 2006 aftershock catalogue — {len(evids)} events</h2>
<p>Click an event id to view its multi-station waveforms (E component; P blue, S red).
ML_tied = ML + 0.41 (local-scale tie; see the data description). Click a column header to sort.</p>
<table id=t><thead><tr><th>evid</th><th>time (UTC)</th><th>ML</th><th>ML_tied</th>
<th>lat</th><th>lon</th><th>depth</th><th>nsta</th><th>gap</th><th>QC</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table>
<script>document.querySelectorAll('th').forEach((h,i)=>h.onclick=()=>{{
const t=document.querySelector('#t tbody');[...t.rows].sort((a,b)=>{{
const x=a.cells[i].innerText,y=b.cells[i].innerText,nx=parseFloat(x),ny=parseFloat(y);
return isNaN(nx)?x.localeCompare(y):nx-ny;}}).forEach(r=>t.appendChild(r));}});</script>"""
    open(f"{a.out}/index.html", "w").write(html)
    print(f"done: {made} previews + index.html -> {a.out}")


if __name__ == "__main__":
    main()
