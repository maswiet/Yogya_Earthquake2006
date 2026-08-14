"""Feasibility assessment: how fine a velocity model could the EQTransformer
catalogue resolve in local earthquake tomography?

Tomographic resolution is set by ray coverage -- ray DENSITY and, decisively,
ray CROSSING-ANGLE diversity -- not by event count alone. With a fixed 12/17-
station geometry over a compact aftershock cloud, more events raise ray density
but the angular diversity is bounded by the network. This script quantifies both
from straight-line source->station rays (a standard first-order proxy; a full
answer needs a LOTOS/checkerboard run, noted below).

For candidate grid sizes it computes, per cell:
  - hit count (rays piercing the cell)  -- density
  - number of distinct 45-deg azimuth sectors with rays  -- crossing diversity
A cell is "resolvable" if hit >= HIT_MIN and it samples >= AZ_MIN sectors.

To isolate what the LARGER catalogue buys, it repeats the whole calculation for
a random 588-event subset (matching Diambama et al. 2019) using the identical
geometry and method, and reports resolvable-fraction vs grid size for both.
"""
import os, re, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HIT_MIN = 20            # rays per cell for a reliable inversion
AZ_MIN = 3             # distinct 45-deg azimuth sectors (of 8) sampled
LAT0, LON0 = -7.92, 110.44
KML = 111.2
def to_km(lat, lon): return (lon-LON0)*KML*np.cos(np.radians(LAT0)), (lat-LAT0)*KML
# target volume (km, relative to LAT0/LON0), covering the aftershock zone
XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX = -18, 18, -18, 18, 0, 20


def stations():
    st = {}
    for line in open(f"{ROOT}/nll/stations_gtsrce.txt"):
        f = line.split()
        if len(f) >= 5 and f[0] == "GTSRCE":
            x, y = to_km(float(f[3]), float(f[4]))
            st[f[1]] = (x, y, 0.0)
    return st


def rays():
    """(event xyz, station xyz) for every P/S pick of every QC event."""
    q = pd.read_csv(f"{ROOT}/full/catalog_quality.csv")
    qpass = set(q[q["pass"]].evid)
    st = stations()
    HYP = sorted(glob.glob(f"{ROOT}/nll/loc_v/yogya_v.2*.grid0.loc.hyp"))
    geo = re.compile(r"Lat\s+(-?[\d.]+)\s+Long\s+(-?[\d.]+)\s+Depth\s+(-?[\d.]+)")
    ev_rays = []          # list of (evx,evy,evz, [station xyz,...])
    for evid, path in enumerate(HYP):
        if evid not in qpass:
            continue
        ex = ey = ez = None; ss = []
        for line in open(path):
            if line.startswith("GEOGRAPHIC"):
                m = geo.search(line)
                if m:
                    ex, ey = to_km(float(m.group(1)), float(m.group(2)))
                    ez = float(m.group(3))
            elif " > " in line and "GAU" in line:
                lf = line.split(" > ")[0].split()
                if lf[4] in ("P", "S") and lf[0] in st:
                    ss.append(st[lf[0]])
        if ex is not None and ss:
            ev_rays.append((ex, ey, ez, ss))
    return ev_rays


def coverage(ev_rays, dx, rng):
    """Return (hit, az_sectors) 3D arrays for grid spacing dx over the target box.
    rng: iterable of event indices to include (for subsampling)."""
    nx = int((XMAX-XMIN)/dx); ny = int((YMAX-YMIN)/dx); nz = int((ZMAX-ZMIN)/dx)
    hit = np.zeros((nx, ny, nz), np.int32)
    azm = np.zeros((nx, ny, nz), np.uint8)      # bitmask of 8 azimuth sectors
    step = dx/2.0
    for i in rng:
        ex, ey, ez, ss = ev_rays[i]
        for sx, sy, sz in ss:
            L = np.hypot(np.hypot(sx-ex, sy-ey), sz-ez)
            n = max(2, int(L/step))
            t = np.linspace(0, 1, n)
            xs = ex+(sx-ex)*t; ys = ey+(sy-ey)*t; zs = ez+(sz-ez)*t
            az = int((np.degrees(np.arctan2(sy-ey, sx-ex)) % 360)//45)
            ix = ((xs-XMIN)/dx).astype(int); iy = ((ys-YMIN)/dx).astype(int)
            iz = ((zs-ZMIN)/dx).astype(int)
            ok = (ix >= 0)&(ix < nx)&(iy >= 0)&(iy < ny)&(iz >= 0)&(iz < nz)
            cells = set(zip(ix[ok], iy[ok], iz[ok]))
            for c in cells:
                hit[c] += 1; azm[c] |= (1 << az)
    az_sectors = np.zeros_like(hit)
    for b in range(8):
        az_sectors += ((azm >> b) & 1).astype(np.int32)
    return hit, az_sectors


def resolvable_fraction(hit, az):
    """fraction of cells inside the aftershock footprint that are resolvable."""
    footprint = hit.sum(axis=2) > 0            # cells ever hit (map view)
    res = (hit >= HIT_MIN) & (az >= AZ_MIN)
    res_map = res.any(axis=2)
    denom = footprint.sum()
    return (res_map & footprint).sum()/max(denom, 1), res


def main():
    ev = rays()
    print(f"QC events with rays: {len(ev)}")
    # fixed random subset of 588 (Diambama size); index-based, deterministic
    idx_full = range(len(ev))
    stride = max(1, len(ev)//588)
    idx_sub = list(range(0, len(ev), stride))[:588]

    grids = [5, 3, 2, 1.5]
    rows = []
    cov_full = {}
    for dx in grids:
        hF, aF = coverage(ev, dx, idx_full)
        hS, aS = coverage(ev, dx, idx_sub)
        fF, resF = resolvable_fraction(hF, aF)
        fS, _ = resolvable_fraction(hS, aS)
        cov_full[dx] = (hF, aF, resF)
        rows.append((dx, fF, fS))
        print(f"  grid {dx:>4} km: resolvable footprint  full={100*fF:4.0f}%  "
              f"588-subset={100*fS:4.0f}%  (hit>= {HIT_MIN}, az>= {AZ_MIN}/8)")

    # figure
    fig = plt.figure(figsize=(16, 9))
    # A: resolvable fraction vs grid
    ax = fig.add_subplot(2, 3, 1)
    r = np.array(rows)
    ax.plot(r[:, 0], 100*r[:, 1], "o-", color="tab:blue", label=f"full ({len(ev)} ev)")
    ax.plot(r[:, 0], 100*r[:, 2], "s--", color="tab:red", label="588-event subset")
    ax.set_xlabel("grid cell size (km)"); ax.set_ylabel("resolvable footprint (%)")
    ax.set_title("A  Resolvable fraction vs grid size", fontsize=10)
    ax.invert_xaxis(); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # B/C: hit-count map at 2 km, two depths
    dx = 2.0; hF, aF, resF = cov_full[dx]
    for k, (zc, ttl) in enumerate([(8, "8 km depth"), (12, "12 km depth")]):
        ax = fig.add_subplot(2, 3, 2+k)
        iz = int((zc-ZMIN)/dx)
        m = hF[:, :, iz].T
        im = ax.imshow(np.log10(m+1), origin="lower",
                       extent=[XMIN, XMAX, YMIN, YMAX], cmap="viridis", aspect="equal")
        ax.contour(np.linspace(XMIN, XMAX, m.shape[1]), np.linspace(YMIN, YMAX, m.shape[0]),
                   (resF[:, :, iz].T).astype(float), levels=[0.5], colors="white", linewidths=1)
        fig.colorbar(im, ax=ax, label="log10 rays/cell")
        ax.set_title(f"B{k}  Ray hits @ {ttl}, 2 km grid\n(white = resolvable)", fontsize=9)
        ax.set_xlabel("E (km)"); ax.set_ylabel("N (km)")

    # D: cross-section hit count (N-S slice through x~0), 2 km
    ax = fig.add_subplot(2, 3, 5)
    ix0 = int((0-XMIN)/dx)
    sec = hF[ix0, :, :].T
    im = ax.imshow(np.log10(sec+1), origin="upper",
                   extent=[YMIN, YMAX, ZMAX, ZMIN], cmap="viridis", aspect="auto")
    fig.colorbar(im, ax=ax, label="log10 rays/cell")
    ax.set_title("D  N-S cross-section ray hits (2 km)", fontsize=9)
    ax.set_xlabel("N (km)"); ax.set_ylabel("depth (km)")

    # E: azimuthal coverage map at 2 km, 10 km depth
    ax = fig.add_subplot(2, 3, 6)
    iz = int((10-ZMIN)/dx)
    im = ax.imshow(aF[:, :, iz].T, origin="lower", extent=[XMIN, XMAX, YMIN, YMAX],
                   cmap="magma", vmin=0, vmax=8, aspect="equal")
    fig.colorbar(im, ax=ax, label="azimuth sectors (of 8)")
    ax.set_title("E  Ray-azimuth diversity @ 10 km, 2 km grid", fontsize=9)
    ax.set_xlabel("E (km)"); ax.set_ylabel("N (km)")

    fig.suptitle("Tomographic resolution feasibility from the EQTransformer catalogue "
                 "(straight-ray coverage proxy)", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = f"{ROOT}/figures/tomo_resolution.png"
    plt.savefig(out, dpi=130)
    print(f"\nthresholds: hit>= {HIT_MIN} rays AND >= {AZ_MIN}/8 azimuth sectors per cell")
    print("straight-ray proxy; a definitive answer needs a LOTOS checkerboard run.")
    print(f"wrote {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
