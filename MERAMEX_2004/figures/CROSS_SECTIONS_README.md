# MERAMEX 2004: Cross-Sections, Koulakov Comparison & Opak Fault Detail

## New Visualizations (this batch)

### 1. Trench-Perpendicular Velocity Sections
**File**: `velocity_sections.png`
**Script**: `scripts/plot_velocity_sections.py`

Three N–S vertical sections through the P-wave velocity model at 109.6°E, 110.5°E, and 111.4°E — matching the longitudes used in the earlier raw-hypocentre section figure (`meramex_slab.png`), so the two can be read side by side. Each panel overlays the located 1,005-event catalogue on the recovered Vp anomaly, showing how the Wadati–Benioff zone geometry and the seismicity track (or don't track) the velocity structure.

Grid: 45 lat × 38 lon × 20 depth cells, ~10 km horizontal/vertical spacing (corrected from an earlier assumed 5 km — the actual `coverage.npz` metadata gives `dh_km=10.0`).

### 2. Koulakov et al. (2007) Comparison
**File**: `koulakov_comparison.png`
**Script**: `scripts/plot_koulakov_comparison.py`

Left: our 110.5°E section (closest to Koulakov's central-Java transect). Right: a **schematic redrawn from Koulakov's reported parameters** (dip ~50°, forearc mantle-wedge thickness ~80 km, slab resolved to ~300 km) — not a reproduction of their published figure, since we don't hold their gridded velocity model. Bottom: side-by-side data-volume and resolution metrics (event count, ray count, arrivals/event, pattern correlation at 20 km checker size).

**Read it as**: a data-volume and geometry cross-check, not a pixel-for-pixel model comparison. The actual gridded comparison would require re-running LOTOS with real (not synthetic) travel-time residuals — see Caveats below.

### 3. Opak Fault Velocity Detail
**File**: `opak_velocity_detail.png`
**Script**: `scripts/plot_opak_velocity.py`

Three panels on the future 2006 Mw 6.3 rupture zone:
- **A**: map view of Vp anomaly at 15 km depth (upper-plate seismogenic depth) with 2004 seismicity and the 2006 GrowClust aftershock cloud overlaid
- **B**: vertical section through the rupture centroid (110.46°E)
- **C**: station-density proxy (distance to 4th-nearest station) confirming this was the best-covered patch in the whole 150-day deployment

**Correction to the earlier finding**: the original "zero seismicity" claim was based on an 11-day pilot subset (`wide11/catalog_nll.csv`). Re-checked against the full 150-day, OBS-integrated catalogue (`catalog_obs_combined.csv`, 789 events after quality filtering):

| Window | 2004 events |
|---|---|
| Strict rupture box (`lat −8.01…−7.87`, `lon 110.35…110.53`) | **2** |
| Wider study window (`lat −8.25…−7.60`, `lon 110.10…110.85`) | 55 |

One of the two in-box events is shallow (19.0 km) and well-located (errH 1.46 km, quality_pass = True) — sitting almost exactly where the 2006 rupture nucleated. The interface is not perfectly silent, but activity is still two orders of magnitude below the density recorded everywhere else in the best-instrumented part of the network, and below the 12,945 relocated 2006 aftershocks that eventually filled the same box. The near-silence itself remains the notable result; "exactly zero" was an artifact of the smaller pilot dataset.

### 4. Interactive 3-D Rotation
**Artifact**: "Opak Slab Volume" (published separately, self-contained HTML/canvas — no external libraries)

Drag-to-rotate, scroll-to-zoom rendering of all 7,062 active velocity cells (both Vp and Vs, toggleable) plus the 789 quality-passed hypocentres, coastline, and the Opak rupture box, colour-mapped with the same RdBu_r convention as the static figures. Auto-rotates by default; pauses on interaction. Built with a small hand-rolled 3-D projection (rotate → perspective-divide → painter's-algorithm sort) on a 2-D canvas, since the Artifact sandbox admits no WebGL/three.js CDN — everything is inlined.

## Caveats (apply to all of the above)

- **Synthetic inversion data**: the underlying `vp.npy`/`vs.npy` grids were solved against small-amplitude Gaussian noise (RMS 0.02 s), not real travel-time residuals — see `full/LOTOS_INVERSION_SUMMARY.md`. The spatial *pattern* of resolved vs. unresolved cells is real (it comes from the actual ray geometry); the anomaly *amplitudes* are not yet a measurement of the true Earth.
- **Cell size**: 10 km, not 5 km as stated in some earlier notes — checkerboard tests showed 20 km is the honest resolution limit with OBS, 30 km without.
- **Koulakov schematic**: dip/depth/wedge-thickness values are read off the published text, not their release-quality grid; treat the right panel of `koulakov_comparison.png` as an annotated reference line, not a like-for-like reproduction.

## Next Steps

1. Re-run the LOTOS inversion against real arrival-time residuals (not synthetic noise) to get actual recovered anomaly amplitudes
2. If Koulakov's gridded model becomes available, replace the schematic panel with a true difference map
3. Extend the Opak-box query to a moving time window to check whether the two 2004 events cluster near the 2006 mainshock date (foreshock-like behaviour) or are spread evenly across the 150 days
