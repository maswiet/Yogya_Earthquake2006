---
title: "Yogya 2006 data-driven 1-D velocity model (Wadati + VELEST)"
type: output
status: active
created: 2026-07-04
updated: 2026-07-04
sources:
  - "[[wiki/sources/yogya-2006-aftershock-edl-dataset]]"
tags:
  - velocity-model
  - velest
  - deliverable
---

## Summary

1-D P/S velocity structure derived **from the EQTransformer arrival times** by the
travel-time method, culminating in a VELEST minimum-1D simultaneous inversion
(velocities + hypocenters + station corrections). Replaces the assumed Central
Java model with a data-driven one.

## Direct travel-time method (195,351 phases from 16,876 NLLoc events)

- **Wadati diagram** (76,008 P-S pairs): **Vp/Vs = 1.735**.
- **Travel-time vs hypocentral distance**: apparent **Vp = 5.93 km/s**,
  **Vs = 3.39 km/s** (ratio 1.75). Slight far-offset curvature ⇒ layered crust.
- Script: `eqt/scripts/parse_arrivals.py` → `eqt/full/ttime_analysis.png`,
  `eqt/full/arrivals.csv`.

## VELEST minimum-1D inversion

- Built VELEST from source (REAL repo, gfortran, arm64). Inputs from 800
  well-located events (gap<160°, RMS<0.30, ≥10 phases). Per-period station codes
  used as separate sites. Vp/Vs fixed 1.735, station corrections on, deep layers
  (≥20 km) damped (poor ray coverage).
- **Convergence:** RMS **0.358 → 0.152 s** over ~22 iterations.
- **Inverted P model** (depth km : Vp km/s): 0:4.94, 4:5.74, 8:6.14, then rising
  to ~7.2 by 12 km. **Vs**: 0:2.77, 4:3.31, 8:3.69. Vp/Vs ≈ 1.73–1.78.
  Well-resolved 0–8 km; ≥12 km poorly constrained (events mostly 5–12 km) → fixed.
- **Station P-corrections** (relative): +1.47 s TF18, +1.34 TF10a, +1.30 TF15b,
  +1.20 TF09a … down to −0.55 TF11b. Pattern matches geology — large positive =
  slow near-surface (Bantul basin, west); negative = fast bedrock (eastern hills).
- Files: `eqt/velest/` (yogya.cmn/.sta/.pha/.mod, yogya.OUT, yogya.finalcnv,
  velest_result.png). Scripts `eqt/scripts/{gen_velest,parse_velest}.py`.

## Coordinate gotcha (documented)

VELEST uses **West-positive longitude**: the origin `olon` must be entered as
**negative** for East (−110.44 for 110.44°E), else station x-coordinates overflow.

## Refined, geologically-constrained model (2026-07-04)

Rebuilt with a priori geology (fixing BOTH shallow and deep layers; inverting only
the well-sampled mid-crust), which also fixes the earlier deep artifact:

- **Shallow FIXED (damp 999):** Merapi sediment — 2.5 km/s (−2 km), 2.9 (0 km),
  4.3 at the **~0.7 km sediment base** (user's constraint). Lateral W/E
  sediment-vs-limestone contrast is carried by station corrections + elevation.
- **Mid-crust INVERTED:** 4.65 (2 km) → 5.49 (4–7 km) → 6.30 (10 km) → 6.39 (13 km).
- **Deep FIXED (damp 999):** 6.55 (16 km), 6.80 (22 km), 7.20 (30 km), 8.0 (40 km)
  — realistic crust, **removes the earlier 7.16 km/s @ 12 km artifact**.
- Vp/Vs ≈ 1.73. Elevation corrections ON (real station elevations 44–260 m).
  Converged RMS **0.142 s** (slightly better than the free-deep run's 0.152).

### Station corrections ↔ Opak-fault geology

`eqt/scripts/interpret_stacorr.py` → `eqt/figures/station_corrections_opak.png`.
Clear at the extremes: **far-west BUM/TF09a/PEL slow (+0.5…+0.8 s = Merapi
sediment)**; **far-east KRI/TF11/KEM fast (−0.6…−1.25 s = limestone)**. Mean W−E
contrast ≈ +0.34 s (west slower), consistent with thick (>700 m) western sediment.
Mid-zone stations near the fault are mixed (transition + relocated-site sites).

## Relocation with the VELEST model (2026-07-04)

Rebuilt NLLoc P & S travel-time grids from the VELEST model and applied the
VELEST station corrections via `LOCDELAY`; relocated all 16,876 events
(`eqt/nll/loc_v/`, `eqt/full/catalog_velest.csv`, scripts `gen_velest_reloc.py`,
`compare_reloc.py`).

- RMS median **0.132 s**; errH median 2.3 km; 8,919 well-constrained.
- **Depths deepen by median +5.4 km** vs the Central-Java model (median 9.4 → 15.0
  km); epicentres shift ~3.5 km.
- **Interpretation / caveat:** the deepening is an **artifact of the VELEST deep
  layers** (Vp≈7.16 km/s already at 12 km — unrealistically fast, and
  poorly constrained since events are mostly 5–12 km deep). A too-fast deep medium
  biases relocated depths downward. The **shallow structure (0–8 km) is reliable;
  the Central-Java-model depths (~9–13 km) are more physically credible** for the
  aftershock zone.
- **Recommended fix:** re-invert VELEST with the deep layers (≥12 km) fixed to
  ~6.3–6.6 km/s (or add larger-offset/blast data), then re-relocate.

### Re-relocation with the REFINED model (2026-07-04) — bias fixed

Rebuilt grids from the refined (fixed-shallow + fixed-deep) model + elevation-aware
station corrections; relocated all 16,876 events. **Depth bias removed:**

| Model | depth median | well-constrained | RMS | errH |
|-------|-------------|------------------|-----|------|
| VELEST v1 (7.16 artifact) | 15.0 km | 8,919 | 0.132 | 2.3 km |
| **VELEST refined** | **10.5 km** | **12,447** | **0.113** | **1.4 km** |
| Central-Java (ref) | 9.4 km | 12,844 | 0.079 | 1.1 km |

Depth vs Central-Java now only +1.0 km (was +5.4). `catalog_velest.csv` = refined;
`catalog_velest_v1.csv` = old biased. Relief map `aftershock_relief_map.png` and
`velest_relocation_compare.png` regenerated with refined depths.

## Next steps

- Refine the VELEST deep layers and re-relocate for unbiased depths.
- HypoDD double-difference relocation for fine fault structure.

## Links

- [[wiki/outputs/yogya-2006-eqt-catalog]]
- [[wiki/syntheses/eqtransformer-yogya-2006-run]]
