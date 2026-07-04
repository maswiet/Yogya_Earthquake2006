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

## Next steps

- Optionally re-run NLLoc/HypoDD with this data-driven model + station corrections
  for improved absolute depths.
- Deeper layers (>12 km) need events/quarry blasts at larger distance to resolve.

## Links

- [[wiki/outputs/yogya-2006-eqt-catalog]]
- [[wiki/syntheses/eqtransformer-yogya-2006-run]]
