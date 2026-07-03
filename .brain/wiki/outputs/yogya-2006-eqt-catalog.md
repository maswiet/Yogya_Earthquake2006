---
title: "Yogya 2006 EQTransformer Aftershock Catalog (v1)"
type: output
status: active
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/yogya-2006-aftershock-edl-dataset]]"
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
  - "[[wiki/sources/eqtransformer-github-repo]]"
tags:
  - catalog
  - deliverable
  - yogya-2006
---

## Summary

First EQTransformer-derived aftershock catalog for the 2006 Yogyakarta sequence,
produced from the temporary 12-station EDL array
([[wiki/sources/yogya-2006-aftershock-edl-dataset]]) using the `original`
(STEAD) EQTransformer weights via SeisBench, associated + located with PyOcto and
a Central Java 1-D velocity model.

## Result (v1)

- **Picks:** 397,292 (197,159 P / 200,133 S) over 960 station-days (days 152–242).
- **Associated events:** 17,150 (≥8 picks incl. ≥3 P+S; ⇒ all ≥4 stations).
- **Final catalog (artifact-filtered):** **16,876 located aftershocks**
  (within 0.4° of array, depth < 35 km). Removed 274 domain-edge / depth-pinned
  spurious associations. 12,534 events use ≥5 stations.
- **Span:** 2006-06-03 → 2006-08-29.
- **Depth:** median 12.9 km (10–90%: 3.5–16.0 km) — crustal.
- **Rate:** mean 192/day, peak 521/day on 2006-06-17; clear Omori-like decay.
- **Space:** tight cluster filling the array (~−7.9, 110.44), elongated NE–SW
  (consistent with the Opak-fault zone).

## Files

- `eqt/full/catalog_eqt.csv` — final catalog (time, lat, lon, depth, n_picks, n_stations).
- `eqt/full/events_full.csv` — all 17,150 associated events (pre-filter).
- `eqt/full/picks_full.csv` — all 397k phase picks.
- `eqt/full/catalog_summary.png` — 4-panel summary figure.
- Scripts: `eqt/scripts/{extract_coords,reloc_coords,preprocess,run_full,associate_full,finalize_catalog}.py`.

## NonLinLoc relocation (v2, 2026-07-04)

All 16,876 events relocated with **NonLinLoc** (built from source, arm64), 3D
travel-time grids from the Central Java 1-D model, oct-tree search, EDT_OT_WT,
Vp/Vs=1.75 (S from P grids). Ran in ~8 min.

- **Quality:** RMS median **0.079 s**; location errors median **errH 1.1 km,
  errZ 1.4 km**; gap median 136°.
- **Well-constrained subset** (gap<180°, errH<5 km, RMS<0.5): **12,844 events**.
- **Depth:** median **9.6 km** (shallower than PyOcto's 12.9).
- **Cross-validation:** median epicenter shift NLLoc↔PyOcto = **1.4 km**
  (depth −1.5 km) — two independent methods agree.
- **Structure:** sharper than PyOcto — a clear **NE–SW fault lineament** with a
  depth gradient (shallow SW → deeper NE), consistent with the debated
  Opak-fault geometry.
- Files: `eqt/full/catalog_nll.csv` (all), `eqt/full/catalog_nll_good.csv`
  (well-constrained), `eqt/full/nll_compare.png`. Pipeline: `eqt/nll/` control
  files + `eqt/scripts/{gen_nll,parse_nll}.py`; NLL build in `eqt/tools/`.

## Caveats / next steps

- NLLoc locations use the default Central Java 1-D model; swapping in the
  reference-study model would sharpen absolute depths. HypoDD double-difference
  relocation is the natural next refinement for fine fault structure.
- **No magnitudes:** raw miniSEED lacks instrument response → no calibrated ML or
  magnitude of completeness. Needs sensor gain/response to add.
- **Comparison pending:** need a reference Yogya 2006 aftershock catalog to make
  the Mousavi/Tottori-style "N× more events" statement quantitative
  ([[wiki/claims/eqtransformer-doubles-detections-tottori]]).

## Links

- [[wiki/syntheses/eqtransformer-yogya-2006-run]]
- [[wiki/entities/yogya-2006-temp-aftershock-network]]
- [[wiki/questions/applying-eqtransformer-to-yogya-2006]]
