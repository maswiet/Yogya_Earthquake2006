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

## HypoDD double-difference relocation (2026-07-08)

Sharpened the fault with catalog-differential-time double-difference relocation
(HypoDD v2.1b built from source, refined VELEST 1-D model).

- Input: 3,500 high-quality events (gap<150, rms<0.22, ≥12 phases) → ph2dt →
  1.45 M differential times → **3,445 events relocated**.
- **Relative precision:** formal errors **ex≈18 m, ey≈17 m, ez≈20 m** (vs ~1.4 km
  single-event errH) — the fault plane collapses to a sharp, steeply-dipping
  structure; across-strike core width visibly narrows.
- Depth median 9.9 km. Files: `full/catalog_hypodd.csv`, `hypodd/hypoDD.reloc`,
  `figures/hypodd_relocation.png`, `figures/hypodd_relief_map.png` (over topo);
  scripts `gen_hypodd.py`, `plot_hypodd.py`.
- **Temporal analysis** (`figures/temporal_aftershocks.png`, `plot_temporal.py`):
  map coloured by time + along-strike-vs-time, depth-vs-time, daily-rate/cumulative.
  Whole fault (~N57°E) activated early and decayed (Omori); activity persisted in
  the central zone; depth range (0–18 km) stable over time; no strong along-strike
  migration front.
- Build note: hypoDD's large static arrays stress the dyld shared-cache mapping
  in this sandbox; a smaller build (MAXEVE=4000, MAXDATA=1.5 M) loads reliably.

### Can HypoDD be staged/batched for all ~10k events? (2026-07-08)

- **No valid batch-then-merge** for a single connected cluster: double-difference
  solves *relative* positions within a linked cluster; arbitrary batches lose
  inter-batch links and each gets its own arbitrary datum → merge is inconsistent.
  The correct way to add events is a **single simultaneous** run.
- **3,500 quality events is scientifically sufficient** (standard practice; adding
  poorly-constrained events adds noise, not resolution).
- **Scaling up demonstrably fails on this hardware:** binaries sized for >8k events
  fail to load (dyld shared-cache mapping); at 5.5–7.8k events they load but the
  LSQR **crashes at runtime** (stack/memory). **Reliable ceiling ≈ 3,500–4,000
  events** here.
- **Proper path for the full ~10k+ catalog:** use **GrowClust** (Trugman & Shearer
  2017) — hierarchical relative relocation, far lower memory, built for large
  datasets — or a machine with more RAM / no sandbox. Recommended follow-up.

## GrowClust — full-catalog double-difference relocation (2026-07-08) ✅

GrowClust (built from source, dttrugman/GrowClust) succeeded where HypoDD hit
hardware walls — small binary, allocatable memory, no dyld/array limits.

- Inputs: **13,716 events** (gap<200, ≥6 phases; ph2dt rebuilt to MEV=17000) →
  368k pairs, **3.56 M catalog differential times** (dt.ct → GrowClust xcordata,
  tdif=tt1−tt2). Refined VELEST 1-D model (extended to 40 km), Vp/Vs=1.735.
- **Ran in ~8 s.** **13,251 events (97%) relocated** in 543 clusters; largest
  cluster **6,917 events** (the main fault). Depth median stable (10.6→10.7 km).
- Across-strike section collapses the diffuse ±7 km cloud into a sharp,
  steeply-dipping fault plane for the FULL catalog (≈4× the 3,445-event HypoDD run).
- Files: `full/catalog_growclust.csv`, `growclust/OUT/out.growclust_cat`,
  `figures/growclust_relocation.png`, `figures/growclust_relief_map.png`; scripts
  `gen_growclust.py`, `plot_growclust.py`.
- **Conclusion:** GrowClust is the correct tool for the full catalog here; it
  relocated all detected quality events simultaneously (not batched), confirming
  the earlier answer that batching is invalid but full simultaneous relocation is
  achievable with a memory-efficient method.

### Tectonic/structural map (2026-07-08)

`figures/growclust_tectonic_map.png` (`plot_tectonic_map.py`): GrowClust catalog
(depth-coloured 0–25 km, paper colour scheme) on grayscale hillshade, with the
four faults **Opak / Oyo / Ngalang / Nglipar** + Muria–Progo lineament digitized
from **Ramdhan et al. (2025)** Fig. 1, plus the four mainshock epicentres
(BMKG/GFZ/USGS/GCMT), cities, and XN stations. Aftershocks delineate the
**Opak–Ngalang fault zone** (cloud between Opak and Ngalang), matching the paper's
source interpretation; shallow-SW → deep-NE gradient. Full Rahardjo et al. (1995)
lithology not reproduced (needs source shapefile); only a generalized Southern
Mountains limestone domain is shaded.

## Caveats / next steps

- NLLoc locations use the default Central Java 1-D model; swapping in the
  reference-study model would sharpen absolute depths. HypoDD double-difference
  relocation is the natural next refinement for fine fault structure.
- ~~No magnitudes~~ → **RESOLVED (2026-07-04):** instrument response provided
  (XN network, L4-3D geophone, 1.7e8 counts/(m/s)). **ML computed** for all events
  via Wood-Anderson simulation — see Magnitudes below.

## Magnitudes (ML) — 2026-07-04

- Response: L4-3D 1 Hz geophone (PAZ f0=1, h=0.707, sens 1.7e8 counts/(m/s)).
  Per pick: remove response → Wood-Anderson simulate → peak |A| in S-window →
  ML = log10(A_mm) + 1.110·log10(R/100) + 0.00189·(R−100) + 3.0 (Hutton & Boore 1987).
  Per event = median over ≥3 stations. 97,692 amplitude readings.
- **16,876 events with ML.** (ML recomputed 2026-07-04 with refined-model
  distances: median **−0.07**, max **3.48**, **Mc ≈ 0.0**, **b = 0.88 ± 0.01**.
  Earlier v1-depth values were median 0.07 / b 0.96.)
- **Omori decay:** daily rate ~300 → ~20/day; fitted **p ≈ 1.05** (c poorly
  constrained — deployment began ~6 days post-mainshock). Peak 382/day on Jun 17.
- Figures: `magnitude_gutenberg_richter.png`, `aftershock_magnitude_map.png`
  (size~ML, colour~depth), `aftershock_rate_decay.png`.
- Files: `eqt/full/catalog_magnitude.csv`, `eqt/full/amplitudes.csv`,
  `eqt/figures/magnitude_gutenberg_richter.png`; scripts
  `eqt/scripts/{build_amplitudes,compute_magnitudes}.py`.
- Caveat: Hutton-Boore (S. California) distance term used as default (no local ML
  scale for Java); the short-period L4-3D may slightly clip the largest events.
- **Comparison pending:** need a reference Yogya 2006 aftershock catalog to make
  the Mousavi/Tottori-style "N× more events" statement quantitative
  ([[wiki/claims/eqtransformer-doubles-detections-tottori]]).

## Links

- [[wiki/syntheses/eqtransformer-yogya-2006-run]]
- [[wiki/entities/yogya-2006-temp-aftershock-network]]
- [[wiki/questions/applying-eqtransformer-to-yogya-2006]]
