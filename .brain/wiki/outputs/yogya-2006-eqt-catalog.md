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

### GrowClust bootstrap uncertainties (nboot=100, 2026-07-08)

Re-ran GrowClust with **nboot=100** bootstrap resamples of the differential-time
data → relative-location uncertainties for 13,251 clustered events:
- **horizontal eh: median 274 m** (90th pct 426 m)
- **vertical ez: median 305 m** (90th pct 471 m)
- **origin-time et: median 47 ms** (90th pct 86 ms)

Sub-500 m relative precision confirms the sharp fault plane is well-resolved.
`catalog_growclust.csv` now carries eh/ez/et; figure
`figures/growclust_uncertainty.png` (`plot_gc_uncertainty.py`). Bootstrap run
took ~70 min (100 re-clusterings of 368k pairs).

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

### Geology-overlay map (2026-07-08) — real lithology

User provided the **Rahardjo et al. (1995) geology shapefile** (WGS84;
`eqt/data/geology_yogyakarta/`). `plot_geology_map.py` clips it to the region and
groups the 24 local map units into the paper's 8 age-classes (e.g. Tmwl=Wonosari
limestone, Tmng/Tms=Semilir/Nglanggeran Miocene volcaniclastics, Teon=Eocene
Nanggulan, KTm1=pre-Tertiary metamorphics, Tpdi=Pendul diorite), overlain on
hillshade with faults + GrowClust events + legend.
`figures/growclust_geology_map.png`.
**Finding:** aftershocks concentrate in the **Lower–Mid Miocene volcaniclastics**
along the boundary with the **Southern Mountains limestone**; shallow (W, near
GFZ/Bantul) → deep (NE).

**Update (2026-07-08):** removed the hand-digitized fault traces (their positions
were unreliable) and added the **GCMT focal mechanism** beachball — 2006 Yogya
Mw 6.4, NP1 323/77/−176, **NP2 232/86/−13** (near-vertical **left-lateral
strike-slip**; Ekström et al. 2012). The NE–SW nodal plane (strike 232≈N52E)
matches the N57E aftershock lineament. Proper fault traces need a verified
fault GIS layer (deferred).

**Layout update (2026-07-08):** moved the geology legend and the title OUTSIDE the
map frame (legend below, 2-col; title above) so they no longer cover hypocentres;
added **three agency beachballs** — GFZ (red ≈229/85/−9), USGS (gold 231/87/3),
GCMT (blue 232/86/−13) — each colour-matched to its epicentre star and offset with
a connector. All three are near-vertical left-lateral strike-slip and mutually
consistent. Scale bar moved to top-left.

## Caveats / next steps

- NLLoc locations use the default Central Java 1-D model; swapping in the
  reference-study model would sharpen absolute depths. HypoDD double-difference
  relocation is the natural next refinement for fine fault structure.
- ~~No magnitudes~~ → **RESOLVED (2026-07-04):** instrument response provided
  (XN network, L4-3D geophone, 1.7e8 counts/(m/s)). **ML computed** for all events
  via Wood-Anderson simulation — see Magnitudes below.

## Magnitudes (ML) — 2026-07-04, **substantially corrected 2026-07-20**

- Response: L4-3D 1 Hz geophone (PAZ f0=1, h=0.707, sens 1.7e8 counts/(m/s)).
  Per pick: remove response → Wood-Anderson simulate → **bandpass 1–20 Hz** →
  peak |A| in S-window →
  ML = log10(A_mm) + 1.110·log10(R/100) + 0.00189·(R−100) + 3.0 (Hutton & Boore 1987).
  Per event = median over ≥3 stations, **minus a per-station ML correction**.
  97,691 amplitude readings.

### ⚠️ Three measurement bugs found & fixed 2026-07-20 (pre-submission)

Triggered by plotting example waveforms of the ML<−1 "tail": the smallest events
turned out to be strong earthquakes (SNR 100–1000), not noise. Root causes:

1. **`obspy.simulate()` taper on 24-h traces.** Deconvolution was run on the
   merged day-long trace; the default 5% cosine taper = 72 min at each end, so
   events near a day boundary were suppressed up to ~200× (**+2.3 ML**). 39% of
   ML<−1 events fell within ±10 min of midnight (28× enrichment, confined to
   that band). Fix: deconvolve a padded per-pick window (`wa_window`, ±60 s).
   Midnight enrichment for ML<−1 → **0.0%** after fix.
2. **50 Hz mains hum** (esp. TF10b, 200 Hz sampling). Peak was measured on an
   unfiltered trace, so hum set the amplitude of weak events (+0.27 ML bias at
   ML<0 for TF10b). Fix: 1–20 Hz bandpass before the peak. TF10b weak-event
   amplitude drops ~15×; strong events untouched (×1.1).
3. **Max-curvature underestimates Mc by ~0.6.** MaxC gives Mc=−0.20 but b then
   climbs monotonically with the cut-off (0.77→0.91) = residual incompleteness.
   Replaced with **b-value stability (MBS, Woessner & Wiemer 2005)**: `mbs_mc()`.

- **Also applied:** per-station ML corrections (`station_ml_corrections.py`,
  alternating-median solve, datum = amplitude-weighted mean zero). Range −0.37
  (TF16) … +0.35 (TF17); tightens per-event station scatter 0.185 → 0.128 (raw),
  final per-event ML_std median **0.230**. **r = −0.03 with VELEST P travel-time
  corrections** — amplitude and travel-time site response are decoupled here.

### Final magnitude results (2026-07-20)

- **16,876 events with ML.** Range **−1.81 … 3.55**, median **−0.22**.
- **Mc = +0.50** (b-stability), **b = 0.89 ± 0.02** (N = 2,258). MaxC value
  (Mc=−0.20, b=0.75) retained for reference only. Old value **Mc≈0.0 / b=0.88
  is superseded and was contaminated by the taper bug + wrong Mc estimator.**
- The single-station ambient noise floor at the median nearest-station distance
  (11.9 km) is ML **−0.11**; network Mc sits **0.6 units above** it — the cost of
  requiring detection at enough stations to associate + locate (≥8 phases, ≥3
  amplitudes). A single-station floor is a lower bound on Mc, not an estimate.
- Sub-Mc events are **real, not spurious**: example ML −1.01 event has SNR 12–34.
  The sample is incomplete below Mc; the detections are genuine.
- **Omori decay:** fitted **p ≈ 1.05**; peak 382/day on Jun 17.

### Quality screen & false-detection rate (`screen_catalog.py`)

- **11,790 / 16,876 pass (69.9%).** Scale-free criteria: ≤25% phases badly
  fitted (|res|>0.5 s), RMS≤0.5 s, ≥8 phases, ≥2 S, gap≤180°.
- **Do NOT threshold on max|res|**: it scales with phase count (median nphs
  10→21, max|res| 0.61→1.41 from ML<0 to ML>1), so it rejects the best-recorded
  events and inverts the pass-rate vs magnitude. Scale-free screen gives the
  physically expected rising pass rate (65%→78%→73%). Dominant reject reason is
  **gap>180° (4,299)** — network geometry, not detection quality.

### Ramdhan et al. 2025 comparison — ⚠️ claim not yet defensible

Ramdhan et al. (Nat. Hazards 121, 2025, s11069-025-07440-8) relocated
**2,141 events**. Our size ratio depends entirely on the completeness threshold:

| threshold | our events (quality-passed) | ratio vs 2,141 |
|---|---|---|
| all detections | 11,790 | 5.51× |
| ML ≥ 0.0 | ~4,600 | 2.15× |
| **ML ≥ +0.50 (our Mc)** | **~1,770** | **0.83×** |

At a defensible common completeness the catalogues are **comparable, not
7.8× larger**. **Blocking need:** Ramdhan's magnitude-of-completeness / FMD —
without it no "N× more events" statement can be made. This is now the critical
path to publication, not an optional extra.
([[wiki/claims/eqtransformer-doubles-detections-tottori]])

- Figures: `magnitude_gutenberg_richter.png`, `aftershock_magnitude_map.png`,
  `aftershock_rate_decay.png`, `detection_noise_floor.png` (4-panel noise-floor /
  FMD / nearest-station diagnostic), `example_waveforms.png` (P/S + SNR for 3
  ML bands), `station_ml_corrections.png`.
- Files: `eqt/full/{catalog_magnitude,amplitudes,catalog_quality}.csv`,
  `eqt/config/station_ml_corrections.json`; scripts
  `eqt/scripts/{build_amplitudes,compute_magnitudes,screen_catalog,station_ml_corrections,plot_noise_floor,plot_example_waveforms}.py`.
- Caveat: Hutton-Boore (S. California) distance term used as default (no local ML
  scale for Java).

### FMD roll-off above ML~1.5 is REAL, not clipping (resolved 2026-07-20)

`check_clipping.py` read raw counts for the 12–15 largest events at their nearest
stations. Digitiser rail = **2^23 = 8,388,608 counts**; the largest event
(ML 3.55) peaks at 1.37e6 = **16% of the rail** (median of the 12 largest: 8%),
with **no flat-topping** (≤2 samples within 1% of peak per S window). Widening the
1–20 Hz measurement band to 0.3–20 or 1–45 Hz changes large-event amplitude by
**<3%** — so no low-frequency content is lost either. **The events do not clip;
max ML 3.55 is a real upper value, not an instrument-limited lower bound.** The
GR deficit (obs/GR-pred 1.0 at Mc → 0.09 at ML 3.0) is a finite-catalogue /
genuine-maximum-magnitude effect. Figure `clipping_check.png`. *(Retracts the
earlier "L4-3D clipping ⇒ max ML is a lower bound" caveat.)*
- Minor: at the 20 Hz upper corner the SMALLEST events (ML~−0.4) lose 47–65%
  amplitude vs a 1–45 Hz band (their corner freq > 20 Hz) — a mild low bias on
  the smallest magnitudes; 20 Hz kept as the standard (noise trade-off).

### TF16 −0.37 deficit is a site/path effect, not an instrument fault (resolved 2026-07-20)

`diagnose_tf16.py`: the deficit **grows with distance** (−0.21 at <10 km →
−0.48 at 20–30 km; a gain error would be flat), is **stable in time** (no drift,
−0.32…−0.49 across weeks 22–34; not a failing sensor), and both horizontals are
healthy (**N/E amplitude ratio 0.89** over 1,696 readings; not a one-channel
wiring/gain fault). ⇒ a **genuine site/path response** — hard limestone east of
the Opak fault, de-amplifying 1–20 Hz S — fully absorbed by the station ML
correction. Explains the r=−0.03 vs VELEST P corrections: amplitude (near-surface
kappa) and P travel-time (deeper velocity) sample different depths. Figure
`tf16_diagnostic.png`. *(Retracts the earlier "possible instrument-gain issue".)*

## Links

- [[wiki/syntheses/eqtransformer-yogya-2006-run]]
- [[wiki/entities/yogya-2006-temp-aftershock-network]]
- [[wiki/questions/applying-eqtransformer-to-yogya-2006]]
