---
title: "Yogya 2006 — SRL paper plan (methodology + data product framing)"
type: output
status: active
created: 2026-07-20
updated: 2026-07-20
sources:
  - "[[wiki/outputs/yogya-2006-eqt-catalog]]"
  - "[[wiki/outputs/yogya-2006-1d-velocity-model]]"
  - "[[wiki/outputs/yogya-2006-pangandaran-triggering]]"
tags:
  - paper
  - srl
  - deliverable
---

## Decision (2026-07-20)

Ramdhan et al. (2025) magnitude/FMD data is hard to obtain, so the "N× more
events than Ramdhan" framing is dropped. The paper is reframed as a
**methodology + open data product** contribution to *Seismological Research
Letters*, with fault structure as independent corroboration rather than a
head-to-head count. This removes the single blocking dependency.

## One-line contribution

> The first EQTransformer-derived, quality-screened, double-difference-relocated
> aftershock catalogue for the 2006 Yogyakarta (Mw 6.4) earthquake — released as
> an open data product — together with a reproducible local-magnitude workflow
> that documents three measurement pitfalls specific to running deep-learning
> pickers on legacy temporary-array data.

## Why an editor bites (the three legs)

1. **Methods, community-useful.** Groups worldwide are now re-running
   SeisBench/EQTransformer on legacy temporary-array archives. We document three
   ML-measurement traps that each shift magnitudes by ≥0.3 unit and are easy to
   miss: (a) `obspy.simulate()`'s 5% cosine taper applied to day-long traces
   silently suppresses near-day-boundary events by up to ~200× (2.3 ML); (b)
   unfiltered peak amplitudes let mains hum (50 Hz) set the magnitude of weak
   events; (c) maximum-curvature underestimates Mc by ~0.6 here — caught by a
   b-value-stability + noise-floor cross-check. Each is shown with the fix and
   the before/after catalogue effect. **This is the novelty that carries the
   paper** — a cautionary + best-practice note with a real dataset behind it.
2. **Open, characterised data product.** 16,876 detections → 11,790
   quality-screened → 13,251 GrowClust-relocated (sub-500 m relative precision),
   with completeness fully characterised (Mc = +0.50 by b-stability, b = 0.89 ±
   0.02) and a physical detection-floor model. Released on Zenodo with a DOI.
   SRL explicitly values reusable data products.
3. **Independent structural corroboration.** 13,251 precisely relocated events
   (543 clusters; main-fault cluster 6,917 events) independently support the
   Opak–Ngalang source zone that Ramdhan et al. inferred from 2,141 events —
   strike N57°E, matching GCMT NP2 (232/86/−13, near-vertical left-lateral),
   with a shallow-SW → deep-NE gradient, aftershocks concentrated in the
   Lower–Mid Miocene volcaniclastics at the Southern Mountains limestone
   boundary. Corroboration from an independent, ~6× denser catalogue is a
   result in its own right — and needs only their published fault interpretation,
   not their magnitude data.

## Positioning vs Ramdhan et al. (2025)

**Complementary, not competing.** They established the fault interpretation; we
provide (i) an independent, denser, openly-released relocated catalogue that
corroborates it, and (ii) the reproducible magnitude/completeness methodology
their study did not detail. Cite them as the reference interpretation we test.
Explicitly avoid any "more events" claim — state that a like-for-like catalogue
comparison awaits a common completeness threshold (their FMD not available).

## Draft abstract (skeleton, ~200 words)

> The 26 May 2006 Yogyakarta earthquake (Mw 6.4) was followed by a temporary
> deployment of 12 short-period (L4-3D) stations recording continuously for ~3
> months. We reprocess this archive with the EQTransformer deep-learning
> phase picker, associate and locate 16,876 events, and relocate 13,251 by
> the GrowClust double-difference method to a median relative precision of
> 274 m horizontal / 305 m vertical. We derive local magnitudes with a
> reproducible Wood-Anderson workflow and characterise catalogue completeness
> (Mc = +0.50, b = 0.89 ± 0.02) against a physical single-station noise floor,
> showing that network completeness sits ~0.6 magnitude units above the
> single-station limit. In building the magnitude catalogue we identify and
> correct three measurement pitfalls that broadly affect deep-learning
> reprocessing of legacy temporary-array data. The relocated hypocentres
> delineate a near-vertical NE–SW (N57°E) structure consistent with the GCMT
> mechanism and with the Opak–Ngalang source zone, independently corroborating
> earlier interpretations from a ~6× denser catalogue. We release the picks,
> magnitudes, quality flags, relocations, and processing code as an open data
> product.

## Section outline → existing assets

1. **Introduction** — 2006 Yogya Mw 6.4; the Opak/Ngalang/Oyo source debate;
   the case for reprocessing legacy temporary-array data with ML pickers.
2. **Data & network** — XN 12-station L4-3D array; `stations_periods.json`
   (station relocations mid-deployment); data hygiene.
3. **Detection, association, location** — EQTransformer/SeisBench → PyOcto →
   NonLinLoc; data-driven VELEST minimum-1D model + station corrections
   ([[wiki/outputs/yogya-2006-1d-velocity-model]]).
4. **Local magnitudes & completeness** — WA workflow; **the three pitfalls**
   (taper / hum / MaxC-Mc) with before/after; station ML corrections;
   noise-floor vs network-Mc. Figs: `magnitude_gutenberg_richter`,
   `detection_noise_floor`, `example_waveforms`, `station_ml_corrections`,
   `clipping_check`, `tf16_diagnostic`.
5. **Quality screening & the catalogue** — scale-free screen (69.9% pass);
   false-detection discussion; the released product.
6. **Relocation & fault structure** — GrowClust 13,251 events, bootstrap
   uncertainties; N57°E, GCMT consistency; geology overlay. Figs:
   `growclust_relocation`, `growclust_uncertainty`, `growclust_geology_map`,
   `temporal_aftershocks`.
7. **(Optional) Pangandaran non-triggering** — short section or supplement
   ([[wiki/outputs/yogya-2006-pangandaran-triggering]]).
8. **Data & code availability** — Zenodo DOI; GitHub scripts.

## Pre-submission checklist

- [ ] **Pick validation vs manual ground truth** — precision/recall + P/S
  residuals on ~100–200 events. Reviewers WILL ask for an ML-picker paper.
- [ ] Zenodo release (picks, magnitudes, quality flags, relocations, model) + DOI.
- [ ] Local-ML caveat: Hutton-Boore distance term (no Java-specific scale) —
  state and bound the effect.
- [ ] Confirm SRL article type (regular article vs data-product/Electronic
  Seismologist framing).
- [ ] Optional: tapered/modified-GR corner magnitude for the real FMD roll-off.

## Open threads carried forward

- Pick validation vs manual picks is the largest remaining reviewer risk.
- Ramdhan FMD/Mc still desirable if it becomes available (would re-enable a
  quantitative catalogue-size comparison), but no longer blocking.
