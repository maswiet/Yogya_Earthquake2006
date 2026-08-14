---
title: "Diambama et al. 2019 — Yogyakarta 2006 tomography (GJI) — OUR GROUP'S prior work"
type: source
status: active
created: 2026-07-20
updated: 2026-07-20
tags:
  - source
  - yogya-2006
  - ngalang-fault
  - priority
  - our-group
---

## Citation

Diambama, A.D., Anggraini, A., Nukman, M., Lühr, B.-G. & Suryanto, W. (2019).
*Velocity structure of the earthquake zone of the M6.3 Yogyakarta earthquake
2006 from a seismic tomography study.* **Geophysical Journal International**
216, 439–452. doi:10.1093/gji/ggy430. Corresponding author: W. Suryanto
(ws@ugm.ac.id), UGM Seismology Research Group.

**This is our own group's prior publication.** File:
`/Users/maswiet/Documents/Manuscript/GJI_Wiwit_.pdf`.

## Why it matters — PRIORITY on the Ngalang Fault

**Our group was the first to associate the 2006 Yogyakarta earthquake source with
the Ngalang Fault, in 2019 — six years before Ramdhan et al. (2025).** Verbatim:

- §5.2: "We assume that this unnamed fault trend is possibly associated with the
  extension of the NE–SW fault zone at Ngalang River (namely Ngalang Fault) as
  shown in Geological Map of Surakarta–Giritontro Quadrangle of Surono et al.
  (1992)."
- Conclusion: "...structural features of NE–SW located to the east of the Opak
  Fault... possibly related to the extension of the fault zone at the Ngalang
  River."
- Fig. 14 caption: "Inferred of the extension of Ngalang Fault."

Ramdhan et al. (2025) later named the Ngalang Fault as the mainshock source from
2,141 relocated events; our group had already proposed that association in 2019.

## Same dataset lineage as the EQTransformer catalogue

- **Same GFZ/German Task Force temporary deployment** (with BMKG + UGM): Mark
  **L4-3D 1 Hz** sensors + **Earth Data Logger (EDL)**, set up from 2006-05-31.
  The tomography used **10 stations** (our EQTransformer catalogue uses 12 from
  the same array).
- Aftershocks from **Anggraini (2013)** PhD (Univ. Potsdam) — manually picked and
  relocated. Aftershocks located **10–15 km east of the Opak Fault**, max depth
  ~20 km, on a NE–SW trend.
- **Manual picks available: 3,769 P + 3,407 S phases; 588 events** used in the
  tomography (single-event location errors ±800 m lateral / ±1200 m vertical).
- Catalogue magnitude range **ML 0.02–3.55** — the **max ML 3.55 exactly matches
  our independent EQTransformer+Wood-Anderson catalogue**, a consistency check on
  our magnitude scale at the high end (and independent support that the largest
  events are not clipped — see [[wiki/outputs/yogya-2006-eqt-catalog]]).

## Fault interpretation (their result)

- Tomography (LOTOS-13, P+S, Vp/Vs) images an "unnamed fault" delineated by a
  velocity anomaly to **5–7 km depth**, interpreted as a **strike-slip fault with
  a reverse component, dipping east, striking NE–SW** — i.e. the Ngalang Fault
  extension.
- Severe damage correlates with **high Vp/Vs** zones = unconsolidated sediment
  (agrees with Walter et al. 2008).
- Uplifted hanging wall on the **west** side of a reverse fault; graben structure
  beneath Yogyakarta = very-low-velocity filled sediment.

## Two things this unlocks for the new paper

1. **Priority narrative.** The EQTransformer study is the CONTINUATION of our
   group's own program: 2019 tomography first delineated + named the Ngalang
   association (resolved only to 5–7 km); the new 13,251-event GrowClust
   catalogue images the SAME fault to sub-500 m over 0–18 km — a genuine
   resolution advance over our own prior work, not a competition with Ramdhan.
2. **Pick validation ground truth.** Anggraini (2013) manual P/S picks on the
   SAME raw archive are the natural benchmark for EQTransformer precision/recall
   + residuals — the top remaining pre-submission task
   ([[wiki/outputs/yogya-2006-srl-paper-plan]]). Need to obtain the Anggraini
   2013 pick files / catalogue from the group.

## Anggraini catalogue obtained + compared (2026-07-20)

`Bantul2006_Aftershock_Catalogue.xlsx` (2 sheets: Absolute_location 590 events;
DD_location 524 HypoDD-relocated). Columns: date/month/year, lon, lat, depth, ML.
**No sub-day origin time and no phase picks** → event-by-event pairing impossible
(nearest-neighbour ML correlation r=0.06 = matches co-located but different
events). Covers **only 3–7 June 2006** (5 days, ~118/day); our catalogue spans
3 Jun–29 Aug, so the comparison window is 3–7 Jun. Script
`eqt/scripts/compare_anggraini.py`, figure `figures/compare_anggraini.png`.

**Findings (distributional):**
- **Detection recovery strong:** 576/590 manual events have a co-located
  same-day counterpart; we add ~324 more (914 QC events vs 590). Spatial and
  depth distributions agree (both delineate the NE–SW Ngalang structure). This
  is our event-level validation in lieu of pick-level precision/recall.
- **Systematic ML offset (now pinned by event-matched picks 2026-08-14): +0.41,
  r=0.95** (see the pick-validation section below). Our absolute ML reads 0.41
  low vs the manual scale; b-value is offset-invariant (b=0.89 stands). The
  earlier distributional estimate here (+0.86 matched-location / +1.25 FMD-tail)
  is SUPERSEDED — it was inflated by our catalogue's much greater small-event
  completeness. Tie to the local scale by +0.41.

## Pick-level validation done (2026-08-14) — the reviewer-critical result

The group provided Anggraini's manual PICK files
([[raw/2026-08-14_anggraini-picks/provenance]]): `phase_300.dat` (588 events,
3776 P + 3414 S, with full h:m:s origin times), `station.dat`, `stat_ft.dat`.
Script `eqt/scripts/validate_picks.py`, figure `figures/pick_validation.png`.
Manual station codes mapped to our TF codes by June-3-7 coordinates
(WON=TF12, PEL=TF13, RAT=TF14, WAN=TF16, PAL=TF19, BUM=TF18 are the exact
co-located core; YOG=TF07, PRA=TF10 movers).

**Results — EQTransformer passes cleanly:**
- **Event recovery: 528/588 (90%)** matched to a manual event within ±5 s
  origin time (independent of location — proper temporal match, not spatial).
- **Pick precision (per-event demeaned): P MAD 0.02 s, S MAD 0.06 s**; 98% of P
  and 83% of S within 0.3 s of the manual pick. 4,761 matched picks across the 6
  co-located core stations, all with ~0.00 median (no mis-mapping). S has a small
  positive tail (our S slightly late — normal for automatic S).
- **ML: event-matched offset +0.41 (mean +0.43, std 0.19), r = 0.95.** This
  SUPERSEDES the earlier distributional estimate (0.86–1.25), which was inflated
  by completeness differences (our catalogue has many more small events). The
  true tie to the local scale is **+0.41**; b-value is offset-invariant so
  b=0.89 stands. Clean constant offset → easy calibration, not a distortion.

**Notes / loose ends:**
- Manual stations **NGL and TRI** (heavily used: 964 + 528 picks) do NOT
  correspond to our processed XN stations and contributed no matched picks —
  Anggraini's network included stations we did not process. Worth confirming
  whether we should add them.
- 60/588 (10%) manual events unmatched — candidates we missed or located as
  time-outliers; worth a spot check.

## Links

- [[wiki/outputs/yogya-2006-eqt-catalog]]
- [[wiki/outputs/yogya-2006-srl-paper-plan]]
- [[wiki/outputs/yogya-2006-1d-velocity-model]]
- [[wiki/entities/yogya-2006-temp-aftershock-network]]
- [[wiki/entities/lotos-tomography-code]] — the LOTOS-13 engine used for this tomography
- [[wiki/sources/zenodo-lotos-koulakov-2021]] — citable LOTOS code release
