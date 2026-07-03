---
title: "Yogya 2006 Temporary Aftershock Network (station table)"
type: entity
status: active
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/yogya-2006-aftershock-edl-dataset]]"
tags:
  - network
  - stations
  - yogya-2006
---

## Summary

12 temporary EDL stations recording aftershocks of the 27 May 2006 Yogyakarta
earthquake. Codes assigned for processing: `TFxx` (from folder `tf30xx`), network
placeholder `YK`. Coordinates = **median of Yogya-region GPS fixes** parsed from
`.msg`/`.gps` (see `eqt/scripts/extract_coords.py`, output
`eqt/config/stations_raw.json`).

## Station coordinates (median GPS)

| Code | Folder | lat | lon | elev (m) | GPS spread | note |
|------|--------|-----|-----|----------|-----------|------|
| TF07 | tf3007 | -8.11676 | 110.56059 | 118 | ~36 km | RELOCATED; has DE test data |
| TF09 | tf3009 | -7.95466 | 110.36673 | 46  | ~9 km  | check for relocation |
| TF10 | tf3010 | -7.92250 | 110.39363 | 71  | ~19 km | RELOCATED |
| TF11 | tf3011 | -7.77729 | 110.64065 | 153 | ~2 km  | ~stable |
| TF12 | tf3012 | -7.88689 | 110.43568 | 120 | ~0.14 km | stable; has 05 test data |
| TF13 | tf3013 | -7.83626 | 110.30570 | 123 | ~0.37 km | stable |
| TF14 | tf3014 | -8.02800 | 110.34204 | 51  | ~0.05 km | stable |
| TF15 | tf3015 | -7.75104 | 110.49077 | 173 | ~13 km | RELOCATED |
| TF16 | tf3016 | -7.91239 | 110.52202 | 212 | ~0.06 km | stable |
| TF17 | tf3017 | -7.95410 | 110.56618 | 201 | ~28 km | RELOCATED; has DE test data |
| TF18 | tf3018 | -7.91790 | 110.19820 | 41  | ~0.05 km | stable |
| TF19 | tf3019 | -7.99533 | 110.48706 | 177 | ~0.1 km | stable; has 04 test data |

- Spread = min-max range of retained fixes (×111 km); large spread ⇒ likely
  mid-campaign relocation → needs per-period coordinates.
- Pilot uses the stable central stations **TF12, TF14, TF16, TF18, TF19**.

## Links

- [[wiki/sources/yogya-2006-aftershock-edl-dataset]]
- [[wiki/syntheses/eqtransformer-yogya-2006-run]]
