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

## Per-period coordinates (relocations)

Reconstructed from per-day GPS (`eqt/scripts/reloc_coords.py` →
`eqt/config/stations_periods.json`). Genuine mid-campaign moves:

| Code | Period (jday) | lat | lon | note |
|------|--------------|-----|-----|------|
| TF07 | 152–160 | -7.79267 | 110.36480 | then moves ~35 km |
| TF07 | 161–198 | -8.11676 | 110.56050 | (station ends ~day 198) |
| TF09 | 153–161 | -7.92593 | 110.28870 | then ~9 km |
| TF09 | 162–235 | -7.95466 | 110.36670 | main site |
| TF10 | 154–163 | -7.74986 | 110.49330 | then ~28 km |
| TF10 | 164–201 | -7.92249 | 110.39360 | main site |
| TF15 | 154–164 | -7.86424 | 110.55880 | then ~15 km |
| TF15 | 165–241 | -7.75104 | 110.49070 | main site |

- **TF11**: ~2 km GPS wobble across day 158/159 — treat as one site.
- **TF17**: single site from day 160 (earlier ~28 km spread was outliers).
- Location must assign each pick the site valid for its day (implemented by
  period-specific station codes, e.g. `TF07a`/`TF07b`, in association).

## Official FDSN metadata (network XN, provided 2026-07-04)

Network **XN** (GFZ:XN2006). Real station codes map to our TF-codes by coordinate
(`eqt/config/stations_xn_meta.txt`):

| XN | our code | lat | lon | start |
|----|----------|-----|-----|-------|
| WON | TF12 | -7.88691 | 110.43568 | 05-30 |
| PEL | TF13 | -7.83628 | 110.30560 | 05-30 |
| RAT | TF14 | -8.02800 | 110.34200 | 05-30 |
| WAN | TF16 | -7.91245 | 110.52201 | 05-30 |
| BUM | TF18 | -7.91793 | 110.19823 | 05-30 |
| PAL | TF19 | -7.99540 | 110.48708 | 05-30 |
| NGL | TF15a | -7.86430 | 110.55864 | 05-30 |
| PRA | TF10a | -7.74983 | 110.49333 | 05-30 |
| IMO | TF10b | -7.92250 | 110.39362 | 05-30 |
| BOG | TF17 | -7.95415 | 110.56620 | 06-01 |
| KRI | TF11b | -7.77720 | 110.64080 | 06-07 |
| KEM | TF07b | -8.11674 | 110.56060 | 06-09 |
| KARA | TF09b | -7.95470 | 110.36675 | 06-10 |

- **Instrument (all stations):** Mark/Sercel **L4-3D** 1 Hz geophone,
  velocity sensitivity **1.7e8 counts/(m/s)**, 100 Hz, channels HHZ/HHN/HHE.
  → enables **ML** via Wood-Anderson simulation (see
  [[wiki/outputs/yogya-2006-eqt-catalog]] magnitudes).
- Note: PRA(TF10a) and IMO(TF10b) are the two sites of the relocated logger
  tf3010; NGL(TF15a) is TF15's first site. Early sites of TF07a/TF09a have no XN
  code (pre-deployment/transit before the logger reached its named site).

## Links

- [[wiki/sources/yogya-2006-aftershock-edl-dataset]]
- [[wiki/syntheses/eqtransformer-yogya-2006-run]]
