---
title: "Yogya 2006 Aftershock EDL Dataset (DATA-GFZ-Gempa-JOgja-tahap-2)"
type: source
status: active
created: 2026-07-03
updated: 2026-07-03
source_author: "GFZ Potsdam temporary aftershock deployment"
source_url: "local: /Volumes/Untitled 1/DATA-GFZ-Gempa-JOgja-tahap-2"
tags:
  - dataset
  - yogya-2006
  - waveforms
  - raw-data
---

## Summary

Raw continuous seismic recordings from the **temporary aftershock deployment**
following the **27 May 2006 (M6.3) Yogyakarta earthquake**, on an external volume
`"/Volumes/Untitled 1/DATA-GFZ-Gempa-JOgja-tahap-2"`. GFZ EarthData Logger (EDL)
data. This is the primary input for the EQTransformer reprocessing effort
([[wiki/syntheses/eqtransformer-yogya-2006-run]]).

## Layout

- One folder per station: `tf3007, tf3009–tf3019` (12 stations; **no tf3008**).
- Inside each: julian-day subfolders (`153`, `154`, …). **Julian-day numbers
  repeat across years** — the folders also contain pre-deployment data.
- Per 30-minute segment: `*.pri0/.pri1/.pri2` = the 3 seismometer components
  (valid **miniSEED**), plus `.gps/.gst/.msg/.pll` status/log files.
- Component mapping (user-confirmed): **p0 = Z, p1 = N, p2 = E**.
- Native sampling: **200 Hz** (resampled to 100 Hz for EQTransformer).
- SEED header station id inside miniSEED = `e30XX` (network placeholder `nn`,
  location `11`, channels `p0/p1/p2`).
- Separate `YogyaData/` folder = regional **IA-network broadband** daily files
  (BHZ/BHN/BHE, e.g. `IA.BAKI..BHZ.D.2006.146`) — NOT the temp array.

## Critical data-hygiene facts

- **Year mixing:** filename encodes year as `e30XX`**`YY`**`MMDDHH…`. Stations
  **tf3007 (04/05/06), tf3012 (05/06), tf3017 (04/05/06), tf3019 (04/06)**
  contain 2004–2005 **pre-deployment test data recorded in Germany** (GPS ≈
  +52.38, +13.06 = GFZ Potsdam). Others are 2006-only. → **Must filter by
  filename year `06` AND Yogya-region GPS box** before processing.
- **Station relocations:** GPS-fix spread shows some loggers moved mid-campaign
  (TF07 ≈36 km, TF17 ≈28 km, TF10 ≈19 km, TF15 ≈13 km, TF09 ≈9 km spread). A
  single median location is wrong for those — handle per-time-period at scale.
  Stable single-site stations: **TF12, TF13, TF14, TF16, TF18, TF19** (spread
  ≤ ~0.4 km; TF11 ≈2 km).
- Deployment window: julian days ~**153–241** (June–Aug 2006), ~90 days/station.

## Links

- [[wiki/entities/yogya-2006-temp-aftershock-network]] — station coordinate table.
- [[wiki/syntheses/eqtransformer-yogya-2006-run]] — processing pipeline + status.
- [[wiki/questions/applying-eqtransformer-to-yogya-2006]]

## Open Questions

- Official FDSN network code for the deployment (using placeholder `YK`) — `needs-review`.
- Sensor type / band code (assumed short-period → `EH?`) — `needs-review`.
- Per-period coordinates for the relocated stations — to resolve at scale.
