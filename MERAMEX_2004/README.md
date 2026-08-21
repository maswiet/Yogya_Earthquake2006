# MERAMEX 2004 — machine-learning reprocessing

Reprocessing of the **MERapi AMphibious EXperiment** (UGM + GFZ + IFM-GEOMAR,
May–October 2004, Central Java) with the same EQTransformer → PyOcto →
NonLinLoc chain that produced the 2006 Yogyakarta aftershock catalog in
[`../eqt`](../eqt).

## The archive

Two read-only NTFS volumes, ~512 GB, **already miniSEED** — no GIPPtools
conversion is needed:

| Volume | Content | DOY (2004) |
|---|---|---|
| `/Volumes/Untitled` | `MERAMEX DATA` — land station-days + 4 OBS recovery folders | 127–260 |
| `/Volumes/Untitled 1` | `MERAMEX_LANJUTAN` | 261–282 |

Layout is `<DOY-folder>/<STATION>/<DDD>/<segments>`, with two nested exceptions
(`150-152/` and `DOY 162_290405/DOY 161_280405/`). **The station code lives in
the directory name only** — the EDL miniSEED header carries the recorder serial
instead.

| Family | Sensor | Rate | Encoding | Header channels |
|---|---|---|---|---|
| EDL (`.PRI0/1/2`) | L4-3D 1 Hz, Trillium T40 | 100 Hz | STEIM1, little-endian | `p0/p1/p2` → Z/N/E |
| SAM (Kiel) | Trillium T3 / 3T | 100 Hz | STEIM1, big-endian | `SPZ/SPN/SPE` |
| SAM (GFZ) | Güralp 3ESP | 100 Hz | STEIM1, big-endian | `BHZ/BHN/BHE` |
| OBH/OBS | HTI hydrophone, DPG | 50 Hz | STEIM2, big-endian | `P0…P3` |

Totals from `config/station_days.csv`: **15,717 station-days** — 14,012 land
(496 GB) and 1,705 OBS (16.5 GB) — across **143 station codes / 169 deployment
periods** (139 EDL, 16 SAM, 9 OBH, 5 OBS). 21 codes were re-occupied at a new
site mid-campaign, so geometry is always resolved *per julian day*.

## Layout

```
config/   INFO.DAT + parsed station table, per-period geometry, station-day index
scripts/  pipeline (see below)
pilot/    pilot picks / associations / catalog
full/     full-campaign products
nll/      NonLinLoc control files, travel-time grids, location output
figures/  summary figures
logs/     run logs
```

## Pipeline

| Step | Script | Env | Output |
|---|---|---|---|
| 1. index the archive | `index_data.py` | base | `config/station_days.csv` |
| 2. per-period geometry | `build_periods.py` | base | `config/stations_periods.json` |
| 3. detect + pick | `run_picks.py` | `eqt` | `*/picks.csv` (streaming, resumable) |
| 4. associate | `associate.py` | `assoc` | `*/events.csv` + assignments |
| 5. NLL inputs | `gen_nll.py` | base | `nll/nll_*.in`, `nll/obs/*.obs` |
| 6. locate | `Vel2Grid`/`Grid2Time`/`NLLoc` | — | `nll/loc/*.hyp` |
| 7. catalog | `parse_nll.py` | base | `*/catalog_nll.csv` |
| 8. summary figure | `plot_pilot.py` | base | `figures/*_summary.png` |
| 9. GMT map | `plot_map_gmt.py` | `gmt` | `figures/*_map_gmt.png` |
| 10. slab sections | `plot_slab_section.py` | `gmt` | `figures/*_slab_section.png` |

`run_pilot.sh` chains steps 4–8 for any working set:
`WORK=full TAG=full bash run_pilot.sh`. Shared archive I/O lives in `mxio.py`.
`plot_event.py` draws a record section for one associated event, and
`plot_opak.py` overlays the 2004 catalogue on the 2006 aftershock cloud.

The two GMT scripts follow the house style of [`../eqt/scripts/plot_relief.py`](../eqt/scripts/plot_relief.py)
— shaded relief, `geo` CPT, turbo depth colours — but keep real bathymetry in the
ocean instead of a flat blue fill, because the trench and forearc are the point.
`load_earth_relief` is called with `resolution="30s", registration="pixel"`:
that is the one tile set whose local cache spans the whole Java margin down to
11°S.

## Runs

| Run | Stations | Days | Station-days | Picks | Events | Median errH |
|---|---|---|---|---|---|---|
| `pilot/` | 43 (Yogya–Merapi) | 155–165 | 452 | 2,577 | 67 | 8.5 km |
| `wide11/` | 108 (all land) | 155–165 | 1,126 | 5,400 | 79 | 6.7 km |
| `full/` | 108 land + 14 OBS | **127–282** | 15,717 | running | — | — |

`wide11` was the control that showed extra land stations buy **accuracy, not
count** (+18 % events but ×2.9 more that pass the quality cut) — see
[`wide11/SUMMARY.md`](wide11/SUMMARY.md). The full run therefore scales with
*time*, and the OBS are the addition that matters for the offshore majority of
the catalogue.

### OBS handling

The 14 OBH/OBS sites are picked separately at a lower threshold (0.2 instead of
0.3): they are noisier, and EQTransformer's STEAD training set contains no
seafloor data. The eight `OH*` sites are **hydrophone-only** — one channel, so
no S phase is physically possible there and any S pick from them must be
dropped. **Clock drift is not yet corrected**: OBS have no GPS underwater, so
the plan is to locate with land stations only, measure each OBS station's
residual against time, fit a linear drift, and relocate. Until that is done, OBS
arrival times should be treated as uncalibrated.

Picking never writes waveforms to disk — the archive is ~500 GB and the internal
disk has ~80 GB free, so `run_picks.py` reads a station-day, picks it, and keeps
only the picks. `done.txt` makes any run resumable.

## Velocity model

Crust = the VELEST minimum-1D model derived from the 2006 Yogyakarta aftershocks
(`../eqt/nll/nll_vel_v.in`); below 40 km the ak135 mantle is appended so
Wadati–Benioff events are not forced into a crustal half-space. NonLinLoc uses
2-D (distance–depth) travel-time grids, which is the right choice for a 1-D
model and tens of stations.

## Reference catalogs to beat

Published MERAMEX local-earthquake catalogs: **292 events** (Wagner et al. 2007,
GJI) and later **505 processed / 344 used** (Koulakov et al. 2007/2009).
