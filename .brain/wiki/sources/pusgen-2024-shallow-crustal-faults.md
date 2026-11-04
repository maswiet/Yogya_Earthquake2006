---
title: "PUSGEN 2024 Shallow-Crustal Active-Fault Database"
type: source
status: active
created: 2026-11-04
tags:
  - fault-database
  - pusgen-2024
  - yogya-2006
---

## Summary

National Indonesian active-fault source model from **PUSGEN 2024** (Pusat Studi
Gempa Nasional / National Center for Earthquake Studies) — the "Peta Sumber dan
Bahaya Gempa Indonesia 2024" shallow-crustal (401) fault set. User-provided as the
**authoritative replacement** for the previously hand-digitized fault traces
(which were removed in the geology map because their positions were unreliable —
see [[wiki/outputs/yogya-2006-eqt-catalog]]).

## Provenance

`/Users/maswiet/Downloads/PUSGEN2024_Shallow_Crustal_Polygon/`:
- `01_..._Trace_and_Polygon_Primary.gpkg` — 3 layers: `trace_frozen_401`
  (fault surface traces, LineString), `polygon_primary_r5` (fault-source surface
  projections, dip-aware + 5 km buffer, Polygon), `primary_overlap_ambiguity`.
- `02_..._Polygon_Vertex_Coordinates.csv`, `03_..._Polygon_Parameters.csv`,
  `04_..._Polygon_Map.pdf`, `05_README_Methodology.pdf`. CRS EPSG:4326.

Extracted to repo: `eqt/data/pusgen2024/{pusgen_traces,pusgen_polygons}.geojson`
(the 4 in-region faults only), via `eqt/scripts/plot_pusgen_map.py` helper.

## In-region faults (Yogya 110.15–110.75°E, −8.15 to −7.63°)

| ID | Fault / segment | Type | L (km) | Dip | Mmax |
|----|-----------------|------|--------|-----|------|
| **JAV67** | **Opak Fault** | **SS-LL** | **47** | **50°** | **6.9** |
| JAV68 | Merapi–Merbabu | SS | 28 | 90° | 6.5 |
| JAV72 | Mataram – Tambakbayan | SS-LL | 14 | 90° | 6.1 |
| JAV73 | Mataram – Dengkeng | SS-LL | 15 | 90° | 6.2 |

Opak trace runs (110.277,−8.025) → (110.556,−7.710), i.e. **NE–SW**.

## Key Points / Findings (overlay on GrowClust catalog)

- **Opak = the 2006 Yogya source** in PUSGEN: SS-LL, **Mmax 6.9**, and notably
  **dip 50°** (not vertical) — consistent with the steeply-dipping plane the
  GrowClust double-difference relocation resolves, and with the GCMT/USGS/GFZ
  near-vertical left-lateral mechanisms.
- **Epicentral offset:** the relocated aftershock cloud sits **~5–15 km E/NE of
  the mapped Opak *surface* trace**, in a parallel NE–SW band. The PUSGEN
  dip-aware *source polygon* (which projects the fault plane down-dip) does cover
  the cloud, so the offset is the expected surface-trace-vs-at-depth separation
  of a dipping/blind structure rather than a contradiction — matches the
  long-standing result that the 2006 rupture was east of the Opak surface trace.
- **Mataram (Tambakbayan/Dengkeng) and Merapi–Merbabu** faults lie N of the
  sequence and are not illuminated by the aftershocks.

## Links

- Figure: `eqt/figures/growclust_pusgen_map.png`
- [[wiki/outputs/yogya-2006-eqt-catalog]] (GrowClust catalog + mechanisms)
- Supersedes the hand-digitized Ramdhan-et-al-style fault traces for mapping.

## Open Questions

- Opak `Dip_dir_parsed` is blank in the parameter table (SS fault) — PUSGEN does
  not assert an E/W dip direction, so the hanging-wall side is not fixed from this
  source alone. A cross-section of GrowClust hypocentres across strike would test
  the dip direction directly (`needs-review`).
