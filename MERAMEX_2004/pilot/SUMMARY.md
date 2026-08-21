# MERAMEX pilot — Stage 1 result

**Scope:** 43 land stations within 65 km of the Yogyakarta–Merapi centre,
DOY 155–165 (3–13 June 2004) = 452 station-days, 15.5 GB read.

## Settings

| Step | Setting |
|---|---|
| Picker | SeisBench `EQTransformer.from_pretrained("original")`, MPS |
| Thresholds | P = S = detection = 0.30 |
| Components | EDL `p0/p1/p2` → Z/N/E; SAM `SPZ/SPN/SPE` and `BHZ/BHN/BHE` |
| Associator | PyOcto, `n_picks = 8`, `n_p_and_s = 3`, depth 0–250 km |
| Association area | lat −11.25…−6.62, lon 109.30…111.59 |
| Velocity model | VELEST minimum-1D (Yogya 2006) crust + ak135 mantle |
| Location | NonLinLoc, OCT search, 2-D travel-time grids, 601×601×254 km at 1 km |

## Numbers

| | |
|---|---|
| Picking throughput | 5.1 s per station-day (MPS), 38.5 min total |
| Picks | 2,577 — 1,351 P (mean prob 0.61), 1,226 S (mean prob 0.49) |
| Picks per station | median 39, range 2–249 over 11 days |
| Picks used by the associator | 1,316 (51 %) |
| Events located | **67** in 11 days = **6.1 / day** |
| Median phases per event | 17 |
| Median RMS | 0.27 s |
| Depth | median 41 km, 10–90 % = 22–107 km |

## What the pilot shows

- The N–S section recovers the **north-dipping Wadati–Benioff zone** without
  being told about it: ~20–30 km depth near −9.1°, ~120 km near −8.2°.
- Record sections (`figures/event_013_deep.png`, `event_007_shallow.png`) show
  impulsive P on the vertical with clean moveout to 90 km, confirming the
  `p0 = Z` convention for this deployment.
- Threshold sensitivity is mild — 50–81 events across `n_picks` 6–10 and
  `n_p_and_s` 2–4 — so the catalog is not a threshold artifact.

## Known limits of this pilot

- Seismicity is dominated by the **offshore subduction zone**, and the array is
  land-only, so the median azimuthal gap is 318°. Only 5 of 67 events fall
  inside the array. **The 14 OBH/OBS stations are the fix** and should be added
  before any depth-sensitive interpretation.
- 43 of 122 land stations were used; the full network will lower the detection
  threshold and tighten geometry.
- No magnitudes yet (needs instrument responses for L4-3D / T40 / T3 / 3ESP).

## Projection to the full campaign

14,012 land station-days at 5.1 s ≈ **20 h of picking**. At the pilot rate alone
(6.1 events/day × ~150 days) the full run yields ~900 events; with the complete
network the rate should be several times higher. Published MERAMEX catalogs for
comparison: 292 events (Wagner et al. 2007), 505 processed / 344 used
(Koulakov et al. 2007).
