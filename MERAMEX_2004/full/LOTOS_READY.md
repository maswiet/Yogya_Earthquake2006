# MERAMEX 2004 — Full Campaign LOTOS Ready

## Catalogue Summary
- **Events**: 1005 (1,005 locations with OBS refinement)
- **Picks**: 93047 (67,853 land + 25,194 OBS)
- **Stations**: 145 total (128 land + 14 OBS + 3 auxiliary)
- **Location precision**: 317 events with errH ≤ 5 km (vs 122 land-only)
- **Median depth**: 37.9 km (range 3.0–148.1 km)

## Ray Coverage
From earlier coverage assessment:
- **Total traced rays**: 21,363+ rays (P + S combined)
- **P-wave rays**: ~12,000 
- **S-wave rays**: ~9,400
- **Cells with ≥5 rays**: 39,869 (5 km × 5 km × 5 km grid, 0–200 km depth)

## Comparison to Koulakov et al. 2007
| Metric | Koulakov | This study | Gain |
|--------|----------|-----------|------|
| Events | 292 | 1,005 | 3.4× |
| Rays | ~13,000 | ~21,400 | 1.6× |
| Rays/event | ~45 | 21.3 | −53% |
| Resolved cells (≥5 rays) | — | 39,869 | +54% with OBS |
| Pattern correlation @ 20 km | — | 0.68–0.71 (0–40 km depth) | +0.07–0.12 vs land-only |

## OBS Integration Summary
- **OBS clock drift**: Measured from combined catalogue
  - OS12: −0.026 s/day (10 event pairs)
  - Others: < 5 pairs each (unreliable)
  - **Decision**: Apply no corrections (zero affected events in final catalogue)
  
- **OBS contribution to precision**: 317 events with errH ≤ 5 km
  - 2.6× improvement over land-only (122 events)
  - Most significant below 20 km depth
  
- **OBS contribution to ray coverage**: 
  - 54% more resolved cells
  - Strongest gain offshore and in mantle wedge (25–60 km)
  - Essential for imaging subduction interface

## Workflow Status
✅ Phase picking: EQTransformer (128 land stations + 14 OBS sites)
✅ Event association: PyOcto (land-driven detection, 1,005 events)
✅ Hypocenter location: NonLinLoc with OBS refinement
✅ OBS clock drift: Measured and assessed (not applied—negligible impact)
✅ Ray tracing: Completed for coverage and resolution assessment
✅ **Catalogue validation**: Ready for tomographic inversion

## Next: LOTOS 3-D Inversion
- Input: catalog_obs_combined.csv (1,005 events)
- Picks: picks_combined.csv (93,047 arrivals)
- Velocity model: VELEST 1-D (16 layers, −3 to 210 km)
- Grid: 5 km × 5 km × 5 km (108.8–112.2°E, 10.0–6.0°S, 0–200 km)
- Checkerboards: 20 km and 30 km cells for resolution testing
- Expected output: 3-D P/S velocity anomalies, Wadati–Benioff zone image
