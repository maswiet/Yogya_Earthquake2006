# OBS Clock Drift Measurement & Correction

## Summary
Measured clock drift for all OBS stations using the combined catalogue (1,005 events from land+OBS picks). 

## Results
- **OS12**: −0.026 s/day drift, offset +311.325 s, 10 event pairs
- **Other stations**: < 5 event pairs each (unreliable)

## Impact Assessment
- OBS picks **do not appear in association assignments** (association is land-driven)
- OBS picks do contribute to **location refinement** via NLLoc travel-time inversion
- **Zero events affected** by OS12 drift correction
- Small drift (±0.026 s/day) with large constant offsets (±300 s) indicates systematic timing offset, not true clock drift

## Conclusion
OS12 drift correction has **no practical impact** on the final catalogue. Decision: proceed with current solution to tomography without applying drift corrections. OBS network already delivers its main benefit—improved location precision through additional ray paths (317 events with errH ≤ 5 km vs. 122 land-only).

## Files Generated
- `full/obs_drift_v2.txt`: measurements from combined catalogue
- `full/picks_obs_drift_corrected.csv`: OS12-corrected picks (not used)
- `full/obs_drift_applied.txt`: correction formula and metadata
