# Can this catalogue improve on Koulakov (2007)? — ray coverage and checkerboard

Run before committing to LOTOS. Everything here uses the same path-length
matrices, so the coverage maps and the resolution tests cannot disagree.

## Data volume

| | Events | P+S arrivals | Arrivals per event |
|---|---|---|---|
| Koulakov et al. 2007 | 292 | ~13,000 | ~45 |
| This study | 1,005 | **22,368** | 22 |

**3.4× the events but only 1.7× the rays.** The extra events are small ones
recorded by fewer stations, so the ray count does not scale with the event
count. 21,363 of the arrivals trace successfully to a first-arriving ray
(11,772 P, 9,591 S); the rest are deeper than the tracer's branch coverage.

## Method

Rays are traced two-point through the same 1-D model used for the locations,
interpolated between layer tops so every arrival is a turning ray. Path lengths
are accumulated on a 5 km × 5 km × 5 km grid (267,000 cells over
108.8–112.2°E, 10.0–6.0°S, 0–200 km).

The synthetic inversion carries the terms that actually compete with velocity
structure:

* one free origin time per event — the dominant trade-off in local-earthquake
  tomography,
* separate P and S slowness fields,
* damping and first-difference smoothing (tuned by scan; damp 1, smooth 3),
* Gaussian noise at the pick uncertainties used in the location, 0.1 s for P and
  0.2 s for S.

A noise-free test with fixed sources would look far better than anything the
data can deliver, and would not be informative.

**These are not LOTOS numbers.** The absolute values depend on the grid,
regularisation and noise assumptions. What is robust is the *comparison between
the three configurations*, which uses identical settings throughout.

## Result — 30 km checkerboard, P wave (correlation / amplitude recovered)

| Depth | 13,000 rays (Koulakov volume) | 21,363 rays (this study, land) | + ocean-bottom stations |
|---|---|---|---|
| 0–5 km | 0.79 / 76 % | 0.81 / 77 % | 0.81 / 75 % |
| 10–15 km | 0.71 / 66 % | 0.76 / 69 % | **0.78 / 72 %** |
| 15–20 km | 0.68 / 60 % | 0.74 / 62 % | **0.78 / 71 %** |
| 20–25 km | 0.67 / 57 % | 0.70 / 59 % | **0.77 / 70 %** |
| 25–30 km | 0.58 / 45 % | 0.65 / 53 % | **0.73 / 63 %** |
| 30–35 km | 0.58 / 46 % | 0.61 / 51 % | **0.73 / 62 %** |
| 40–45 km | 0.64 / 56 % | 0.67 / 61 % | **0.72 / 67 %** |
| 50–55 km | 0.60 / 49 % | 0.63 / 56 % | **0.74 / 68 %** |
| 55–60 km | 0.37 / 28 % | 0.42 / 31 % | 0.56 / 45 % |

Cells carrying ≥ 5 rays: **19,652 → 25,888 → 39,869**.

## What this says

**The extra events alone are not the upgrade.** Land-only, this catalogue gains
about +0.03 to +0.07 in pattern correlation and 25–30 % more resolved cells over
Koulakov's data volume. Real, but not a step change — because ray coverage is
limited by the *geometry* of sources and receivers, not by how many rays follow
the same paths. Tripling the event count while halving the arrivals per event
mostly re-samples the volume that was already sampled.

**The ocean-bottom stations are the upgrade.** Adding them lifts correlation by
0.07–0.12 through 15–60 km and raises the resolved cell count by **54 %**, with
the largest gains exactly where the land-only geometry is weakest: offshore and
below 25 km, i.e. the subduction interface and the mantle wedge. At a 20 km
checker the land-only solution falls below the conventional 0.7 threshold from
15 km down; with OBS it holds 0.66–0.71 to 40 km.

This is robust to how optimistic the OBS assumption is. The table above pairs
every event with each OBS site within 250 km; a conservative 150 km cutoff still
gives +31 % cells and 0.68–0.81 correlation through 0–55 km.

## Consequence for the plan

The OBS clock-drift correction moves from "nice to have" to the **critical
path**. It is worth more to the tomography than everything the extra 713 events
bought. Sequence:

1. Measure each OBS station's residual against the land-only catalogue as a
   function of time, fit the linear drift, correct the picks.
2. Screen OS12 on the same residuals (see `../full/OBS_QC.md`).
3. Re-associate and re-locate with OBS included.
4. Then run LOTOS on the combined set, with 30 km checkers for comparability
   with Koulakov 2007 and 20 km to show what the new geometry buys.

## Caveats

* The synthetic OBS rays assume the stations detect every catalogue event inside
  the cutoff. Real detection will be lower, so treat the OBS column as the
  ceiling and the 150 km variant as the floor.
* The 15 km and 10 km checkers were also run: on the 5 km grid they give
  0.48–0.78 and 0.44–0.61 in the crust. 15 km is marginal in the top 10 km only;
  10 km is not resolved anywhere. 20 km is the honest limit to claim, and only
  with OBS.
* S-wave recovery runs 0.05–0.10 below P throughout, as expected from the larger
  pick uncertainty.

## Figures

* `../figures/tomo_coverage.png` — rays per cell at 5, 15, 25, 40, 60, 100 km
* `../figures/checker_20km.png` — input vs recovered, 20 km checkers
* `../figures/tomo_resolution_curves.png` — the three configurations vs depth
