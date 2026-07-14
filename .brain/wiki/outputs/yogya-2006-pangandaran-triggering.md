---
title: "Pangandaran Mw7.7 (17 Jul 2006) triggering test on Yogya aftershocks"
type: output
status: active
created: 2026-07-14
updated: 2026-07-14
sources:
  - "[[wiki/outputs/yogya-2006-eqt-catalog]]"
tags:
  - stress-triggering
  - aftershocks
  - completeness
---

## Summary

Tested whether the **17 Jul 2006 Pangandaran Mw7.7 tsunami earthquake** (USGS
`usp000ensm`; 08:19:27 UTC; −9.284°/107.419°; 20 km; shallow subduction thrust
NP1 290/10/102) perturbed the ongoing **Yogyakarta 2006 aftershock** rate.
**Conclusion: no detectable static or dynamic triggering.** The apparent
post-17-Jul drop in small-event rate is a **catalog-completeness artifact**, not a
real change in seismicity.

Script: `eqt/scripts/analyze_pangandaran.py` → `eqt/figures/pangandaran_rate.png`,
`eqt/full/pangandaran_summary.txt`. Catalog: `full/catalog_magnitude.csv`
(gap<180). Pangandaran falls at **t = 51.4 d** after the Yogya mainshock
(26 May 2006 22:53:58 UTC).

## Key Points

- **Magnitude-dependent β is the decisive diagnostic.** β-statistic (Matthews &
  Reasenberg; obs vs pre-event local rate) computed vs a sliding ML threshold:
  - +3 d window: β = **−6.4 at ML≥0.0** → rises to **≈0 at ML≥0.8–1.0** and flat above.
  - +7 d window: β = **−11 at ML≥0.0** → −1.5 to −2 above Mc.
  - A *real* rate change would be magnitude-independent. A change that only
    affects the smallest events and **vanishes above Mc is completeness**, not
    physics.
- **Robust (ML≥1.0) rate is unchanged.** Pre 10 d ≈ 4.6/d vs post 3 d ≈ 6.0/d
  (ratio 1.30, β≈+0.1) — i.e. slightly *higher*, statistically flat. Only the
  ML≥0.3 rate "drops" (ratio ≈0.79), consistent with M7.7 coda-masking + any
  network/completeness change from mid-July on.
- **No immediate spike.** 6-hour binned rate ±7 d around 08:19 UTC shows no step
  or burst at t=0 (the one +1.8 d bin excursion appears in *both* magnitude bands
  = an ordinary local aftershock burst, not a coincident trigger).
- **Omori baseline** (curve_fit, full window excluding the Pangandaran window,
  bounded): K≈3.3e4, c≈26.5 d, p≈1.60. The sequence is vigorously decaying; the
  rate rides the Omori curve straight through 17 Jul.

## Stress mechanisms (order of magnitude, at R≈364 km epicentral / ~250 km to rupture)

- **Static Coulomb** ~ M0/(μR³) ≈ **3×10⁻¹² bar** — 7+ orders below the ~0.1–1
  bar nominal triggering threshold. Ruled out (static stress falls ~1/r³).
- **Dynamic** ~ μ·PGV/Vs with PGV~1 cm/s (long-period-rich slow tsunami
  earthquake) ≈ **~0.86 bar** — *above* the nominal dynamic threshold, yet **no
  rate increase is observed**. So even a non-trivial dynamic transient did not
  measurably modulate an already-high, fast-decaying aftershock rate.

## Interpretation

Pangandaran did **not** measurably affect the Yogya aftershock population.
Static triggering is impossible at this distance; dynamic shaking was
order-0.1–1 bar but produced no detectable rate change — plausibly because the
Opak-fault aftershock volume was already responding vigorously to the mainshock
stress and/or was not in a near-failure state receptive to a remote transient.
The one seemingly-suggestive signal (small-event deficit) is an observational
artifact.

## Caveats / Open Questions

- Completeness/network changes after mid-July are inferred from the
  magnitude-dependence, not independently documented (station up-time log would
  confirm). `needs-review`.
- Analysis is bounded by the temp-array recording window (through ~late Aug 2006).
- PGV at Yogya is an order-of-magnitude estimate; no on-scale regional record was
  used.

## Links

- Catalog: [[wiki/outputs/yogya-2006-eqt-catalog]]
- Velocity model: [[wiki/outputs/yogya-2006-1d-velocity-model]]
- Figure: `eqt/figures/pangandaran_rate.png`
