---
title: "2000 Tottori Earthquake Sequence (Japan)"
type: entity
status: seed
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
tags:
  - earthquake
  - dataset
  - japan
---

## Summary

The **2000 Tottori earthquake** sequence in Japan is the real-data validation
case for EQTransformer. Mousavi et al. applied the model to **5 weeks of
continuous seismic data** recorded during the sequence.

## Facts

- Used as an independent field test (not part of training).
- Result: EQTransformer detected and located **~2× more earthquakes** than the
  reference catalog, using **less than 1/3** of the available seismic stations
  ([[wiki/claims/eqtransformer-doubles-detections-tottori]]).

## Links

- [[wiki/entities/eqtransformer-model]]
- [[wiki/sources/earthquake-transformer-mousavi-2020]]

## Open Questions

- Mainshock magnitude/date and reference-catalog details not captured from the
  primary source — `needs-review` (general knowledge: the 2000 Western Tottori
  event was ~Mw 6.6, 6 Oct 2000, but confirm before citing).
