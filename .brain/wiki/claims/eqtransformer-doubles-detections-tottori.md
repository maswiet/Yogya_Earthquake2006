---
title: "EQTransformer detected ~2× more earthquakes with <1/3 of stations (Tottori 2000)"
type: claim
status: active
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
tags:
  - result
  - benchmark
---

## Claim

Applied to 5 weeks of continuous data from the 2000 Tottori (Japan) sequence,
EQTransformer **detected and located about twice as many earthquakes** as the
reference, while using **fewer than one third** of the available seismic stations.

## Evidence

- Verbatim abstract, [[wiki/sources/earthquake-transformer-mousavi-2020]]:
  "...we were able to detect and locate two times more earthquakes using only a
  portion (less than 1/3) of seismic stations."
- Confidence: **high** for the headline figures (stated in the abstract). The
  underlying per-station counts and location-quality metrics were not captured
  from the paper body — `needs-review` for those specifics.

## Significance

- Demonstrates that a deep-learning picker can extract a substantially more
  complete catalog from **sparser** networks — directly relevant to regions with
  limited station coverage.

## Links

- [[wiki/entities/eqtransformer-model]]
- [[wiki/entities/2000-tottori-earthquake-sequence]]
- [[wiki/concepts/earthquake-signal-detection]]
