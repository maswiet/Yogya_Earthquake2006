---
title: "Can EQTransformer improve the Yogyakarta 2006 aftershock catalog?"
type: question
status: active
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
tags:
  - open-question
  - method-choice
---

## Question

Is EQTransformer ([[wiki/entities/eqtransformer-model]]) a suitable tool for
detecting and phase-picking aftershocks of the **2006 Yogyakarta earthquake**
from continuous waveform data — potentially producing a denser/more complete
catalog than manual or classical processing?

## Why it matters

- The Tottori 2000 result shows ~2× more detected events from < 1/3 of stations
  ([[wiki/claims/eqtransformer-doubles-detections-tottori]]). If it generalizes,
  a sparse Yogya-region network could yield a much richer aftershock catalog.

## What we'd need to confirm

- Availability and format of continuous waveform data for the Yogya 2006 period
  and network (stations, sampling rate, components).
- Whether the globally-trained model generalizes to Indonesian
  crustal/instrument conditions, or needs fine-tuning/transfer learning.
- Downstream location workflow (velocity model, associator) to turn picks into
  hypocenters.

## Status

Open — no project data assessed yet. First seismic-ML source in the vault.

## Links

- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/concepts/earthquake-signal-detection]]
- [[wiki/concepts/deep-learning-seismic-processing]]
