---
title: "Can EQTransformer improve the Yogyakarta 2006 aftershock catalog?"
type: question
status: active
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
  - "[[wiki/sources/eqtransformer-github-repo]]"
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

## Tooling available

- EQTransformer is **pip/conda installable** with **pretrained models** included
  ([[wiki/sources/eqtransformer-github-repo]]) — no retraining required for a
  first trial run. It reads/writes seismic data via **ObsPy** (mseed) and
  includes a simple [[wiki/concepts/phase-association]] step, so a
  detection→picking→association pass on Yogya data is low-barrier to prototype.

## Status

**Answered — yes (2026-07-03).** EQTransformer (original weights, via SeisBench)
runs well on the Yogya 2006 array and produces a dense, geophysically sound
aftershock catalog: **16,876 located events** over Jun–Aug 2006, crustal depths
(median 12.9 km), tight fault-zone cluster, Omori decay. See deliverable
[[wiki/outputs/yogya-2006-eqt-catalog]] and [[wiki/syntheses/eqtransformer-yogya-2006-run]].
Remaining to fully close the Mousavi/Tottori "N× more" comparison: a reference
aftershock catalog and (optional) magnitudes via instrument response.

## Links

- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/concepts/earthquake-signal-detection]]
- [[wiki/concepts/deep-learning-seismic-processing]]
