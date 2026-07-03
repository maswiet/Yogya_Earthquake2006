---
title: "Seismic Phase Picking"
type: concept
status: seed
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
tags:
  - seismology
  - phase-picking
---

## Summary

Seismic **phase picking** is identifying the precise arrival times of seismic
wave phases — primarily the **P wave** (first, compressional) and **S wave**
(later, shear) — in a recorded waveform. Accurate picks are the input to
earthquake location, magnitude, and tomography workflows.

## Key Points

- The **S-phase pick** is the harder task: S arrivals emerge within the P-wave
  coda and are lower-amplitude relative to noise. EQTransformer highlights S
  picking as a challenge it handles well ([[wiki/sources/earthquake-transformer-mousavi-2020]]).
- Traditionally done by human analysts or classical detectors (e.g. STA/LTA,
  autoregressive/AIC pickers); deep-learning pickers now approach human
  precision.
- Picking quality directly limits **hypocenter location** accuracy — the core
  need for aftershock/microearthquake studies like the Yogya 2006 work.

## Links

- [[wiki/concepts/earthquake-signal-detection]] — detection and picking are
  performed jointly by EQTransformer.
- [[wiki/concepts/deep-learning-seismic-processing]]
- [[wiki/questions/applying-eqtransformer-to-yogya-2006]]

## Open Questions

- What P/S pick precision (in seconds) does EQTransformer report vs. human
  analysts? Not captured from primary source yet — `needs-review`.
