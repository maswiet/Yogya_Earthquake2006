---
title: "Earthquake Signal Detection"
type: concept
status: seed
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
tags:
  - seismology
  - earthquake-detection
---

## Summary

Earthquake **signal detection** is deciding whether a segment of continuous
seismic data contains an earthquake signal versus noise. It is the first step in
building an earthquake catalog from raw recordings, and is hardest under **high
background noise** and for **small (micro) earthquakes**.

## Key Points

- Higher detection **sensitivity** surfaces more and smaller events, growing the
  catalog — but also raises false-positive risk that must be managed
  ([[wiki/sources/earthquake-transformer-mousavi-2020]]).
- Detection and **phase picking** are related tasks; EQTransformer performs them
  in tandem so information shared between them improves both.
- More complete detection with fewer stations lowers the network density needed
  to characterize a sequence (Tottori 2000: ~2× events from < 1/3 of stations).

## Links

- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/concepts/deep-learning-seismic-processing]]
- [[wiki/claims/eqtransformer-doubles-detections-tottori]]

## Open Questions

- How does detection sensitivity trade against false positives at the low-SNR
  tail? Not captured from primary source — `needs-review`.
