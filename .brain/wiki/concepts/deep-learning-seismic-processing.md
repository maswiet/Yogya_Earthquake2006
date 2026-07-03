---
title: "Deep-Learning Seismic Processing (Attention & Multi-Task)"
type: concept
status: seed
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
tags:
  - deep-learning
  - seismology
  - attention
---

## Summary

Using neural networks — increasingly **attention/transformer**-based — to process
raw seismic waveforms for detection and phase picking. Two ideas central to
EQTransformer:

1. **Hierarchical attention mechanism** — lets the model weight informative parts
   of the waveform (the phase arrivals) within the full signal.
2. **Multi-task tandem learning** — training detection and P/S picking jointly so
   each task's information reinforces the others, beating single-task models.

## Key Points

- Reported to outperform both prior **deep-learning** pickers and **traditional**
  (non-ML) detection/picking algorithms ([[wiki/sources/earthquake-transformer-mousavi-2020]]).
- Trained on **global** seismic data, aiming for a general-purpose model rather
  than a network-specific one.
- Practical payoff: comparable-to-human picks with far higher throughput, enabling
  denser catalogs from sparse networks.

## Links

- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/concepts/earthquake-signal-detection]]
- [[wiki/entities/eqtransformer-model]]

## Open Questions

- Exact architecture (CNN/LSTM/transformer layer counts, decoder branches) and
  training dataset name/size were not captured from a primary source —
  `needs-review`; read the full paper/methods to fill in.
