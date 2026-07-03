---
title: "Earthquake Transformer (Mousavi et al., 2020)"
type: source
status: active
created: 2026-07-03
updated: 2026-07-03
source_author: "S. Mostafa Mousavi; William L. Ellsworth; Weiqiang Zhu; Lindsay Y. Chuang; Gregory C. Beroza"
source_created: 2020-01-01
source_url: "https://www.nature.com/articles/s41467-020-17591-w"
doi: "10.1038/s41467-020-17591-w"
raw: "[[raw/2026-07-03_earthquake-transformer-mousavi-2020/source]]"
tags:
  - deep-learning
  - seismology
  - phase-picking
  - earthquake-detection
---

## Summary

Introduces **EQTransformer**, a global deep-learning model that performs
earthquake signal **detection** and **P/S seismic phase picking**
*simultaneously*. Doing both related tasks in tandem — sharing information
between the discrete phase arrivals and the full waveform through a
**hierarchical attention mechanism** — improves performance on each task
individually. Published in *Nature Communications* 11:3952 (2020).

Full citation: Mousavi, S.M., Ellsworth, W.L., Zhu, W., Chuang, L.Y., Beroza,
G.C. (2020). *Nature Communications* 11, 3952. DOI 10.1038/s41467-020-17591-w.

## Key Points

- Single model, three tasks in tandem: event detection + P-phase pick + S-phase
  pick. Multi-task tandem learning improves each individual task.
- Uses a hierarchical **attention mechanism** over earthquake waveforms.
- Reported to **outperform** prior deep-learning and traditional
  phase-picking/detection algorithms.
- Validation case: **5 weeks of continuous data** from the **2000 Tottori
  (Japan)** aftershock sequence — detected & located **~2× more earthquakes**
  using **< 1/3 of the seismic stations**.
- P and S picks are close to human-analyst manual picks; higher sensitivity
  surfaces more and **smaller** events (relevant to microearthquake monitoring).
- Trained on global seismic data (per project docs); layer-level architecture and
  training-set name/size not captured from a primary source — see needs-review.

## Extracted Concepts

- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/concepts/earthquake-signal-detection]]
- [[wiki/concepts/deep-learning-seismic-processing]]

## Extracted Entities

- [[wiki/entities/eqtransformer-model]]
- [[wiki/entities/2000-tottori-earthquake-sequence]]

## Claims

- [[wiki/claims/eqtransformer-doubles-detections-tottori]]

## Conflicts / Updates

- None yet (first source in vault).

## Open Questions

- [[wiki/questions/applying-eqtransformer-to-yogya-2006]]
