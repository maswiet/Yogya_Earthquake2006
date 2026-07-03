---
title: "EQTransformer (model)"
type: entity
status: active
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
tags:
  - model
  - deep-learning
  - tool
---

## Summary

**EQTransformer** is the global deep-learning model introduced by Mousavi et al.
(2020) for simultaneous earthquake **detection** and **P/S phase picking**, built
around a hierarchical attention mechanism. It is distributed as an open-source
Python package (docs: https://eqtransformer.readthedocs.io/).

## Facts

- Origin paper: [[wiki/sources/earthquake-transformer-mousavi-2020]] (Nature
  Communications 11:3952, 2020; DOI 10.1038/s41467-020-17591-w).
- Authors / group: S. Mostafa Mousavi, William L. Ellsworth, Weiqiang Zhu,
  Lindsay Y. Chuang, Gregory C. Beroza (Stanford-affiliated seismology group).
- Outputs: earthquake-detection probability + P-arrival pick + S-arrival pick.
- Trained on global seismic data; intended as a general-purpose picker.
- Validated on the [[wiki/entities/2000-tottori-earthquake-sequence]].

## Relevance to this project

- Candidate tool for building/densifying an aftershock catalog for the
  **Yogyakarta 2006** earthquake from continuous waveform data — see
  [[wiki/questions/applying-eqtransformer-to-yogya-2006]].

## Links

- [[wiki/concepts/deep-learning-seismic-processing]]
- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/claims/eqtransformer-doubles-detections-tottori]]
