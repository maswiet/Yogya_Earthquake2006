---
title: "Phase Association"
type: concept
status: seed
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/eqtransformer-github-repo]]"
tags:
  - seismology
  - workflow
---

## Summary

**Phase association** is grouping individual P/S phase picks from many stations
that belong to the *same* earthquake, so the combined arrivals can be located as
one event. It is the bridge between per-station picking and a located catalog.

## Key Points

- Sits downstream of detection + picking and upstream of hypocenter location.
- EQTransformer's package ships a **simple** built-in associator
  ([[wiki/sources/eqtransformer-github-repo]]); larger studies often pair pickers
  with dedicated associators (e.g. GaMMA, REAL, PhaseLink) — `needs-review`,
  those alternatives are general-knowledge, not from an ingested source.

## Links

- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/concepts/earthquake-signal-detection]]
- [[wiki/entities/eqtransformer-model]]
- [[wiki/questions/applying-eqtransformer-to-yogya-2006]] — associator + location
  are the steps needed to turn Yogya picks into a catalog.
