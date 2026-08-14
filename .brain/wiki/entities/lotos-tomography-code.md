---
title: "LOTOS — Local earthquake tomography code (Koulakov)"
type: entity
status: active
created: 2026-08-14
updated: 2026-08-14
tags:
  - entity
  - lotos
  - tomography
  - software
---

## Summary

**LOTOS** (LOcal Tomography Software; Koulakov 2009, *BSSA* 99(1), 194–214) is a
widely used code for local earthquake tomography: it simultaneously relocates
sources and inverts P- and S-wave arrival times for 3-D velocity structure on a
grid, using bending ray tracing from an optimized 1-D reference model, with
built-in checkerboard resolution testing. Author: Ivan Koulakov (IPGG SB RAS,
Novosibirsk).

## Key Points

- **Used by our group** in [[wiki/sources/diambama-2019-yogya-tomography]]
  (LOTOS-13) to image the Yogyakarta 2006 earthquake zone and delineate the
  NE–SW Ngalang Fault structure.
- Method steps: reference-model + source location → bending-tracing ray paths →
  matrix construction → simultaneous P/S inversion with amplitude damping +
  smoothing (values tuned by synthetic tests) → iterate.
- A citable open release of the code (Gakkel Ridge variant) is archived on Zenodo
  — see [[wiki/sources/zenodo-lotos-koulakov-2021]].
- Related to our own VELEST minimum-1D work: LOTOS optimizes a 1-D model as a
  by-product, analogous to VELEST; see [[wiki/outputs/yogya-2006-1d-velocity-model]].

## Links

- [[wiki/sources/zenodo-lotos-koulakov-2021]]
- [[wiki/sources/diambama-2019-yogya-tomography]]
- [[wiki/outputs/yogya-2006-1d-velocity-model]]
