---
title: "LOTOS code release (Koulakov 2021, Zenodo) — Gakkel Ridge tomography"
type: source
status: active
created: 2026-08-14
updated: 2026-08-14
source_author: "Koulakov, Ivan (IPGG SB RAS)"
source_created: 2021-08-31
source_url: "https://zenodo.org/records/5338981"
raw: "[[raw/2026-08-14_zenodo-lotos-koulakov/provenance]]"
tags:
  - source
  - lotos
  - tomography
  - software
---

## Summary

Open (CC BY 4.0) Zenodo release of the **LOTOS local-earthquake-tomography code**
by its author, Ivan Koulakov (DOI 10.5281/zenodo.5338981, 2021-08-31). This
particular release is the **Gakkel Ridge (Arctic, ~85°E)** application — a Windows
build (VS2010 + Intel Visual Fortran) with the full program listing, initial data,
parameters, and README, adapted for a mobile network on ice floes (adds a
constant-velocity water layer). File: `lotos_gakkel_release.zip` (29.6 MB, not
downloaded — metadata captured only).

## Key Points

- LOTOS (Koulakov 2009) is the tomography engine used by our group in
  [[wiki/sources/diambama-2019-yogya-tomography]] (LOTOS-13) for the Yogyakarta
  2006 velocity model — so this is a citable, versioned copy of that code from
  the method's author.
- This release is **NOT Yogyakarta data**: it is the Arctic Gakkel-Ridge
  application. Only the code/method is transferable; the water-layer adaptation
  is specific to the marine setting.
- Workflow LOTOS implements: 1-D reference model + source relocation → bending
  ray tracing → simultaneous P/S velocity inversion on a grid, with checkerboard
  resolution tests — the same procedure described in Diambama et al. 2019.

## Extracted Entities

- [[wiki/entities/lotos-tomography-code]]

## Conflicts / Updates

- None. Complements the Diambama 2019 source (their methods reference LOTOS-13).

## Open Questions

- Which LOTOS version did Diambama et al. 2019 actually run (LOTOS-13 vs this
  later release), and is a re-run with a current LOTOS worthwhile for the SRL
  paper's velocity-model section? `needs-review`.
