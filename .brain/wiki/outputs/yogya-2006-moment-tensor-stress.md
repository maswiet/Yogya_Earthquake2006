---
title: "Aftershock moment tensors & stress direction — Yogya 2006"
type: output
status: active
created: 2026-11-05
updated: 2026-11-05
sources:
  - "[[wiki/syntheses/eqtransformer-yogya-2006-run]]"
tags:
  - yogya-2006
  - moment-tensor
  - stress-inversion
  - bayesisola
  - result
---

## Goal

Compute aftershock focal mechanisms (full-waveform moment tensors) with a usable
confidence level, then derive the statistical stress direction and overlay SHmax
on topography. User asked to start from ML>1.5 for a larger population.

## Method

- **Moment tensors:** BayesISOLA (Vackář et al. 2017) deviatoric full-waveform
  inversion. Axitra Green's functions on the VELEST 1-D model
  ([[wiki/outputs/yogya-2006-1d-velocity-model]]); L4-3D response (1 Hz, 0.707).
  Band 0.8–2.0 Hz (ML≥2.5), 1.0–2.5 (≥2.0), 1.2–3.0 (else). Per-event SAC windows
  from raw EDL (PRE/POST 90 s). Grid depth_unc=12 km, max_points=300.
  Key speed fix: Axitra `xl` patched from 2000 km (regional default) to 8×max_dist
  (~256 km) for our ~32 km network → tens-of-× faster.
- **Stress inversion:** Michael (1984, 1987) linear inversion with bootstrap
  (3000) + **random nodal-plane selection** per event (SATSI/MSATSI-style);
  SHmax via Lund & Townend (2007). Self-validated: imposed SHmax 20°/R 0.5 →
  recovered 19.4°/0.49. Scripts: `eqt/scripts/stress_inversion.py` (numpy,
  `--selftest`), `eqt/scripts/plot_stress_map.py` (pygmt topo overlay).

## Results

- **Batch:** 219 events ML≥1.5 (gap<180), 217 inverted OK, 2 trivial errors.
  Quality yield: **VR≥0.5 → 27**, VR≥0.4 → 53, VR≥0.3 → 97.
- **VR is magnitude-independent** — best single solution is an **ML 1.82** event
  (VR 0.66). SNR/location/path control quality, not size. (Vindicates ML>1.5.)
- **Robust, plane-independent observable (from 27 VR≥0.5 MT):**
  - Extension (Shmin) **horizontal ~N112°E ± 20°**, plunge ~16° (T-axes cluster;
    consistent at VR≥0.4, n=53: N107°E ± 21°).
  - **SHmax ~N22°E ± 20°** (NNE–SSW), ⟂ horizontal extension.
  - Regime **strike-slip-dominant transtension**: 17 SS, 9 normal, 1 reverse.
- **Formal full stress tensor is NOT resolved:** bootstrap cones ~70° on all
  three axes, SHmax histogram ~flat, **R ≈ 0.30** (σ1≈σ2 near-degenerate),
  mean misfit 23°. Honest reason: short-period L4-3D + 1-D model give per-event
  MT uncertainty too large, and random-nodal-plane bootstrap for a
  strike-slip-heavy set amplifies σ1/σ2 trade-off. P/T axes are intrinsic to each
  MT (plane-independent) → the extension direction survives; the tensor does not.

## Tectonic reading

WNW–ESE horizontal extension / NNE SHmax in transtension is consistent with the
**Opak fault** (NE–SW, oblique left-lateral; PUSGEN dip 50°,
[[wiki/sources/pusgen-2024-shallow-crustal-faults]]) as the controlling structure;
the aftershock cloud (beachballs) sits ~5–15 km E of the Opak surface trace.
Mw from MT runs ~0.9–1.0 unit above our ML (ML 2.84 → Mw 2.98), matching the known
ML offset vs the Anggraini manual scale.

## Files

- Catalog: `eqt/mt/results/mt_catalog.csv` (219 rows), `stress_result.json`.
- Figures: `eqt/figures/stress_map.png` (SHmax on topography + faults + beachballs),
  `eqt/figures/stress_stereonet.png` (σ-axis bootstrap clouds + P/T axes + R hist).
- Scripts: `eqt/scripts/{mt_extract,mt_run_one,mt_batch,mt_static_inputs,
  stress_inversion,plot_stress_map}.py`; patched clone `eqt/tools/BayesISOLA`.

## Follow-ups

- Reporting SHmax N22°E±20° / transtension is defensible; the full tensor is not.
- To tighten: add first-motion polarities (HASH) — dense short-period local array
  suits P-polarity mechanisms better than full-waveform MT here.
- Across-strike cross-section to tie beachball P/T axes to Opak dip.
