---
title: "Aftershock moment tensors & stress direction — Yogya 2006"
type: output
status: active
created: 2026-11-05
updated: 2026-11-10
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
- **CRITICAL QC (2026-11-09):** the catalog `nsta` is stations *extracted*, not
  *inverted*. BayesISOLA auto-selects stations by distance vs. frequency; with the
  high-fmax bands given to small events, often only the 1–3 nearest stations are
  inverted. **VR anti-correlates with stations actually used (r = −0.49)**: a
  "VR 0.77" small event (e.g. 20060609T083059) fit **one** station (TF19, 5 km) —
  trivially high VR, azimuthally unconstrained mechanism. ML↔nsta_used = +0.89
  (bigger events get lower fmax → more stations inverted → lower but meaningful VR).
  **Fix:** gate on `nsta_used` (parsed from each `out/log.txt`), NOT VR alone.
  `nsta_used` added as a catalog column.
- **Primary quality gate: nsta_used ≥ 5, VR ≥ 0.3, DC ≥ 50 → n = 36.**
- **Robust, plane-independent observable (36 multi-station MT):**
  - Extension (Shmin) **horizontal ~N106°E ± 31°**, plunge ~12–16° (T-axes cluster).
  - **SHmax ~N16°E** (NNE–SSW), ⟂ horizontal extension.
  - Regime **transtension**; the reliable multi-station set leans **normal-dominant**
    (σ1 plunge ~54°, 16/36 strike-slip). (Old VR≥0.5 gate looked more strike-slip,
    σ1 plunge 35° — biased by 1-station fits.)
  - **The extension direction is STABLE across all gates** (N100–112°E): VR≥0.5
    n=27 → N112±20; used≥5 VR≥0.3 n=36 → N106±31; used≥6 n=23 → N100±28. Only the
    spread tightens/loosens; the headline holds. The ±20° of the old gate was
    artificially tight.
- **Formal full stress tensor is NOT resolved:** bootstrap cones ~68° on all
  three axes, SHmax histogram ~flat, **R ≈ 0.28** (σ1≈σ2 near-degenerate),
  mean misfit 34°. P/T axes are intrinsic to each MT (plane-independent) → the
  extension direction survives; the full tensor does not.

## First-motion cross-check (HASH-style) — INDEPENDENT CONFIRMATION (2026-11-10)

Built a prototype first-motion focal-mechanism pipeline to test the MT extension
direction with a method that does **not** depend on the problematic 1-D waveform
fit. Scripts: `eqt/scripts/hash_core.py` (ray tracer + polarity + radiation-sign
grid search + P/T + quality), `hash_one.py` (one event + beachball/first-motion
figure), `hash_batch.py`, `hash_summary.py`.

- **Three ingredients, two built from scratch** (NLL stores neither): azimuth from
  NLL `.hyp` `SAzim`; **takeoff angle** from an own flat-layer ray tracer on
  `crustal.dat` (direct + head waves) — validated vs NLL `TTpred` to 0.01–0.14 s
  at 6–30 km; **P first-motion polarity** measured on the raw SAC verticals at the
  NLL pick (EQT records no `FM`), SNR-weighted.
- **Catalog:** `eqt/mt/results/hash_catalog.csv`, **n = 202** events with ≥6 usable
  polarities (quality A=90, B=98, C=14). Grid search s/d/r at 5°.
- **Per-event mechanisms are underdetermined** (7–11 polarities, network az-gaps;
  STDR often <0.5; the "best" nodal planes are tie-degenerate). So we stack the
  **acceptable-set-mean P/T axes** (plane-independent), not nodal planes.
- **Stacked extension (T-axis) = N102°E ± 27° (R=0.65), SHmax ≈ N12°E**, P axes
  steep → normal-dominant.
- **Agreement with MT (N106°E ± 31°, SHmax N16°E) within ~4°**, from two fully
  independent methods. First motions corroborate the WNW–ESE extension / NNE SHmax
  robustly, and with a slightly tighter spread.
- Figures: `eqt/figures/hash_summary.png` (rose + P/T stereonet + comparison),
  `eqt/figures/hash_firstmotion_20060820T005344.png` (9-station worked example,
  misfit 0/7).
- Caveat: polarity reliability is limited by EQT pick accuracy and emergent onsets;
  weights down-weight low-SNR picks but a manual/ML polarity re-pick would improve
  individual mechanisms. The stacked direction is robust regardless.

## Tectonic reading

WNW–ESE horizontal extension / NNE SHmax in transtension is consistent with the
**Opak fault** (NE–SW, oblique left-lateral; PUSGEN dip 50°,
[[wiki/sources/pusgen-2024-shallow-crustal-faults]]) as the controlling structure;
the aftershock cloud (beachballs) sits ~5–15 km E of the Opak surface trace.
Mw from MT runs ~0.9–1.0 unit above our ML (ML 2.84 → Mw 2.98), matching the known
ML offset vs the Anggraini manual scale.

## Files

- Catalog: `eqt/mt/results/mt_catalog.csv` (219 rows, incl. `nsta_used`),
  `stress_result.json`.
- Figures: `eqt/figures/stress_map.png` (SHmax on topography + faults + beachballs),
  `eqt/figures/stress_stereonet.png` (σ-axis bootstrap clouds + P/T axes + R hist),
  `eqt/figures/mt_waveform_fit_*_20060820T005344.png` (obs vs synthetic, 9-station
  illustration event ML2.84 VR0.57), `mt_beachball_20060820T005344.png`.
- Scripts: `eqt/scripts/{mt_extract,mt_run_one,mt_batch,mt_static_inputs,
  stress_inversion,plot_stress_map,mt_waveform_plot}.py`; patched clone
  `eqt/tools/BayesISOLA`.

## Follow-ups

- Reporting SHmax N12–16°E / WNW–ESE extension / transtension is defensible and now
  **confirmed by two independent methods** (MT + first motions); the full stress
  tensor is still not resolved.
- DONE: first-motion (HASH-style) cross-check — N102°E ± 27° extension, agrees with
  MT within ~4°. Next refinement would be ML/manual polarity re-picking to tighten
  per-event mechanisms and a formal stress inversion on the first-motion mechanisms.
- Across-strike cross-section to tie beachball P/T axes to Opak dip.
- `updated: 2026-11-10`.
