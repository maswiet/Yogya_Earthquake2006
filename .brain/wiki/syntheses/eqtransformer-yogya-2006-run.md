---
title: "Running EQTransformer on the Yogya 2006 aftershock array"
type: synthesis
status: active
created: 2026-07-03
updated: 2026-07-03
sources:
  - "[[wiki/sources/earthquake-transformer-mousavi-2020]]"
  - "[[wiki/sources/eqtransformer-github-repo]]"
  - "[[wiki/sources/yogya-2006-aftershock-edl-dataset]]"
tags:
  - pipeline
  - yogya-2006
  - workflow
  - active-work
---

## Goal

Reprocess the 2006 Yogyakarta aftershock recordings with the EQTransformer model
to detect/locate **more aftershocks** than the reference catalog — the
Mousavi-style result shown on the 2000 Tottori sequence
([[wiki/claims/eqtransformer-doubles-detections-tottori]]).

## Working directory

`/Users/maswiet/Work/Yogya_Earthquake2006/eqt/` — `scripts/`, `config/`,
`logs/`, `pilot/`. Data on external volume (see
[[wiki/sources/yogya-2006-aftershock-edl-dataset]]).

## Environment (arm64 Mac)

- The original `EQTransformer` PyPI package pins TF 2.5 / Keras 2.3 and does not
  install cleanly on arm64 / Python 3.12. **Decision: run the same EQTransformer
  model + original STEAD weights via SeisBench** (PyTorch), which installs
  cleanly and is maintained. Results are equivalent.
- conda env **`eqt`**: `seisbench 0.11.7`, `obspy 1.5`, `torch`, numpy 2.
  `EQTransformer.from_pretrained("original")`.
- conda env **`assoc`** (separate, numpy<2): `pyocto` for association + location.
  Kept separate because PyOcto needs numpy<2 and would break SeisBench.

## Pipeline

1. **extract_coords.py** → per-station median GPS (Yogya box filter) →
   `config/stations_raw.json`. DONE.
2. **preprocess.py** → merge 30-min `.pri0/1/2` (2006 only) → relabel p0/1/2 to
   `EHZ/EHN/EHE`, merge gaps, resample 200→100 Hz → per-station-day miniSEED
   under `pilot/mseed/<STA>/`. DONE + validated.
3. **run_eqt.py** → EQTransformer (`original`) detection + P/S picks → picks CSV.
   Thresholds P=S=det=0.3. ~2.5 min/station-day on CPU.
4. **Association + location** (PyOcto, 1-D velocity model) → events. TODO.
5. **Compare** event count / completeness vs reference catalog. TODO.

## Decisions locked (user)

- Components: **p0=Z, p1=N, p2=E**.
- **Pilot first**: stable stations **TF12/14/16/18/19**, days **154–159**
  (early June, peak aftershocks), then scale to 12 stations × ~90 days.
- Take pipeline **through location + compare**.

## Status (2026-07-03)

- Pilot preprocessing: 30 station-day mseed files written.
- Smoke test (TF14/155): **259 P + 271 S picks, 290 detections** in 152 s — picker
  works; pick density consistent with active aftershock sequence (association
  will reject noise).
- Pilot picking (30 files) running.

## Needed external inputs (for step 5)

- **1-D velocity model** for the Yogya/Central Java region (P & S). `needs-review`
  — ask user or adopt a published Central Java model.
- **Reference aftershock catalog** for the same period to compare against
  (event count + magnitude of completeness). `needs-review` — ask user.

## Open Questions / risks

- CPU-only throughput: full run ≈ 45 h single-thread — consider MPS/parallelism.
- Relocated stations (TF07/10/15/17, maybe TF09) need per-period coordinates.
- Threshold choice affects event count; keep probabilities to re-filter.
- Fair comparison requires matching region/time/magnitude cuts to the reference.

## Links

- [[wiki/entities/eqtransformer-model]]
- [[wiki/entities/yogya-2006-temp-aftershock-network]]
- [[wiki/concepts/phase-association]]
- [[wiki/questions/applying-eqtransformer-to-yogya-2006]]
