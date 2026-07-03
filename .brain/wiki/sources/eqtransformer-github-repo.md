---
title: "EQTransformer GitHub Repository (smousavi05/EQTransformer)"
type: source
status: active
created: 2026-07-03
updated: 2026-07-03
source_author: "S. Mostafa Mousavi (smousavi05)"
source_url: "https://github.com/smousavi05/EQTransformer"
raw: "[[raw/2026-07-03_eqtransformer-github/source]]"
tags:
  - software
  - deep-learning
  - seismology
  - tool
---

## Summary

Official open-source (MIT) Python package implementing **EQTransformer** — the
AI earthquake detector and P/S phase picker from
[[wiki/sources/earthquake-transformer-mousavi-2020]]. Provides an end-to-end
workflow: download continuous data → detect events → pick P/S arrivals →
(simple) phase association, using **pretrained models** so it can be run without
retraining.

## Key Points

- **Install:** `pip install EQTransformer` or `conda install -c smousavi05 eqtransformer`. License **MIT**.
- **Stack:** Python 3 + **TensorFlow** (>=2.5.0), **ObsPy** (seismic I/O), Pandas, Jupyter. Apple-M1 path via `tensorflow-macos` (Python 3.10).
- **Modules:** continuous-data download/preprocess; `mseed_predictor` for mseed streams; detection + picking inference; training/testing; phase association.
- **Two pretrained models** ship in `ModelsAndSampleData/`:
  - *Original* (`EqT_model.h5`) — tuned to minimize **false negatives** (catch more).
  - *Conservative* — tuned to minimize **false positives**; recommended for **travel-time tomography** and **template matching**.
- Outputs include prediction probabilities and uncertainty estimates.
- Docs: https://rebrand.ly/EQT-documentations · Examples: https://rebrand.ly/EQT-examples · Google Colab notebook included.

## Extracted Concepts

- [[wiki/concepts/phase-association]] (new)
- [[wiki/concepts/earthquake-signal-detection]] (false-negative vs false-positive model choice)
- [[wiki/concepts/seismic-phase-picking]]
- [[wiki/concepts/deep-learning-seismic-processing]]

## Extracted Entities

- [[wiki/entities/eqtransformer-model]] (updated with software/package facts)

## Claims

- [[wiki/claims/eqtransformer-doubles-detections-tottori]]

## Conflicts / Updates

- Complements [[wiki/sources/earthquake-transformer-mousavi-2020]] — the paper
  describes the method; this repo is the runnable implementation. No conflicts.

## Open Questions

- [[wiki/questions/applying-eqtransformer-to-yogya-2006]] — the pip-installable
  tool + pretrained models lower the barrier to a Yogya 2006 trial run.
