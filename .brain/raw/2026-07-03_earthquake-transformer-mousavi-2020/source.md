# Raw Capture — Earthquake Transformer (Mousavi et al., 2020)

## Provenance

- Requested URL: https://www.nature.com/articles/s41467-020-17591-w
- Final access: Nature URL redirects through an IdP cookie-auth handshake
  (`idp.nature.com/authorize`) that automated fetch cannot follow. Bibliographic
  data and the verbatim abstract were captured from open mirrors:
  - https://ideas.repec.org/a/nat/natcom/v11y2020i1d10.1038_s41467-020-17591-w.html (abstract, citation)
  - https://eqtransformer.readthedocs.io/ (architecture note)
  - Web search (title/authors/citation cross-check)
- Retrieval date: 2026-07-03
- Access limitations: Full PDF body (methods, figures, tables) NOT captured;
  only the abstract is verbatim. Architecture/training-data specifics beyond the
  abstract are marked needs-review downstream.

## Bibliographic Record

- Title: Earthquake transformer—an attentive deep-learning model for simultaneous
  earthquake detection and phase picking
- Authors: S. Mostafa Mousavi; William L. Ellsworth; Weiqiang Zhu;
  Lindsay Y. Chuang; Gregory C. Beroza
- Journal: Nature Communications
- Volume 11, Issue 1, Article number 3952 (pages 1–12)
- Year: 2020
- DOI: 10.1038/s41467-020-17591-w

## Abstract (verbatim)

"Earthquake signal detection and seismic phase picking are challenging tasks in
the processing of noisy data and the monitoring of microearthquakes. Here we
present a global deep-learning model for simultaneous earthquake detection and
phase picking. Performing these two related tasks in tandem improves model
performance in each individual task by combining information in phases and in the
full waveform of earthquake signals by using a hierarchical attention mechanism.
We show that our model outperforms previous deep-learning and traditional
phase-picking and detection algorithms. Applying our model to 5 weeks of
continuous data recorded during 2000 Tottori earthquakes in Japan, we were able
to detect and locate two times more earthquakes using only a portion (less than
1/3) of seismic stations. Our model picks P and S phases with precision close to
manual picks by human analysts; however, its high efficiency and higher
sensitivity can result in detecting and characterizing more and smaller events."

## Additional captured note

- The model (widely known as "EQTransformer") uses a hierarchical architecture
  with an attention mechanism, designed for earthquake signals, and was trained
  on global seismic data (per eqtransformer.readthedocs.io). Layer-level details
  and training-set name/size were not captured from a primary source.
