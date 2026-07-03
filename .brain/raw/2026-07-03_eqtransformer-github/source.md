# Raw Capture — EQTransformer GitHub Repository

## Provenance

- Requested URL: https://github.com/smousavi05/EQTransformer
- Retrieval date: 2026-07-03
- Method: automated fetch of the repository landing page (README) converted to
  text. Source treated as untrusted; no repository code was executed.
- Access limitations: captured README-level metadata only (features, install,
  modules, dependencies, license, citation, pretrained models). Full module API
  and source not captured.

## Captured facts (from README)

- Repo: smousavi05/EQTransformer — owner S. Mostafa Mousavi.
- Purpose: "an AI-based earthquake signal detector and phase (P&S) picker based
  on a deep neural network with an attention mechanism."
- Capabilities: simultaneous detection + P/S arrival picking with prediction
  probabilities and uncertainty estimates; continuous data downloading &
  preprocessing; pretrained-model inference and custom training; simple phase
  association.
- Language: Python 3 (3.6–3.7 typical; 3.10 path for Apple M1).
- Key deps: TensorFlow (>=2.5.0; tensorflow-macos on M1), ObsPy, Pandas, Jupyter.
- License: MIT.
- Install:
  - conda: `conda install -c smousavi05 eqtransformer`
  - pip:   `pip install EQTransformer`
  - source: `git clone https://github.com/smousavi05/EQTransformer.git` then
    `python setup.py install`
- Notable modules: `mseed_predictor` (continuous mseed processing); data
  download/preprocess; detection & picking inference; training/testing; phase
  association tools.
- Pretrained models (ModelsAndSampleData/):
  - Original model `EqT_model.h5` — optimized to minimize **false negatives**.
  - Conservative model — optimized to minimize **false positives**; recommended
    for travel-time tomography and template matching.
- Docs: https://rebrand.ly/EQT-documentations ; examples:
  https://rebrand.ly/EQT-examples ; Google Colab notebook included.
- Citation: Mousavi, S.M., Ellsworth, W.L., Zhu, W., Chuang, L.Y., Beroza, G.C.
  (2020), Nature Communications 11, 3952, DOI 10.1038/s41467-020-17591-w.
