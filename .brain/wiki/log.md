<!-- afk-research:managed v1 -->
# Log

Append-only chronological record for the `.brain` second brain. Use headings in this format:

```markdown
## [YYYY-MM-DD] type | Title
```

Allowed types: `setup`, `ingest`, `query`, `lint`, `maintenance`, `export`, `import`, `schema`.

## [2026-07-03] ingest | Earthquake Transformer (Mousavi et al., 2020)

- **Trigger:** `/ingest https://www.nature.com/articles/s41467-020-17591-w`
- **Source access:** Nature URL redirected through an IdP cookie handshake that
  automated fetch cannot follow; verbatim abstract + citation captured from open
  mirrors (ideas.repec.org, eqtransformer.readthedocs.io) and cross-checked via
  web search. Paper body (methods/figures) NOT captured.
- **Files changed (created):**
  - `raw/2026-07-03_earthquake-transformer-mousavi-2020/source.md`
  - `wiki/sources/earthquake-transformer-mousavi-2020.md`
  - `wiki/concepts/seismic-phase-picking.md`
  - `wiki/concepts/earthquake-signal-detection.md`
  - `wiki/concepts/deep-learning-seismic-processing.md`
  - `wiki/entities/eqtransformer-model.md`
  - `wiki/entities/2000-tottori-earthquake-sequence.md`
  - `wiki/claims/eqtransformer-doubles-detections-tottori.md`
  - `wiki/questions/applying-eqtransformer-to-yogya-2006.md`
- **Files changed (updated):** `wiki/index.md`, `wiki/log.md`
- **Key result:** First seismic-ML source ingested. EQTransformer = global
  deep-learning model for simultaneous detection + P/S phase picking; Tottori
  2000 field test found ~2× more earthquakes using <1/3 of stations.
- **Follow-ups:** (1) Read full paper to capture architecture/training-set
  details and per-station Tottori metrics (`needs-review`). (2) Assess Yogya 2006
  waveform data availability and EQTransformer generalization —
  [[wiki/questions/applying-eqtransformer-to-yogya-2006]].
