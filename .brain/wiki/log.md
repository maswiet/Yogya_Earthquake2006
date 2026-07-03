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

## [2026-07-03] ingest | EQTransformer GitHub Repository (smousavi05/EQTransformer)

- **Trigger:** `/ingest https://github.com/smousavi05/EQTransformer`
- **Source access:** README-level metadata fetched from the repo landing page.
  No repository code executed. Full module API/source not captured.
- **Files changed (created):**
  - `raw/2026-07-03_eqtransformer-github/source.md`
  - `wiki/sources/eqtransformer-github-repo.md`
  - `wiki/concepts/phase-association.md`
- **Files changed (updated):**
  - `wiki/entities/eqtransformer-model.md` (added software/package facts + source link)
  - `wiki/concepts/earthquake-signal-detection.md` (two-model false-neg/false-pos tradeoff; Conflicts/Updates)
  - `wiki/questions/applying-eqtransformer-to-yogya-2006.md` (tooling-available note)
  - `wiki/index.md`, `wiki/log.md`
- **Key result:** Runnable implementation of the already-ingested model.
  MIT-licensed, `pip install EQTransformer`, TensorFlow+ObsPy, ships two
  pretrained models (Original = min false negatives; Conservative = min false
  positives). Confirms a low-barrier path to a Yogya 2006 trial run.
- **Follow-ups:** Locate Yogya 2006 continuous waveform data + choose a location
  workflow (velocity model + associator) to turn picks into hypocenters.

## [2026-07-03] ingest | Yogya 2006 aftershock EDL dataset + EQT pipeline kickoff

- **Trigger:** User pointed to raw dataset `/Volumes/Untitled 1/DATA-GFZ-Gempa-
  JOgja-tahap-2` and asked to run EQTransformer to find more aftershocks
  (Mousavi/Tottori-style).
- **Files created:** `wiki/sources/yogya-2006-aftershock-edl-dataset.md`,
  `wiki/entities/yogya-2006-temp-aftershock-network.md`,
  `wiki/syntheses/eqtransformer-yogya-2006-run.md`.
- **Files updated:** `wiki/questions/applying-eqtransformer-to-yogya-2006.md`
  (status → in progress), `wiki/index.md`.
- **Key results / durable facts:**
  - 12-station temp EDL array (tf3007, tf3009–tf3019; no tf3008). Components
    p0=Z/p1=N/p2=E, 200 Hz, miniSEED in 30-min `.pri0/1/2` segments.
  - Data hygiene: folders mix 2004–05 Germany pre-deployment data; filter by
    filename year `06` + Yogya GPS box. Several stations relocated mid-campaign
    (TF07/10/15/17, maybe TF09) → need per-period coords.
  - Environment: original EQTransformer won't install on arm64; running the same
    model + original STEAD weights via **SeisBench** (env `eqt`); **PyOcto** in
    separate env `assoc` for association/location.
  - Pipeline built + validated; pilot (TF12/14/16/18/19, days 154–159) picking
    running. Smoke test TF14/155 = 259 P + 271 S picks.
- **Follow-ups:** Need 1-D velocity model + reference aftershock catalog for the
  locate-and-compare step (asking user).

## [2026-07-03] maintenance | EQT Yogya pilot validated end-to-end

- **Trigger:** Pilot run of the full pipeline on 5 stable stations × days 154–159.
- **Files changed:** `wiki/syntheses/eqtransformer-yogya-2006-run.md` (status +
  scale-up plan); scripts under `eqt/scripts/` (run_eqt, associate_locate,
  plot_pilot); outputs under `eqt/pilot/`.
- **Key result:** Pipeline works. MPS ~19× faster than CPU. Pilot: 36,939 picks →
  2,538 located events (808 with ≥4 stations); epicenters cluster in the array,
  depths 5–15 km (median 11), Omori-like decay. Velocity model = published
  Central Java 1-D (default, swappable). Reference-catalog comparison deferred
  (user decides later).
- **Follow-ups:** Build streaming preprocess→pick driver; run all 12 stations ×
  ~90 days on MPS (~2 h); relocation-aware coords for location; obtain reference
  catalog for the '2× more events' comparison.
