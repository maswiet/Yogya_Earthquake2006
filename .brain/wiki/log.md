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

## [2026-07-03] export | Yogya 2006 EQT aftershock catalog v1 produced

- **Trigger:** Full run completed; association + finalization.
- **Files changed:** `wiki/outputs/yogya-2006-eqt-catalog.md` (new);
  `wiki/questions/applying-eqtransformer-to-yogya-2006.md` (answered),
  `wiki/index.md`. Deliverables under `eqt/full/` (catalog_eqt.csv,
  events_full.csv, picks_full.csv, catalog_summary.png).
- **Key result:** 960 station-days → 397k picks → 17,150 associated →
  **16,876 located aftershocks** (artifact-filtered), Jun 3–Aug 29 2006, median
  depth 12.9 km, Omori decay, tight fault-zone cluster. Pipeline fully validated.
- **Follow-ups:** reference-catalog comparison (deferred by user); magnitudes need
  instrument response; optional relocation with NonLinLoc/HypoDD + reference
  velocity model.

## [2026-07-04] maintenance | NonLinLoc relocation of the EQT catalog

- **Trigger:** User asked to proceed to NonLinLoc relocation.
- **Files changed:** `wiki/outputs/yogya-2006-eqt-catalog.md` (NLLoc section).
  Deliverables `eqt/full/catalog_nll*.csv`, `eqt/full/nll_compare.png`;
  `eqt/nll/` control files + grids; scripts `gen_nll.py`, `parse_nll.py`.
- **Key result:** Built NonLinLoc from source (arm64); 3D grids from Central Java
  model; relocated all 16,876 events in ~8 min. RMS median 0.079 s, errH 1.1 km;
  12,844 well-constrained. Resolves a NE–SW fault lineament; agrees with PyOcto
  to 1.4 km median. Gotchas: needed 3D grids (2D setup gave 3 s RMS) and
  LOCGRID SAVE (NO_SAVE suppressed .hyp); pandas-3 datetime[us] scaling in parse.
- **Follow-ups:** reference-catalog comparison; magnitudes; optional HypoDD.

## [2026-07-04] export | Data-driven 1-D velocity (Wadati + VELEST)

- **Trigger:** User asked to derive 1-D velocity from distance–arrival-time via
  travel-time method (e.g. VELEST).
- **Files:** `wiki/outputs/yogya-2006-1d-velocity-model.md` (new); `wiki/index.md`.
  Deliverables `eqt/full/{arrivals.csv,ttime_analysis.png}`, `eqt/velest/*`;
  scripts `parse_arrivals.py`, `gen_velest.py`, `parse_velest.py`.
- **Key result:** Wadati Vp/Vs=1.735; apparent Vp=5.93/Vs=3.39 km/s. VELEST
  minimum-1D (built from source, arm64) converged RMS 0.358→0.152 s: Vp 4.9→6.1
  km/s over 0–8 km, Vp/Vs~1.74, geologically sensible station corrections
  (TF18 +1.47 s basin … TF11b −0.55 s). Gotcha: VELEST longitude is West-positive
  (enter East as negative olon).
- **Follow-ups:** relocate with data-driven model + station corrections; resolve
  >12 km layers needs larger-offset data.

## [2026-07-04] export | ML magnitudes + relief map + VELEST relocation

- **Trigger:** User provided XN station metadata + L4-3D instrument response;
  asked for improved relief map (labels) and magnitudes.
- **Files:** `wiki/entities/yogya-2006-temp-aftershock-network.md` (XN names +
  response), `wiki/outputs/yogya-2006-eqt-catalog.md` (Magnitudes),
  `wiki/outputs/yogya-2006-1d-velocity-model.md` (VELEST reloc). Figures under
  `eqt/figures/` (all committed to GitHub): aftershock_distribution,
  aftershock_relief_map, velest_1d_model, velest_relocation_compare,
  magnitude_gutenberg_richter. Scripts: build_amplitudes, compute_magnitudes,
  plot_relief, gen_velest_reloc, compare_reloc.
- **Key results:** VELEST-model relocation (depths bias deep +5.4 km from deep-
  layer artifact); PyGMT relief map (light relief, blue sea, XN labels);
  **ML for 16,876 events (median 0.07, max 3.65); Gutenberg-Richter Mc=0.20,
  b=0.96±0.01.**
- **Follow-ups:** local ML calibration for Java; refine VELEST deep layers;
  reference-catalog comparison; HypoDD.

## [2026-07-04] maintenance | Refined VELEST model (geology) + de-biased relocation

- **Trigger:** User asked to fix VELEST shallow + deep layers using Opak-fault
  geology (W = Merapi sediment >700 m slow; E = limestone fast).
- **Files:** `wiki/outputs/yogya-2006-1d-velocity-model.md`. Figures (committed):
  velest_1d_model, station_corrections_opak, aftershock_relief_map,
  velest_relocation_compare. Scripts: gen_velest (fixed shallow sediment @0.7 km +
  fixed deep + elevations + iuseelev), interpret_stacorr.
- **Key result:** Refined model — shallow sediment fixed (2.5→4.3 to 0.7 km), deep
  fixed realistic (6.55→8.0, artifact gone), mid-crust inverted; RMS 0.142.
  Station corrections show W-sediment slow (+0.5..0.8 s BUM/PEL) vs E-limestone
  fast (−0.6..−1.25 s KRI/TF11/KEM). **Re-relocation de-biased depths: 15.0 →
  10.5 km**, RMS 0.132→0.113, well-constrained 8.9k→12.4k.
- **Follow-ups:** HypoDD; local ML scale.

## [2026-07-14] query | Pangandaran Mw7.7 triggering test on Yogya aftershocks

- **Trigger:** User asked whether the 17 Jul 2006 Pangandaran Mw7.7 tsunami
  earthquake perturbed the Yogya aftershock population (deep stress-triggering
  analysis). USGS event usp000ensm.
- **Files:** `wiki/outputs/yogya-2006-pangandaran-triggering.md` (new),
  `wiki/index.md`. Code (committed to repo): `eqt/scripts/analyze_pangandaran.py`,
  `eqt/figures/pangandaran_rate.png`.
- **Key result:** NO detectable static or dynamic triggering. β-vs-ML-threshold
  shows the apparent post-17-Jul small-event deficit vanishes above Mc
  (β≈−6 at ML≥0 → ≈0 at ML≥1.0) = completeness/coda-masking artifact, not real.
  Robust ML≥1.0 rate unchanged (β≈+0.1). Static Coulomb ~3e-12 bar (negligible);
  dynamic ~0.86 bar (above threshold yet no observed rate change).
- **Follow-ups:** confirm mid-July station up-time (completeness) from log;
  regional PGV record for Yogya if available.
