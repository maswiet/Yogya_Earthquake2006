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

## [2026-07-20] maintenance | Magnitude pipeline: 3 measurement bugs fixed + station corrections (pre-SRL)

- **Trigger:** Preparing an SRL submission; user asked to see example waveforms
  of the sub-ML-1 "tail". Plotting them exposed that the smallest catalog events
  were strong earthquakes (SNR 100–1000), not noise — starting a bug hunt.
- **Files (repo):** `eqt/scripts/build_amplitudes.py` (per-pick windowed
  deconvolution + 1–20 Hz bandpass + volume auto-detect),
  `compute_magnitudes.py` (station corrections + `mbs_mc` b-stability Mc),
  `screen_catalog.py` (new), `station_ml_corrections.py` (new),
  `plot_noise_floor.py` (new), `plot_example_waveforms.py` (new),
  `plot_decay.py` (Mc from mbs). New data: `eqt/full/catalog_quality.csv`,
  `eqt/config/station_ml_corrections.json`. Figures regenerated.
- **Files (brain):** `wiki/outputs/yogya-2006-eqt-catalog.md` (Magnitudes section
  rewritten), `wiki/outputs/yogya-2006-pangandaran-triggering.md` (re-run note),
  `wiki/index.md`, `wiki/log.md`.
- **Key results:**
  - **Bug 1 — simulate() taper on 24-h traces:** 5% cosine taper = 72 min/end;
    events near a day boundary suppressed up to ~200× (+2.3 ML). 39% of ML<−1
    events were within ±10 min of midnight (28× enrichment). Fixed by
    deconvolving a padded per-pick window; midnight enrichment → 0.0%.
  - **Bug 2 — 50 Hz mains hum** (TF10b, 200 Hz): unfiltered peak set weak-event
    amplitudes (+0.27 ML @ ML<0). Fixed by 1–20 Hz bandpass.
  - **Bug 3 — max-curvature underestimates Mc by ~0.6.** b climbs 0.77→0.91 with
    cut-off ⇒ residual incompleteness. Replaced with b-stability (MBS).
  - Per-station ML corrections (−0.37 TF16 … +0.35 TF17); scatter 0.185→0.128.
    **r=−0.03 vs VELEST P corrections** (amplitude vs travel-time site response
    decoupled).
  - **Final: Mc=+0.50, b=0.89±0.02 (N=2258); ML −1.81..3.55.** Old Mc≈0/b=0.88
    superseded. Quality screen: 11,790/16,876 pass (69.9%). Sub-Mc events are
    real (SNR 12–34) — sample incomplete, detections genuine.
  - **Ramdhan 2025 comparison NOT yet defensible:** ratio is 5.51× (all) but
    0.83× at our Mc=+0.50 vs their 2,141. Need their FMD/Mc — critical path.
- **Follow-ups:** obtain Ramdhan Mc/catalog; investigate FMD roll-off >ML1.5
  (L4-3D clipping ⇒ max ML is a lower bound); resolve TF16 −0.37 gain anomaly;
  apply station corrections consistently in any future ML recompute.

## [2026-07-20] query | Resolved FMD roll-off (no clipping) + TF16 anomaly (site effect)

- **Trigger:** Two open magnitude threads before SRL, after Ramdhan data proved
  hard to obtain: (a) does the L4-3D clip the largest events, and (b) is TF16's
  −0.37 deficit an instrument fault?
- **Files (repo):** `eqt/scripts/check_clipping.py` (new),
  `eqt/scripts/diagnose_tf16.py` (new); figures `clipping_check.png`,
  `tf16_diagnostic.png`.
- **Files (brain):** `wiki/outputs/yogya-2006-eqt-catalog.md` (two caveats
  retracted/resolved), `wiki/index.md`, `wiki/log.md`.
- **Key results:**
  - **No clipping.** Digitiser rail = 2^23 counts; largest event (ML 3.55) peaks
    at 16% of rail (median of 12 largest = 8%), no flat-topping. Band-widening
    changes large-event amplitude <3%. ⇒ **max ML 3.55 is real**, and the GR
    roll-off (obs/GR 1.0→0.09 from Mc to ML 3.0) is a finite-catalogue / real
    max-magnitude effect. Retracts the "clipping ⇒ lower bound" caveat.
  - **TF16 is a site/path effect, not a gain fault.** Deficit grows with distance
    (−0.21<10 km → −0.48 at 20–30 km), stable in time, both horizontals healthy
    (N/E 0.89). Hard limestone E of Opak; decoupled from VELEST P-corr (different
    sampling depth). Absorbed by station correction. Retracts "instrument-gain
    issue".
  - Minor: 20 Hz upper corner mildly under-measures the smallest events
    (ML~−0.4 lose 47–65% vs 1–45 Hz); kept as standard.
- **Follow-ups:** Ramdhan FMD/Mc remains the only blocker for the "N× events"
  claim; consider a modified-GR (tapered) fit to report a corner magnitude.

## [2026-07-20] planning | Reframe SRL paper: methodology + open data product

- **Trigger:** Ramdhan et al. 2025 FMD/Mc hard to obtain; user chose to reframe
  away from the "N× more events" claim rather than block on it.
- **Files (brain):** `wiki/outputs/yogya-2006-srl-paper-plan.md` (new),
  `wiki/index.md`, `wiki/log.md`.
- **Key result:** New framing = 3 legs — (1) reproducible ML workflow documenting
  three deep-learning-reprocessing pitfalls (taper on day-long traces, mains hum,
  MaxC-underestimates-Mc); (2) open, quality-screened, GrowClust-relocated
  catalogue with characterised completeness (Mc=+0.50, b=0.89); (3) independent
  structural corroboration of the Opak–Ngalang source (N57°E ≈ GCMT NP2) from a
  ~6× denser catalogue. Ramdhan positioned as complementary reference, not
  competitor — needs only their fault interpretation, not their magnitudes.
  Draft abstract + section→asset outline + pre-submission checklist recorded.
- **Follow-ups:** pick validation vs manual ground truth (top reviewer risk);
  Zenodo DOI release; confirm SRL article type; optional tapered-GR corner mag.

## [2026-07-20] ingest | Diambama et al. 2019 (GJI) — our group's Ngalang Fault priority

- **Trigger:** User provided their group's own 2019 GJI paper
  (`/Users/maswiet/Documents/Manuscript/GJI_Wiwit_.pdf`) noting "our group first
  proposed the Ngalang fault."
- **Files (brain):** `wiki/sources/diambama-2019-yogya-tomography.md` (new),
  `wiki/outputs/yogya-2006-srl-paper-plan.md` (priority narrative),
  `wiki/index.md`, `wiki/log.md`.
- **Key result:** Diambama, Anggraini, Nukman, Lühr & Suryanto (2019, GJI 216,
  439–452) explicitly associated the 2006 source with the **Ngalang Fault**
  (§5.2, Conclusion, Fig. 14) — **6 years before Ramdhan et al. 2025**. Same GFZ
  L4-3D/EDL deployment, same aftershock dataset (Anggraini 2013). Tomography
  resolved the fault to only 5–7 km; the new EQTransformer catalogue images the
  same structure to sub-500 m over 0–18 km. Paper reframed as the resolution
  culmination of our own program, with Ramdhan as independent corroboration.
  Two unlocks: (1) priority is ours; (2) **Anggraini 2013 manual picks (3,769 P
  + 3,407 S, 588 events) = in-house pick-validation ground truth.** Max ML 3.55
  matches our catalogue exactly (magnitude cross-check).
- **Follow-ups:** obtain Anggraini 2013 pick files for EQTransformer
  precision/recall validation; add Diambama 2019 + Anggraini 2013 to references.

## [2026-07-20] query | Compare EQT catalogue vs Anggraini 2013 manual catalogue

- **Trigger:** User provided the Anggraini dissertation catalogue
  (`~/Downloads/Bantul2006_Aftershock_Catalogue.xlsx`).
- **Files:** `eqt/scripts/compare_anggraini.py` + `figures/compare_anggraini.png`
  (repo); `wiki/sources/diambama-2019-yogya-tomography.md`,
  `wiki/outputs/{yogya-2006-eqt-catalog,yogya-2006-srl-paper-plan}.md`,
  `wiki/index.md`, `wiki/log.md` (brain).
- **Key result:** Catalogue = 590 events (Absolute) / 524 (DD), **3–7 Jun 2006
  only**, columns date+lon+lat+depth+ML — **no origin times, no picks** (so no
  pick-level precision/recall; event pairing impossible, NN ML r=0.06). In the
  3–7 Jun window: **strong detection recovery** (576/590 co-located same-day) +
  spatial/depth agreement (same NE–SW Ngalang structure) = event-level
  validation. **⚠️ Systematic ML offset: our absolute ML ~0.9–1.25 units BELOW
  Anggraini's manual scale** (matched-median +0.86; FMD-tail +1.25). b-value is
  offset-invariant so b=0.89 stands; absolute Mc/magnitudes need a local tie.
- **Follow-ups:** determine Anggraini's magnitude method (ML vs duration Md) to
  explain the offset; obtain phase-pick files for true pick validation; decide
  recalibrate-vs-report for the absolute scale.

## [2026-08-14] ingest | Zenodo LOTOS code release (Koulakov 2021)

- **Trigger:** User ran `/ingest https://zenodo.org/records/5338981`.
- **Files (brain):** `raw/2026-08-14_zenodo-lotos-koulakov/provenance.md` (new),
  `wiki/sources/zenodo-lotos-koulakov-2021.md` (new),
  `wiki/entities/lotos-tomography-code.md` (new),
  `wiki/sources/diambama-2019-yogya-tomography.md` (backlinks),
  `wiki/index.md`, `wiki/log.md`.
- **Key result:** Zenodo record = open (CC BY 4.0) LOTOS tomography code release by
  Ivan Koulakov (DOI 10.5281/zenodo.5338981, 2021), the **Gakkel Ridge (Arctic)**
  variant with a water layer — NOT Yogyakarta data. Relevant because LOTOS
  (Koulakov 2009) is the SAME engine (LOTOS-13) our group used in Diambama et al.
  2019. Captured metadata/provenance only (29.6 MB zip not downloaded); created a
  source page + a LOTOS entity, cross-linked to the Diambama source and the 1-D
  velocity-model output.
- **Follow-ups (`needs-review`):** confirm which LOTOS version Diambama 2019 ran;
  decide whether a re-run with a current LOTOS is worthwhile for the SRL paper's
  velocity-model section.

## [2026-08-14] query | Pick-level validation vs Anggraini 2013 manual picks

- **Trigger:** User provided Anggraini's manual pick files (phase_300.dat +
  station.dat + stat_ft.dat) — the pick-level ground truth.
- **Files:** `eqt/scripts/validate_picks.py` + `figures/pick_validation.png`
  (repo); `raw/2026-08-14_anggraini-picks/` (immutable copies + provenance),
  `wiki/sources/diambama-2019-yogya-tomography.md`,
  `wiki/outputs/{yogya-2006-eqt-catalog,yogya-2006-srl-paper-plan}.md`,
  `wiki/index.md`, `wiki/log.md` (brain).
- **Key result — EQTransformer passes cleanly:** phase_300.dat = 588 events,
  3776 P + 3414 S, with full h:m:s origin times (which the xlsx lacked).
  Matched events by origin time (±5 s): **528/588 (90%) recovery**. Pick
  precision (per-event demeaned): **P MAD 0.02 s, S MAD 0.06 s** (98% P / 83% S
  within 0.3 s) over 4,761 matched picks on 6 co-located core stations
  (WON/TF12, PEL/TF13, RAT/TF14, WAN/TF16, PAL/TF19, BUM/TF18). **Event-matched
  ML: offset +0.41, r=0.95** — clean constant tie to the local scale; SUPERSEDES
  the earlier distributional 0.86–1.25 (inflated by completeness). b-value
  offset-invariant, unchanged.
- **Follow-ups:** manual stations NGL/TRI (964+528 picks) are outside our XN set
  — decide whether to add; spot-check the 10% unmatched events; decide whether to
  report ML raw and/or locally-tied (+0.41) in the paper.

## [2026-08-14] query | Tomographic resolution feasibility from the EQT catalogue

- **Trigger:** User asked how fine a velocity model the new catalogue could
  resolve if used for tomography.
- **Files:** `eqt/scripts/tomo_resolution.py` + `figures/tomo_resolution.png`
  (repo). Straight source->station ray coverage proxy (hit count + 45-deg
  azimuth-sector diversity per cell); 588-event subset (Diambama size) run with
  identical geometry to isolate the catalogue-size gain.
- **Key result:** 142,032 rays (71k P + 71k S) from 11,790 QC events, 17 station
  positions = **20x Diambama 2019's 7,176 rays**. Resolvable footprint
  (hit>=20 & >=3/8 azimuth sectors): **5 km 87%, 3 km 69%, 2 km 62%, 1.5 km 56%**
  (full) vs 58/48/39/36% for the 588-event subset. Azimuth diversity reaches
  **7-8/8 sectors in the source-zone core**, dropping to 1-3 at the margins and
  below ~15 km. **Verdict: ~2-3 km resolution achievable in the aftershock core
  (4-14 km depth), 5 km across nearly the whole footprint** -- vs Diambama's
  effective 5-10 km (5 km "center only"). The 20x ray gain buys ~1.5-2x
  resolvable fraction, NOT 20x: the fixed 12-station geometry (azimuth diversity)
  is the ceiling, not ray count.
- **Caveats/follow-ups:** straight-ray proxy + heuristic thresholds; a definitive
  answer needs a LOTOS/checkerboard run. Resolution is best INSIDE the source
  volume; margins and >15 km stay coarse regardless of catalogue size.

## [2026-08-14] planning | Seismica Data Report + event-browser data product

- **Trigger:** User asked whether the work is publication-worthy, targeting a
  Seismica data-based Report, with a GitHub clickable event table -> per-event
  waveform image, reusable for ML.
- **Files:** `eqt/scripts/build_event_browser.py` (new POC), `.gitignore`
  (ignore figures/events/), `wiki/outputs/yogya-2006-srl-paper-plan.md`,
  `wiki/log.md` (brain).
- **Key result:** Assessed as worth publishing — Seismica "Report" (data-focused)
  is a strong fit. Grounded in Seismica guidelines: <=10k words/10 items; data
  MUST be DOI-citable (Zenodo/GFZ), GitHub NOT archival. Architecture: Zenodo =
  catalogue CSVs + ML-ready SeisBench/HDF5 waveform archive (~2.8 GB downsampled)
  + preview PNGs; GitHub = code + interactive event browser (sortable HTML table
  linking to previews). POC browser built (works; empty previews this run because
  the raw EDL volume unmounted mid-session — pipeline itself proven earlier).
  "Most events" is a valid Data Report descriptor given honest completeness
  reporting + the pick validation QA.
- **Follow-ups:** remount raw volume, run `build_event_browser.py --all` (hours);
  build the SeisBench/HDF5 ML archive; create the Zenodo deposit (Seismica
  community) for the DOI; wire GitHub Pages browser to Zenodo-hosted images.

## [2026-08-14] export | SeisBench ML waveform archive builder

- **Trigger:** User approved building the SeisBench archive (the ML-ready data
  product for the Seismica report); raw EDL volume remounted.
- **Files:** `eqt/scripts/build_seisbench.py` (new), `.gitignore` (ignore
  eqt/seisbench/), `wiki/log.md`. Output `eqt/seisbench/{metadata.csv,
  waveforms.hdf5}` (gitignored -> Zenodo).
- **Key result:** Writes one trace per (event, station): 3-component (ZNE) 60 s
  window at 100 Hz around the P pick, P/S arrival-sample labels, full event +
  station metadata (ML + the +0.41 local-scale tie, gap, QC flag, 80/10/10
  train/dev/test split). Raw counts (STEAD/EQT convention). Verified round-trip
  via seisbench.data.WaveformDataset: P@sample 1500, S-P gaps physical, splits
  correct, waveforms (3,6000) float32. Event-browser POC also confirmed working
  with real waveforms once the volume was remounted. Full run launched
  (945 station-days, ~71k traces expected).
- **Follow-ups:** on completion, spot-check a few waveforms + labels; create the
  Zenodo deposit (Seismica community) with metadata.csv + waveforms.hdf5 +
  catalogue CSVs + preview PNGs; run `build_event_browser.py --all`; draft the
  Seismica Data Report manuscript.

## [2026-08-14] planning | SRL Data Mine vs Seismica venue assessment

- **Trigger:** User leaning toward SRL Data Mine instead of Seismica.
- **Files (brain):** `wiki/outputs/yogya-2006-srl-paper-plan.md`, `wiki/log.md`.
- **Key result:** Fetched SRL Data Mine author info. Max 6,000 words/10 figs,
  fixed 8-section template, judged on data access + reusability, NO new
  methods/hypothesis testing. Two honest concerns for our project: (1) it targets
  *recent* deployments — ours is a 2006 legacy deployment reprocessed into a new
  derived catalogue (confirm eligibility with the editor); (2) our pitfalls +
  pick-validation must ride as data-quality notes, not headline methods (Seismica
  would let them breathe). Both venues need the Zenodo DOI. Recorded the Data Mine
  8-section -> asset map. Decision pending user.
- **Follow-ups:** user to pick venue (+ optionally email Data Mine editor re
  legacy-data eligibility); then draft to the chosen template; Zenodo deposit is
  the shared next production step.

## [2026-08-14] export | SeisBench archive COMPLETE + verified (97,660 traces)

- **Follow-up to the SeisBench builder entry above.** Full run finished:
  **97,660 traces** (waveforms.hdf5 6.6 GB + metadata.csv 23 MB) from all 16,876
  events, 17 stations. P labels on all traces; S on 94,051 (96%); 71,107 QC-passed
  traces. Splits train 78,075 / dev 9,633 / test 9,952. Verified via
  seisbench.data.WaveformDataset: S-P all positive (median 2.22 s), P-SNR ~5-9.5
  at the P label, labels visually aligned with arrivals across ML +3.55..-0.25
  (`eqt/figures/seisbench_check.png`). The ML-ready data product is ready for the
  Zenodo deposit. 6.6 GB at 100 Hz fits Zenodo (50 GB); optional downsample if
  desired.
- **Follow-ups:** Zenodo deposit (Seismica/SRL Data Mine — catalogue CSVs +
  waveforms.hdf5 + metadata.csv + preview PNGs); run build_event_browser.py --all;
  publication-clean the figures.

## [2026-08-14] export | Assemble Zenodo deposit folder

- **Trigger:** User asked to run the full event browser and assemble the Zenodo
  deposit; intended target /Volumes/Untitled 1 (large free space).
- **Gotcha:** /Volumes/Untitled 1 is **NTFS, mounted read-only** on macOS — cannot
  write there. Assembled instead at **~/Yogya2006_Zenodo** on the internal disk
  (98 GB free; ~8 GB deposit). User can upload to Zenodo from there (browser/API
  upload does not need the external drive) or copy to a Mac-formatted volume.
- **Files:** `eqt/scripts/build_gallery.py` (new — previews + browser sourced
  from the SeisBench HDF5, NOT the raw volume: fast, drive-independent, matches
  deposited data), `manuscript/zenodo_deposit_README.md` (tracked copy of the
  deposit README), `wiki/log.md`.
- **Deposit contents (~/Yogya2006_Zenodo/):** catalog/ (magnitude, quality,
  growclust CSVs + stations.csv), waveforms_ml/ (waveforms.hdf5 6.6 GB +
  metadata.csv), event_browser/ (index.html + per-event previews), README.md.
- **Design note:** the event gallery is generated from waveforms.hdf5, so
  previews are guaranteed consistent with the archived data and no raw-volume
  re-read (~hours) is needed. Full run launched (16,876 previews, ~60 min).
- **Follow-ups:** user creates the Zenodo record (Seismica/SRL community) + DOI,
  fills the README placeholders (raw-data repo/DOI, citation), completes author
  list + references in the manuscript, publication-cleans the figures.

## [2026-08-14] export | Zenodo deposit COMPLETE (7.4 GB, verified)

- **Follow-up:** Gallery finished — **16,876 event previews (873 MB)** generated
  from the SeisBench HDF5 + index.html. Full deposit at ~/Yogya2006_Zenodo =
  **7.4 GB**: catalog/ (4 CSVs), waveforms_ml/ (6.6 GB hdf5 + metadata, 97,660
  rows), event_browser/ (16,876 previews + sortable index), README.md. Integrity
  verified: all previews present, median 51 KB (real waveforms), labels aligned.
  **The full open data product is assembled and ready to upload to Zenodo.**
- **All that remains is user-side:** create the Zenodo record + DOI; fill README/
  manuscript placeholders (raw-data repo/DOI, authors, references); confirm Data
  Mine legacy-data eligibility with the editor; publication-clean the figures.
