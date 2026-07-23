<!-- afk-research:managed v1 -->
# Index

This is the content map for the `.brain` second brain. Read this file before answering from or editing the wiki. Update it on every ingest and every durable query that changes the vault.

Last updated: 2026-07-20

## Start Here

- [[AGENTS|AGENTS.md]] - Operating schema for the LLM wiki agent.
- [[CLAUDE|CLAUDE.md]] - Claude-compatible pointer to the schema.
- [[wiki/log|log.md]] - Append-only chronological activity log.

## Folder Conventions

- `raw/` - Immutable source captures and provenance records.
- `wiki/` - All maintained knowledge and operating files.
- `wiki/index.md` - This content map.
- `wiki/log.md` - Append-only activity log.
- `wiki/sources/` - Source summaries and extracted takeaways.
- `wiki/concepts/` - Durable ideas and patterns.
- `wiki/entities/` - People, organizations, projects, datasets, and named things.
- `wiki/claims/` - Evidence-bearing claims worth tracking.
- `wiki/questions/` - Research questions and durable answers.
- `wiki/syntheses/` - Multi-source analysis and evolving theses.
- `wiki/outputs/` - Exportable artifacts and examples.
- `wiki/templates/` - Reusable page templates and workflow checklists.
- `wiki/inbox/` - Unprocessed material waiting for ingest.
- `wiki/scratch/` - Temporary agent work notes.
- `wiki/archive/` - Superseded or inactive material.

## Sources

- [[wiki/sources/earthquake-transformer-mousavi-2020]] - EQTransformer: deep-learning model for simultaneous earthquake detection + P/S phase picking (Nat. Commun. 2020).
- [[wiki/sources/eqtransformer-github-repo]] - Official MIT-licensed Python package (TensorFlow/ObsPy) implementing EQTransformer with pretrained models.
- [[wiki/sources/yogya-2006-aftershock-edl-dataset]] - Raw GFZ EDL aftershock recordings (external volume); layout + critical data-hygiene facts.

## Concepts

- [[wiki/concepts/seismic-phase-picking]] - Identifying P/S wave arrival times; S picking is the hard case.
- [[wiki/concepts/earthquake-signal-detection]] - Deciding signal-vs-noise; sensitivity vs. false positives.
- [[wiki/concepts/deep-learning-seismic-processing]] - Attention + multi-task tandem learning on waveforms.
- [[wiki/concepts/phase-association]] - Grouping per-station picks into single events; bridge to location.

## Entities

- [[wiki/entities/eqtransformer-model]] - The EQTransformer model / open-source package.
- [[wiki/entities/2000-tottori-earthquake-sequence]] - Japan sequence used as EQTransformer's field-test dataset.
- [[wiki/entities/yogya-2006-temp-aftershock-network]] - 12-station temp EDL array + coordinate table.

## Syntheses

- [[wiki/syntheses/eqtransformer-yogya-2006-run]] - Pipeline, decisions, env, and live status for running EQTransformer on the Yogya 2006 data.

## Claims

- [[wiki/claims/eqtransformer-doubles-detections-tottori]] - ~2× more earthquakes detected with <1/3 of stations (Tottori 2000).

## Questions

- [[wiki/questions/applying-eqtransformer-to-yogya-2006]] - Can EQTransformer improve the Yogyakarta 2006 aftershock catalog?

## Outputs

- [[wiki/outputs/yogya-2006-eqt-catalog]] - EQTransformer Yogya 2006 aftershock catalog (16,876 located; 11,790 quality-passed). Magnitudes corrected 2026-07-20: Mc=+0.50, b=0.89±0.02.
- [[wiki/outputs/yogya-2006-1d-velocity-model]] - Data-driven 1-D velocity (Wadati Vp/Vs=1.735; VELEST minimum-1D + station corrections).
- [[wiki/outputs/yogya-2006-pangandaran-triggering]] - Pangandaran Mw7.7 (17 Jul 2006) triggering test: no detectable triggering; small-event deficit is a completeness artifact.
- [[wiki/outputs/yogya-2006-srl-paper-plan]] - SRL paper plan (2026-07-20): reframed as methodology + open data product; "N× more events" claim dropped; pick validation is the main remaining task.

## Templates

- [[wiki/templates/source-page]] - Template for source summary pages.
- [[wiki/templates/concept-page]] - Template for concept pages.
- [[wiki/templates/entity-page]] - Template for entity pages.
- [[wiki/templates/ingest-checklist]] - Checklist for future ingest work.

## Open Threads

- Assess whether EQTransformer generalizes to Yogyakarta 2006 aftershock data — see [[wiki/questions/applying-eqtransformer-to-yogya-2006]].
- **SRL paper (reframed 2026-07-20):** methodology + open data product, not a count comparison — see [[wiki/outputs/yogya-2006-srl-paper-plan]]. Ramdhan FMD no longer blocking. **Main remaining task: pick validation vs manual ground truth** (biggest reviewer risk for an ML-picker paper). Then Zenodo DOI release.
- ~~Open magnitude issues: FMD roll-off / TF16 anomaly~~ **RESOLVED 2026-07-20**: FMD roll-off is real (events reach only 16% of the digitiser rail, no clipping ⇒ max ML 3.55 is real, not a lower bound); TF16 −0.37 is a genuine site/path effect (distance-dependent, time-stable, both components healthy), absorbed by the station correction. See [[wiki/outputs/yogya-2006-eqt-catalog]].
- `needs-review`: EQTransformer's exact architecture, training-set name/size, and per-station Tottori metrics were not captured from a primary source — fill in from the full paper/methods.
- Keep all second-brain work inside `.brain`, with content organized under only `raw/` and `wiki/`.
- Brain operations that change `.brain` must commit and push those `.brain` changes to the repository remote before reporting completion.
- Brain imports are non-destructive: imported knowledge enriches the current vault, while schema and operational-file collisions are preserved as provenance unless explicitly approved.
- Brain setup should also ensure generated export/import artifacts under `.outputs/` are ignored by Git.
