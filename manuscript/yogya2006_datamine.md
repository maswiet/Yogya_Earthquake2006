# A deep-learning-derived aftershock catalogue and machine-learning-ready waveform archive for the 2006 Yogyakarta (M<sub>w</sub> 6.4) earthquake

*Draft for SRL Data Mine — target: 6,000 words, ≤10 figures. Written to the Data Mine 8-section template. Figures/numbers current as of the working catalogue; the Zenodo DOI is a placeholder pending deposit.*

**Authors (placeholder):** W. Suryanto¹\*, A. Anggraini¹, [others], B.-G. Lühr²
¹ Seismology Research Group, Department of Physics, Universitas Gadjah Mada, Yogyakarta, Indonesia
² GFZ German Research Centre for Geosciences, Potsdam, Germany
\* Corresponding author: ws@ugm.ac.id

---

## 1. Abstract

The 26 May 2006 Yogyakarta earthquake (M<sub>w</sub> 6.4) was one of Indonesia's
most damaging crustal events, and its aftershock sequence was recorded by a
temporary 12-station short-period network operated by the German Task Force for
Earthquakes (GFZ) with Indonesian partners. We reprocess this legacy continuous
archive with a deep-learning phase picker (EQTransformer) to produce an open,
quality-screened aftershock catalogue of **16,876 located events** spanning
3 June–29 August 2006, of which **11,790 pass a scale-free quality screen** and
**13,251 are relocated** by double-difference (GrowClust) to a median relative
precision of 274 m (horizontal) / 305 m (vertical). Local magnitudes span
M<sub>L</sub> −1.8 to 3.6 (magnitude of completeness M<sub>c</sub> ≈ 0.5,
b = 0.89 ± 0.02); an event-matched comparison against the group's earlier manual
catalogue ties the automatic magnitude scale to the local scale by a constant
+0.41 unit (r = 0.95). Automatic picks match co-located manual picks with a
median absolute residual of 0.02 s (P) and 0.06 s (S), and 90 % of the manual
events are recovered. We release the catalogue, per-event quality metrics,
relocations, and a **machine-learning-ready waveform archive** (SeisBench/HDF5
format: per event-station three-component windows with labelled P and S arrivals)
via a DOI-citable repository, together with an interactive event browser. The
dataset supports reuse in machine-learning phase-picking, local-earthquake
tomography, statistical seismology, and hazard studies.

**Keywords:** Yogyakarta earthquake; aftershocks; EQTransformer; earthquake
catalogue; machine learning; open data.

---

## 2. Introduction (Experiment Motivation)

The 26 May 2006 Yogyakarta earthquake (M<sub>w</sub> 6.4; 22:53 UTC / 05:54 local
27 May) struck the densely populated Bantul region south of Yogyakarta, Central
Java, causing ~5,700 fatalities and economic losses near USD 3.1 billion
(International Recovery Platform & UGM, 2009). The event ruptured a crustal fault
in the transition between the Yogyakarta depression and the Southern Mountains,
but the causative structure has been debated for nearly two decades: the surface
Opak Fault, versus an unnamed NE–SW structure ~10–15 km to its east that our
group first associated with the Ngalang Fault from local-earthquake tomography
(Diambama et al., 2019), later independently identified as the source by Ramdhan
et al. (2025).

The aftershock sequence was recorded by a temporary deployment (Section 3) whose
continuous data have, to date, been analysed only through manually picked
catalogues of a few hundred events (Anggraini, 2013; Diambama et al., 2019).
Deep-learning phase pickers now make it feasible to reprocess such legacy
archives at a completeness and event count far beyond what manual analysis
allowed, and — critically — to release the results as reusable, labelled data
products. This contribution documents such a reprocessed catalogue and its
accompanying machine-learning-ready waveform archive. Consistent with the Data
Mine scope, we describe the dataset, its production, and its quality and
reusability; we do not test tectonic hypotheses or introduce new methods.

---

## 3. Instrument Deployment and Details (Experiment Design)

Following the mainshock, the German Task Force for Earthquakes (GFZ Potsdam), with
the Indonesian Agency for Meteorology, Climatology and Geophysics (BMKG) and
Universitas Gadjah Mada, deployed a temporary network of **12 short-period
stations** (network code XN) around the aftershock zone, operating from
31 May 2006. Each station comprised a **Mark L4-3D three-component 1 Hz
geophone** (velocity sensitivity 1.7 × 10⁸ counts m⁻¹ s) recorded on an **Earth
Data Logger (EDL)** at 100–200 Hz, with GPS timing. Average inter-station spacing
was ~10 km (≈4 km near the Opak Fault, ≈16 km at the margins). A subset of
stations was relocated during the deployment; we track the time-dependent station
positions and use the position active at each event time. Station coordinates are
provided with the data release.

The reprocessing here uses the continuous vertical- and horizontal-component
records for days 154–241 of 2006 (3 June–29 August), covering 944 station-days.

---

## 4. Overall Data Quality and Availability

### 4.1 Catalogue production

Continuous records were processed with **EQTransformer** (Mousavi et al., 2020)
via SeisBench (Woollam et al., 2022) using the pretrained `original` (STEAD)
weights to detect events and pick P and S arrivals. Picks were associated into
events with PyOcto and located with **NonLinLoc** (Lomax et al., 2000) in a
data-driven 1-D velocity model derived for this sequence with **VELEST**
(Kissling et al., 1994; V<sub>p</sub>/V<sub>s</sub> = 1.735), with station
corrections that map the sediment-versus-limestone contrast across the Opak
Fault. The procedure yields **16,876 located events** over 3 June–29 August 2006.

### 4.2 Quality screening

We provide a per-event quality table with scale-free metrics (fraction of badly
fitted phases, RMS, phase and S-pick counts, azimuthal gap) and a boolean pass
flag; **11,790 events (69.9 %) pass**. We deliberately avoid thresholding the
single largest residual, which scales with phase count and would reject the
best-recorded events; the screen therefore admits a higher fraction of larger,
better-recorded events (rejection falls from ~44 % below M<sub>L</sub> 0 to ~20 %
at M<sub>L</sub> 0.5–1). The dominant rejection cause is azimuthal gap > 180°
(a network-geometry limit), not detection quality. Users may re-screen from the
provided metrics.

### 4.3 Relocation

The full quality catalogue was relocated by double-difference with **GrowClust**
(Trugman & Shearer, 2017) using 3.56 million catalogue differential times:
**13,251 events (97 %) relocated** in 543 clusters (largest cluster 6,917
events). Bootstrap resampling (n = 100) gives median relative uncertainties of
**274 m horizontal, 305 m vertical, 47 ms origin time**. Relocated hypocentres,
cluster IDs, and per-event uncertainties are included.

### 4.4 Local magnitudes and completeness

Local magnitudes were computed from Wood-Anderson-simulated S-wave amplitudes
(Hutton & Boore, 1987) after per-pick instrument-response removal, a 1–20 Hz
band-pass (to reject 50 Hz mains noise present at some sites), and per-station
site corrections. Magnitudes span **M<sub>L</sub> −1.8 to 3.6** (median −0.22).
Completeness assessed by b-value stability (Woessner & Wiemer, 2005) is
**M<sub>c</sub> ≈ 0.5**, with **b = 0.89 ± 0.02** (Aki–Utsu). The maximum-
curvature estimate (M<sub>c</sub> = −0.2) is reported for reference but
underestimates completeness for this network. *(Figure: frequency–magnitude
distribution and detection-limit analysis.)*

Two site/processing effects are documented as data-quality notes and are already
corrected in the released magnitudes: (i) instrument-simulation on day-long
traces can suppress amplitudes of events near a UTC day boundary — handled by
per-pick windowed deconvolution; (ii) a 50 Hz mains tone at some stations
(notably one hard-rock site) is removed by the band-pass. One station shows a
distance-dependent, time-stable amplitude deficit consistent with a hard-rock
site response east of the Opak Fault, absorbed by its station correction.

### 4.5 Validation against manual picks

Because the same sequence was manually picked for earlier studies (Anggraini,
2013), we validate the automatic catalogue directly. Over the common 3–7 June
window, **90 % of manual events (528/588) are recovered** within 5 s origin-time.
Automatic picks match co-located manual picks with **median absolute residuals of
0.02 s (P) and 0.06 s (S)** (per-event demeaned; 4,761 matched picks). Event-
matched magnitudes correlate at **r = 0.95** with a constant offset of **+0.41**
(our M<sub>L</sub> reads 0.41 low relative to the manual scale); we provide both
the raw M<sub>L</sub> and the local-scale-tied value. The b-value is invariant to
this offset. *(Figure: pick-residual histograms + magnitude cross-plot.)*

### 4.6 Data availability and formats

All products are archived at [Zenodo DOI — placeholder] under CC-BY-4.0:

- **Event catalogue** (CSV): origin time, hypocentre, M<sub>L</sub> (raw and
  local-tied), quality metrics + pass flag, GrowClust relocation + uncertainties.
- **Machine-learning waveform archive** (SeisBench/HDF5 + metadata CSV): one trace
  per event-station — three-component (Z,N,E) 60 s windows at 100 Hz around the P
  pick, with labelled P and S arrival samples, event and station metadata, and a
  train/dev/test split. Directly loadable with `seisbench.data.WaveformDataset`.
- **Per-event waveform previews** (PNG) and an **interactive event browser**
  (sortable HTML table linking each event to its multi-station waveforms), served
  from the code repository.

---

## 5. Initial Observations

The catalogue delineates a compact, NE–SW-elongated (≈N57°E) aftershock zone
between the Opak and Ngalang faults, ~10–15 km east of the surface Opak trace,
with hypocentres concentrated at **2.8–14.9 km depth (median 10.5 km)**. After
double-difference relocation the diffuse cloud collapses onto a steeply dipping
plane consistent with the GCMT mechanism (near-vertical NE–SW left-lateral
strike-slip) and with the Ngalang Fault association of Diambama et al. (2019).
The daily event rate decays in an Omori-like manner (p ≈ 1.05) from a peak of
382 events/day on 17 June. *(Figures: epicentre/relief map; depth cross-section;
along-strike-vs-time; frequency–magnitude distribution.)*

---

## 6. Initial Results (optional)

**Reusability for tomography.** A straight-ray coverage analysis indicates the
catalogue's ~142,000 ray paths (≈20× the earlier manual dataset) support a
velocity-model resolution of ~2–3 km in the source-volume core (4–14 km depth)
and ~5 km across the aftershock footprint, versus the 5–10 km of prior tomography
— a resolution improvement bounded by the fixed station geometry rather than
event count. A definitive resolution test requires a full inversion and is left
to future work. *(Figure: ray-coverage / azimuthal-diversity maps.)*

**Reusability for machine learning.** The labelled waveform archive (§4.6) is
provided in a community-standard format for training and benchmarking
deep-learning phase pickers on tropical, short-period, temporary-network data —
a setting under-represented in existing training sets.

---

## 7. Summary

We release an open, deep-learning-derived aftershock catalogue for the 2006
Yogyakarta earthquake — 16,876 located events (11,790 quality-screened; 13,251
double-difference relocated to sub-500 m relative precision) — validated
event-for-event against manual picks (90 % recovery; P/S residual medians
0.02 / 0.06 s; magnitude r = 0.95). Alongside the catalogue we provide a
machine-learning-ready, labelled waveform archive and an interactive event
browser. The dataset substantially extends the completeness and event count
available for this important sequence and is designed for reuse in
machine-learning, tomography, statistical-seismology, and hazard applications.

---

## 8. Data and Resources

The event catalogue, quality metrics, relocations, local magnitudes, and the
SeisBench/HDF5 waveform archive are archived at [Zenodo DOI — placeholder]
(CC-BY-4.0). Processing code, the interactive event browser, and figure scripts
are at [GitHub URL]. EQTransformer is from Mousavi et al. (2020) via SeisBench
(Woollam et al., 2022); NonLinLoc from Lomax et al. (2000); VELEST from Kissling
et al. (1994); GrowClust from Trugman & Shearer (2017). The raw continuous data
were collected by the GFZ German Task Force for Earthquakes with BMKG and UGM.
[Note per Data Mine policy: state the raw-data repository/DOI and any embargo.]

---

## References (to complete)

Anggraini, A. (2013). *The 26 May 2006 Yogyakarta earthquake, aftershocks and
interactions.* PhD dissertation, University of Potsdam.

Diambama, A. D., Anggraini, A., Nukman, M., Lühr, B.-G., & Suryanto, W. (2019).
Velocity structure of the earthquake zone of the M6.3 Yogyakarta earthquake 2006
from a seismic tomography study. *Geophys. J. Int.*, 216, 439–452.

Hutton, L. K., & Boore, D. M. (1987). The M<sub>L</sub> scale in southern
California. *Bull. Seismol. Soc. Am.*, 77, 2074–2094.

Kissling, E., et al. (1994). Initial reference models in local earthquake
tomography. *J. Geophys. Res.*, 99, 19635–19646.

Lomax, A., et al. (2000). Probabilistic earthquake location in 3D and layered
models (NonLinLoc). In *Advances in Seismic Event Location*.

Mousavi, S. M., et al. (2020). Earthquake transformer — an attentive
deep-learning model for simultaneous earthquake detection and phase picking.
*Nat. Commun.*, 11, 3952.

Ramdhan, M., et al. (2025). Aftershock sequence of the Yogyakarta earthquake 2006
(M<sub>w</sub> ~6.4)… *Nat. Hazards*, 121. doi:10.1007/s11069-025-07440-8.

Trugman, D. T., & Shearer, P. M. (2017). GrowClust: a hierarchical clustering
algorithm for relative earthquake relocation. *Seismol. Res. Lett.*, 88, 379–391.

Woessner, J., & Wiemer, S. (2005). Assessing the quality of earthquake
catalogues: estimating the magnitude of completeness and its uncertainty. *Bull.
Seismol. Soc. Am.*, 95, 684–698.

Woollam, J., et al. (2022). SeisBench — a toolbox for machine learning in
seismology. *Seismol. Res. Lett.*, 93, 1695–1709.
