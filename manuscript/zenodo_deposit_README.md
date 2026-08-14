# Yogyakarta 2006 (Mw 6.4) aftershock catalogue and ML-ready waveform archive

An open, deep-learning-derived aftershock catalogue for the 26 May 2006
Yogyakarta earthquake, with a machine-learning-ready labelled waveform archive.
Produced from the temporary GFZ/BMKG/UGM short-period network (network code XN)
by reprocessing the continuous archive with EQTransformer.

**Companion article:** [SRL Data Mine / Seismica — citation + DOI to be added]
**License:** CC-BY-4.0
**Contact:** W. Suryanto, Universitas Gadjah Mada (ws@ugm.ac.id)

---

## Contents

```
Yogya2006_Zenodo/
├── README.md                     ← this file
├── catalog/
│   ├── catalog_magnitude.csv     16,876 located events: time, hypocentre, ML, quality metrics
│   ├── catalog_quality.csv       per-event quality metrics + pass flag (11,790 pass)
│   ├── catalog_growclust.csv     13,251 GrowClust double-difference relocations + uncertainties
│   └── stations.csv              station codes + coordinates
├── waveforms_ml/
│   ├── waveforms.hdf5            SeisBench waveform archive (97,660 traces, 6.6 GB)
│   └── metadata.csv             per-trace metadata + P/S labels
└── event_browser/
    ├── index.html               sortable event table linking to per-event previews
    └── previews/ev*.png         per-event multi-station waveform previews
```

## 1. Event catalogue (`catalog/`)

- **16,876 located events**, 3 June–29 August 2006. Located with NonLinLoc in a
  data-driven VELEST 1-D model (Vp/Vs = 1.735) with station corrections.
- **11,790 pass** a scale-free quality screen (`catalog_quality.csv`: fraction of
  badly fitted phases, RMS, phase/S counts, azimuthal gap, `pass` flag). Users may
  re-screen from the provided metrics.
- **13,251 relocated** by GrowClust double-difference (`catalog_growclust.csv`),
  median relative precision 274 m horizontal / 305 m vertical / 47 ms
  (bootstrap n=100).

### Magnitudes — important note on the scale

Local magnitudes (`ML` in `catalog_magnitude.csv`) are Wood-Anderson M_L
(Hutton & Boore 1987) from S-wave amplitudes after response removal, 1–20 Hz
band-pass, and per-station site corrections. Range M_L −1.8 to 3.6;
completeness M_c ≈ 0.5 (b-value stability); b = 0.89 ± 0.02.

An event-matched comparison against the manual catalogue of Anggraini (2013)
gives a **constant offset of +0.41** (r = 0.95): our M_L reads 0.41 unit **below**
the local manual scale. To compare with the local/regional catalogues, use
**M_L(tied) = M_L + 0.41**. The b-value is invariant to this offset.

## 2. Machine-learning waveform archive (`waveforms_ml/`)

SeisBench-format dataset (Woollam et al. 2022): **97,660 traces**, one per
(event, station).

- **Waveforms:** three-component (Z, N, E), 60 s at 100 Hz (6000 samples),
  starting 15 s before the P pick. Raw counts (L4-3D response NOT removed).
- **Labels:** `trace_p_arrival_sample` (all traces, at sample 1500),
  `trace_s_arrival_sample` (94,051 traces).
- **Metadata:** origin time, hypocentre, `source_magnitude` (raw M_L) and
  `source_magnitude_local_tie` (M_L + 0.41), azimuthal gap, quality flag,
  station coordinates, epicentral distance, and a `split` column
  (train / dev / test = 78,075 / 9,633 / 9,952).

Load with SeisBench:

```python
import seisbench.data as sbd
data = sbd.WaveformDataset("waveforms_ml")   # folder with metadata.csv + waveforms.hdf5
wave = data.get_waveforms(0)                 # (3, 6000) array, Z,N,E
p_sample = data.metadata.trace_p_arrival_sample.iloc[0]
```

Intended for training and benchmarking deep-learning phase pickers on tropical,
short-period, temporary-network data. Filter with `source_quality_passed` for the
71,107 quality-screened traces.

## 3. Event browser (`event_browser/`)

Open `index.html` in a browser: a sortable table of all events (time, location,
depth, M_L, M_L_tied, station count, gap, quality) where each event id links to a
preview of its multi-station waveforms with the P/S picks. Previews are generated
from the archived waveforms, so they match `waveforms.hdf5` exactly.

## Provenance

Raw continuous data: temporary deployment by the German Task Force for
Earthquakes (GFZ Potsdam) with BMKG and Universitas Gadjah Mada, from 31 May 2006
(12 Mark L4-3D 1 Hz sensors on Earth Data Loggers, network code XN). The same
archive underlies the manual catalogues of Anggraini (2013) and Diambama et al.
(2019). [State the raw-data repository / network DOI and any embargo here.]

## How to cite

[Article citation + this dataset's Zenodo DOI — to be added on deposit.]

## Key references

- Mousavi et al. (2020), *Nat. Commun.* 11, 3952 — EQTransformer.
- Woollam et al. (2022), *SRL* 93, 1695–1709 — SeisBench.
- Trugman & Shearer (2017), *SRL* 88, 379–391 — GrowClust.
- Hutton & Boore (1987), *BSSA* 77, 2074–2094 — M_L scale.
- Diambama et al. (2019), *GJI* 216, 439–452 — Ngalang Fault, same deployment.
- Anggraini (2013), PhD dissertation, Univ. Potsdam — manual catalogue.
