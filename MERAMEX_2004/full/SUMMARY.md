# MERAMEX 2004 — full campaign result

**6 May – 7 Oct 2004, 15,271 station-days, 93,047 picks, 1,005 located events.**

## Picking

| | station-days | picks | P | S |
|---|---|---|---|---|
| Land (EDL + SAM), 128 sites | 13,948 | 67,853 | 37,302 | 30,551 |
| OBH/OBS, 14 sites | 1,323 | 25,194 | 19,321 | 5,873 |

Mean pick probability 0.61 (P) / 0.52 (S) on land; 4.9 picks per station-day.
The eight hydrophone-only OH sites contribute 3,029 P picks after the
`--duplicate-1c` fix (see `OBS_QC.md`).

## Catalogue (land-only association)

**1,005 events.** For comparison, the published MERAMEX catalogues are 292
events (Wagner et al. 2007) and 505 processed / 344 used (Koulakov et al. 2007),
so this is **3.4× the first and 2× the second** — and the projection made after
the 11-day control run (~1,000) was right.

| | |
|---|---|
| Median RMS | 0.29 s |
| Median azimuthal gap | 306° |
| Median errH / errZ | 7.2 / 9.1 km |
| Events passing the quality cut | 235 |
| Events with gap ≤ 180° | 72 |
| Events with errH ≤ 3 km | 122 |
| Median depth | 24.1 km (5–95 %: 1.7–88.0 km) |
| Events inside the array (lat > −8.2°) | 241 |
| Crustal events on land (z < 25 km, lat > −8.2°) | 111 |

OBS picks are **not** in this catalogue. Their clocks are uncalibrated (no GPS
underwater, no drift file on the drives), so the land-only solution is the
reference against which OBS drift will be measured.

## Wadati–Benioff zone

Fitting events 30–160 km deep in ±0.5° longitude bands:

| Profile | n | Dip | Scatter about the fit |
|---|---|---|---|
| 109.6°E | 48 | 35.7° | 10.8 km |
| 110.5°E | 154 | 25.3° | 28.7 km |
| 111.4°E | 33 | 33.1° | 20.2 km |

The 110.5°E number should not be read as a dip: that band contains the dense
shallow crustal population beneath the array as well as the slab, and a single
straight line through both is meaningless. The 109.6°E fit, with 10.8 km
scatter, is the only one clean enough to quote. Separating the two populations
is a task for the relocated catalogue.

## Relative relocation

`run_reloc.sh` on 398 selected events (catalog differential times, no
cross-correlation yet):

* **HypoDD** — 145 events relocated in 14 clusters, largest 66 events.
  Median errH **0.55 km**, errZ **0.45 km**, from 7.2 / 9.1 km before.
* **GrowClust** — 398 events, 26 multi-event clusters, largest 16.

### Crustal lineaments that emerge

| cid | n | Centre | Depth | Long × short | Ratio | Strike |
|---|---|---|---|---|---|---|
| 4 | 18 | −7.514 / 109.875 | 12 km | 17.8 × 5.0 km | 3.6 | 149° |
| 5 | 11 | −7.632 / 111.243 | 7 km | 21.7 × 9.8 km | 2.2 | 89° |
| 6 | 9 | −7.319 / 109.504 | 9 km | 30.7 × 3.2 km | **9.5** | 72° |
| 11 | 5 | −7.232 / 109.912 | 4 km | 3.8 × 0.4 km | **9.8** | 151° |
| 13 | 4 | −6.719 / 110.721 | 10 km | 5.1 × 0.3 km | **16.0** | 62° |

GrowClust, which is stricter, resolves several of the same features to
sub-kilometre width: 2.0 × 0.7 km at −7.516/109.873, 5.7 × 0.2 km striking 64°
on the north coast, 6.8 × 0.2 km striking 45° near Pemalang.

## The Opak / 2006-rupture question

Counting **crustal** events (z < 25 km) within 20 km of each centre, over the
full 152 days:

| Target | Crustal events within 20 km |
|---|---|
| Dieng | 13 |
| Lawu | 15 |
| **2006 rupture zone** | **4** |
| **Merapi** | **1** |

This revises the 11-day answer, which found zero and could not distinguish
"quiet" from "not enough time". The rupture zone is **not dead** — it produced
four shallow crustal events in five months, one of them inside the 2006
aftershock box itself (22 Aug 2004, −7.994/110.419, 19 km, gap 137°, errH
1.5 km). But it runs at roughly a quarter the rate of Dieng or Lawu, two areas
with comparable station coverage.

Merapi at one event is the quietest place in the network, two years before its
2006 eruption.

A 95 % Poisson upper bound on the rate inside the 2006 aftershock box is
**13 events/yr** — the first quantitative limit on pre-seismic activity there,
and the number that a locked-asperity interpretation has to be consistent with.

## Figures

* `figures/full_map_gmt.png` — catalogue over relief and Java-margin bathymetry
* `figures/full_slab_section.png` — three true-scale trench-perpendicular sections
* `figures/full_reloc.png` — NonLinLoc vs HypoDD at identical scale
