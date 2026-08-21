# MERAMEX — full land network, same 11 days

Control experiment against [`../pilot/SUMMARY.md`](../pilot/SUMMARY.md): identical
days (DOY 155–165), identical settings, but **all 108 land stations** instead of
the 43 nearest Yogyakarta–Merapi.

## Cost

1,126 station-days, 41 GB read, 61 min picking on MPS (5.2 s per station-day).

## What more stations bought

| | 43 stations | 108 stations | change |
|---|---|---|---|
| Picks | 2,577 | 5,400 | ×2.1 |
| **Events located** | 67 | **79** | **+18 %** |
| Events passing quality cut | 7 | **20** | **×2.9** |
| Median azimuthal gap | 318° | 307° | −11° |
| Events with gap ≤ 250° | 1 | 9 | ×9 |
| Median horizontal error | 8.5 km | 6.7 km | −21 % |
| Median vertical error | 11.7 km | 9.5 km | −19 % |
| Events with errH ≤ 5 km | 16 | 28 | ×1.8 |
| Crustal events (z < 25 km) | 10 | 19 | ×1.9 |

All 67 events from the smaller set were recovered — nothing was lost.

**The extra stations buy accuracy, not count.** The 74 added sites sit a median
168 km from the seismicity centroid (versus 111 km for the 43), because most of
them are the northern Semarang–Rembang–Karimunjawa belt. They widen the aperture,
so 11 of the 12 new events are at the *edges* of the network — East Java
(110.9–112.7°E), the west (108.4–109.4°E), the north coast — not new small events
in the centre. Detection in the array centre was already saturated.

**Implication for the full campaign:** the event count will scale with *time*,
not with adding stations. Expect roughly 7 events/day × ~150 days ≈ **1,000
events** — about 3× the 292 of Wagner et al. (2007) and 2× the 505 processed by
Koulakov et al. (2007), but with far better-constrained hypocentres. The 10–30×
figure quoted after the pilot was too optimistic and is withdrawn.

## The Opak / 2006 rupture question

| | |
|---|---|
| Stations within 10 / 20 / 30 km of the 2006 rupture centroid | 4 / 8 / 18 |
| 2004 events within 20 km | **0** |
| 2004 events inside the 2006 aftershock box | **0** |
| 2004 events within 50 km of Merapi | **0** |
| Well-located shallow crustal events elsewhere on land | 7 (errH 1.9–5.0 km, gap down to 134°) |

The positive controls matter: the network demonstrably locates shallow crustal
earthquakes on land to ~2 km when they happen. It found none at the place with
its own densest coverage.

Independent check that this is not a picker/associator blind spot: of 783 P–S
pairs, 130 have S−P < 3 s (source within ~24 km of that station), but **none is
seen by three or more stations at once** — they are single-site noise, not
missed earthquakes. The S−P distribution peaks at 12.3 s ≈ 100 km, i.e. the
seismicity really is distant.

**Caveat that still stands:** 11 days. At a background rate of a few events per
year in that box the expected count is ≈ 0.1, so zero is not yet evidence of
quiescence. The 150-day run converts this into a measurable upper bound.

Figure: `../figures/opak_2004_vs_2006.png`
