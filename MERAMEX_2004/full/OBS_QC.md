# OBS / OBH quality control

First full OBS pass: **22,165 picks over 1,323 station-days**, threshold 0.20
(land runs at 0.30 — the seafloor sites are noisier and EQTransformer's STEAD
training set contains no ocean-bottom data).

## Which sites produced anything

| Folder | INFO.DAT | Channels | Picks | P | S | P/S | picks/day |
|---|---|---|---|---|---|---|---|
| OS02 | OB2 | 3 | 875 | 381 | 494 | 0.8 | 7.6 |
| OS06 | OB6 | 3 | 2,138 | 939 | 1,199 | 0.8 | 18.6 |
| OS07 | OB7 | 3 | 975 | 502 | 473 | 1.1 | 14.3 |
| OS08 | OB8 | 3 | 2,704 | 1,377 | 1,327 | 1.0 | 21.6 |
| OS09 | OB9 | 3 | 2,717 | 1,357 | 1,360 | 1.0 | 21.6 |
| **OS12** | **O12** | 3 | **12,756** | **11,736** | **1,020** | **11.5** | **109** |
| OH01/03/04/05/10/11/13/14 | OB1, OB3-5, O10, O11, O13, O14 | 1 | 0 | — | — | — | — |

## OH sites: one channel

Eight of the fourteen are hydrophone-only. The 3-component model skipped them
silently (`len(st) >= 3`), so the first pass got nothing from them at all.
`run_picks.py --duplicate-1c` copies the single trace onto all three components,
which makes the model run: OH05 and OH13 give 3–4 P per day, OH01 still nothing.
**S picks from these sites are discarded** — a hydrophone measures pressure and
cannot record a shear wave. Their entries were cleared from `done_obs.txt` for a
re-pick.

## OS12 is an outlier and is not yet trusted

109 picks/day against 7.6–21.6 for the other five, and — more diagnostic — a
**P/S ratio of 11.5** where every other site sits near 1.0. The model is firing
"P" on impulsive transients that have no shear arrival behind them.

Checked and ruled out: a fixed instrumental tick. Inter-pick intervals are
broadly distributed (median 150 s, 5–95 % 14–507 s) with no dominant period, and
the rate is steady at ~97 P/day across all 117 days rather than episodic. So it
reads as a persistently noisy site — bottom current, instrument knock, or a
mis-wired vertical — rather than a periodic artefact.

OS12 alone contributes 58 % of all OBS picks. Left unscreened it would dominate
any association involving the offshore array. **Decision deferred to the
association stage**: the honest test is whether OS12 picks line up with events
located independently by the land network, which needs the full land catalogue.
Until then its picks are kept in the file but should not be trusted as arrivals.

## Coincidence between OBS sites is low

Only 0.5 % of OBS picks fall inside a 20 s window shared by three or more OBS
stations. That is not by itself alarming — the sites are 50–100 km apart and
most events will be recorded by the land array plus one or two OBS — but it does
mean the offshore array cannot be validated against itself. It has to be
validated against the land catalogue.

## Still outstanding: clock drift

OBS have no GPS underwater. No drift-correction file was found on either drive,
so the arrival times must be treated as uncalibrated until the residuals of each
OBS station are measured against land-only locations and fitted with a linear
drift. See the OBS section of `../README.md`.
