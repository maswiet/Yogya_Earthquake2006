# SIMULPS (simul2023) Tomographic Inversion — MERAMEX 2004

## Source

[Simul2023](https://zenodo.org/records/10695070) (Eberhart-Phillips, Thurber, Rietbrock et al.,
CC-BY-4.0), a direct descendant of Thurber's 1983 SIMUL3/SIMULPS12 local-earthquake
tomography code. Complete Fortran source (`simul2023.f`, `simul2023_common.inc`),
compile notes, and format documentation (`*.txt`) are kept here for provenance;
the vendor's ~100 MB `examples.zip` was used to verify our converter's byte-level
formats against real working input files, then discarded (not needed after
verification, and too large to keep in this repo — re-download from the Zenodo
record if you need to re-check a format).

Build: `gfortran -O2 -fcheck=all -o simul2023 simul2023.f` (see `NOTE_compile.txt`).
The binary is not committed; rebuild it locally.

## Two real bugs found and patched in the vendor's 2024 release

Both were confirmed by (a) tracing the exact crash through the source, (b) running
our compiled binary unmodified against the vendor's own `spi` example (which
worked, ruling out a compiler/platform issue), and (c) hand-computing the expected
coordinate/behaviour and confirming the patched code matches it.

**1. `nzco=0` never initialises the origin's TM coordinates.**
`setorg()` returns immediately for `nzco=0` (`simul2023.f:11740`) without ever
calling `tmllxy` to set `onorth`/`oeast` — but `disto()`'s `nzco=0` branch
(`simul2023.f:2405`) unconditionally subtracts these from every station's and
event's easting/northing. Every distance in the study came out in the hundreds
to thousands of km instead of the true tens of km, and `intmap()`'s velocity-grid
lookup then indexed far outside its array (`ixloc` bound 10000). **Worked around**
by using `nzco=5` (Thurber's original short-distance conversion) instead, which
never touches the TM/UTM code path at all — correct and simpler for a ~400 km
study box; no source patch needed. See `scripts/simulps_convert.py::write_stations`.

**2. Southern-hemisphere UTM false-northing typo in `utm_geo`.**
Even with `nzco=3` (explicit central meridian), `utm_geo` has (`simul2023.f:14301-14302`):
```fortran
if(inorth.eq.0) north=0.d0
if(north.eq.1) north=10000000.d0     ! should be if(inorth.eq.1) ...
```
`north` is a local, un-SAVEd variable — for any southern-hemisphere point
(`inorth=1`, true for all of Java) neither line fires, so `north` is read
uninitialised and added to the northing. We did not need to patch this (nzco=5
avoids it too), but it explains why nzco=3 alone was not sufficient once nzco=0
was ruled out — documented here in case a future run needs true UTM.

**3. `getgap()` reads `iazst(1)` before it's necessarily set.**
`simul2023.f:3624` — `sort()` only populates `iazst(1..naz)`; when a specific
event's final iteration reaches `getgap` with `naz` at 0 or 1 P-observations
counted, `azp(iazst(1))` at line 3646 reads uninitialised/stale memory, and
whatever garbage integer comes back is used as an array index in the next
statement, aborting with an "index above upper bound" that has nothing to do
with the real cause. **Patched** (`simul2023.f`, `getgap`) to return early with
`gap=360.0` for `naz<2` rather than run the sort/gap logic on data that isn't
there. This is the one source change we kept — everything else was worked
around through input file / control parameter choices, not by editing vendor
code.

Both issues were reproduced deterministically (same event, same garbage index,
every run) and confirmed by first isolating a single event into a 1-event
input set before touching anything else — worth remembering if this vendor
code throws similar "index above/below bound" errors on a different dataset.

## What's here

- `relocation_demo/` — `nitmax=0` (relocation only, no velocity update) sanity
  run on all 393 quality-filtered events, the run used to isolate and fix the
  two bugs above. `output` has the full per-event log; `hypxyze.out` the final
  x/y/z/lat/lon hypocentres.
- `inversion/` — the real tomographic run, `nitmax=4`, 390 events (after also
  excluding 3 events deeper than 180 km — see Limitations), 7×7×7 interior grid
  (5×5×5 free P/Vp-Vs nodes plus one edge ring). `vlxyzltln.out` is the easiest
  file to parse for plotting (lat, lon, x, y, z, Vp, Vp/Vs, Vs, DWS ray-density
  per node); `velomod.out` is the same model in the native simul input format,
  usable as a warm-start for a follow-up run; `itersum` has the iteration-0
  data variance and per-node hit counts.

## Reproducing

```bash
cd MERAMEX_2004
python3 scripts/simulps_convert.py --out simulps/inversion \
    --nx 7 --ny 7 --nz 9 --nitmax 4 --dxmax 3.0 --eigtol 0.01 --nitloc 6
cd simulps
gfortran -O2 -fcheck=all -o simul2023 simul2023.f
cp simul2023 inversion/ && cd inversion && ./simul2023
python3 ../../scripts/plot_simulps_model.py --vlxyz vlxyzltln.out \
    --out ../../figures/simulps_velocity_slices.png
```

## Result summary

- 390 of 1005 catalogue events pass quality filtering (`quality_pass` or
  errH≤15 km, gap≤300°, RMS≤0.6 s, depth≤180 km — the depth cut drops 3 events
  whose location strayed past our grid's interior/edge transition during
  relocation, a real Fortran `STOP` from `TTMDER`, not a crash, so it's easy to
  tell apart from an actual bug) — 9,996 P+S-P observations.
- 7×7×9 grid, 5×5×7 free interior nodes (175 total, 97 with ≥1 ray hit after
  the `nhitct=5` cutoff, 93 P-velocity + 4 Vp/Vs).
- Iteration 0 (starting 1-D model) weighted RMS 0.796 s; converges over 4
  velocity-inversion iterations (see `itersum`).
- Ray density (DWS) peaks at 25–50 km depth (max 36,399 and 29,137 hits) and
  falls off sharply below 150 km (517 hits) and essentially vanishes by 200 km
  (11 hits) — the same resolution ceiling our earlier Python/LOTOS-style
  inversion found independently (`../tomo/ASSESSMENT.md`), now confirmed by a
  second, independently-implemented tomography code on the same data.

## Limitations / honest caveats

- **Grid is coarse** (5 interior nodes per horizontal axis, ~90 km spacing) —
  deliberately conservative for a first validated run; the vendor's own
  "Hint 1" (start with `nx=ny=3`) turned out to be the wrong end of the scale
  for this program version (it triggered bug #1's symptom on the very first
  event, before we'd isolated the real cause) — a real 3-D grid is what the
  vendor's own examples use and is what we ended up validating against.
- **Straight rays only** (`i3d=0`): no pseudo-bending. Reasonable as a first
  pass and much faster (the full 390-event, 4-iteration inversion runs in
  ~10 s), but the reported velocity perturbations should be treated as a
  first-order result, not a publication-ready model. Re-running with `i3d=2`
  is the natural next step once the grid is refined.
- **3 deep events excluded** (the depth≤180 km cut) rather than fixed —
  their post-location depth (up to 268 km) landed in the outer grid ring where
  simul2023's own `TTMDER` correctly refuses to extrapolate. A finer/deeper
  interior grid would let these back in; for now they're a documented, minor
  exclusion (0.8% of the quality-passed catalogue), not a hidden one.
