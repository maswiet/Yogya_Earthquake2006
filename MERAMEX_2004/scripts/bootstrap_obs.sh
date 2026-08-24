#!/bin/bash
# Bootstrap OBS integration: associate + relocate with OBS picks included

set -uo pipefail
cd "$(dirname "$0")/.."

echo "=== BOOTSTRAP STEP 1: Combine land + OBS picks ==="
python3 - <<'PY'
import pandas as pd
land = pd.read_csv("full/picks_land.csv")
obs = pd.read_csv("full/picks_obs.csv")
obs['kind'] = 'OBS'
combined = pd.concat([land, obs], ignore_index=True)
print(f"land: {len(land)} | OBS: {len(obs)} | combined: {len(combined)}")
combined.to_csv("full/picks_combined.csv", index=False)
PY

echo ""
echo "=== STEP 2: Associate ==="
python3 scripts/associate.py --picks full/picks_combined.csv --out full/events_combined.csv --tmpdir full/assoc_obs_tmp --n_picks 8 --n_p_and_s 3 --zmax 250 2>&1 | tail -3

echo ""
echo "=== STEP 3: Generate NLL control files ==="
python3 scripts/gen_nll.py --events full/events_combined.csv --assignments full/events_combined_assignments.csv --outdir nll --tag obs --nxy 601 --ndist 801 --nz 254 2>&1 | head -2

B=/Users/maswiet/Work/Yogya_Earthquake2006/eqt/tools/NonLinLoc/src/bin
echo "Computing travel-time grids..."
$B/Vel2Grid nll/nll_vel_obs.in > /dev/null 2>&1 &
$B/Vel2Grid nll/nll_vel_obs_S.in > /dev/null 2>&1 &
wait
$B/Grid2Time nll/nll_time_obs_P.in > /dev/null 2>&1 &
$B/Grid2Time nll/nll_time_obs_S.in > /dev/null 2>&1 &
wait

echo ""
echo "=== STEP 4: Locate ==="
time $B/NLLoc nll/nll_loc_obs.in > /dev/null 2>&1

echo ""
echo "=== STEP 5: Parse locations ==="
python3 scripts/parse_nll.py --hyp nll/loc/obs.sum.grid0.loc.hyp --out full/catalog_obs_combined.csv --max-gap 300 2>&1 | tail -2

echo ""
echo "BOOTSTRAP DONE: catalog_obs_combined.csv ready"
