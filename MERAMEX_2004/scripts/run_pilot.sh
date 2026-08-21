#!/usr/bin/env bash
# End-to-end MERAMEX chain: picks -> events -> NonLinLoc -> catalog -> figure.
# Assumes run_picks.py has already produced $ROOT/$WORK/picks.csv.
#
#   WORK=wide11 TAG=wide11 TITLE="..." bash run_pilot.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$ROOT/.." && pwd)"
NLLBIN="$REPO/eqt/tools/NonLinLoc/src/bin"
TAG="${TAG:-pilot}"
WORK="${WORK:-pilot}"
NXY="${NXY:-601}"
NDIST="${NDIST:-801}"
NZ="${NZ:-254}"
MAXGAP="${MAXGAP:-300}"
TITLE="${TITLE:-MERAMEX 2004 $WORK}"

source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh

echo "=== 1. associate (pyocto, env assoc) ==="
conda activate assoc
python "$HERE/associate.py" \
  --picks   "$ROOT/$WORK/picks.csv" \
  --out     "$ROOT/$WORK/events.csv" \
  --tmpdir  "$ROOT/$WORK/octo_tmp" \
  --n_picks "${N_PICKS:-8}" --n_p_and_s "${N_PS:-3}" --zmax "${ZMAX:-200}"
conda deactivate

echo "=== 2. NonLinLoc control files ==="
python "$HERE/gen_nll.py" \
  --events      "$ROOT/$WORK/events.csv" \
  --assignments "$ROOT/$WORK/events_assignments.csv" \
  --outdir      "$ROOT/nll" --tag "$TAG" --nxy "$NXY" --ndist "$NDIST" --nz "$NZ"

echo "=== 3. travel-time grids ==="
cd "$ROOT/nll"
"$NLLBIN/Vel2Grid"  "nll_vel_${TAG}.in"    > "$ROOT/logs/vel_P.log"  2>&1
"$NLLBIN/Vel2Grid"  "nll_vel_${TAG}_S.in"  > "$ROOT/logs/vel_S.log"  2>&1
"$NLLBIN/Grid2Time" "nll_time_${TAG}_P.in" > "$ROOT/logs/gt_P.log"   2>&1
"$NLLBIN/Grid2Time" "nll_time_${TAG}_S.in" > "$ROOT/logs/gt_S.log"   2>&1
ls time/ | wc -l | xargs echo "  travel-time grid files:"

echo "=== 4. locate ==="
"$NLLBIN/NLLoc" "nll_loc_${TAG}.in" > "$ROOT/logs/nlloc.log" 2>&1
tail -3 "$ROOT/logs/nlloc.log"

echo "=== 5. catalog + figure ==="
cd "$HERE"
python parse_nll.py --hyp "$ROOT/nll/loc/${TAG}.sum.grid0.loc.hyp" \
                    --out "$ROOT/$WORK/catalog_nll.csv" --max-gap "$MAXGAP"
python plot_pilot.py --catalog "$ROOT/$WORK/catalog_nll.csv" \
                     --stations "$ROOT/$WORK/events_stations.csv" \
                     --out "$ROOT/figures/${WORK}_summary.png" --title "$TITLE"
echo "CHAIN DONE ($WORK)"
