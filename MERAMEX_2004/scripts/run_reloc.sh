#!/usr/bin/env bash
# Relative relocation chain: NLLoc locations -> ph2dt -> HypoDD -> GrowClust.
#
#   TAG=full bash run_reloc.sh
#
# Optional overrides: GAPMAX NPHMIN RMSMAX ZMAX MAXSEP MINLNK MINOBS
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$ROOT/.." && pwd)"
PH2DT="$REPO/eqt/tools/HypoDD/src/ph2dt/ph2dt"
HYPODD="$REPO/eqt/tools/HypoDD/src/hypoDD/hypoDD"
GROWCLUST="$REPO/eqt/tools/GrowClust/SRC/growclust"
TAG="${TAG:-wide11}"
OUTDIR="${OUTDIR:-$ROOT/$TAG}"

echo "=== 1. build relocation inputs from nll/loc/${TAG}.* ==="
python "$HERE/gen_reloc.py" --tag "$TAG" --outroot "$ROOT" \
  --gapmax "${GAPMAX:-300}" --nphmin "${NPHMIN:-8}" --rmsmax "${RMSMAX:-0.6}" \
  --zmax "${ZMAX:-200}" --maxsep "${MAXSEP:-20}" --minlnk "${MINLNK:-6}" \
  --minobs "${MINOBS:-6}"

echo "=== 2. ph2dt: catalog differential times ==="
cd "$ROOT/hypodd"
"$PH2DT" ph2dt.inp > "$ROOT/logs/ph2dt_${TAG}.log" 2>&1
grep -E "^ (outliers|associated|P-phase|S-phase|outlier)" "$ROOT/logs/ph2dt_${TAG}.log" | head -8 || true
wc -l dt.ct event.sel | sed 's/^/  /'

echo "=== 3. HypoDD ==="
"$HYPODD" hypoDD.inp > "$ROOT/logs/hypodd_${TAG}.log" 2>&1 || {
  tail -20 "$ROOT/logs/hypodd_${TAG}.log"; exit 1; }
tail -4 "$ROOT/logs/hypodd_${TAG}.log"
python "$HERE/parse_reloc.py" --hypodd hypoDD.reloc --out "$OUTDIR/catalog_hypodd.csv"

echo "=== 4. GrowClust ==="
python "$HERE/ct2xcor.py" --ct "$ROOT/hypodd/dt.ct" \
  --out "$ROOT/growclust/IN/xcordata.txt"
cd "$ROOT/growclust"
"$GROWCLUST" growclust.inp > "$ROOT/logs/growclust_${TAG}.log" 2>&1 || {
  tail -20 "$ROOT/logs/growclust_${TAG}.log"; exit 1; }
tail -4 "$ROOT/logs/growclust_${TAG}.log"
python "$HERE/parse_reloc.py" --growclust OUT/out.growclust_cat \
  --out "$OUTDIR/catalog_growclust.csv"

echo "RELOC DONE ($TAG)"
