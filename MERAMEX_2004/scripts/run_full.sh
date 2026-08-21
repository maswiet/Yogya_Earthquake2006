#!/usr/bin/env bash
# Unattended full-campaign picking. run_picks.py is resumable through its
# progress file, so this wrapper just keeps restarting it until the work list is
# empty: a long MPS run can be killed by the OS, and a bounded --max-items keeps
# any accumulated GPU memory from mattering.
#
#   bash run_full.sh land        # EDL + SAM, DOY 127-282
#   bash run_full.sh obs         # OBH/OBS at a lower threshold
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
WHAT="${1:-land}"
CHUNK="${CHUNK:-1500}"
MAXTRIES="${MAXTRIES:-60}"

case "$WHAT" in
  land) KINDS="EDL,SAM"; THR=0.3; EXTRA="" ;;
  obs)  KINDS="OBS";     THR=0.2; EXTRA="--duplicate-1c" ;;
  *) echo "usage: run_full.sh [land|obs]" >&2; exit 2 ;;
esac

OUT="$ROOT/full/picks_${WHAT}.csv"
PROG="$ROOT/full/done_${WHAT}.txt"
LOG="$ROOT/logs/picking_${WHAT}.log"
mkdir -p "$ROOT/full" "$ROOT/logs"

for try in $(seq 1 "$MAXTRIES"); do
  echo "=== pass $try  ($(date '+%F %T')) ===" | tee -a "$LOG"
  python "$HERE/run_picks.py" --out "$OUT" --progress "$PROG" \
      --kinds "$KINDS" --dmin 127 --dmax 282 --device mps \
      --pthr "$THR" --sthr "$THR" --dthr "$THR" $EXTRA \
      --max-items "$CHUNK" >> "$LOG" 2>&1
  rc=$?
  left=$(grep -c "still to do" /dev/null 2>/dev/null; \
         tail -200 "$LOG" | grep -oE "[0-9]+ still to do" | tail -1 | grep -oE "^[0-9]+")
  done_n=$(wc -l < "$PROG" 2>/dev/null || echo 0)
  echo "  pass $try rc=$rc | done=$done_n | remaining at pass start=${left:-?}" | tee -a "$LOG"
  if [ "${left:-1}" = "0" ]; then
    echo "ALL DONE ($WHAT): $done_n station-days" | tee -a "$LOG"
    break
  fi
done
