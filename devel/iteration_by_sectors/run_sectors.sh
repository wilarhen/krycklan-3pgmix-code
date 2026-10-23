#!/usr/bin/env bash
## ---------------------------------------------------------------------------
## run_sectors.sh
##
## Runs polygon_grids_v02.r over every sector, one R process per sector.
##
##   ./run_sectors.sh                 all 100 sectors, 10 per side
##   ./run_sectors.sh 1 100           explicit range
##   ./run_sectors.sh 62 62 TRUE      one sector, saving the raw RData
##
## Each sector gets its own log under logs/. Completed sectors are skipped on
## re-run, so an interrupted run resumes where it stopped.
## ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT=${SCRIPT:-polygon_grids_v02.r}
CELLS_PER_SIDE=${CELLS_PER_SIDE:-10}
FROM=${1:-1}
TO=${2:-$((CELLS_PER_SIDE * CELLS_PER_SIDE))}
SAVE_RDATA=${3:-FALSE}

LOGDIR=logs
DONEDIR=logs/.done
mkdir -p "$LOGDIR" "$DONEDIR"

echo "sectors $FROM to $TO | $CELLS_PER_SIDE per side | save_rdata=$SAVE_RDATA"
echo "logs in $LOGDIR/"
START=$(date +%s)
ok=0; skipped=0; failed=0; empty=0

for i in $(seq "$FROM" "$TO"); do
  stamp=$(printf '%03d' "$i")

  if [ -f "$DONEDIR/sector_${stamp}.done" ]; then
    skipped=$((skipped+1))
    continue
  fi

  echo "--- sector $i / $TO  ($(date +%H:%M:%S)) ---"
  t0=$(date +%s)

  if Rscript "$SCRIPT" "$i" "$CELLS_PER_SIDE" "$SAVE_RDATA" \
       > "$LOGDIR/sector_${stamp}.log" 2>&1; then
    t1=$(date +%s)
    if grep -q "produced no simulated cells" "$LOGDIR/sector_${stamp}.log"; then
      echo "    empty (outside catchment), $((t1-t0))s"
      empty=$((empty+1))
    else
      echo "    ok, $((t1-t0))s"
      ok=$((ok+1))
    fi
    touch "$DONEDIR/sector_${stamp}.done"
  else
    echo "    FAILED, see $LOGDIR/sector_${stamp}.log"
    tail -5 "$LOGDIR/sector_${stamp}.log" | sed 's/^/      /'
    failed=$((failed+1))
  fi
done

END=$(date +%s)
echo
echo "=========================================================="
printf 'done in %d min | ok %d | empty %d | skipped %d | failed %d\n' \
       $(( (END-START)/60 )) "$ok" "$empty" "$skipped" "$failed"
echo "=========================================================="

## completeness across all sectors, from the per-sector logs
echo
echo "cells simulated per sector:"
grep -h "cells simulated" "$LOGDIR"/sector_*.log 2>/dev/null | sed 's/^/  /'

## collect every failure log written by FIX 8
FAILCSV=$(find . -name 'failures_*.csv' 2>/dev/null | wc -l)
if [ "$FAILCSV" -gt 0 ]; then
  echo
  echo "model failure logs written: $FAILCSV"
  echo "combine with:"
  echo "  Rscript -e \"f<-list.files(pattern='failures_.*csv',recursive=TRUE,full.names=TRUE);"
  echo "               d<-do.call(rbind,lapply(f,read.csv));"
  echo "               print(sort(table(d\\\$msg),decreasing=TRUE)); nrow(d)\""
fi

## to re-run everything from scratch:  rm -rf logs/.done
