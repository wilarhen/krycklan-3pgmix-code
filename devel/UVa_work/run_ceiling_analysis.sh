#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_ceiling_analysis.sh
#
# Runs the three analyses that answer the "this is just an input evaluation"
# objection. All three are read-only with respect to the cube; they only need
# the plot-level extraction table.
#
#   1. edit ceiling_config.R  (PLOT_DATA, COL, PREDICTORS)
#   2. bash run_ceiling_analysis.sh
#
# Options:
#   bash run_ceiling_analysis.sh 23     run only script 23
#   bash run_ceiling_analysis.sh check  validate the config, run nothing
# ---------------------------------------------------------------------------
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
mkdir -p logs

command -v Rscript >/dev/null || { echo "Rscript not on PATH"; exit 1; }
[ -f ceiling_config.R ] || { echo "ceiling_config.R not found in $HERE"; exit 1; }

# ---- config check ----------------------------------------------------------
echo "checking ceiling_config.R ..."
Rscript -e '
  source("ceiling_config.R")
  d <- load_plots()
  p <- to_plot_level(d)
  cat("plot-years:", nrow(d), " plots:", nrow(p), "\n")
  print(table(p$species))
  cat("columns resolved OK\n")
' || { echo; echo "config check FAILED -- fix ceiling_config.R and rerun"; exit 1; }

[ "${1:-}" = "check" ] && { echo "check only, stopping."; exit 0; }

# ---- run -------------------------------------------------------------------
SCRIPTS=(23_ceiling_multivariate.R 24_ceiling_inference.R 25_light_saturation.R)
if [ $# -gt 0 ]; then
    SCRIPTS=()
    for a in "$@"; do
        m=$(ls "${a}"_*.R 2>/dev/null | head -1)
        [ -n "$m" ] && SCRIPTS+=("$m") || echo "no script matching '$a'"
    done
fi

fail=0
for s in "${SCRIPTS[@]}"; do
    log="logs/${s%.R}.log"
    echo
    echo "=============================================================="
    echo "running $s   -> $log"
    echo "=============================================================="
    if Rscript "$s" 2>&1 | tee "$log"; then :; else fail=1; echo "FAILED: $s"; fi
done

echo
echo "outputs in the OUT_DIR set in ceiling_config.R; logs in $HERE/logs"
exit $fail
