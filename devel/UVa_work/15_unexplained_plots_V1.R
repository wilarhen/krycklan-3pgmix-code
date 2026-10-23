#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 15_unexplained_plots_V1.R
##
## SELF-CONTAINED.
##
## Script 14 attributed 4 of the 9 missing plots to absent input values plus
## the grid offset. Five remain: 39, 74, 89, 226, 485. They have complete
## values in all four structural layers and were still not simulated.
##
## Remaining candidates, in order of likelihood:
##   (a) outside the catchment polygon used for masking
##   (b) missing soil or climate inputs
##   (c) run_3PG failed and the error was discarded by tryCatch
##
## (c) cannot be tested from the outputs, only inferred by elimination. If (a)
## and (b) are both cleared, add error logging before the rerun.
##
## RUN
##   Rscript 15_unexplained_plots_V1.R 2>&1 | tee unexplained.txt
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra); library(readxl); library(dplyr)
})

bp  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
unexplained <- c(39, 74, 89, 226, 485)
missing_all <- c(39, 74, 89, 226, 237, 313, 473, 485, 516)

rule <- function(x) cat('\n', strrep('=', 70), '\n', x, '\n', strrep('=', 70), '\n', sep='')

## ---- plots -----------------------------------------------------------------
cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))
pv  <- vect(cbind(cc$E, cc$N), type = 'points', crs = 'EPSG:3006')
sel <- cc$Plot_ID %in% unexplained
cat('plots with coordinates:', nrow(cc), '| unexplained:', sum(sel), '\n')

## ---- (a) catchment mask ----------------------------------------------------
rule('(a) ARE THEY INSIDE THE CATCHMENT POLYGON?')

shp <- list.files(file.path(bp, 'krycklan_borders'), '\\.shp$',
                  full.names = TRUE, recursive = TRUE)
if (!length(shp)) cat('no shapefile found under krycklan_borders\n')
for (s in shp) {
  kb <- tryCatch(project(vect(s), 'EPSG:3006'), error = function(e) NULL)
  if (is.null(kb)) { cat(basename(s), ': unreadable\n'); next }
  inside <- !is.na(relate(pv, kb, 'intersects')[, 1]) &
             relate(pv, kb, 'intersects')[, 1]
  cat(sprintf('%-34s all plots inside %2d/%2d | unexplained inside %d/%d | missing inside %d/%d\n',
              basename(s), sum(inside), nrow(cc),
              sum(inside[sel]), sum(sel),
              sum(inside[cc$Plot_ID %in% missing_all]), length(missing_all)))
  if (any(!inside[sel]))
    cat('    OUTSIDE:', paste(cc$Plot_ID[sel & !inside], collapse = ', '), '\n')
}

## ---- (b) soil coverage -----------------------------------------------------
rule('(b) DO THE SOIL LAYERS COVER THEM?')

soil_dirs <- c('ESDB-Raster-Library-1k-GeoTIFF-20240507/TEXT',
               'ESDB-Raster-Library-1k-GeoTIFF-20240507/AWC_TOP')
for (sd in soil_dirs) {
  f <- list.files(file.path(bp, sd), '\\.tif$', full.names = TRUE)
  if (!length(f)) { cat(sd, ': no tif found\n'); next }
  r <- rast(f[1])[[1]]
  v <- tryCatch(extract(r, project(pv, crs(r)), ID = FALSE)[, 1],
                error = function(e) rep(NA_real_, nrow(cc)))
  cat(sprintf('%-12s NA at unexplained %d/%d | NA at all plots %d/%d | values: %s\n',
              basename(sd), sum(!is.finite(v[sel])), sum(sel),
              sum(!is.finite(v)), nrow(cc),
              paste(sort(unique(v[is.finite(v)])), collapse = ', ')))
}

## ---- (c) elimination -------------------------------------------------------
rule('(c) VERDICT')

cat('If the unexplained plots are inside the polygon and have soil values,\n')
cat('the remaining candidate is model failure discarded by tryCatch in\n')
cat('polygon_grids_v01.r. That cannot be tested from the saved outputs.\n\n')
cat('ACTION FOR THE RERUN: replace\n')
cat('    tryCatch(..., error = function(e) NULL)\n')
cat('with a version that records the cell id and the error message, e.g.\n\n')
cat('    tryCatch(run_3PG(...),\n')
cat('             error = function(e) {\n')
cat('               fails[[length(fails)+1]] <<- data.frame(\n')
cat('                 cell = nm, msg = conditionMessage(e))\n')
cat('               NULL\n')
cat('             })\n\n')
cat('and write the failure table alongside the cube. The fraction of cells\n')
cat('that simulated successfully is a number reviewers will ask for, and it\n')
cat('is needed for the Results section regardless.\n')
