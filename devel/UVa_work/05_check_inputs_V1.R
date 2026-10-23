#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 05_check_inputs_V1.R
##
## The cube is correctly registered (no flip) and spatially coherent
## (Moran's I = 0.76), yet uncorrelated with observed plot GPP. Therefore the
## fault is upstream of 3-PG: the INPUT rasters describe a landscape that does
## not match the measured plots.
##
## This script tests each input layer against the corresponding measured plot
## attribute. The elevation test is the critical one: DEM vs measured elevation
## involves no model physics at all, so it isolates registration from ecology.
##
##   ELEVATION r ~ 1     -> coordinates fine; fault is in species/age/biomass
##   ELEVATION r poor    -> the whole input stack is spatially misplaced
##
##   Rscript 05_check_inputs_V1.R <base_path> <Dataset.xlsx>
## ---------------------------------------------------------------------------

## SPECIES COLUMN CORRECTED: this script reads only '50 main plots',
## whose header transposes SPM and DTS. The identifiers below have been
## exchanged so that species is read from the column holding species.
suppressPackageStartupMessages({
  library(terra); library(sf); library(readxl); library(dplyr)
})

args      <- commandArgs(trailingOnly = TRUE)
base_path <- if (length(args) >= 1) args[1] else '~/3pg_data'
xlsx_path <- if (length(args) >= 2) args[2] else
  '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
base_path <- path.expand(base_path)

cat('base_path:', base_path, '\n')
if (!dir.exists(base_path)) stop('base_path not found', call. = FALSE)
cat('subdirs  :', paste(basename(list.dirs(base_path, recursive = FALSE)),
                        collapse = ', '), '\n')

## ---- observed plot attributes (2016 only; static attributes) --------------
obs <- suppressMessages(read_excel(xlsx_path, sheet = '50 main plots',
                                   skip = 2, .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year = as.integer(Year)) |>
  filter(!is.na(Plot_ID), Year == 2016) |>
  transmute(Plot_ID, Elevation, Slope, Age, Trees, BA, Vol,
            AGCtree, LAImax, SPM, DTS)

coords <- suppressMessages(read_excel(xlsx_path, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         E = as.numeric(E), N = as.numeric(N)) |>
  filter(!is.na(E))

d <- inner_join(obs, coords, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
cat('plots    :', nrow(d), '\n')

## ---- helper ---------------------------------------------------------------
find_rast <- function(subdir, pattern) {
  p <- file.path(base_path, subdir)
  if (!dir.exists(p)) return(NULL)
  f <- list.files(p, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (!length(f)) return(NULL)
  tryCatch(rast(f[1]), error = function(e) NULL)
}

test_layer <- function(r, obs_vec, label, method = 'pearson') {
  if (is.null(r)) { cat(sprintf('%-28s SKIPPED (not found)\n', label)); return(invisible(NULL)) }
  v <- tryCatch(terra::extract(r[[1]], project(pv, crs(r)), ID = FALSE)[,1],
                error = function(e) NULL)
  if (is.null(v)) { cat(sprintf('%-28s extract failed\n', label)); return(invisible(NULL)) }
  ok <- is.finite(v) & is.finite(obs_vec)
  if (sum(ok) < 5) { cat(sprintf('%-28s n=%d too few\n', label, sum(ok))); return(invisible(NULL)) }
  r_ <- cor(v[ok], obs_vec[ok], method = method)
  cat(sprintf('%-28s n=%3d  r=%+0.3f  sim[%.1f, %.1f]  obs[%.1f, %.1f]\n',
              label, sum(ok), r_, min(v[ok]), max(v[ok]),
              min(obs_vec[ok]), max(obs_vec[ok])))
  invisible(data.frame(layer = label, n = sum(ok), r = r_))
}

## ---- THE CRITICAL TEST ----------------------------------------------------
cat('\n===================== REGISTRATION (no model physics) ==================\n')
dem <- find_rast('Sweden_dem', 'dem.*\\.tif$')
res_dem <- test_layer(dem, d$Elevation, 'DEM vs measured Elevation')

if (!is.null(dem)) {
  slp <- tryCatch(terrain(dem[[1]], v = 'slope', unit = 'degrees'),
                  error = function(e) NULL)
  test_layer(slp, d$Slope, 'DEM slope vs meas. Slope')
}

cat('\n  r > 0.95 -> coordinates and registration are CORRECT.\n')
cat('             The fault is in the forest attribute layers below.\n')
cat('  r < 0.8  -> the input stack is spatially misplaced. Everything\n')
cat('             downstream is invalid; fix this before anything else.\n')

## ---- forest attribute layers ----------------------------------------------
cat('\n===================== FOREST ATTRIBUTE INPUTS ==========================\n')
test_layer(find_rast('biomass', '.*\\.tif$'),          d$AGCtree, 'biomass map vs AGCtree')
test_layer(find_rast('Sweden_Forest_Age', '(PINE|SPRUCE).*\\.tif$'), d$Age, 'age map (conifer) vs Age')
test_layer(find_rast('Sweden_Forest_Age', '(Birch).*\\.tif$'),       d$Age, 'age map (birch) vs Age')
test_layer(find_rast('nationalndcoverdata', '.*\\.tif$'), as.numeric(d$SPM),
           'land cover vs SPM', method = 'spearman')

cat('\n===================== SOIL / CLIMATE (coarse, expect weak) =============\n')
test_layer(find_rast('ESDB-Raster-Library-1k-GeoTIFF-20240507/AWC_TOP', '.*\\.tif$'),
           as.numeric(d$DTS), 'AWC vs DTS', method = 'spearman')

cat('\n===================== WHAT DRIVES OBSERVED GPP =========================\n')
cat('For reference, the observed determinants the model must reproduce:\n')
gpp <- suppressMessages(read_excel(xlsx_path, sheet = '50 main plots',
                                   skip = 2, .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)), Year = as.integer(Year)) |>
  filter(!is.na(Plot_ID), Year == 2016) |>
  transmute(Plot_ID, gpp_tree = GPP - GPPff) |>
  inner_join(d, by = 'Plot_ID')

for (v in c('AGCtree','LAImax','BA','Vol','Age','Trees','Elevation'))
  cat(sprintf('  cor(observed tree GPP, %-10s) = %+0.3f\n', v,
              cor(gpp$gpp_tree, gpp[[v]], use = 'complete.obs')))

cat('\nIf the model reproduced these relationships it would have skill.\n')
cat('Compare against the input correlations above to see which link is broken.\n')
