#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 14_missing_plots_V1.R
##
## SELF-CONTAINED. Assumes nothing from any previous session.
##
## Script 13 showed that resampling onto the offset cube grid (origin 4, -4)
## instead of the structural grid (origin 0, 0) lost information, and that the
## age layers gained 4 NA values at the plots in the process. Separately, the
## first validation run found 9 plots entirely absent from the simulated
## domain.
##
## THIS SCRIPT TESTS whether those two facts are the same fact:
##   (a) which plots have no age value in the SOURCE product
##   (b) which plots lost their age value to resampling onto the offset grid
##   (c) whether (a) + (b) accounts for the 9 missing plots
##   (d) how many plots changed land cover class through resampling
##       (correlation is meaningless for class codes, so count changes instead)
##
## RUN
##   Rscript 14_missing_plots_V1.R 2>&1 | tee missing_plots.txt
## ---------------------------------------------------------------------------

## SPECIES COLUMN CORRECTED: this script reads only '50 main plots',
## whose header transposes SPM and DTS. The identifiers below have been
## exchanged so that species is read from the column holding species.
suppressPackageStartupMessages({
  library(terra); library(readxl); library(dplyr)
})

## ================================================================= CONFIG ===
bp   <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
cube_file <- file.path(bp, 'out_3PG_netcdf/GPP.nc')
xl   <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
out  <- 'validation_out'; dir.create(out, showWarnings = FALSE)

f_bio    <- file.path(bp, 'biomass/BIO.tif')
f_agepin <- file.path(bp, 'Sweden_Forest_Age/PINE_AGE_2025.tif')
f_agespr <- file.path(bp, 'Sweden_Forest_Age/SPRUCE_AGE_2025.tif')
f_spp    <- file.path(bp, 'nationalndcoverdata/nmd2018bas_ogeneraliserad_v1_0.tif')

## plots reported missing by the first validation run
missing_known <- c(39, 74, 89, 226, 237, 313, 473, 485, 516)

rule <- function(x) cat('\n', strrep('=', 74), '\n', x, '\n', strrep('=', 74), '\n', sep='')

## ============================================================== 1. SETUP ====
rule('1. SETUP')

obs <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2,
                                   .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year    = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year == 2016) |>
  transmute(Plot_ID, Age, AGCtree, sp = ifelse(SPM == 1, 'Pine', 'Spruce'))

cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))

d  <- inner_join(obs, cc, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
cat('plots:', nrow(d), '\n')

cube <- rast(cube_file)
grid_cube <- cube[[1]]
e   <- ext(cube)
win <- ext(e$xmin - 500, e$xmax + 500, e$ymin - 500, e$ymax + 500)
cat('cube grid origin', origin(grid_cube), '| res', res(grid_cube), '\n')

## helper: extract from source and from the offset grid
pair <- function(path, method) {
  rc <- crop(rast(path)[[1]], win)
  list(src = extract(rc, pv, ID = FALSE)[, 1],
       cub = extract(project(rc, grid_cube, method = method), pv, ID = FALSE)[, 1])
}

## ================================================ 2. WHERE PLOTS ARE LOST ===
rule('2. WHERE THE PLOTS ARE LOST')

layers <- list(
  biomass    = pair(f_bio,    'bilinear'),
  age_pine   = pair(f_agepin, 'bilinear'),
  age_spruce = pair(f_agespr, 'bilinear'),
  species    = pair(f_spp,    'near'))

for (nm in names(layers)) {
  L <- layers[[nm]]
  no_src  <- d$Plot_ID[!is.finite(L$src)]
  lost_al <- d$Plot_ID[is.finite(L$src) & !is.finite(L$cub)]
  cat(sprintf('\n%s\n', toupper(nm)))
  cat('  missing at source      :', length(no_src),
      if (length(no_src))  paste0(' [', paste(sort(no_src), collapse=', '), ']') else '', '\n')
  cat('  lost through alignment :', length(lost_al),
      if (length(lost_al)) paste0(' [', paste(sort(lost_al), collapse=', '), ']') else '', '\n')
}

## the model needs ALL inputs present, so the union is what gets dropped
present_src <- Reduce(`&`, lapply(layers, function(L) is.finite(L$src)))
present_cub <- Reduce(`&`, lapply(layers, function(L) is.finite(L$cub)))

drop_src <- sort(d$Plot_ID[!present_src])
drop_cub <- sort(d$Plot_ID[!present_cub])

cat('\n--- UNION ACROSS ALL LAYERS (what reduce(intersect) would drop) ---\n')
cat('dropped if built on the SOURCE grid  :', length(drop_src),
    paste0('[', paste(drop_src, collapse=', '), ']'), '\n')
cat('dropped if built on the CUBE grid    :', length(drop_cub),
    paste0('[', paste(drop_cub, collapse=', '), ']'), '\n')
cat('avoidable by fixing the grid         :', length(setdiff(drop_cub, drop_src)),
    paste0('[', paste(sort(setdiff(drop_cub, drop_src)), collapse=', '), ']'), '\n')

## ============================================ 3. AGAINST THE KNOWN LIST =====
rule('3. COMPARISON WITH THE PLOTS ACTUALLY MISSING FROM THE CUBE')

cat('reported missing by validation run :', paste(sort(missing_known), collapse=', '), '\n')
cat('predicted missing (cube grid)      :', paste(drop_cub, collapse=', '), '\n\n')
cat('  correctly predicted :', length(intersect(drop_cub, missing_known)), '/',
    length(missing_known), '\n')
cat('  predicted but present:', paste(setdiff(drop_cub, missing_known), collapse=', '), '\n')
cat('  missing but unexplained:', paste(setdiff(missing_known, drop_cub), collapse=', '), '\n')

if (setequal(drop_cub, missing_known)) {
  cat('\n  EXACT MATCH. The plot loss is fully explained by missing input values\n')
  cat('  plus resampling onto the offset grid.\n')
} else {
  cat('\n  Not an exact match. The remainder is explained by something else,\n')
  cat('  most likely the catchment mask or a sector boundary.\n')
}

## ============================================== 4. SPECIES CLASS CHANGES ====
rule('4. LAND COVER CLASS CHANGES THROUGH RESAMPLING')

S  <- layers$species
ok <- is.finite(S$src) & is.finite(S$cub)
changed <- ok & (S$src != S$cub)
cat('plots with a changed land cover class:', sum(changed), '/', sum(ok), '\n')
if (any(changed)) {
  cat('\nplot-level changes:\n')
  print(data.frame(Plot_ID = d$Plot_ID[changed],
                   observed_sp = d$sp[changed],
                   nmd_source = S$src[changed],
                   nmd_cube   = S$cub[changed]), row.names = FALSE)
  cat('\nNMD 2018 codes: 111 pine, 112 spruce, 113 mixed conifer,\n')
  cat('114 mixed, 115+ broadleaf, 118+ non-forest classes.\n')  # VERIFY codes
  cat('A change across the forest / non-forest boundary is the serious case.\n')
}

cat('\ncontingency table (source rows, cube columns):\n')
print(table(source = S$src[ok], cube = S$cub[ok]))

## ==================================================== 5. SUMMARY ============
rule('5. SUMMARY')

summ <- data.frame(
  quantity = c('plots', 'dropped on source grid', 'dropped on cube grid',
               'avoidable by fixing grid', 'species class changed'),
  value = c(nrow(d), length(drop_src), length(drop_cub),
            length(setdiff(drop_cub, drop_src)), sum(changed)))
print(summ, row.names = FALSE)

write.csv(data.frame(Plot_ID = d$Plot_ID,
                     dropped_source = !present_src,
                     dropped_cube   = !present_cub,
                     species_changed = changed),
          file.path(out, 'plot_loss_attribution.csv'), row.names = FALSE)
cat('\nwrote', file.path(out, 'plot_loss_attribution.csv'), '\n')
