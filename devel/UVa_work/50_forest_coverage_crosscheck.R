## ---------------------------------------------------------------------------
## 50_forest_coverage_crosscheck.R
##
## Why 4 349 of the 4 803 ha of mapped forest in the Krycklan parcel carry
## simulated values. Crosses every mapped-forest
## cell (NMD 111-117, 121-127, cell centre inside the largest polygon of
## krycklan_borders.shp) with the corrected cube's own coverage (gpp.tif,
## band 1, non-NA) and with the inputs the pipeline requires
## (polygon_grids_v03.r: a cell is simulated only if, for at least one species
## of its NMD class, that species' age raster has a value, and biomass, site
## and climate inputs are complete). Read-only on all inputs.
##
## Classes written: simulated; no stand age for any species of the cell; the
## residue (age present but not simulated), with checks on biomass, DEM and
## zero biomass. Also tests whether the cells without an age are the biomass
## product's structural zeros (value 0 where mean height < 3 m).
##
## Outputs (in the directory given as the first argument, default results/):
##   forest_coverage_crosscheck.csv           counts and hectares per class
##   forest_coverage_crosscheck_by_class.csv  the same by NMD class
## First run 17 Oct 2026 on iufor: 480 357 forest cells (4 803.57 ha);
## 434 938 simulated (4 349.38 ha); 35 815 without an age (358.15 ha), of which
## 8 144 have biomass 0; residue 9 604 (96.04 ha), not explained.
##
##   Rscript 50_forest_coverage_crosscheck.R [out_dir]   (on iufor; B = ~/datos/3PG)
## ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(terra))
B <- '/home/wilar/datos/3PG'
out <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else 'results'
g <- rast(file.path(B, 'outputs_3PG', 'cube', 'gpp.tif'))[[1]]
parts <- disagg(vect(file.path(B, 'krycklan_borders', 'krycklan_borders.shp')))
kry <- project(parts[which.max(expanse(parts, unit = 'ha'))], crs(g))
on <- function(r, m = 'near') { r <- crop(r, project(as.polygons(ext(g), crs = crs(g)), crs(r)), snap = 'out')
  if (!compareGeom(r, g, stopOnError = FALSE)) r <- resample(project(r, crs(g), method = m), g, method = m); r }
inside <- values(!is.na(mask(rast(g, vals = 1), kry)), mat = FALSE) %in% TRUE
nmd <- values(on(rast(file.path(B, 'nationalndcoverdata', 'nmd2018bas_ogeneraliserad_v1_0.tif'))), mat = FALSE)
A <- file.path(B, 'Sweden_Forest_Age')
ap <- values(on(rast(file.path(A, 'PINE_AGE_2025.tif'))), mat = FALSE)
as <- values(on(rast(file.path(A, 'SPRUCE_AGE_2025.tif'))), mat = FALSE)
ab <- values(on(rast(file.path(A, 'Birch_ALS80.tif'))), mat = FALSE)
bio <- values(on(rast(file.path(B, 'biomass', 'BIO.tif'))), mat = FALSE)
dem <- values(on(rast(file.path(B, 'Sweden_dem', 'dem.tif')), 'average'), mat = FALSE)
sim <- !is.na(values(g, mat = FALSE))
lk <- read.csv(file.path(B, 'nationalndcoverdata', 'spp_codes_from_land.csv'))
for (k in 1:3) lk[[paste0('Layer_', k)]] <- trimws(lk[[paste0('Layer_', k)]])
spk <- function(k) { x <- lk[[paste0('Layer_', k)]][match(nmd, lk$code)]; x[x %in% c('NA', '')] <- NA; x }
age_for <- function(s) ifelse(s == 'Pinus sylvestris', ap, ifelse(s == 'Picea abies', as, ifelse(s == 'Betula', ab, NA)))
slot_age <- Reduce(`|`, lapply(1:3, function(k) { s <- spk(k); !is.na(s) & is.finite(age_for(s)) }))
any_age <- is.finite(ap) | is.finite(as) | is.finite(ab)
forest <- inside & nmd %in% c(111:117, 121:127)
u <- forest & !sim; noage <- u & !slot_age; res <- u & slot_age; zero <- is.finite(bio) & bio == 0
row <- function(item, m) data.frame(item = item, cells = sum(m), ha = sum(m) / 100)
S <- rbind(row('mapped forest', forest), row('simulated', forest & sim), row('not simulated', u),
           row('not simulated: no stand age for any species of the cell', noage),
           row('  of which no value in any age raster', u & !any_age),
           row('  of which biomass 0 (structural zero below 3 m)', noage & zero),
           row('  of which biomass > 0', noage & is.finite(bio) & bio > 0),
           row('not simulated: residue, stand age present', res),
           row('  residue lacking biomass', res & !is.finite(bio)), row('  residue lacking DEM', res & !is.finite(dem)),
           row('  residue with biomass 0', res & zero),
           row('check: simulated forest without a stand age', forest & sim & !slot_age),
           row('check: simulated non-forest inside the parcel', inside & sim & !forest),
           row('all forest cells with biomass 0', forest & zero))
print(S, row.names = FALSE)
cls <- sort(unique(nmd[forest]))
C <- do.call(rbind, lapply(cls, function(c) data.frame(nmd = c, forest = sum(forest & nmd == c), simulated = sum(forest & sim & nmd == c),
  no_age = sum(noage & nmd == c), residue = sum(res & nmd == c))))
print(C, row.names = FALSE)
write.csv(S, file.path(out, 'forest_coverage_crosscheck.csv'), row.names = FALSE)
write.csv(C, file.path(out, 'forest_coverage_crosscheck_by_class.csv'), row.names = FALSE)
