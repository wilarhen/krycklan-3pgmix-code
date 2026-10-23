#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 13_alignment_audit_V1.R
##
## SELF-CONTAINED. Assumes nothing from any previous session.
##
## QUESTION
##   The cube grid has origin (4, -4). The three structural products (biomass,
##   stand age, species) share origin (0, 0) at 10 m. The model grid is
##   therefore offset by 4 m from the grid its information-bearing inputs are
##   natively on, so every cell was interpolated from four source cells instead
##   of being read directly.
##
##   This script measures what that cost, and whether adopting a structural
##   product as the reference grid removes it.
##
## WHAT IT REPORTS
##   1. Grid geometry of every layer and of the cube
##   2. Loss incurred on the ACTUAL cube grid, at the plot locations
##   3. Loss on the PROPOSED structural grid (expected: none)
##   4. Effect of bilinear vs nearest resampling for age and species
##   5. Propagation: does aligned biomass still track measured biomass?
##
## RUN
##   Rscript 13_alignment_audit_V1.R 2>&1 | tee alignment_audit.txt
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
f_dem    <- file.path(bp, 'Sweden_dem/dem.tif')

rule <- function(x) cat('\n', strrep('=', 74), '\n', x, '\n', strrep('=', 74), '\n', sep='')

## ============================================================== 1. PLOTS ====
rule('1. VALIDATION PLOTS')

obs <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2,
                                   .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year    = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year == 2016) |>
  transmute(Plot_ID, AGCtree, Age, LAImax,
            sp = ifelse(SPM == 1, 'Pine', 'Spruce'))

cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))

d  <- inner_join(obs, cc, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
cat('plots:', nrow(d), '| pine:', sum(d$sp == 'Pine'),
    '| spruce:', sum(d$sp == 'Spruce'), '\n')

## ==================================================== 2. GRID GEOMETRY ======
rule('2. GRID GEOMETRY OF EVERY LAYER')

cube <- rast(cube_file)
lyrs <- list(cube = cube[[1]],
             biomass = rast(f_bio)[[1]],
             age_pine = rast(f_agepin)[[1]],
             age_spruce = rast(f_agespr)[[1]],
             species = rast(f_spp)[[1]],
             dem = rast(f_dem)[[1]])

geo <- do.call(rbind, lapply(names(lyrs), function(n) {
  r <- lyrs[[n]]
  data.frame(layer = n,
             origin_x = origin(r)[1], origin_y = origin(r)[2],
             res_x = res(r)[1], res_y = res(r)[2],
             crs = crs(r, describe = TRUE)$code,
             stringsAsFactors = FALSE)
}))
print(geo, row.names = FALSE, digits = 6)

cat('\nOFFSET of the cube grid from the structural grid:\n')
cat(sprintf('  x: %.2f m   y: %.2f m\n',
            origin(cube)[1] - origin(lyrs$biomass)[1],
            origin(cube)[2] - origin(lyrs$biomass)[2]))
cat('  A non-zero offset means every structural value fed to the model was\n')
cat('  interpolated from four source cells rather than read directly.\n')

## ======================================== 3. BUILD THE TWO TARGET GRIDS =====
rule('3. TARGET GRIDS')

## (a) the grid the cube actually used
grid_cube <- cube[[1]]

## (b) the proposed grid: same coverage, but on the structural origin (0,0)
e  <- ext(cube)
e2 <- ext(floor(e$xmin/10)*10, ceiling(e$xmax/10)*10,
          floor(e$ymin/10)*10, ceiling(e$ymax/10)*10)
grid_struct <- rast(e2, resolution = 10, crs = 'EPSG:3006')

cat('grid_cube   origin', origin(grid_cube),   ' res', res(grid_cube),   '\n')
cat('grid_struct origin', origin(grid_struct), ' res', res(grid_struct), '\n')

## crop window with margin, in the source CRS
win <- ext(e$xmin - 500, e$xmax + 500, e$ymin - 500, e$ymax + 500)

## ============================================== 4. MEASURE THE LOSS =========
rule('4. RESAMPLING LOSS AT THE PLOT LOCATIONS')

measure <- function(path, label, target, target_label, method) {
  r  <- rast(path)[[1]]
  rc <- crop(r, win)
  a  <- project(rc, target, method = method)
  v0 <- extract(rc, pv, ID = FALSE)[, 1]      # source, no resampling
  v1 <- extract(a,  pv, ID = FALSE)[, 1]      # after alignment
  ok <- is.finite(v0) & is.finite(v1)
  if (sum(ok) < 5) return(NULL)
  rr <- cor(v0[ok], v1[ok])
  data.frame(layer = label, grid = target_label, method = method,
             n = sum(ok), r = rr, var_retained = 100 * rr^2,
             mean_abs_diff = mean(abs(v0[ok] - v1[ok])),
             max_abs_diff = max(abs(v0[ok] - v1[ok])),
             na_gained = sum(is.finite(v0) & !is.finite(v1)),
             stringsAsFactors = FALSE)
}

jobs <- list(
  list(f_bio,    'biomass',    'bilinear'),
  list(f_agepin, 'age_pine',   'bilinear'),
  list(f_agepin, 'age_pine',   'near'),
  list(f_agespr, 'age_spruce', 'bilinear'),
  list(f_agespr, 'age_spruce', 'near'),
  list(f_spp,    'species',    'near'))

res_all <- list()
for (j in jobs) {
  res_all[[length(res_all)+1]] <- measure(j[[1]], j[[2]], grid_cube,   'cube (offset)',   j[[3]])
  res_all[[length(res_all)+1]] <- measure(j[[1]], j[[2]], grid_struct, 'structural (0,0)', j[[3]])
}
tab <- bind_rows(res_all)
print(as.data.frame(tab), row.names = FALSE, digits = 4)

cat('\nREAD THIS AS\n')
cat('  grid = structural: r should be 1.000 (no resampling needed)\n')
cat('  grid = cube      : r < 1 quantifies what the 4 m offset cost\n')
cat('  bilinear vs near on age: bilinear invents ages at stand boundaries\n')

## ================================= 5. DOES IT REACH THE VALIDATION? =========
rule('5. PROPAGATION TO THE MEASURED STAND ATTRIBUTES')

prop_check <- function(path, label, obs_vec, obs_label, method = 'bilinear') {
  r  <- rast(path)[[1]]
  rc <- crop(r, win)
  v_src  <- extract(rc, pv, ID = FALSE)[, 1]
  v_cube <- extract(project(rc, grid_cube,   method = method), pv, ID = FALSE)[, 1]
  v_str  <- extract(project(rc, grid_struct, method = method), pv, ID = FALSE)[, 1]
  f <- function(v) { ok <- is.finite(v) & is.finite(obs_vec)
                     if (sum(ok) < 5) NA_real_ else cor(v[ok], obs_vec[ok]) }
  cat(sprintf('%-12s vs %-10s  source r=%+.4f | cube grid r=%+.4f | struct grid r=%+.4f\n',
              label, obs_label, f(v_src), f(v_cube), f(v_str)))
}

cat('The quantity that matters: does alignment weaken the relationship between\n')
cat('an input layer and the stand attribute it is meant to represent?\n\n')
prop_check(f_bio,    'biomass',    d$AGCtree / 0.47, 'AGCtree')
prop_check(f_agepin, 'age_pine',   d$Age,            'Age')
prop_check(f_agespr, 'age_spruce', d$Age,            'Age')

## =============================================== 6. VERDICT ================
rule('6. VERDICT')

tc <- tab[tab$grid == 'cube (offset)' & tab$method == 'bilinear', ]
ts <- tab[tab$grid == 'structural (0,0)' & tab$method == 'bilinear', ]
if (nrow(tc)) {
  cat('Variance retained on the cube grid (bilinear):\n')
  for (i in seq_len(nrow(tc)))
    cat(sprintf('  %-12s %.1f %%\n', tc$layer[i], tc$var_retained[i]))
}
if (nrow(ts) && all(ts$r > 0.9999)) {
  cat('\nOn the structural grid all layers are reproduced exactly.\n')
  cat('CONCLUSION: adopting a structural product as the reference grid removes\n')
  cat('the loss entirely. The elevation product aggregates onto that grid by an\n')
  cat('integer factor (2 m -> 10 m), so nothing is lost there either.\n')
} else {
  cat('\nStructural grid did NOT reproduce exactly. Inspect the table above.\n')
}

write.csv(tab, file.path(out, 'alignment_audit.csv'), row.names = FALSE)
write.csv(geo, file.path(out, 'grid_geometry.csv'), row.names = FALSE)
cat('\nwrote', file.path(out, 'alignment_audit.csv'), 'and grid_geometry.csv\n')

## ---------------------------------------------------------------------------
## IF THE LOSS IS CONFIRMED, THE FIX IS ONE LINE IN polygon_grids_v01.r
##
##   ## before
##   kryck_elev <- agg_and_crop.(dem, 10, kryck_bord)      # 2 m DEM -> offset grid
##   ref <- kryck_elev
##
##   ## after
##   ref <- crop(rast('.../biomass/BIO.tif'), kryck_bord)  # structural grid (0,0)
##   kryck_elev <- align_to_ref(dem, ref, method = 'bilinear')
##
## And in internals.r, align_to_ref currently resamples twice
## (project onto a default grid, then resample onto ref). One step suffices:
##
##   align_to_ref <- function(r, ref, method = 'bilinear') {
##     if (terra::compareGeom(r, ref, stopOnError = FALSE)) return(r)
##     terra::project(r, ref, method = method)
##   }
##
## Use method = 'near' for age and species. Bilinear interpolation of a stand
## age surface across a clearcut boundary produces ages that exist nowhere.
## ---------------------------------------------------------------------------
