#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 06_input_output_link_V1.R
##
## Inputs are sound (BIO vs AGCtree r=0.87, ages r=0.73-0.79, DEM r=0.96) but
## the cube is uncorrelated with observed GPP. So the pipeline loses the signal
## internally. This asks the decisive question:
##
##   Does the simulated GPP still depend on the biomass and age that were fed
##   into it?
##
##   NO  -> the model washed out its initial state. Likely causes: 15-18 years
##          of spin-up from a 2001 start with uniform climate, constant FR=0.5
##          and near-uniform 1 km soil, plus the mean()-across-species summary.
##   YES -> the inputs survive but the physics maps them wrongly.
##
##   Rscript 06_input_output_link_V1.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra); library(readxl); library(dplyr); library(tidyr)
})

bp   <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
cube_file <- file.path(bp, 'out_3PG_netcdf/GPP.nc')
xl   <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
cf   <- 0.47

## ---- observed -------------------------------------------------------------
o <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2,
                                 .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year == 2017) |>
  transmute(Plot_ID,
            obs_gpp  = GPP - GPPff,
            obs_biom = AGCtree / cf,
            obs_age  = Age,
            obs_lai  = LAImax,
            obs_ba   = BA)

cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E))

d  <- inner_join(o, cc, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
cat('plots:', nrow(d), '\n')

## ---- simulated GPP, 2017 annual sum ---------------------------------------
cube <- rast(cube_file)
lt   <- time(cube)
i17  <- which(format(lt, '%Y') == '2017')
sim  <- terra::extract(cube[[i17]], pv, ID = FALSE)
d$sim_gpp <- rowSums(sim, na.rm = FALSE) * cf * 100

## ---- the model's OWN inputs at the same points ----------------------------
gv <- function(f) terra::extract(rast(f)[[1]], project(pv, crs(rast(f))),
                                 ID = FALSE)[, 1]
d$in_biom   <- gv(file.path(bp, 'biomass/BIO.tif'))
d$in_pine   <- gv(file.path(bp, 'Sweden_Forest_Age/PINE_AGE_2025.tif'))
d$in_spruce <- gv(file.path(bp, 'Sweden_Forest_Age/SPRUCE_AGE_2025.tif'))

v <- c('sim_gpp','obs_gpp','in_biom','obs_biom','in_pine','in_spruce','obs_age',
       'obs_lai','obs_ba')
m <- cor(d[, v], use = 'pairwise.complete.obs')

cat('\n================ CORRELATION MATRIX ================\n')
print(round(m, 3))

cat('\n================ THE THREE THAT MATTER =============\n')
cat(sprintf('  input biomass -> observed biomass   r = %+0.3f   (inputs are good)\n',
            m['in_biom','obs_biom']))
cat(sprintf('  observed biomass -> observed GPP    r = %+0.3f   (real ecology)\n',
            m['obs_biom','obs_gpp']))
cat(sprintf('  input biomass -> SIMULATED GPP      r = %+0.3f   <-- THE ANSWER\n',
            m['in_biom','sim_gpp']))

cat('\n  If the last one is near zero, 3-PG discarded the initial biomass it was\n')
cat('  given. With a 2001 start, uniform climate, FR fixed at 0.5 and 1 km soil,\n')
cat('  every stand converges to the same trajectory long before 2017. That is\n')
cat('  exactly the range compression seen in validation (sim 142-925 vs\n')
cat('  observed 122-1712).\n')

cat('\n================ HOW UNIFORM ARE THE DRIVERS? ======\n')
for (nm in c('ESDB-Raster-Library-1k-GeoTIFF-20240507/AWC_TOP/AWC_TOP.tif',
             'ESDB-Raster-Library-1k-GeoTIFF-20240507/TEXT/TEXT.tif')) {
  f <- file.path(bp, nm)
  if (file.exists(f)) {
    u <- unique(gv(f))
    cat(sprintf('  %-12s distinct values across the 50 plots: %s\n',
                basename(nm), paste(sort(u[!is.na(u)]), collapse = ', ')))
  }
}
cat('  A single value means that driver contributes NO spatial variation.\n')

write.csv(d, 'validation_out/input_output_link.csv', row.names = FALSE)
cat('\nwrote validation_out/input_output_link.csv\n')
