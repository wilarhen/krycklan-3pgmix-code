#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 16_preflight_corrected_V1.R    ("Step A")
##
## SELF-CONTAINED. Assumes nothing from any previous session.
##
## PURPOSE
##   Decision point before committing compute to the full cube. Runs 3-PG at
##   the 57 plot locations with every correction applied, and answers three
##   questions:
##
##     Q1  Does simulated LAI land inside the measured range (0.28 to 6.75)?
##     Q2  Does r(GPP, APAR) stay near 0.88, i.e. is the height fix holding?
##     Q3  Does pine still reach its empirical ceiling (r = 0.634)?
##
##   If all three hold, proceed to one sector then the full cube. If any fails,
##   stop here rather than after days of compute.
##
## CORRECTIONS APPLIED
##   - foliage from Marklund fractions, NOT the branch column
##   - height_model = 2 with the Forrester et al. 2021 family (33, 37, 39)
##   - stems_n from a species mean (~1000) rather than a constant 2000
##   - inputs read DIRECTLY from source rasters, so no alignment and no grid
##     offset (point mode bypasses the resampling that cost the cube)
##   - errors logged with plot id and message, never discarded
##
## FACTORIAL
##   inputs : mapped (BIO, age rasters, NMD)  vs observed (plot measurements)
##   stems  : species mean                    vs observed
##   temps  : published parameters            vs boreal override
##
## RUN
##   Rscript 16_preflight_corrected_V1.R 2>&1 | tee preflight.txt
## ---------------------------------------------------------------------------

## SPECIES COLUMN CORRECTED: this script reads only '50 main plots',
## whose header transposes SPM and DTS. The identifiers below have been
## exchanged so that species is read from the column holding species.
suppressPackageStartupMessages({
  library(terra); library(readxl); library(dplyr); library(tidyr)
  library(r3PG); library(ggplot2)
})

## ================================================================= CONFIG ===
bp  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
rd  <- file.path(bp, '3PG_10m/x731866.017196386_y7132899.71658545.RData')
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
## Writes to results/ for the alpha-variant run, so that the
## preflight outputs in validation_out/, which the ceiling chain, 26 and 27
## read, are not overwritten.
out <- 'results'; dir.create(out, showWarnings = FALSE)

f_bio    <- file.path(bp, 'biomass/BIO.tif')
f_agepin <- file.path(bp, 'Sweden_Forest_Age/PINE_AGE_2025.tif')
f_agespr <- file.path(bp, 'Sweden_Forest_Age/SPRUCE_AGE_2025.tif')
f_spp    <- file.path(bp, 'nationalndcoverdata/nmd2018bas_ogeneraliserad_v1_0.tif')
f_dem    <- file.path(bp, 'Sweden_dem/dem.tif')

cf        <- 0.47          # dry matter -> carbon
from_date <- '2016-01'
to_date   <- '2018-12'

## Foliage as a fraction of ABOVEGROUND biomass (Marklund 1988). % VERIFY
fol_frac  <- c('Pinus sylvestris' = 0.045, 'Picea abies' = 0.085)
## Root as a fraction of aboveground, from the partitioning table
root_frac <- c('Pinus sylvestris' = 0.25/0.76, 'Picea abies' = 0.25/0.75)
## Parameter sets: all Forrester et al. 2021, all assume height_model = 2
pset      <- c('Pinus sylvestris' = 'Pinus sylvestris 39',
               'Picea abies'      = 'Picea abies 37')
## Boreal temperature override
bor <- list('Pinus sylvestris' = c(Tmin = -5, Topt = 13, Tmax = 30, alphaCx = 0.050),
            'Picea abies'      = c(Tmin = -3, Topt = 14, Tmax = 30, alphaCx = 0.050))

## The 'boreal' override above changes the temperature terms AND alphaCx in
## one step, so a bias change under it cannot be attributed to either. These
## separate them. 0.070 is the value reported for Switzerland by the WSL
## group (Trotsiuk et al., supplementary material) and is tested here as a
## sensitivity, not adopted: no calibration is performed anywhere in this
## analysis.
bor_temp <- list('Pinus sylvestris' = c(Tmin = -5, Topt = 13, Tmax = 30),
                 'Picea abies'      = c(Tmin = -3, Topt = 14, Tmax = 30))
a050 <- list('Pinus sylvestris' = c(alphaCx = 0.050),
             'Picea abies'      = c(alphaCx = 0.050))
## a070: Picea abies posterior median 0.0696 (95 % interval 0.0569-0.0832),
## Trotsiuk et al. 2020, Glob. Change Biol. 26:2463-2476, Table S3. A second
## published parameterisation for spruce. Pine keeps its published 0.0257.
a070 <- list('Picea abies'      = c(alphaCx = 0.070))

overrides <- list(boreal = bor, temp = bor_temp, a050 = a050, a070 = a070)

## measured reference ranges, for the pass/fail checks
LAI_LO <- 0.28; LAI_HI <- 6.75
CEIL_PINE <- 0.634; CEIL_SPRUCE <- 0.268

rule <- function(x) cat('\n', strrep('=', 74), '\n', x, '\n', strrep('=', 74), '\n', sep='')

## ============================================================ 1. CLIMATE ====
rule('1. CLIMATE')
nm <- load(rd); xx <- get(nm[1])
cl <- xx[[1]] |> filter(group == 'climate') |>
  select(date, variable, value) |> distinct(date, variable, .keep_all = TRUE) |>
  pivot_wider(names_from = variable, values_from = value) |> arrange(date)
rm(xx); invisible(gc())
d_climate <- cl |>
  transmute(year = as.integer(format(date, '%Y')),
            month = as.integer(format(date, '%m')),
            tmp_min, tmp_max, tmp_ave, prcp,
            srad = solar_rad, frost_days) |> as.data.frame()
cat('months:', nrow(d_climate), '|', min(d_climate$year), '-', max(d_climate$year), '\n')
stopifnot(all(2016:2018 %in% d_climate$year))

## ============================================================== 2. PLOTS ====
rule('2. PLOTS AND INPUTS (read directly from source rasters)')

obs_all <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2,
                                       .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year    = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year %in% 2016:2018)

gpp_obs <- obs_all |> transmute(Plot_ID, Year, obs_gpp = GPP - GPPff, obs_lai = LAImax)

d <- obs_all |> filter(Year == 2016) |>
  transmute(Plot_ID, Age, Trees, AGCtree, LAImax, SPM,
            sp_obs = ifelse(SPM == 1, 'Pinus sylvestris', 'Picea abies'))

cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))

d  <- inner_join(d, cc, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')

gv <- function(f) { r <- rast(f)[[1]]
                    terra::extract(r, project(pv, crs(r)), ID = FALSE)[, 1] }
d$in_bio    <- gv(f_bio)
d$in_pine   <- gv(f_agepin)
d$in_spruce <- gv(f_agespr)
d$in_nmd    <- gv(f_spp)
d$altitude  <- gv(f_dem)
d$latitude  <- crds(project(pv, 'EPSG:4326'))[, 2]

## NMD 2018: 111 pine, 112 spruce, 113 mixed conifer. % VERIFY codes
d$sp_map  <- ifelse(d$in_nmd == 111, 'Pinus sylvestris',
             ifelse(d$in_nmd == 112, 'Picea abies', d$sp_obs))
d$age_map <- ifelse(d$sp_map == 'Pinus sylvestris', d$in_pine, d$in_spruce)

cat('plots:', nrow(d), '| pine:', sum(d$sp_obs == 'Pinus sylvestris'),
    '| spruce:', sum(d$sp_obs == 'Picea abies'), '\n')
cat('species map agreement with observed:', round(100*mean(d$sp_map == d$sp_obs), 1), '%\n')
cat('plots with complete mapped inputs:',
    sum(is.finite(d$in_bio) & is.finite(d$age_map)), '/', nrow(d), '\n')

## species mean stem density, the replacement for the constant 2000
stems_mean <- tapply(d$Trees, d$sp_obs, function(x) mean(x, na.rm = TRUE))
cat('species mean stem density:',
    paste(sprintf('%s = %.0f', names(stems_mean), stems_mean), collapse = ' | '), '\n')

## ============================================================ 3. RUNNER =====
rule('3. RUNNING THE FACTORIAL')

fails <- list()

run_one <- function(i, v) {
  sp  <- if (v$inputs == 'observed') d$sp_obs[i] else d$sp_map[i]
  age <- if (v$inputs == 'observed') d$Age[i]    else d$age_map[i]
  agb <- if (v$inputs == 'observed') d$AGCtree[i]/cf else d$in_bio[i]
  stn <- if (v$stems  == 'observed') d$Trees[i]  else unname(stems_mean[sp])

  if (!is.finite(agb) || agb <= 0 || !is.finite(age) || age <= 0 ||
      !is.finite(stn) || stn <= 0 || !sp %in% names(fol_frac)) {
    fails[[length(fails)+1]] <<- data.frame(Plot_ID = d$Plot_ID[i],
      variant = v$name, msg = 'incomplete or unsupported input')
    return(NULL)
  }

  fol <- agb * fol_frac[[sp]]
  sprow <- data.frame(species = sp,
                      planted = sprintf('%04d-01', 2016 - round(age)),
                      fertility = 0.5, stems_n = stn,
                      biom_stem = agb * (1 - fol_frac[[sp]]),
                      biom_root = agb * root_frac[[sp]],
                      biom_foliage = fol, stringsAsFactors = FALSE)
  site <- data.frame(latitude = d$latitude[i], altitude = d$altitude[i],
                     soil_class = 2, asw_i = 150, asw_min = 0, asw_max = 200,
                     from = from_date, to = to_date)

  par <- get_parameters('parameters', sp_names = pset[[sp]]); names(par)[2] <- sp
  ov <- overrides[[v$temps]]          # NULL when temps == 'published'
  if (!is.null(ov)) {
    o <- ov[[sp]]
    for (k in names(o)) if (k %in% par$parameter) par[par$parameter == k, sp] <- o[[k]]
  }
  sz <- get_parameters('sizeDist', sp_names = pset[[sp]]); names(sz)[2] <- sp

  tryCatch(
    run_3PG(site = site, species = sprow, climate = d_climate, thinning = NULL,
            parameters = prepare_parameters(par, sp_names = sp),
            size_dist  = prepare_sizeDist(sz,  sp_names = sp),
            settings = list(light_model = 2, transp_model = 2, phys_model = 2,
                            height_model = 2,          ## <-- THE FIX
                            correct_bias = 0, calculate_d13c = 0),
            check_input = TRUE, df_out = TRUE),
    error = function(e) {
      fails[[length(fails)+1]] <<- data.frame(Plot_ID = d$Plot_ID[i],
        variant = v$name, msg = conditionMessage(e))
      NULL })
}

pk <- function(o, var, yr = '2017', fun = mean) {
  z <- o[o$variable == var & format(o$date, '%Y') == yr, ]
  if (!nrow(z)) return(NA_real_)
  fun(z$value, na.rm = TRUE)
}


## ---------------------------------------------------------------------------
## 36_alphacx_zero_bias.R  (generated from 16_preflight_corrected_V1.R)
##
## Diagnostic only. Sweeps alphaCx with everything else at the corrected
## configuration (height form 2, pine set 39, spruce set 37, published
## temperatures, observed inputs and stems) and finds, per species, the value
## that sets the mean bias against measured tree GPP to zero on plot means.
## The paper keeps the published parameters; this says what production level
## those parameters imply relative to measurement.
## ---------------------------------------------------------------------------

ACX <- c(0.0257, 0.027, 0.030, 0.035, 0.040, 0.045, 0.0465, 0.050,
         0.055, 0.060, 0.0635, 0.065, 0.070)
overrides <- lapply(ACX, function(a)
  list('Pinus sylvestris' = c(alphaCx = a), 'Picea abies' = c(alphaCx = a)))
names(overrides) <- sprintf('a%0.4f', ACX)
variants <- lapply(names(overrides), function(nm)
  list(name = nm, inputs = 'observed', stems = 'observed', temps = nm))

res <- list(); fails <- list()
for (v in variants) {
  rr <- list()
  for (i in seq_len(nrow(d))) {
    o <- run_one(i, v); if (is.null(o)) next
    ann <- o[o$variable == 'gpp', ] |>
      mutate(Year = as.integer(format(date, '%Y'))) |>
      group_by(Year) |> summarise(sim_gpp = sum(value) * cf * 100, .groups = 'drop')
    rr[[length(rr)+1]] <- cbind(Plot_ID = d$Plot_ID[i], variant = v$name,
                                sp = d$sp_obs[i], ann)
  }
  if (length(rr)) res[[v$name]] <- bind_rows(rr)
}
r <- bind_rows(res) |> inner_join(gpp_obs, by = c('Plot_ID','Year'))

stat <- function(z) data.frame(n = nrow(z),
  sim = mean(z$sim_gpp), obs = mean(z$obs_gpp),
  bias = mean(z$sim_gpp - z$obs_gpp),
  RMSE = sqrt(mean((z$sim_gpp - z$obs_gpp)^2)),
  r = suppressWarnings(cor(z$sim_gpp, z$obs_gpp, use = 'complete.obs')))

out <- do.call(rbind, lapply(names(res), function(nm) {
  a <- aggregate(cbind(sim_gpp, obs_gpp) ~ Plot_ID + sp, r[r$variant == nm, ], mean)
  do.call(rbind, lapply(sort(unique(a$sp)), function(s)
    cbind(alphaCx = as.numeric(sub('^a', '', nm)), species = s, stat(a[a$sp == s, ]))))
}))
out <- out[order(out$species, out$alphaCx), ]
print(out, row.names = FALSE, digits = 4)
write.csv(out, 'results/alphacx_zero_bias.csv', row.names = FALSE)

cat('\nValue of alphaCx giving zero mean bias, by linear interpolation on the grid:\n')
for (s in unique(out$species)) {
  z <- out[out$species == s, ]
  f <- approx(z$bias, z$alphaCx, xout = 0)$y
  cat(sprintf('  %-18s alphaCx = %.4f\n', s, f))
}
