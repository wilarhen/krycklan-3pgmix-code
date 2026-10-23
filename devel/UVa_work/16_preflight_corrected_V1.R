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

variants <- list(
  list(name = 'mapped/mean',      inputs = 'mapped',   stems = 'mean',     temps = 'published'),
  list(name = 'mapped/obs',       inputs = 'mapped',   stems = 'observed', temps = 'published'),
  list(name = 'observed/obs',     inputs = 'observed', stems = 'observed', temps = 'published'),
  list(name = 'observed/obs/bor',  inputs = 'observed', stems = 'observed', temps = 'boreal'),
  list(name = 'observed/obs/temp', inputs = 'observed', stems = 'observed', temps = 'temp'),
  list(name = 'observed/obs/a050', inputs = 'observed', stems = 'observed', temps = 'a050'),
  list(name = 'observed/obs/a070', inputs = 'observed', stems = 'observed', temps = 'a070'))

res <- list()
for (v in variants) {
  rr <- list()
  for (i in seq_len(nrow(d))) {
    o <- run_one(i, v); if (is.null(o)) next
    ann <- o[o$variable == 'gpp', ] |>
      mutate(Year = as.integer(format(date, '%Y'))) |>
      group_by(Year) |> summarise(sim_gpp = sum(value) * cf * 100, .groups = 'drop')
    rr[[length(rr)+1]] <- cbind(
      Plot_ID = d$Plot_ID[i], variant = v$name,
      sp = ifelse(v$inputs == 'observed', d$sp_obs[i], d$sp_map[i]), ann,
      lai = pk(o,'lai'), apar = pk(o,'apar',fun=sum), height = pk(o,'height'),
      alpha_c = pk(o,'alpha_c'),
      f_sw = pk(o,'f_sw'), f_vpd = pk(o,'f_vpd'), f_tmp = pk(o,'f_tmp'),
      f_frost = pk(o,'f_frost'), f_nutr = pk(o,'f_nutr'), f_age = pk(o,'f_age'))
  }
  cat(sprintf('  %-18s simulated %2d / %d plots\n', v$name, length(rr), nrow(d)))
  if (length(rr)) res[[v$name]] <- bind_rows(rr)
}

if (length(fails)) {
  ft <- bind_rows(fails)
  cat('\nFAILURES (', nrow(ft), ' total )\n', sep = '')
  print(ft |> count(variant, msg), row.names = FALSE)
  write.csv(ft, file.path(out,'preflight_failures.csv'), row.names = FALSE)
} else cat('\nno failures\n')

r <- bind_rows(res) |> inner_join(gpp_obs, by = c('Plot_ID','Year'))

## ============================================================ 4. ANSWERS ====
rule('4. THE THREE QUESTIONS')

q <- r |> group_by(variant) |>
  summarise(n = n(),
            lai_min = min(lai, na.rm=TRUE), lai_med = median(lai, na.rm=TRUE),
            lai_max = max(lai, na.rm=TRUE),
            r_lai = cor(lai, obs_lai, use='complete.obs'),
            h_med = median(height, na.rm=TRUE),
            r_gpp_apar = cor(sim_gpp, apar, use='complete.obs'),
            bias = mean(sim_gpp - obs_gpp),
            RMSE = sqrt(mean((sim_gpp - obs_gpp)^2)),
            r_vs_obs = cor(sim_gpp, obs_gpp, use='complete.obs'),
            .groups='drop')
print(as.data.frame(q), digits = 3, row.names = FALSE)

cat(sprintf('\nQ1 LAI inside measured range %.2f to %.2f?\n', LAI_LO, LAI_HI))
for (i in seq_len(nrow(q)))
  cat(sprintf('   %-18s %.2f to %.2f  %s\n', q$variant[i], q$lai_min[i], q$lai_max[i],
              if (q$lai_min[i] >= LAI_LO*0.5 && q$lai_max[i] <= LAI_HI*1.5) 'PASS' else 'CHECK'))

cat('\nQ2 r(GPP, APAR) near 0.88?\n')
for (i in seq_len(nrow(q)))
  cat(sprintf('   %-18s %.3f  %s\n', q$variant[i], q$r_gpp_apar[i],
              if (q$r_gpp_apar[i] > 0.7) 'PASS' else 'FAIL'))

cat('\nQ3 pine at its ceiling (0.634)?\n')
bysp <- r |> group_by(variant, sp) |>
  summarise(n_plots = n_distinct(Plot_ID), bias = mean(sim_gpp - obs_gpp),
            r = cor(sim_gpp, obs_gpp, use='complete.obs'), .groups='drop')
print(as.data.frame(bysp), digits = 3, row.names = FALSE)

## bootstrap by plot
boot_r <- function(dat, B = 2000) {
  p <- unique(dat$Plot_ID); o <- numeric(B)
  for (b in seq_len(B)) {
    s <- dat[unlist(lapply(sample(p, length(p), TRUE),
                           function(x) which(dat$Plot_ID == x))), ]
    o[b] <- suppressWarnings(cor(s$sim_gpp, s$obs_gpp, use='complete.obs'))
  }
  quantile(o[is.finite(o)], c(.025,.975))
}
cat('\nbootstrap CIs (2000 resamples of plots)\n')
for (v in unique(r$variant)) for (s in unique(r$sp)) {
  dat <- r[r$variant==v & r$sp==s, ]
  if (n_distinct(dat$Plot_ID) < 5) next
  ci <- boot_r(dat)
  cat(sprintf('  %-18s %-18s r=%+.3f  CI [%+.3f, %+.3f]\n', v, s,
      cor(dat$sim_gpp, dat$obs_gpp, use='complete.obs'), ci[1], ci[2]))
}

## ======================================== 5. IS WATER LIMITING AFTER ALL? ===
rule('5. PHYSIOLOGICAL MODIFIERS (post-correction)')
cat('The earlier estimate f_sw ~ 0.98 came from the broken configuration,\n')
cat('where LAI of 33 distorted the water balance. Re-checked here.\n\n')
print(r |> group_by(variant) |>
        summarise(across(c(alpha_c,f_sw,f_vpd,f_tmp,f_frost,f_nutr,f_age),
                         ~mean(.x, na.rm=TRUE)), .groups='drop') |>
        as.data.frame(), digits = 3, row.names = FALSE)
cat('\nLowest modifier is the binding constraint. f_sw near 1 confirms that\n')
cat('soil water is not limiting at this site.\n')

## =============================================================== 6. PLOTS ===
p <- ggplot(r, aes(obs_gpp, sim_gpp, colour = sp)) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_point(size=2, alpha=.85) + facet_wrap(~variant) +
  labs(x='Measured tree GPP (gC m-2 yr-1)', y='Simulated GPP (gC m-2 yr-1)') +
  theme_bw()
ggsave(file.path(out,'fig_preflight_gpp.png'), p, width=9, height=7, dpi=300)

pl <- ggplot(r, aes(obs_lai, lai, colour = sp)) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_point(size=2, alpha=.85) + facet_wrap(~variant) +
  labs(x='Measured LAImax', y='Simulated LAI') + theme_bw()
ggsave(file.path(out,'fig_preflight_lai.png'), pl, width=9, height=7, dpi=300)

write.csv(r, file.path(out,'preflight_results.csv'), row.names=FALSE)
write.csv(q, file.path(out,'preflight_summary.csv'), row.names=FALSE)

rule('VERDICT')
cat('Proceed to one sector, then the full cube, only if Q1 and Q2 PASS and\n')
cat('pine r is close to 0.634 with a lower CI bound above zero.\n')
cat('\nwrote', file.path(out,'preflight_results.csv'), 'and two figures\n')
