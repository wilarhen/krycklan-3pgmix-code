#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 08_pointmode_57plots_V1.R
##
## A0: run 3-PG at the Krycklan flux plots instead of 1.2M pixels.
##
## DIAGNOSIS BEING TESTED
##   polygon_grids_v01.r lines 166-184 map the biomass CSV columns
##     Stem_prop / Branch_prop / Stump_and_root_prop
##   onto
##     biom_stem / biom_foliage / biom_root
##   so BRANCH WOOD is fed in as FOLIAGE. Branch/Stem = 0.25/0.50 = 0.50 for
##   Picea, and the saved cube shows biom_foliage/biom_stem = 0.49. Confirmed.
##   Result: LAI 2.5-33.6 (measured LAImax is 0.5-5.6), canopy far above light
##   saturation, APAR maxed everywhere, GPP pinned near 1.0 tDM/ha/month
##   regardless of biomass -> R2 = 0.015 against the plots.
##
## FACTORIAL
##   foliage : 'branch' (reproduces the bug)  vs 'marklund' (the fix)
##   inputs  : 'mapped' (BIO/age rasters/NMD) vs 'observed' (plot measurements)
##
##   'branch' + 'mapped' should reproduce the cube. That is the control: if it
##   does not reproduce LAI ~ 15-30 and a flat GPP, the diagnosis is wrong.
##
## Simulation window is 2016-01 to 2018-12, initialised from the contemporaneous
## stand state. No 20-year spin-up, so the initial state cannot be washed out
## and we test only the state -> production conversion.
##
##   Rscript 08_pointmode_57plots_V1.R
## ---------------------------------------------------------------------------

## SPECIES COLUMN CORRECTED: this script reads only '50 main plots',
## whose header transposes SPM and DTS. The identifiers below have been
## exchanged so that species is read from the column holding species.
suppressPackageStartupMessages({
  library(terra); library(sf); library(readxl)
  library(dplyr); library(tidyr); library(r3PG); library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG ----
bp   <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
rd   <- file.path(bp, '3PG_10m/x731866.017196386_y7132899.71658545.RData')
xl   <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
out  <- 'validation_out'; dir.create(out, showWarnings = FALSE)

cf        <- 0.47      # carbon fraction, tC -> tDM
from_date <- '2016-01'
to_date   <- '2018-12'

## Foliage as a fraction of ABOVEGROUND biomass (Marklund 1988, approx).
fol_marklund <- c('Pinus sylvestris' = 0.045, 'Picea abies' = 0.085)
## Proportions from biom_stocks_from_Hans_Peterson.csv (of TOTAL tree biomass)
prop_tab <- data.frame(
  sp     = c('Pinus sylvestris', 'Picea abies'),
  stem   = c(0.60, 0.50),
  branch = c(0.16, 0.25),
  root   = c(0.25, 0.25))

sink(file.path(out, 'pointmode_log.txt'), split = TRUE)
cat('run:', format(Sys.time()), '\n')

## ------------------------------------------------------------- 1. CLIMATE ---
cat('\n== climate ==\n')
nm <- load(rd); xx <- get(nm[1])
cl <- xx[[1]] |> filter(group == 'climate') |>
  select(date, variable, value) |>
  distinct(date, variable, .keep_all = TRUE) |>
  pivot_wider(names_from = variable, values_from = value) |>
  arrange(date)
rm(xx); invisible(gc())

d_climate <- cl |>
  transmute(year  = as.integer(format(date, '%Y')),
            month = as.integer(format(date, '%m')),
            tmp_min, tmp_max, tmp_ave, prcp,
            srad = solar_rad, frost_days) |>
  as.data.frame()
cat('months:', nrow(d_climate), '|',
    min(d_climate$year), '-', max(d_climate$year), '\n')
stopifnot(all(2016:2018 %in% d_climate$year))

## ------------------------------------------------------------ 2. OBSERVED ---
cat('\n== plots ==\n')
obs <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2,
                                   .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year    = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year %in% 2016:2018)

## static attributes from 2016; GPP kept per year
static <- obs |> filter(Year == 2016) |>
  transmute(Plot_ID, Age, Trees, AGCtree, LAImax, SPM,
            sp_obs = ifelse(SPM == 1, 'Pinus sylvestris', 'Picea abies'))
gpp_obs <- obs |> transmute(Plot_ID, Year, obs_gpp = GPP - GPPff,
                            obs_lai = LAImax)

coords <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))

d <- inner_join(static, coords, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
cat('plots:', nrow(d), '| pine:', sum(d$SPM == 1), 'spruce:', sum(d$SPM == 2), '\n')

## ------------------------------------------------------ 3. MAPPED INPUTS ----
cat('\n== mapped inputs at plots ==\n')
gv <- function(f) {
  r <- rast(f)
  terra::extract(r[[1]], project(pv, crs(r)), ID = FALSE)[, 1]
}
d$in_bio    <- gv(file.path(bp, 'biomass/BIO.tif'))
d$in_pine   <- gv(file.path(bp, 'Sweden_Forest_Age/PINE_AGE_2025.tif'))
d$in_spruce <- gv(file.path(bp, 'Sweden_Forest_Age/SPRUCE_AGE_2025.tif'))
d$in_nmd    <- gv(file.path(bp, 'nationalndcoverdata/nmd2018bas_ogeneraliserad_v1_0.tif'))
d$altitude  <- gv(file.path(bp, 'Sweden_dem/dem.tif'))

## NMD 2018: 111 pine, 112 spruce, 113 mixed conifer, 114+ broadleaf/mixed
d$sp_map <- ifelse(d$in_nmd == 111, 'Pinus sylvestris',
            ifelse(d$in_nmd == 112, 'Picea abies',
                   ifelse(d$SPM == 1, 'Pinus sylvestris', 'Picea abies')))
d$age_map <- ifelse(d$sp_map == 'Pinus sylvestris', d$in_pine, d$in_spruce)

d$latitude <- crds(project(pv, 'EPSG:4326'))[, 2]
cat('latitude range:', round(range(d$latitude), 3), '\n')
cat('altitude range:', round(range(d$altitude, na.rm = TRUE), 1), '\n')
cat('species agreement map vs observed:',
    round(100 * mean(d$sp_map == d$sp_obs), 1), '%\n')

## soil: ESDB takes ~2 values here and f_sw was ~0.98, so soil is not limiting.
d$soil_class <- 2; d$asw_i <- 150; d$asw_min <- 0; d$asw_max <- 200

## ------------------------------------------------------------ 4. VARIANTS ---
## Build the d_species row for one plot under one variant.
make_species <- function(i, foliage, inputs) {
  sp  <- if (inputs == 'observed') d$sp_obs[i] else d$sp_map[i]
  age <- if (inputs == 'observed') d$Age[i]    else d$age_map[i]
  agb <- if (inputs == 'observed') d$AGCtree[i] / cf else d$in_bio[i]
  stn <- if (inputs == 'observed') d$Trees[i]  else 2000

  if (!is.finite(agb) || agb <= 0 || !is.finite(age) || age <= 0) return(NULL)
  p <- prop_tab[prop_tab$sp == sp, ]

  if (foliage == 'branch') {
    ## THE BUG: proportions of total biomass applied to BIO, branch -> foliage
    b_stem <- agb * p$stem
    b_fol  <- agb * p$branch
    b_root <- agb * p$root
  } else {
    ## FIX: 3-PG stem pool includes branches; foliage from Marklund
    b_fol  <- agb * fol_marklund[[sp]]
    b_stem <- agb * (p$stem + p$branch) / (p$stem + p$branch) * (1 - fol_marklund[[sp]])
    b_root <- agb * p$root / (p$stem + p$branch)
  }
  planted <- sprintf('%04d-01', 2016 - round(age))
  data.frame(species = sp, planted = planted, fertility = 0.5,
             stems_n = stn, biom_stem = b_stem,
             biom_root = b_root, biom_foliage = b_fol,
             stringsAsFactors = FALSE)
}

run_one <- function(i, foliage, inputs) {
  sprow <- make_species(i, foliage, inputs)
  if (is.null(sprow)) return(NULL)
  site <- data.frame(latitude = d$latitude[i], altitude = d$altitude[i],
                     soil_class = d$soil_class[i], asw_i = d$asw_i[i],
                     asw_min = d$asw_min[i], asw_max = d$asw_max[i],
                     from = from_date, to = to_date)
  sn <- sprow$species
  pn <- if (sn == 'Pinus sylvestris') 'Pinus sylvestris 30' else 'Picea abies'
  par <- get_parameters('parameters', sp_names = pn)
  names(par)[2] <- sn
  sz  <- get_parameters('sizeDist',   sp_names = pn)
  names(sz)[2]  <- sn

  tryCatch(
    run_3PG(site = site, species = sprow, climate = d_climate,
            thinning = NULL,
            parameters = prepare_parameters(par, sp_names = sn),
            size_dist  = prepare_sizeDist(sz,  sp_names = sn),
            settings = list(light_model = 2, transp_model = 2, phys_model = 2,
                            height_model = 1, correct_bias = 0,
                            calculate_d13c = 0),
            check_input = TRUE, df_out = TRUE),
    error = function(e) { attr(e, 'plot') <- d$Plot_ID[i]; e })
}

summarise_run <- function(o) {
  if (is.null(o) || inherits(o, 'error')) return(NULL)
  g <- o[o$variable == 'gpp', ]
  l <- o[o$variable == 'lai', ]
  data.frame(Year = as.integer(format(g$date, '%Y')), gpp = g$value) |>
    group_by(Year) |> summarise(sim_gpp = sum(gpp) * cf * 100, .groups = 'drop') |>
    mutate(sim_lai = mean(l$value, na.rm = TRUE),
           sim_lai_max = max(l$value, na.rm = TRUE))
}

## ------------------------------------------------------------- 5. EXECUTE ---
grid <- expand.grid(foliage = c('branch','marklund'),
                    inputs  = c('mapped','observed'),
                    stringsAsFactors = FALSE)

all_res <- list(); fails <- 0
for (k in seq_len(nrow(grid))) {
  fo <- grid$foliage[k]; ip <- grid$inputs[k]
  cat('\n== variant:', fo, '/', ip, '==\n')
  rows <- list()
  for (i in seq_len(nrow(d))) {
    o <- run_one(i, fo, ip)
    if (inherits(o, 'error')) { fails <- fails + 1; next }
    s <- summarise_run(o); if (is.null(s)) next
    rows[[length(rows) + 1]] <- cbind(Plot_ID = d$Plot_ID[i], s)
  }
  if (!length(rows)) { cat('  all runs failed\n'); next }
  r <- bind_rows(rows) |> mutate(foliage = fo, inputs = ip)
  all_res[[k]] <- r
  cat('  plots simulated:', n_distinct(r$Plot_ID), '\n')
  cat('  simulated LAI  : ', paste(round(range(r$sim_lai, na.rm = TRUE), 2),
                                   collapse = ' - '), '\n')
}
cat('\nfailed runs:', fails, '\n')

res <- bind_rows(all_res) |>
  inner_join(gpp_obs, by = c('Plot_ID','Year')) |>
  inner_join(select(d, Plot_ID, in_bio, AGCtree), by = 'Plot_ID')

## -------------------------------------------------------------- 6. VERDICT --
cat('\n=================== RESULTS ===================\n')
cat('measured LAImax range : ',
    paste(round(range(gpp_obs$obs_lai, na.rm = TRUE), 2), collapse = ' - '), '\n')
cat('measured tree GPP     : ',
    paste(round(range(gpp_obs$obs_gpp, na.rm = TRUE)), collapse = ' - '),
    ' gC/m2/yr\n\n')

tab <- res |> group_by(foliage, inputs) |>
  summarise(n = n(),
            lai_lo = min(sim_lai, na.rm = TRUE),
            lai_hi = max(sim_lai, na.rm = TRUE),
            gpp_mean = mean(sim_gpp), gpp_lo = min(sim_gpp), gpp_hi = max(sim_gpp),
            bias = mean(sim_gpp - obs_gpp),
            RMSE = sqrt(mean((sim_gpp - obs_gpp)^2)),
            R2_vs_obs  = cor(sim_gpp, obs_gpp)^2,
            r_vs_obs   = cor(sim_gpp, obs_gpp),
            r_vs_biom  = cor(sim_gpp, AGCtree),
            .groups = 'drop')
print(as.data.frame(tab), digits = 3, row.names = FALSE)

cat('\nINTERPRETATION\n')
cat('  branch/mapped   should show LAI >> 5.6 and r_vs_obs ~ 0 (reproduces cube)\n')
cat('  marklund/*      should show LAI within 0.5-5.6 and r_vs_biom clearly > 0\n')
cat('  If marklund fixes LAI but r_vs_obs stays ~0, something else remains.\n')

p <- ggplot(res, aes(obs_gpp, sim_gpp)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(colour = factor(Year)), size = 2, alpha = .85) +
  facet_grid(foliage ~ inputs) +
  labs(x = 'Observed tree GPP (gC m-2 yr-1)',
       y = 'Simulated GPP (gC m-2 yr-1)', colour = 'Year') + theme_bw()
ggsave(file.path(out, 'fig_pointmode_factorial.png'), p,
       width = 8, height = 7, dpi = 300)

pl <- ggplot(res, aes(obs_lai, sim_lai)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2, alpha = .85) +
  facet_grid(foliage ~ inputs, scales = 'free_y') +
  labs(x = 'Measured LAImax', y = 'Simulated LAI') + theme_bw()
ggsave(file.path(out, 'fig_pointmode_lai.png'), pl, width = 8, height = 7, dpi = 300)

write.csv(res, file.path(out, 'pointmode_results.csv'), row.names = FALSE)
write.csv(tab, file.path(out, 'pointmode_summary.csv'), row.names = FALSE)
cat('\nwrote', file.path(out, 'pointmode_results.csv'), 'and two figures\n')
sink()
