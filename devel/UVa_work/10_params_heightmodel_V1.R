#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 10_params_heightmodel_V1.R
##
## TWO FAULTS FOUND IN THE PARAMETER SETUP
##
## (1) HEIGHT EQUATION MISMATCH
##     'Picea abies 37' has nHB = 17.85, which only makes sense in the
##     EXPONENTIAL height form (height_model = 2):
##        1.3 + aH*exp(-nHB/dbh) = 1.3 + 37.73*exp(-17.85/14.7) = 12.5 m
##     Run under the POWER form (height_model = 1) it gives
##        aH * dbh^nHB = 37.73 * 14.7^17.85 = 2.6e22 m
##     which matches the observed median of 2.89e22. Meanwhile
##     'Pinus sylvestris 30' (aH 4.589, nHB 0.474) is a POWER-form set.
##     The two sets assume different equations and were being mixed.
##     FIX: height_model = 2 with parsets 37 + 39 (both Forrester et al. 2021).
##
## (2) TEMPERATURE RESPONSE FROM THE WRONG BIOME
##     'Picea abies 37' (Switzerland) has Topt = 24.95 C and Tmin = 2.47 C.
##     Krycklan's warmest month is ~13-15 C, so f_tmp ~ 0.33 all year, and
##     34 of 50 plots are spruce. There is NO boreal Picea set in r3PG.
##
## VARIANTS
##   A_control  hm=1, Pinus 30 + Picea 37   -> reproduces the overflow
##   B_hm2      hm=2, Pinus 39 + Picea 37   -> consistent published parameters
##   C_boreal   B + boreal Tmin/Topt/Tmax/alphaCx
##   D_finnish  C but pine from Landsberg et al. 2005 Finland (parset 50)
##
##   Rscript 10_params_heightmodel_V1.R
## ---------------------------------------------------------------------------

## SPECIES COLUMN CORRECTED: this script reads only '50 main plots',
## whose header transposes SPM and DTS. The identifiers below have been
## exchanged so that species is read from the column holding species.
suppressPackageStartupMessages({
  library(terra); library(readxl); library(dplyr); library(tidyr)
  library(r3PG); library(ggplot2)
})

bp  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
rd  <- file.path(bp, '3PG_10m/x731866.017196386_y7132899.71658545.RData')
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
out <- 'validation_out'; dir.create(out, showWarnings = FALSE)
cf  <- 0.47
fol_marklund <- c('Pinus sylvestris' = 0.045, 'Picea abies' = 0.085)
root_frac    <- c('Pinus sylvestris' = 0.25/0.76, 'Picea abies' = 0.25/0.75)

sink(file.path(out, 'params_log.txt'), split = TRUE)
cat('run:', format(Sys.time()), '\n')

## ---- climate ---------------------------------------------------------------
nm <- load(rd); xx <- get(nm[1])
cl <- xx[[1]] |> filter(group == 'climate') |>
  select(date, variable, value) |> distinct(date, variable, .keep_all = TRUE) |>
  pivot_wider(names_from = variable, values_from = value) |> arrange(date)
rm(xx); invisible(gc())
d_climate <- cl |> transmute(year = as.integer(format(date, '%Y')),
                             month = as.integer(format(date, '%m')),
                             tmp_min, tmp_max, tmp_ave, prcp,
                             srad = solar_rad, frost_days) |> as.data.frame()
cat('\nsite temperature context (monthly tmp_ave means):\n')
print(round(tapply(d_climate$tmp_ave, d_climate$month, mean), 1))
cat('warmest month mean:', round(max(tapply(d_climate$tmp_ave, d_climate$month, mean)), 1), 'C\n')

## ---- plots -----------------------------------------------------------------
obs_all <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2,
                                       .name_repair = 'unique')) |> as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year %in% 2016:2018)
gpp_obs <- obs_all |> transmute(Plot_ID, Year, obs_gpp = GPP - GPPff, obs_lai = LAImax)
d <- obs_all |> filter(Year == 2016) |>
  transmute(Plot_ID, Age, Trees, AGCtree, LAImax, SPM,
            sp = ifelse(SPM == 1, 'Pinus sylvestris', 'Picea abies'))
cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))
d <- inner_join(d, cc, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
d$latitude <- crds(project(pv, 'EPSG:4326'))[, 2]
dm <- rast(file.path(bp, 'Sweden_dem/dem.tif'))
d$altitude <- terra::extract(dm[[1]], project(pv, crs(dm)), ID = FALSE)[, 1]
cat('\nplots:', nrow(d), '| pine:', sum(d$SPM == 1), 'spruce:', sum(d$SPM == 2), '\n')

## ---- variants --------------------------------------------------------------
bor_spruce <- c(Tmin = -3, Topt = 14, Tmax = 30, alphaCx = 0.050)
bor_pine   <- c(Tmin = -5, Topt = 13, Tmax = 30, alphaCx = 0.050)

variants <- list(
  A_control = list(hm = 1, pine = 'Pinus sylvestris 30', spruce = 'Picea abies 37',
                   ov = NULL),
  B_hm2     = list(hm = 2, pine = 'Pinus sylvestris 39', spruce = 'Picea abies 37',
                   ov = NULL),
  C_boreal  = list(hm = 2, pine = 'Pinus sylvestris 39', spruce = 'Picea abies 37',
                   ov = list('Pinus sylvestris' = bor_pine,
                             'Picea abies'      = bor_spruce)),
  D_finnish = list(hm = 2, pine = 'Pinus sylvestris 50', spruce = 'Picea abies 37',
                   ov = list('Picea abies' = bor_spruce)))

build_par <- function(sp, pset, ov) {
  p <- get_parameters('parameters', sp_names = pset)
  names(p)[2] <- sp
  if (!is.null(ov) && sp %in% names(ov)) {
    o <- ov[[sp]]
    for (k in names(o)) if (k %in% p$parameter) p[p$parameter == k, sp] <- o[[k]]
  }
  p
}

run_one <- function(i, v) {
  sp  <- d$sp[i]; agb <- d$AGCtree[i] / cf; age <- d$Age[i]
  if (!is.finite(agb) || agb <= 0 || !is.finite(age) || age <= 0) return(NULL)
  fol <- agb * fol_marklund[[sp]]
  sprow <- data.frame(species = sp,
                      planted = sprintf('%04d-01', 2016 - round(age)),
                      fertility = 0.5, stems_n = d$Trees[i],
                      biom_stem = agb * (1 - fol_marklund[[sp]]),
                      biom_root = agb * root_frac[[sp]],
                      biom_foliage = fol, stringsAsFactors = FALSE)
  site <- data.frame(latitude = d$latitude[i], altitude = d$altitude[i],
                     soil_class = 2, asw_i = 150, asw_min = 0, asw_max = 200,
                     from = '2016-01', to = '2018-12')
  pset <- if (sp == 'Pinus sylvestris') v$pine else v$spruce
  par  <- build_par(sp, pset, v$ov)
  sz   <- get_parameters('sizeDist', sp_names = pset); names(sz)[2] <- sp
  tryCatch(run_3PG(site = site, species = sprow, climate = d_climate,
                   thinning = NULL,
                   parameters = prepare_parameters(par, sp_names = sp),
                   size_dist  = prepare_sizeDist(sz,  sp_names = sp),
                   settings = list(light_model = 2, transp_model = 2,
                                   phys_model = 2, height_model = v$hm,
                                   correct_bias = 0, calculate_d13c = 0),
                   check_input = TRUE, df_out = TRUE),
           error = function(e) NULL)
}

pk <- function(o, v, yr = '2017', fun = mean) {
  z <- o[o$variable == v & format(o$date, '%Y') == yr, ]
  if (!nrow(z)) return(NA_real_)
  fun(z$value, na.rm = TRUE)
}

res <- list()
for (vn in names(variants)) {
  v <- variants[[vn]]
  cat('\n== variant', vn, '| height_model', v$hm, '|', v$pine, '+', v$spruce,
      '|', if (is.null(v$ov)) 'published' else 'boreal override', '==\n')
  rr <- list(); nf <- 0
  for (i in seq_len(nrow(d))) {
    o <- run_one(i, v); if (is.null(o)) { nf <- nf + 1; next }
    ann <- o[o$variable == 'gpp', ] |>
      mutate(Year = as.integer(format(date, '%Y'))) |>
      group_by(Year) |> summarise(sim_gpp = sum(value) * cf * 100, .groups = 'drop')
    rr[[length(rr) + 1]] <- cbind(
      Plot_ID = d$Plot_ID[i], variant = vn, sp = d$sp[i], ann,
      lai = pk(o, 'lai'), apar = pk(o, 'apar', fun = sum),
      height = pk(o, 'height'), dbh = pk(o, 'dbh'),
      alpha_c = pk(o, 'alpha_c'),
      f_tmp = pk(o, 'f_tmp'), f_frost = pk(o, 'f_frost'),
      f_age = pk(o, 'f_age'), f_nutr = pk(o, 'f_nutr'),
      f_vpd = pk(o, 'f_vpd'), f_sw = pk(o, 'f_sw'))
  }
  cat('  simulated:', length(rr), '| failed:', nf, '\n')
  if (length(rr)) res[[vn]] <- bind_rows(rr)
}

r <- bind_rows(res) |> inner_join(gpp_obs, by = c('Plot_ID','Year'))

## ---- verdict ---------------------------------------------------------------
cat('\n================= HEIGHT SANITY =================\n')
print(r |> group_by(variant) |>
        summarise(n = n_distinct(Plot_ID),
                  h_med = median(height, na.rm = TRUE),
                  h_max = max(height, na.rm = TRUE),
                  dbh_med = median(dbh, na.rm = TRUE), .groups = 'drop') |>
        as.data.frame(), digits = 4, row.names = FALSE)
cat('  realistic Krycklan height is roughly 5-28 m\n')

cat('\n================= MODIFIERS (2017 means) =================\n')
print(r |> group_by(variant) |>
        summarise(across(c(alpha_c, f_tmp, f_frost, f_age, f_nutr, f_vpd, f_sw),
                         ~ mean(.x, na.rm = TRUE)), .groups = 'drop') |>
        as.data.frame(), digits = 3, row.names = FALSE)

cat('\n================= SKILL =================\n')
tab <- r |> group_by(variant) |>
  summarise(n = n(),
            lai_med = median(lai, na.rm = TRUE),
            gpp_mean = mean(sim_gpp), obs_mean = mean(obs_gpp),
            bias = mean(sim_gpp - obs_gpp),
            RMSE = sqrt(mean((sim_gpp - obs_gpp)^2)),
            r_gpp_apar = cor(sim_gpp, apar, use = 'complete.obs'),
            r_gpp_lai  = cor(sim_gpp, lai,  use = 'complete.obs'),
            r_vs_obs   = cor(sim_gpp, obs_gpp, use = 'complete.obs'),
            R2_vs_obs  = cor(sim_gpp, obs_gpp, use = 'complete.obs')^2,
            .groups = 'drop')
print(as.data.frame(tab), digits = 3, row.names = FALSE)

cat('\nmeasured LAImax median:', round(median(gpp_obs$obs_lai, na.rm = TRUE), 2), '\n')
cat('measured tree GPP mean:', round(mean(gpp_obs$obs_gpp, na.rm = TRUE)), 'gC/m2/yr\n')
cat('\nCEILING: observed GPP correlates with observed biomass at r = 0.474,\n')
cat('so R2 ~ 0.22 is the realistic maximum against these plots.\n')

cat('\nWHAT TO LOOK FOR\n')
cat('  A_control  height ~1e22 and r_gpp_apar ~ 0  (confirms the mismatch)\n')
cat('  B_hm2      height 5-30 m. If r_gpp_apar jumps, the height equation was\n')
cat('             the whole problem.\n')
cat('  C_boreal   f_tmp should rise well above 0.33 and bias shrink.\n')
cat('  D_finnish  tests whether a genuinely boreal pine set beats the override.\n')

## by species: spruce is 34 of 50 plots and carries the Swiss parameters
cat('\n================= BY SPECIES =================\n')
print(r |> group_by(variant, sp) |>
        summarise(n = n(), bias = mean(sim_gpp - obs_gpp),
                  r = cor(sim_gpp, obs_gpp, use = 'complete.obs'),
                  f_tmp = mean(f_tmp, na.rm = TRUE), .groups = 'drop') |>
        as.data.frame(), digits = 3, row.names = FALSE)

p <- ggplot(r, aes(obs_gpp, sim_gpp, colour = sp)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2, alpha = .85) + facet_wrap(~variant) +
  labs(x = 'Observed tree GPP (gC m-2 yr-1)',
       y = 'Simulated GPP (gC m-2 yr-1)') + theme_bw()
ggsave(file.path(out, 'fig_params_variants.png'), p, width = 9, height = 7, dpi = 300)

write.csv(r, file.path(out, 'params_results.csv'), row.names = FALSE)
write.csv(tab, file.path(out, 'params_summary.csv'), row.names = FALSE)
cat('\nwrote', file.path(out, 'params_results.csv'), '\n')
sink()
