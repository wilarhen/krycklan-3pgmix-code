#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 09_light_model_V1.R
##
## After the foliage fix, simulated LAI matches measurements (r = 0.875,
## median 2.92 vs 2.74) but simulated GPP is INDEPENDENT of simulated LAI
## (r ~ 0). Under Beer's law that is impossible, so absorbed radiation is not
## responding to leaf area.
##
## Prime suspect: light_model = 2 (3-PGmix layered canopy), which derives APAR
## from height / crown_length / crown_width / canopy_vol_frac rather than from
## LAI directly. The saved cube showed canopy_cover = 0.0 everywhere.
##
## This runs the plots under light_model 1 (Beer's law) and 2 (3-PGmix) and
## reports whether APAR tracks LAI in each.
##
##   Rscript 09_light_model_V1.R
## ---------------------------------------------------------------------------

## SPECIES COLUMN CORRECTED: this script reads only '50 main plots',
## whose header transposes SPM and DTS. The identifiers below have been
## exchanged so that species is read from the column holding species.
suppressPackageStartupMessages({
  library(terra); library(readxl); library(dplyr); library(tidyr); library(r3PG)
})

bp  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
rd  <- file.path(bp, '3PG_10m/x731866.017196386_y7132899.71658545.RData')
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
out <- 'validation_out'; dir.create(out, showWarnings = FALSE)
cf  <- 0.47
fol_marklund <- c('Pinus sylvestris' = 0.045, 'Picea abies' = 0.085)
root_frac    <- c('Pinus sylvestris' = 0.25/0.76, 'Picea abies' = 0.25/0.75)

sink(file.path(out, 'light_model_log.txt'), split = TRUE)

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

## ---- plots -----------------------------------------------------------------
obs <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2,
                                   .name_repair = 'unique')) |> as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year == 2016) |>
  transmute(Plot_ID, Age, Trees, AGCtree, LAImax, SPM,
            sp = ifelse(SPM == 1, 'Pinus sylvestris', 'Picea abies'))
cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
  filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))
d <- inner_join(obs, cc, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
d$latitude <- crds(project(pv, 'EPSG:4326'))[, 2]
dm <- rast(file.path(bp, 'Sweden_dem/dem.tif'))
d$altitude <- terra::extract(dm[[1]], project(pv, crs(dm)), ID = FALSE)[, 1]
cat('plots:', nrow(d), '\n')

## ---- run one plot under a given light model --------------------------------
run_lm <- function(i, lm) {
  sp  <- d$sp[i]; agb <- d$AGCtree[i] / cf; age <- d$Age[i]
  if (!is.finite(agb) || agb <= 0) return(NULL)
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
  pn  <- if (sp == 'Pinus sylvestris') 'Pinus sylvestris 30' else 'Picea abies'
  par <- get_parameters('parameters', sp_names = pn); names(par)[2] <- sp
  sz  <- get_parameters('sizeDist',   sp_names = pn); names(sz)[2]  <- sp
  tryCatch(run_3PG(site = site, species = sprow, climate = d_climate,
                   thinning = NULL,
                   parameters = prepare_parameters(par, sp_names = sp),
                   size_dist  = prepare_sizeDist(sz,  sp_names = sp),
                   settings = list(light_model = lm, transp_model = 2,
                                   phys_model = 2, height_model = 1,
                                   correct_bias = 0, calculate_d13c = 0),
                   check_input = TRUE, df_out = TRUE),
           error = function(e) NULL)
}

pick <- function(o, v, fun = mean) {
  z <- o[o$variable == v & format(o$date, '%Y') == '2017', ]
  if (!nrow(z)) return(NA_real_)
  fun(z$value, na.rm = TRUE)
}

res <- list()
for (lm in c(1, 2)) {
  cat('\n== light_model =', lm, '==\n')
  rr <- list()
  for (i in seq_len(nrow(d))) {
    o <- run_lm(i, lm); if (is.null(o)) next
    rr[[length(rr) + 1]] <- data.frame(
      Plot_ID = d$Plot_ID[i], light_model = lm,
      lai          = pick(o, 'lai'),
      apar         = pick(o, 'apar', sum),
      canopy_cover = pick(o, 'canopy_cover'),
      canopy_vol_frac = pick(o, 'canopy_vol_frac'),
      height       = pick(o, 'height'),
      crown_width  = pick(o, 'crown_width'),
      crown_length = pick(o, 'crown_length'),
      alpha_c      = pick(o, 'alpha_c'),
      gpp          = pick(o, 'gpp', sum) * cf * 100,
      obs_gpp      = NA_real_, obs_lai = d$LAImax[i])
    }
  if (!length(rr)) { cat('  all failed\n'); next }
  res[[length(res) + 1]] <- bind_rows(rr)
  cat('  n =', nrow(res[[length(res)]]), '\n')
}
r <- bind_rows(res)

cat('\n============ DOES APAR RESPOND TO LAI? ============\n')
for (lm in unique(r$light_model)) {
  s <- r[r$light_model == lm, ]
  cat(sprintf('\nlight_model = %d   (n = %d)\n', lm, nrow(s)))
  cat(sprintf('  r(apar, lai)      = %+0.3f   <-- must be strongly positive\n',
              cor(s$apar, s$lai, use = 'complete.obs')))
  cat(sprintf('  r(gpp,  lai)      = %+0.3f\n',
              cor(s$gpp,  s$lai, use = 'complete.obs')))
  cat(sprintf('  r(gpp,  apar)     = %+0.3f\n',
              cor(s$gpp,  s$apar, use = 'complete.obs')))
  cat(sprintf('  r(lai,  obs_lai)  = %+0.3f\n',
              cor(s$lai,  s$obs_lai, use = 'complete.obs')))
  for (v in c('lai','apar','canopy_cover','canopy_vol_frac','height',
              'crown_width','crown_length','alpha_c','gpp')) {
    z <- s[[v]]
    cat(sprintf('  %-16s med %9.3f  range %9.3f - %9.3f\n', v,
                median(z, na.rm = TRUE), min(z, na.rm = TRUE), max(z, na.rm = TRUE)))
  }
}

cat('\nINTERPRETATION\n')
cat('  light_model 1 gives r(apar,lai) >> 0 while 2 does not  -> 3-PGmix layered\n')
cat('    canopy is the fault; switch to light_model = 1 or fix crown geometry.\n')
cat('  BOTH near zero -> APAR is not LAI-driven at all; inspect canopy_cover and\n')
cat('    canopy_vol_frac above (0 means the canopy is treated as absent).\n')

write.csv(r, file.path(out, 'light_model_test.csv'), row.names = FALSE)
cat('\nwrote', file.path(out, 'light_model_test.csv'), '\n')
sink()
