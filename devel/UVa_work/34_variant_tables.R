## ---------------------------------------------------------------------------
## 34_variant_tables.R
##
## Builds the alpha-variant table and the corrected Table 5. Aggregates existing model output; runs no
## model.
##
##   Part 1  results/alpha_variants.csv, from results/preflight_results.csv
##           (script 16, seven variants, run with out = 'results').
##           Separates the alphaCx change from the temperature change that the
##           'bor' override applies together.
##   Part 2  results/table5_corrected.csv, from validation_out/params_results.csv
##           (script 10, rerun 1 Oct after the species correction). Rebuilds
##           every cell of manuscript Table 5 on one unit, plot means, per
##           species. The published table mixed plot-years and plot means and
##           pooled median height and r(GPP, APAR) over species.
##
## Plot means (2016-2018), bootstrap of plots, B = 2000, percentile intervals,
## seed as in ceiling_config.R.
##
## Checks (stop on failure):
##   - observed/obs and observed/obs/bor here equal r, bias and RMSE in
##     validation_out/plotlevel/plotlevel_stats.csv;
##   - script 10's B_hm2 and C_boreal equal script 16's observed/obs and
##     observed/obs/bor.
##
##   Rscript 34_variant_tables.R      (from UVa_work/)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(r3PG))
NB  <- 2000
OUT <- 'results'
TOL <- 1e-9
sink(file.path(OUT, 'variant_tables_log.txt'), split = TRUE)

## ---- plot-mean statistics ---------------------------------------------------
pm_of <- function(d) aggregate(cbind(sim_gpp, obs_gpp, apar, height) ~ Plot_ID,
                               d, mean)
boot_ci <- function(a) {
  set.seed(20260824)        # per call: identical data give identical intervals
  o <- replicate(NB, { b <- a[sample.int(nrow(a), replace = TRUE), ]
                       suppressWarnings(cor(b$sim_gpp, b$obs_gpp)) })
  quantile(o[is.finite(o)], c(.025, .975))
}
stats_of <- function(d) {
  a  <- pm_of(d)
  ci <- if (sd(a$sim_gpp) > 0) boot_ci(a) else c(NA, NA)
  data.frame(n_plots = nrow(a),
             median_h = median(a$height),
             r_gpp_apar = suppressWarnings(cor(a$sim_gpp, a$apar)),
             bias = mean(a$sim_gpp - a$obs_gpp),
             RMSE = sqrt(mean((a$sim_gpp - a$obs_gpp)^2)),
             r = suppressWarnings(cor(a$sim_gpp, a$obs_gpp)),
             r_lo = ci[[1]], r_hi = ci[[2]])
}

## ---- parameter values per variant ------------------------------------------
base_par <- function(set) {
  p <- get_parameters('parameters', sp_names = set)
  v <- setNames(p[[2]], p$parameter)
  v[c('alphaCx', 'Tmin', 'Topt', 'Tmax')]
}
overlay <- function(set, ov) { v <- base_par(set); v[names(ov)] <- ov; v }

## ===========================================================================
## PART 1
## ===========================================================================
cat('==== PART 1: alphaCx and temperature separated (script 16) ====\n')
set.seed(20260824)
pf <- read.csv(file.path(OUT, 'preflight_results.csv'), stringsAsFactors = FALSE)

set16 <- c('Pinus sylvestris' = 'Pinus sylvestris 39', 'Picea abies' = 'Picea abies 37')
ov16 <- list(
  'observed/obs'      = list(),
  'observed/obs/bor'  = list('Pinus sylvestris' = c(Tmin = -5, Topt = 13, Tmax = 30, alphaCx = 0.050),
                             'Picea abies'      = c(Tmin = -3, Topt = 14, Tmax = 30, alphaCx = 0.050)),
  'observed/obs/temp' = list('Pinus sylvestris' = c(Tmin = -5, Topt = 13, Tmax = 30),
                             'Picea abies'      = c(Tmin = -3, Topt = 14, Tmax = 30)),
  'observed/obs/a050' = list('Pinus sylvestris' = c(alphaCx = 0.050),
                             'Picea abies'      = c(alphaCx = 0.050)),
  'observed/obs/a070' = list('Picea abies'      = c(alphaCx = 0.070)))
role <- c('observed/obs' = 'published (Forrester 2021)',
          'observed/obs/bor' = 'override: temperature + alphaCx 0.050',
          'observed/obs/temp' = 'override: temperature only',
          'observed/obs/a050' = 'excursion: alphaCx 0.050 both species',
          'observed/obs/a070' = 'published spruce alternative (Trotsiuk 2020); pine published')

rows <- list()
for (v in names(ov16)) for (s in names(set16)) {
  d <- pf[pf$variant == v & pf$sp == s, ]
  if (!nrow(d)) stop('no rows for ', v, ' / ', s)
  pv <- overlay(set16[[s]], ov16[[v]][[s]])
  rows[[length(rows) + 1]] <- cbind(
    data.frame(variant = v, role = role[[v]], species = s, parameter_set = set16[[s]],
               alphaCx = pv[['alphaCx']], Tmin = pv[['Tmin']], Topt = pv[['Topt']],
               Tmax = pv[['Tmax']]),
    stats_of(d))
}
A4 <- do.call(rbind, rows)
print(transform(A4[, c('variant', 'species', 'alphaCx', 'Tmin', 'Topt', 'Tmax',
                       'n_plots', 'bias', 'RMSE', 'r', 'r_lo', 'r_hi')],
                bias = round(bias, 1), RMSE = round(RMSE, 1), r = round(r, 3),
                r_lo = round(r_lo, 3), r_hi = round(r_hi, 3)), row.names = FALSE)

## check against script 26
pl <- read.csv('validation_out/plotlevel/plotlevel_stats.csv', stringsAsFactors = FALSE)
cat('\nCheck against validation_out/plotlevel/plotlevel_stats.csv:\n')
for (v in c('observed/obs', 'observed/obs/bor')) for (s in names(set16)) {
  a <- A4[A4$variant == v & A4$species == s, ]
  b <- pl[pl$variant == v & pl$species == s, ]
  dd <- c(r = a$r - b$r_corrected, bias = a$bias - b$bias, RMSE = a$RMSE - b$RMSE)
  ok <- all(abs(dd) < TOL)
  cat(sprintf('  %-17s %-17s r %.6f vs %.6f  bias %.3f vs %.3f  RMSE %.3f vs %.3f  %s\n',
              v, s, a$r, b$r_corrected, a$bias, b$bias, a$RMSE, b$RMSE,
              if (ok) 'PASS' else 'FAIL'))
  if (!ok) stop('Part 1 check failed for ', v, ' / ', s)
}
write.csv(A4, file.path(OUT, 'alpha_variants.csv'), row.names = FALSE)

## ===========================================================================
## PART 2
## ===========================================================================
cat('\n==== PART 2: Table 5 on plot means, per species (script 10) ====\n')
set.seed(20260824)
p10 <- read.csv('validation_out/params_results.csv', stringsAsFactors = FALSE)

bor_p <- c(Tmin = -5, Topt = 13, Tmax = 30, alphaCx = 0.050)
bor_s <- c(Tmin = -3, Topt = 14, Tmax = 30, alphaCx = 0.050)
def <- list(
  A_control = list(label = 'A. Original', hm = 1,
                   set = c('Pinus sylvestris' = 'Pinus sylvestris 30', 'Picea abies' = 'Picea abies 37'),
                   ov = list()),
  B_hm2     = list(label = 'B. Height form 2; pine set 30 -> 39', hm = 2,
                   set = c('Pinus sylvestris' = 'Pinus sylvestris 39', 'Picea abies' = 'Picea abies 37'),
                   ov = list()),
  C_boreal  = list(label = 'C. B + boreal Tmin/Topt/Tmax and alphaCx 0.050', hm = 2,
                   set = c('Pinus sylvestris' = 'Pinus sylvestris 39', 'Picea abies' = 'Picea abies 37'),
                   ov = list('Pinus sylvestris' = bor_p, 'Picea abies' = bor_s)))

## published Table 5 (main.tex lines 1189-1196), grouped by soil
pub <- data.frame(
  variant = rep(c('A_control', 'B_hm2', 'C_boreal'), each = 2),
  species = rep(c('Pinus sylvestris', 'Picea abies'), 3),
  pub_n = rep(c(14, 33), 3),
  pub_median_h = c(3.2e22, 3.2e22, 15.7, 15.7, 15.7, 15.7),
  pub_r_gpp_apar = c(0.008, 0.008, 0.927, 0.927, 0.928, 0.928),
  pub_bias = c(-118, -345, -455, -365, 5, 389),
  pub_RMSE = c(NA, NA, 488, 463, 260, 488),
  pub_r = c(0.464, 0.073, 0.702, 0.213, 0.650, 0.085),
  pub_lo = c(0.057, -0.124, 0.451, -0.101, 0.376, -0.307),
  pub_hi = c(0.725, 0.230, 0.883, 0.416, 0.873, 0.373),
  pub_source = c(rep('h, r_GPP,APAR, bias: script 10, plot-years, h and r_GPP,APAR pooled; r, CI: script 12, plot-years', 2),
                 rep('h: script 16, plot-years, pooled; r_GPP,APAR: script 26, plot means, pooled; bias, RMSE, r, CI: script 26, plot means', 4)))

rows <- list()
for (v in names(def)) for (s in c('Pinus sylvestris', 'Picea abies')) {
  d  <- p10[p10$variant == v & p10$sp == s, ]
  pv <- overlay(def[[v]]$set[[s]], def[[v]]$ov[[s]])
  rows[[length(rows) + 1]] <- cbind(
    data.frame(variant = v, label = def[[v]]$label, species = s,
               height_form = def[[v]]$hm, parameter_set = def[[v]]$set[[s]],
               alphaCx = pv[['alphaCx']], Tmin = pv[['Tmin']], Topt = pv[['Topt']],
               Tmax = pv[['Tmax']]),
    stats_of(d))
}
T5 <- merge(do.call(rbind, rows), pub, by = c('variant', 'species'), sort = FALSE)
print(transform(T5[, c('variant', 'species', 'parameter_set', 'alphaCx', 'n_plots',
                       'median_h', 'r_gpp_apar', 'bias', 'RMSE', 'r', 'r_lo', 'r_hi')],
                median_h = signif(median_h, 3), r_gpp_apar = round(r_gpp_apar, 3),
                bias = round(bias, 1), RMSE = round(RMSE, 1), r = round(r, 3),
                r_lo = round(r_lo, 3), r_hi = round(r_hi, 3)), row.names = FALSE)

## check: script 10 B and C equal script 16 observed/obs and /bor
cat('\nCheck: script 10 B_hm2 / C_boreal against script 16 observed/obs / bor:\n')
for (pair in list(c('B_hm2', 'observed/obs'), c('C_boreal', 'observed/obs/bor')))
  for (s in c('Pinus sylvestris', 'Picea abies')) {
    a <- T5[T5$variant == pair[1] & T5$species == s, ]
    b <- A4[A4$variant == pair[2] & A4$species == s, ]
    ok <- abs(a$r - b$r) < TOL && abs(a$bias - b$bias) < TOL
    cat(sprintf('  %-9s = %-17s %-17s r %.6f vs %.6f  bias %.3f vs %.3f  %s\n',
                pair[1], pair[2], s, a$r, b$r, a$bias, b$bias,
                if (ok) 'PASS' else 'FAIL'))
    if (!ok) stop('Part 2 equivalence check failed for ', pair[1], ' / ', s)
  }

write.csv(T5, file.path(OUT, 'table5_corrected.csv'), row.names = FALSE)
cat('\nwrote results/alpha_variants.csv and results/table5_corrected.csv\n')
sink()
