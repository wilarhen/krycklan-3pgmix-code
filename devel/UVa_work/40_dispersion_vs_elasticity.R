## ---------------------------------------------------------------------------
## 40_dispersion_vs_elasticity.R
##
## Does simulated tree GPP vary less between spruce plots than between pine
## plots by as much as the elasticity difference predicts? No model run.
##
## Per species, at the published alphaCx (pine 0.0257, spruce 0.0270) and at
## the zero-bias alphaCx (pine 0.0625, interpolated between grid points 0.06
## and 0.0635; spruce 0.0465): SD and CV of simulated and of measured tree
## GPP. Then spruce-to-pine ratios of
##   SD(sim)/SD(obs), CV(sim), CV(obs),
##   elasticity E, and the light signal E x CV(LAI) (script 25's cv_signal).
##
## Sources: published setting, per-plot values in
## validation_out/preflight_results.csv (observed/obs); zero-bias setting,
## results/alphacx_slope_from_sweep.csv (script 39, derived from the sweep
## summaries, so no intervals exist there); E and cv_signal,
## validation_out/ceiling/light_saturation.csv (k from K_EXT, sets 39 and 37).
## Intervals: plot bootstrap within species, B = 2000, seed 20260824; the E
## and light-signal intervals add the Forrester et al. (2021) k posterior as
## in 37_elasticity_interval.R (split normal, same seed and draw order).
##
##   Rscript 40_dispersion_vs_elasticity.R      (from UVa_work/)
## Output: results/dispersion_vs_elasticity.csv
## ---------------------------------------------------------------------------

NB <- 2000; SEED <- 20260824
OUT <- 'results'
sink(file.path(OUT, 'dispersion_vs_elasticity_log.txt'), split = TRUE)
q  <- function(v) quantile(v, c(.025, .975))
SP <- c('Pinus sylvestris', 'Picea abies')

g  <- read.csv('validation_out/preflight_results.csv', stringsAsFactors = FALSE)
g  <- g[g$variant == 'observed/obs', ]
pm <- aggregate(cbind(sim_gpp, obs_gpp, obs_lai) ~ Plot_ID + sp, g, mean)
sw <- read.csv('results/alphacx_slope_from_sweep.csv', stringsAsFactors = FALSE)
ls <- read.csv('validation_out/ceiling/light_saturation.csv', stringsAsFactors = FALSE)

## ---- per species ----------------------------------------------------------------
pubA <- c('Pinus sylvestris' = 0.0257, 'Picea abies' = 0.0270)
zbA  <- c('Pinus sylvestris' = 0.0625, 'Picea abies' = 0.0465)
zb   <- function(s, v) { z <- sw[sw$species == s, ]; approx(z$alphaCx, z[[v]], zbA[[s]])$y }
per <- do.call(rbind, lapply(SP, function(s) {
  d <- pm[pm$sp == s, ]
  rbind(
    data.frame(species = s, setting = 'published', alphaCx = pubA[[s]], n = nrow(d),
               sd_sim = sd(d$sim_gpp), sd_obs = sd(d$obs_gpp),
               mean_sim = mean(d$sim_gpp), mean_obs = mean(d$obs_gpp)),
    data.frame(species = s, setting = 'zero_bias', alphaCx = zbA[[s]], n = nrow(d),
               sd_sim = zb(s, 'sd_sim'), sd_obs = zb(s, 'sd_obs'),
               mean_sim = zb(s, 'mean_sim'), mean_obs = zb(s, 'mean_obs')))
}))
per$sd_ratio_sim_obs <- per$sd_sim / per$sd_obs
per$cv_sim <- per$sd_sim / per$mean_sim
per$cv_obs <- per$sd_obs / per$mean_obs
per$E <- ls$elasticity[match(per$species, ls$species)]
per$cv_lai <- ls$lai_cv[match(per$species, ls$species)]
per$cv_signal <- ls$cv_signal[match(per$species, ls$species)]
cat('== per species ==\n'); print(per, row.names = FALSE, digits = 4)

## ---- bootstrap ----------------------------------------------------------------------
KP <- data.frame(species = SP, med = c(0.3247, 0.6378), lo = c(0.3012, 0.5934), hi = c(0.4558, 0.6496))
rsplit <- function(n, med, s_lo, s_hi) { z <- abs(rnorm(n)); up <- runif(n) < 0.5
  pmax(ifelse(up, med + s_hi * z, med - s_lo * z), 1e-6) }
Ef <- function(k, L) k * L * exp(-k * L) / (1 - exp(-k * L))
B <- list()
for (i in 1:2) {
  s <- SP[i]; d <- pm[pm$sp == s, ]; n <- nrow(d)
  set.seed(SEED)
  idx <- replicate(NB, sample.int(n, replace = TRUE), simplify = FALSE)
  kd  <- rsplit(NB, KP$med[i], (KP$med[i] - KP$lo[i]) / qnorm(.975), (KP$hi[i] - KP$med[i]) / qnorm(.975))
  B[[s]] <- data.frame(
    E      = mapply(function(k, j) mean(Ef(k, d$obs_lai[j])), kd, idx),
    signal = mapply(function(k, j) mean(Ef(k, d$obs_lai[j])) * sd(d$obs_lai[j]) / mean(d$obs_lai[j]), kd, idx),
    cv_sim = vapply(idx, function(j) sd(d$sim_gpp[j]) / mean(d$sim_gpp[j]), 0),
    cv_obs = vapply(idx, function(j) sd(d$obs_gpp[j]) / mean(d$obs_gpp[j]), 0),
    sd_rat = vapply(idx, function(j) sd(d$sim_gpp[j]) / sd(d$obs_gpp[j]), 0))
}
rb <- function(v) B[['Picea abies']][[v]] / B[['Pinus sylvestris']][[v]]

## ---- spruce-to-pine ratios --------------------------------------------------------
pv <- function(st, v) per[[v]][per$species == 'Picea abies' & per$setting == st] /
                      per[[v]][per$species == 'Pinus sylvestris' & per$setting == st]
row <- function(quantity, setting, est, b = NULL) data.frame(quantity = quantity, setting = setting,
  ratio_spruce_to_pine = est, lo = if (is.null(b)) NA else q(b)[[1]], hi = if (is.null(b)) NA else q(b)[[2]])
R <- rbind(
  row('SD(sim)/SD(obs)', 'published', pv('published', 'sd_ratio_sim_obs'), rb('sd_rat')),
  row('SD(sim)/SD(obs)', 'zero_bias', pv('zero_bias', 'sd_ratio_sim_obs')),
  row('CV(sim)', 'published', pv('published', 'cv_sim'), rb('cv_sim')),
  row('CV(sim)', 'zero_bias', pv('zero_bias', 'cv_sim')),
  row('CV(obs)', 'measured', pv('published', 'cv_obs'), rb('cv_obs')),
  row('elasticity E', 'k posterior + plot bootstrap', pv('published', 'E'), rb('E')),
  row('E x CV(LAI)', 'k posterior + plot bootstrap', pv('published', 'cv_signal'), rb('signal')))
R$share_cv_sim_published_below <- NA
R$share_cv_sim_published_below[R$quantity == 'elasticity E'] <- mean(rb('cv_sim') < rb('E'))
R$share_cv_sim_published_below[R$quantity == 'E x CV(LAI)']  <- mean(rb('cv_sim') < rb('signal'))
cat('\n== spruce-to-pine ratios ==\n'); print(R, row.names = FALSE, digits = 4)

## one file: per-species rows, then the ratios
per$record <- 'per_species'; R$record <- 'ratio_spruce_to_pine'
all <- merge(per, R, all = TRUE, sort = FALSE)
all <- all[, c('record', setdiff(names(all), 'record'))]
write.csv(all, file.path(OUT, 'dispersion_vs_elasticity.csv'), row.names = FALSE, na = '')
cat('\nwrote', file.path(OUT, 'dispersion_vs_elasticity.csv'), '\n')
sink()
