## ---------------------------------------------------------------------------
## 37_elasticity_interval.R
##
## Interval on the elasticity of absorbed radiation with respect to leaf area,
##   E = k L exp(-k L) / (1 - exp(-k L)),
## evaluated per plot on measured maximum LAI and averaged within species, as
## in 25_light_saturation.R. Two sources of uncertainty are propagated:
##
##   k  Forrester et al. (2021) posterior for the sets used in the corrected
##      configuration: Pinus sylvestris set 39, 0.3247 (0.3012 to 0.4558);
##      Picea abies set 37, 0.6378 (0.5934 to 0.6496). The intervals are
##      posterior quantile ranges, not standard errors, and are asymmetric.
##      Each posterior is represented by a split normal with its median at the
##      reported value and separate scales below and above, chosen so that its
##      2.5 and 97.5 percentiles equal the reported bounds exactly. One k is
##      drawn per species per replicate, independently between species,
##      truncated at zero (never reached at these values).
##   L  bootstrap of plots within species (measured maximum LAI, plot means).
##
## For comparison, the same plot resamples give the measured slope of
## log(tree GPP) on log(LAI) and the ceiling r(tree GPP, LAI).
## Plot means, corrected species, 47 plots of observed/obs, B = 2000,
## seed 20260824, percentile intervals.
##
##   Rscript 37_elasticity_interval.R      (from UVa_work/)
## ---------------------------------------------------------------------------

NB <- 2000; SEED <- 20260824
OUT <- 'results'
dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, 'elasticity_interval_log.txt'), split = TRUE)

g <- read.csv('validation_out/preflight_results.csv', stringsAsFactors = FALSE)
g <- g[g$variant == 'observed/obs', ]
a <- aggregate(cbind(obs_gpp, obs_lai) ~ Plot_ID + sp, g, mean)
SP <- c('Pinus sylvestris', 'Picea abies')

KP <- data.frame(species = SP, set = c('Pinus sylvestris 39', 'Picea abies 37'),
                 med = c(0.3247, 0.6378), lo = c(0.3012, 0.5934), hi = c(0.4558, 0.6496))
KP$s_lo <- (KP$med - KP$lo) / qnorm(.975)
KP$s_hi <- (KP$hi - KP$med) / qnorm(.975)

## split normal: P(below median) = 0.5, each half a half-normal with its scale
rsplit <- function(n, med, s_lo, s_hi) {
  z <- abs(rnorm(n)); up <- runif(n) < 0.5
  pmax(ifelse(up, med + s_hi * z, med - s_lo * z), 1e-6)
}
Ef <- function(k, L) k * L * exp(-k * L) / (1 - exp(-k * L))

## script 24's ceiling interval is copied for reference only; the paper reports
## the script 26 interval (plotlevel_stats.csv). See validation_out/ceiling/README.md
ci24 <- read.csv('validation_out/ceiling/ceiling_inference.csv', stringsAsFactors = FALSE)

res <- list(); draws <- list()
for (i in seq_along(SP)) {
  s <- SP[i]; d <- a[a$sp == s, ]; n <- nrow(d)
  ## plot resamples first, so the slope interval is the one obtained with this
  ## seed alone; then the k draws
  set.seed(SEED)
  idx <- replicate(NB, sample.int(n, replace = TRUE), simplify = FALSE)
  kd <- rsplit(NB, KP$med[i], KP$s_lo[i], KP$s_hi[i])
  E_both <- mapply(function(k, j) mean(Ef(k, d$obs_lai[j])), kd, idx)
  E_k    <- vapply(kd, function(k) mean(Ef(k, d$obs_lai)), numeric(1))
  E_L    <- vapply(idx, function(j) mean(Ef(KP$med[i], d$obs_lai[j])), numeric(1))
  slope  <- vapply(idx, function(j) coef(lm(log(d$obs_gpp[j]) ~ log(d$obs_lai[j])))[[2]], numeric(1))
  c24    <- ci24[ci24$species == s, ]
  q <- function(v) quantile(v, c(.025, .975))
  E0 <- mean(Ef(KP$med[i], d$obs_lai))
  sl0 <- coef(lm(log(d$obs_gpp) ~ log(d$obs_lai)))[[2]]
  c0 <- cor(d$obs_gpp, d$obs_lai)
  cat(sprintf('\n%s (%s), n = %d plots\n', s, KP$set[i], n))
  cat(sprintf('  k draws: median %.4f, 2.5%% %.4f, 97.5%% %.4f (target %.4f, %.4f, %.4f)\n',
              median(kd), q(kd)[1], q(kd)[2], KP$med[i], KP$lo[i], KP$hi[i]))
  cat(sprintf('  E at k bounds, measured LAI: %.4f (k %.4f) to %.4f (k %.4f)\n',
              mean(Ef(KP$lo[i], d$obs_lai)), KP$lo[i], mean(Ef(KP$hi[i], d$obs_lai)), KP$hi[i]))
  cat(sprintf('  E point %.4f\n  E, k and LAI   [%.4f, %.4f]\n  E, k only      [%.4f, %.4f]\n  E, LAI only    [%.4f, %.4f]\n',
              E0, q(E_both)[1], q(E_both)[2], q(E_k)[1], q(E_k)[2], q(E_L)[1], q(E_L)[2]))
  cat(sprintf('  measured slope %.4f [%.4f, %.4f]; contains E point: %s\n',
              sl0, q(slope)[1], q(slope)[2], E0 >= q(slope)[1] & E0 <= q(slope)[2]))
  cat(sprintf('  ceiling %.4f [%.4f, %.4f] (script 24, ceiling_inference.csv; recomputed r %.4f)\n',
              c24$r_ceiling, c24$ceil_lo, c24$ceil_hi, c0))
  res[[s]] <- data.frame(species = s, parameter_set = KP$set[i], n_plots = n,
    k = KP$med[i], k_lo = KP$lo[i], k_hi = KP$hi[i],
    E = E0, E_lo = q(E_both)[[1]], E_hi = q(E_both)[[2]],
    E_k_only_lo = q(E_k)[[1]], E_k_only_hi = q(E_k)[[2]],
    E_lai_only_lo = q(E_L)[[1]], E_lai_only_hi = q(E_L)[[2]],
    slope = sl0, slope_lo = q(slope)[[1]], slope_hi = q(slope)[[2]],
    slope_contains_E = E0 >= q(slope)[[1]] & E0 <= q(slope)[[2]],
    ceiling = c24$r_ceiling, ceiling_lo = c24$ceil_lo, ceiling_hi = c24$ceil_hi)
  draws[[s]] <- E_both
}

dE <- draws[['Pinus sylvestris']] - draws[['Picea abies']]
d0 <- res[['Pinus sylvestris']]$E - res[['Picea abies']]$E
qd <- quantile(dE, c(.025, .975))
cat(sprintf('\nE(pine) - E(spruce) = %.4f  95%% interval [%.4f, %.4f]  excludes zero: %s  share <= 0: %.4f\n',
            d0, qd[1], qd[2], qd[1] > 0 | qd[2] < 0, mean(dE <= 0)))

R <- do.call(rbind, res); rownames(R) <- NULL
R$E_diff_pine_minus_spruce <- d0; R$E_diff_lo <- qd[[1]]; R$E_diff_hi <- qd[[2]]
write.csv(R, file.path(OUT, 'elasticity_interval.csv'), row.names = FALSE)
cat('\nwrote', file.path(OUT, 'elasticity_interval.csv'), '\n')
sink()
