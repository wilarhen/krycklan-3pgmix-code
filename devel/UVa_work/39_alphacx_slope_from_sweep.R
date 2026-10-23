## ---------------------------------------------------------------------------
## 39_alphacx_slope_from_sweep.R
##
## Slope of simulated on measured tree GPP, and SD of simulated and measured
## GPP, at each alphaCx of the sweep in results/alphacx_zero_bias.csv
## (script 36). No model run. The sweep stores, per species and alphaCx, the
## plot-mean summaries n, mean sim, mean obs, bias, RMSE and r, not the
## per-plot series, so:
##
##   SD_obs  from the measured plot means (validation_out/preflight_results.csv,
##           observed/obs; identical plots and values in every sweep row)
##   var(sim - obs) = RMSE^2 - bias^2                 (RMSE is a mean over plots)
##   var(sim - obs) = var_sim + var_obs - 2 r sd_sim sd_obs
##           -> quadratic in sd_sim; both roots kept, the selected one is the
##              root that reproduces figure 5 at the published alphaCx
##   slope = r * SD_sim / SD_obs
##
## Variances are on the n divisor inside the derivation (as RMSE is) and
## reported as sample SDs (n - 1). Bootstrap intervals cannot be computed from
## summaries; they need the per-plot series.
##
##   Rscript 39_alphacx_slope_from_sweep.R      (from UVa_work/)
## Output: results/alphacx_slope_from_sweep.csv
## ---------------------------------------------------------------------------

sw <- read.csv('results/alphacx_zero_bias.csv', stringsAsFactors = FALSE)
g  <- read.csv('validation_out/preflight_results.csv', stringsAsFactors = FALSE)
g  <- g[g$variant == 'observed/obs', ]
pm <- aggregate(cbind(sim_gpp, obs_gpp) ~ Plot_ID + sp, g, mean)
f5 <- read.csv('results/fig5_annotations.csv', stringsAsFactors = FALSE)

pv <- function(x) mean((x - mean(x))^2)               # n-divisor variance
rows <- list()
for (i in seq_len(nrow(sw))) {
  s <- sw$species[i]; n <- sw$n[i]
  o <- pm$obs_gpp[pm$sp == s]
  stopifnot(length(o) == n, abs(mean(o) - sw$obs[i]) < 1e-6)
  Vo <- pv(o); Vd <- sw$RMSE[i]^2 - sw$bias[i]^2; r <- sw$r[i]
  disc <- r^2 * Vo - Vo + Vd
  x <- r * sqrt(Vo) + c(1, -1) * sqrt(disc)          # sd_sim roots, n divisor
  rows[[i]] <- data.frame(species = s, alphaCx = sw$alphaCx[i], n = n, r = r,
    bias = sw$bias[i], mean_sim = sw$sim[i], mean_obs = sw$obs[i],
    sd_obs = sqrt(Vo * n / (n - 1)),
    sd_sim_root_plus = x[1] * sqrt(n / (n - 1)), sd_sim_root_minus = x[2] * sqrt(n / (n - 1)),
    slope_root_plus = r * x[1] / sqrt(Vo), slope_root_minus = r * x[2] / sqrt(Vo))
}
D <- do.call(rbind, rows)

## root selection: the published alphaCx must reproduce figure 5 and the
## preflight sim SD
pubA <- c('Pinus sylvestris' = 0.0257, 'Picea abies' = 0.027)
D$root <- NA_character_
for (s in names(pubA)) {
  p  <- D[D$species == s & abs(D$alphaCx - pubA[[s]]) < 1e-9, ]
  sd_true <- sd(pm$sim_gpp[pm$sp == s]); sl_true <- f5$slope_sim_on_obs[f5$species == s]
  pick <- if (abs(p$sd_sim_root_plus - sd_true) < abs(p$sd_sim_root_minus - sd_true)) 'plus' else 'minus'
  cat(sprintf('%s published alphaCx %.4f: sd_sim roots %.4f / %.4f, preflight %.4f -> root %s; slope %.6f vs figure 5 %.6f\n',
      s, pubA[[s]], p$sd_sim_root_plus, p$sd_sim_root_minus, sd_true, pick,
      if (pick == 'plus') p$slope_root_plus else p$slope_root_minus, sl_true))
  D$root[D$species == s] <- pick
}
D$sd_sim <- ifelse(D$root == 'plus', D$sd_sim_root_plus, D$sd_sim_root_minus)
D$slope  <- ifelse(D$root == 'plus', D$slope_root_plus, D$slope_root_minus)
## the other root: is it ever the plausible one? (negative or ~0 means no)
D$other_root_sd <- ifelse(D$root == 'plus', D$sd_sim_root_minus, D$sd_sim_root_plus)

print(D[, c('species', 'alphaCx', 'n', 'r', 'bias', 'sd_obs', 'sd_sim', 'slope', 'other_root_sd')],
      row.names = FALSE, digits = 4)

## zero-bias settings: spruce 0.0465 is a grid point; pine 0.0625 lies between
## 0.06 and 0.0635 and is interpolated linearly in alphaCx
zb <- function(s, a) {
  z <- D[D$species == s, ]
  data.frame(species = s, alphaCx = a,
             sd_sim = approx(z$alphaCx, z$sd_sim, a)$y, slope = approx(z$alphaCx, z$slope, a)$y,
             bias = approx(z$alphaCx, z$bias, a)$y, sd_obs = z$sd_obs[1],
             r = approx(z$alphaCx, z$r, a)$y)
}
Z <- rbind(zb('Pinus sylvestris', 0.0625), zb('Picea abies', 0.0465))
cat('\nzero-bias alphaCx (pine interpolated between 0.06 and 0.0635):\n')
print(Z, row.names = FALSE, digits = 4)

write.csv(D, 'results/alphacx_slope_from_sweep.csv', row.names = FALSE)
