## ---------------------------------------------------------------------------
## 41_tier2_support.R
##
## Values quoted in the manuscript that no other script writes. No model run; arithmetic on existing output.
##
##   results/pooled_vs_within_species.csv
##       r(GPP, APAR) and r(GPP, LAI), pooled across species and within each,
##       plot means, for script 10's variants A, B and C
##       (validation_out/params_results.csv).
##   results/lai_agreement.csv
##       simulated against measured maximum LAI under variant A (foliage
##       corrected, original height form and sets) and B, plot means, 47 plots.
##   results/alphacx_zero_bias_interpolated.csv
##       RMSE, bias and r at the zero-bias alphaCx (pine 0.0625 interpolated
##       between grid points 0.06 and 0.0635; spruce 0.0465 on the grid), and
##       how far simulated GPP departs from proportionality to alphaCx across
##       the sweep (results/alphacx_zero_bias.csv, script 36).
##   results/species_map_agreement.csv
##       NMD species at the plot (preflight mapped/obs) against field species
##       (observed/obs).
##   results/forest_floor_check.csv
##       r(forest-floor GPP, measured LAI) by species on the 47 analysed plots
##       and on the 41 plots covered by the original cube.
##
## Plot bootstrap B = 2000, percentile intervals, seed 20260824 reset per
## interval.
##
##   Rscript 41_tier2_support.R      (from UVa_work/)
## ---------------------------------------------------------------------------

suppressMessages(library(readxl))
NB <- 2000; SEED <- 20260824; OUT <- 'results'
sink(file.path(OUT, 'tier2_support_log.txt'), split = TRUE)
SP <- c('Pinus sylvestris', 'Picea abies')
q  <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
boot_r <- function(x, y) { set.seed(SEED)
  q(replicate(NB, { i <- sample.int(length(x), replace = TRUE); cor(x[i], y[i]) })) }

## ---- pooled against within-species correlations ------------------------------
p10 <- read.csv('validation_out/params_results.csv', stringsAsFactors = FALSE)
rows <- list()
for (v in c('A_control', 'B_hm2', 'C_boreal')) {
  a <- aggregate(cbind(sim_gpp, apar, lai) ~ Plot_ID + sp, p10[p10$variant == v, ], mean)
  for (g in c('pooled', SP)) {
    d <- if (g == 'pooled') a else a[a$sp == g, ]
    rows[[length(rows) + 1]] <- data.frame(variant = v, group = g, n_plots = nrow(d),
      r_gpp_apar = cor(d$sim_gpp, d$apar), r_gpp_lai = cor(d$sim_gpp, d$lai))
  }
}
P <- do.call(rbind, rows)
cat('== pooled against within species (plot means) ==\n'); print(P, row.names = FALSE, digits = 3)
d26 <- read.csv('validation_out/plotlevel/plotlevel_diagnostics.csv', stringsAsFactors = FALSE)
cat(sprintf('check: B pooled r(GPP,APAR) %.4f vs script 26 observed/obs r_apar_pl %.4f\n',
            P$r_gpp_apar[P$variant == 'B_hm2' & P$group == 'pooled'],
            d26$r_apar_pl[d26$variant == 'observed/obs']))
write.csv(P, file.path(OUT, 'pooled_vs_within_species.csv'), row.names = FALSE)

## ---- LAI agreement ---------------------------------------------------------------
rows <- list()
for (v in c('A_control', 'B_hm2')) {
  a <- aggregate(cbind(lai, obs_lai) ~ Plot_ID, p10[p10$variant == v, ], mean)
  rows[[v]] <- data.frame(variant = v, n_plots = nrow(a),
    sim_median = median(a$lai), sim_min = min(a$lai), sim_max = max(a$lai),
    obs_median = median(a$obs_lai), obs_min = min(a$obs_lai), obs_max = max(a$obs_lai),
    r = cor(a$lai, a$obs_lai))
}
LA <- do.call(rbind, rows)
cat('\n== LAI agreement (plot means) ==\n'); print(LA, row.names = FALSE, digits = 3)
write.csv(LA, file.path(OUT, 'lai_agreement.csv'), row.names = FALSE)

## ---- zero-bias alphaCx -------------------------------------------------------------
sw <- read.csv('results/alphacx_zero_bias.csv', stringsAsFactors = FALSE)
pubA <- c('Pinus sylvestris' = 0.0257, 'Picea abies' = 0.0270)
zbA  <- c('Pinus sylvestris' = 0.0625, 'Picea abies' = 0.0465)
Z <- do.call(rbind, lapply(SP, function(s) {
  z <- sw[sw$species == s, ]; ip <- function(v, a) approx(z$alphaCx, z[[v]], a)$y
  per_unit <- z$sim / z$alphaCx; base <- ip('sim', pubA[[s]]) / pubA[[s]]
  data.frame(species = s, alphaCx_published = pubA[[s]], alphaCx_zero_bias = zbA[[s]],
    bias_published = ip('bias', pubA[[s]]), bias_zero = ip('bias', zbA[[s]]),
    RMSE_published = ip('RMSE', pubA[[s]]), RMSE_zero = ip('RMSE', zbA[[s]]),
    RMSE_reduction = 1 - ip('RMSE', zbA[[s]]) / ip('RMSE', pubA[[s]]),
    r_published = ip('r', pubA[[s]]), r_zero = ip('r', zbA[[s]]),
    max_departure_from_proportional = max(abs(per_unit / base - 1)))
}))
cat('\n== zero-bias alphaCx ==\n'); print(Z, row.names = FALSE, digits = 4)
write.csv(Z, file.path(OUT, 'alphacx_zero_bias_interpolated.csv'), row.names = FALSE)

## ---- species map against field ------------------------------------------------------
pf <- read.csv('validation_out/preflight_results.csv', stringsAsFactors = FALSE)
m <- unique(pf[pf$variant == 'mapped/obs', c('Plot_ID', 'sp')])
o <- unique(pf[pf$variant == 'observed/obs', c('Plot_ID', 'sp')])
j <- merge(m, o, by = 'Plot_ID', suffixes = c('_map', '_field'))
MA <- data.frame(n_plots = nrow(j), agree = sum(j$sp_map == j$sp_field),
                 share = mean(j$sp_map == j$sp_field))
cat('\n== species map against field ==\n'); print(MA, row.names = FALSE)
write.csv(MA, file.path(OUT, 'species_map_agreement.csv'), row.names = FALSE)

## ---- forest-floor GPP against measured LAI --------------------------------------------
x <- suppressMessages(read_excel('/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx',
                                 sheet = '50 main plots', skip = 2))
x <- x[!is.na(x$Plot_ID), ]
ff <- aggregate(cbind(gppff = as.numeric(GPPff), lai = as.numeric(LAImax)) ~ Plot_ID, x, mean)
ff <- merge(ff, o, by = 'Plot_ID')                       # 47 analysed plots, field species
vp <- read.csv('validation_out/validation_pairs.csv', stringsAsFactors = FALSE)
cube41 <- unique(vp$Plot_ID[!is.na(vp$obs_tree)])
rows <- list()
for (set in c('47 analysed plots', '41 plots covered by the original cube')) {
  d0 <- if (set == '47 analysed plots') ff else ff[ff$Plot_ID %in% cube41, ]
  for (s in SP) {
    d <- d0[d0$sp == s, ]; ci <- boot_r(d$gppff, d$lai)
    rows[[length(rows) + 1]] <- data.frame(plot_set = set, species = s, n_plots = nrow(d),
      r = cor(d$gppff, d$lai), lo = ci[[1]], hi = ci[[2]])
  }
}
FF <- do.call(rbind, rows)
cat('\n== forest-floor GPP against measured LAI ==\n'); print(FF, row.names = FALSE, digits = 3)
write.csv(FF, file.path(OUT, 'forest_floor_check.csv'), row.names = FALSE)
sink()
