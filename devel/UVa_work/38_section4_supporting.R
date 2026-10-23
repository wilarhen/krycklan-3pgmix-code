## ---------------------------------------------------------------------------
## 38_section4_supporting.R
##
## Supporting numbers for the Discussion of the third question that are not
## produced by any other script:
##
##   results/ceiling_difference.csv   ceiling(pine) - ceiling(spruce), plot
##                                    bootstrap within each species
##   results/ceiling_power.csv        plots per species to detect the
##                                    difference the elasticity model predicts
##   results/apar_vs_lai.csv          r(GPP, Beer-Lambert APAR) against
##                                    r(GPP, LAI), pooled and per species
##   results/wetness_by_species.csv   SLU soil moisture map, slope and
##                                    elevation at the plots, by species
##
## Plot means, corrected species, 47 plots of observed/obs, B = 2000,
## seed 20260824, percentile intervals. k from ceiling_config.R K_EXT (sets
## 39 and 37). The wetness section needs the T7 mounted.
##
##   Rscript 38_section4_supporting.R      (from UVa_work/)
## ---------------------------------------------------------------------------

suppressMessages(library(readxl))
NB <- 2000; SEED <- 20260824
OUT <- 'results'
source('ceiling_config.R')          # K_EXT
sink(file.path(OUT, 'section4_supporting_log.txt'), split = TRUE)

g <- read.csv('validation_out/preflight_results.csv', stringsAsFactors = FALSE)
g <- g[g$variant == 'observed/obs', ]
a <- aggregate(cbind(obs_gpp, obs_lai) ~ Plot_ID + sp, g, mean)
P <- a[a$sp == 'Pinus sylvestris', ]; S <- a[a$sp == 'Picea abies', ]
q <- function(v) quantile(v, c(.025, .975))

## ---- 1. difference between the two ceilings ----------------------------------
obs <- cor(P$obs_gpp, P$obs_lai) - cor(S$obs_gpp, S$obs_lai)
set.seed(SEED)
dd <- replicate(NB, { i <- sample.int(nrow(P), replace = TRUE)
                      j <- sample.int(nrow(S), replace = TRUE)
                      cor(P$obs_gpp[i], P$obs_lai[i]) - cor(S$obs_gpp[j], S$obs_lai[j]) })
ls <- read.csv('validation_out/ceiling/light_saturation.csv', stringsAsFactors = FALSE)
pred_sp <- ls$r_ceiling_pred[ls$species == 'Picea abies']
pine_c  <- cor(P$obs_gpp, P$obs_lai)
D <- data.frame(ceiling_pine = pine_c, ceiling_spruce = cor(S$obs_gpp, S$obs_lai),
                n_pine = nrow(P), n_spruce = nrow(S),
                diff = obs, diff_lo = q(dd)[[1]], diff_hi = q(dd)[[2]],
                predicted_spruce_ceiling = pred_sp, predicted_diff = pine_c - pred_sp,
                predicted_diff_excluded = (pine_c - pred_sp) < q(dd)[[1]] | (pine_c - pred_sp) > q(dd)[[2]],
                share_ge_predicted = mean(dd >= pine_c - pred_sp))
cat('== ceiling difference ==\n'); print(D, row.names = FALSE)
write.csv(D, file.path(OUT, 'ceiling_difference.csv'), row.names = FALSE)

## ---- 2. plots needed -------------------------------------------------------------
## Fisher z. 'one_sample' is the formula in 24_ceiling_inference.R (lines 89-92),
## which produced the 643 / 572 for the model-ceiling gap; 'two_sample' doubles
## the variance term because the two ceilings come from independent plots.
z  <- function(r) 0.5 * log((1 + r) / (1 - r))
kk <- (qnorm(0.975) + qnorm(0.80))^2
dz <- abs(z(pine_c) - z(pred_sp))
W <- data.frame(r_pine = pine_c, r_spruce_predicted = pred_sp, delta_z = dz,
                power = 0.80, alpha = 0.05,
                n_one_sample_formula = ceiling(kk / dz^2 + 3),
                n_per_species_two_sample = ceiling(2 * kk / dz^2 + 3),
                n_available_pine = nrow(P), n_available_spruce = nrow(S))
cat('\n== plots needed ==\n'); print(W, row.names = FALSE)
write.csv(W, file.path(OUT, 'ceiling_power.csv'), row.names = FALSE)

## ---- 3. APAR against LAI ------------------------------------------------------------
a$apar <- 1 - exp(-K_EXT[a$sp] * a$obs_lai)
st <- function(d) c(cor(d$obs_gpp, d$obs_lai), cor(d$obs_gpp, d$apar))
rows <- list()
for (lab in c('pooled', 'Pinus sylvestris', 'Picea abies')) {
  d <- if (lab == 'pooled') a else a[a$sp == lab, ]
  o <- st(d); set.seed(SEED)
  b <- t(replicate(NB, st(d[sample.int(nrow(d), replace = TRUE), ])))
  df <- b[, 2] - b[, 1]
  rows[[lab]] <- data.frame(group = lab, n = nrow(d),
    k_pine = K_EXT[['Pinus sylvestris']], k_spruce = K_EXT[['Picea abies']],
    r_lai = o[1], r_lai_lo = q(b[, 1])[[1]], r_lai_hi = q(b[, 1])[[2]],
    r_apar = o[2], r_apar_lo = q(b[, 2])[[1]], r_apar_hi = q(b[, 2])[[2]],
    diff = o[2] - o[1], diff_lo = q(df)[[1]], diff_hi = q(df)[[2]], share_diff_gt0 = mean(df > 0))
}
A <- do.call(rbind, rows); rownames(A) <- NULL
cat('\n== APAR against LAI ==\n'); print(A, row.names = FALSE, digits = 4)
write.csv(A, file.path(OUT, 'apar_vs_lai.csv'), row.names = FALSE)

## ---- 4. wetness and topography by species ----------------------------------------
r <- terra::rast('/run/media/wilar/T7 Shield/slu_data_bases/3PG/SLU_soil_map/SLUMarkfuktighetskarta.tif')
p <- read.csv('validation_out/itc_vs_measured.csv')          # coordinates only
p$sp <- g$sp[match(p$Plot_ID, g$Plot_ID)]; p <- p[!is.na(p$sp), ]
v <- terra::project(terra::vect(p, geom = c('E', 'N'), crs = 'EPSG:3006'), terra::crs(r))
p$wet <- terra::extract(r, v)[, 2]
x <- suppressMessages(read_excel('/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx',
                                 sheet = '50 main plots', skip = 2))
x <- x[!is.na(x$Plot_ID) & !duplicated(x$Plot_ID), c('Plot_ID', 'Slope', 'Elevation')]
p <- merge(p, x, by = 'Plot_ID')
wt <- wilcox.test(wet ~ sp, data = p)$p.value
fl <- fligner.test(Slope ~ factor(sp), data = p)$p.value
H <- do.call(rbind, lapply(c('Pinus sylvestris', 'Picea abies'), function(s) {
  d <- p[p$sp == s, ]
  data.frame(species = s, n = nrow(d), wetness_mean = mean(d$wet), wetness_median = median(d$wet),
             wetness_sd = sd(d$wet), slope_mean = mean(d$Slope), slope_sd = sd(d$Slope),
             elevation_mean = mean(d$Elevation), elevation_sd = sd(d$Elevation))
}))
H$wetness_wilcoxon_p <- wt; H$slope_spread_fligner_p <- fl
cat('\n== wetness and topography by species ==\n'); print(H, row.names = FALSE, digits = 4)
write.csv(H, file.path(OUT, 'wetness_by_species.csv'), row.names = FALSE)
sink()
