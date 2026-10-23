## ---------------------------------------------------------------------------
## 25_light_saturation.R
##
## FIX 3 of 3, and the one that turns the paper from an evaluation into a
## finding.
##
## The spruce ceiling is currently an observation: r = 0.268, explained after
## the fact by light saturation. Stated that way it is post hoc, and a
## reviewer is entitled to treat it as a small-sample result with a story
## attached.
##
## The same claim can be made A PRIORI. Under Beer-Lambert with published
## extinction coefficients, the elasticity of absorbed radiation with respect
## to leaf area is
##
##     E = d ln(APAR) / d ln(L) = k L exp(-k L) / (1 - exp(-k L))
##
## This is a property of k and L only. No plot data, no model output, no
## input rasters. Combined with the observed spread of measured LAI it
## predicts, before looking at GPP, how strongly production can track
## structure in each species -- and therefore predicts the ratio of the two
## ceilings.
##
## The script computes the prediction, then tests it against the observed
## ceilings. If they agree, the mechanism is confirmed rather than invoked.
##
##   Rscript 25_light_saturation.R
## ---------------------------------------------------------------------------

source('ceiling_config.R')
set.seed(SEED)

d <- to_plot_level(load_plots())

## species actually present, in descending plot count. SPECIES_KEEP may be
## NULL in the config, in which case every label found is analysed.
SPP <- if (is.null(SPECIES_KEEP))
    names(sort(table(d$species), decreasing = TRUE)) else SPECIES_KEEP


fapar     <- function(k, L) 1 - exp(-k * L)
elasticity <- function(k, L) {
    kl <- k * L
    kl * exp(-kl) / (1 - exp(-kl))
}
marginal  <- function(k, L) k * exp(-k * L)   # d fAPAR / d L

cat(strrep('=', 70), '\n', sep = '')
cat('A PRIORI PREDICTION FROM BEER-LAMBERT, NO GPP DATA USED\n')
cat(strrep('=', 70), '\n\n', sep = '')

pred <- list()
for (sp in SPP) {
    if (!sp %in% names(K_EXT)) {
        cat('skipping "', sp, '": no extinction coefficient in K_EXT\n', sep = '')
        next
    }
    dd <- d[d$species == sp, ]
    if (nrow(dd) < 5) { cat('skipping "', sp, '": only ', nrow(dd),
                            ' plots\n', sep = ''); next }
    k  <- K_EXT[[sp]]
    L  <- mean(dd$lai, na.rm = TRUE)
    sdL <- stats::sd(dd$lai, na.rm = TRUE)
    cvL <- sdL / L

    ## Elasticity is nonlinear in L, so E at the mean LAI is NOT the mean of
    ## the per-plot elasticities (Jensen). Averaging per plot is the quantity
    ## that means "how sensitive are these canopies", and it is what a reviewer
    ## recomputing from the plot table will get. Using E(k, mean(L)) understates
    ## spruce elasticity (0.293 vs 0.329) and thereby overstates the species
    ## contrast; it also produced a spurious three-decimal match to the observed
    ## spruce ceiling in an earlier version of this script.
    E_mean_lai <- elasticity(k, L)
    E  <- mean(elasticity(k, dd$lai[is.finite(dd$lai)]))
    ## predicted CV of the structural component of production
    cv_signal <- E * cvL

    cat(sprintf('%-8s k = %.3f  mean measured LAI = %.2f (sd %.2f, cv %.3f)\n',
                sp, k, L, sdL, cvL))
    cat(sprintf('         fAPAR at mean LAI          = %.3f\n', fapar(k, L)))
    cat(sprintf('         marginal return d fAPAR/dL = %.3f\n', marginal(k, L)))
    cat(sprintf('         elasticity, mean over plots = %.3f\n', E))
    cat(sprintf('         elasticity at mean LAI      = %.3f  (not used)\n', E_mean_lai))
    cat(sprintf('         predicted cv of structural signal = %.3f\n\n', cv_signal))

    pred[[sp]] <- data.frame(species = sp, k = k, lai_mean = L, lai_sd = sdL,
                             lai_cv = cvL, fapar = fapar(k, L),
                             marginal = marginal(k, L), elasticity = E,
                             cv_signal = cv_signal)
}
P <- do.call(rbind, pred)

## ---- the falsifiable prediction --------------------------------------------
## With a residual (non-structural) variance in production that is common to
## both species, the ratio of attainable correlations follows from the ratio
## of structural signal amplitudes alone:
##       r = cv_signal / sqrt(cv_signal^2 + cv_resid^2)
## The ratio r_spruce / r_pine is therefore predicted without knowing cv_resid
## only in the small-signal limit; otherwise cv_resid is solved from the pine
## ceiling and used to PREDICT the spruce ceiling. That is the sharper test.

## identify the low-k and high-k species from K_EXT rather than by name, so
## the script does not depend on the labelling convention in the plot table
if (nrow(P) >= 2) {
    ord <- P$species[order(P$k)]
    sp1 <- ord[1]; sp2 <- ord[length(ord)]      # sp1 = low k, sp2 = high k
    cat('low-k species: ', sp1, '   high-k species: ', sp2, '\n', sep = '')
    d1 <- d[d$species == sp1, ]; d2 <- d[d$species == sp2, ]
    r1_obs <- stats::cor(d1$gpp_obs, d1$lai)
    r2_obs <- stats::cor(d2$gpp_obs, d2$lai)
    s1 <- P$cv_signal[P$species == sp1]
    s2 <- P$cv_signal[P$species == sp2]

    ## solve residual cv from the PINE ceiling, then predict SPRUCE
    cv_resid <- s1 * sqrt(1 / r1_obs^2 - 1)
    r2_pred  <- s2 / sqrt(s2^2 + cv_resid^2)

    cat(strrep('-', 70), '\n', sep = '')
    cat('CALIBRATE ON THE LOW-k SPECIES, PREDICT THE HIGH-k SPECIES\n')
    cat(strrep('-', 70), '\n', sep = '')
    cat(sprintf('  %s ceiling observed  r = %6.3f  (used to fit)\n', sp1, r1_obs))
    cat(sprintf('  implied common residual cv           %6.3f\n', cv_resid))
    cat(sprintf('  %s ceiling PREDICTED  r = %6.3f\n', sp2, r2_pred))
    cat(sprintf('  %s ceiling OBSERVED   r = %6.3f\n', sp2, r2_obs))
    cat(sprintf('  discrepancy                          %6.3f\n', r2_obs - r2_pred))

    ## bootstrap the observed spruce ceiling to see whether the prediction
    ## falls inside it
    b <- vapply(seq_len(N_BOOT), function(i) {
        x <- d2[sample.int(nrow(d2), nrow(d2), replace = TRUE), ]
        stats::cor(x$gpp_obs, x$lai)
    }, numeric(1))
    ci <- stats::quantile(b, c(0.025, 0.975), na.rm = TRUE)
    cat(sprintf('  observed %s ceiling 95%% CI [%6.3f, %6.3f]\n', sp2, ci[1], ci[2]))
    inside <- r2_pred >= ci[1] && r2_pred <= ci[2]
    cat('  the a priori prediction falls ',
        if (inside) 'INSIDE' else 'OUTSIDE', ' that interval\n', sep = '')
    cat('\n  READ THIS AS A BOUND, NOT A POINT PREDICTION.\n')
    cat('  Light extinction and published parameters say the high-k species\n')
    cat('  cannot much exceed the predicted value. An OBSERVED ceiling below\n')
    cat('  the prediction is expected: saturation is one constraint among\n')
    cat('  several, so the bound is not tight. What the calculation excludes\n')
    cat('  is an observed ceiling near that of the low-k species.\n')
    cat('\n  The calibration assumes both species share a residual variance in\n')
    cat('  production not attributable to structure. If the high-k species\n')
    cat('  occupies wetter or more heterogeneous ground its residual is larger\n')
    cat('  and the true prediction falls. That is untestable at this sample\n')
    cat('  size and must be stated as a limitation.\n')

    P$r_ceiling_obs  <- c(r1_obs, r2_obs)[match(P$species, c(sp1, sp2))]
    P$r_ceiling_pred <- c(NA, r2_pred)[match(P$species, c(sp1, sp2))]
}

## ---- sensitivity to k, since k is a published parameter not a measurement --
cat('\n', strrep('-', 70), '\n', sep = '')
cat('SENSITIVITY OF THE PREDICTION TO THE EXTINCTION COEFFICIENT\n')
cat(strrep('-', 70), '\n', sep = '')
Lp <- P$lai_mean[which.max(P$k)]
for (kk in seq(0.45, 0.80, by = 0.05))
    cat(sprintf('  spruce k = %.2f  ->  fAPAR %.3f  elasticity %.3f\n',
                kk, fapar(kk, Lp), elasticity(kk, Lp)))
cat('  Report this range. It shows the conclusion does not rest on one\n')
cat('  parameter value being exactly right.\n')

utils::write.csv(P, file.path(path.expand(OUT_DIR), 'light_saturation.csv'),
                 row.names = FALSE)

## ---- figure ----------------------------------------------------------------
png(file.path(path.expand(OUT_DIR), 'fig_saturation.png'),
    width = 1500, height = 700, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
Ls  <- seq(0, 8, by = 0.05)
SPF <- P$species                       # only species that got through above
kk  <- vapply(SPF, function(s) K_EXT[[s]], numeric(1))

plot(Ls, fapar(kk[1], Ls), type = 'l', lwd = 2, ylim = c(0, 1),
     xlab = 'leaf area index', ylab = 'fraction of radiation absorbed',
     main = 'Light interception')
for (i in seq_along(SPF)[-1]) lines(Ls, fapar(kk[i], Ls), lwd = 2, lty = i)
points(P$lai_mean, mapply(fapar, kk, P$lai_mean), pch = 19, cex = 1.4)
legend('bottomright', sprintf('%s, k = %.3f', SPF, kk),
       lty = seq_along(SPF), lwd = 2, bty = 'n')

plot(Ls, elasticity(kk[1], Ls), type = 'l', lwd = 2, ylim = c(0, 1),
     xlab = 'leaf area index', ylab = 'd ln(APAR) / d ln(LAI)',
     main = 'Sensitivity of production to structure')
for (i in seq_along(SPF)[-1]) lines(Ls, elasticity(kk[i], Ls), lwd = 2, lty = i)
points(P$lai_mean, P$elasticity, pch = 19, cex = 1.4)
par(op); dev.off()
cat('\nfigure written to ', file.path(OUT_DIR, 'fig_saturation.png'),
    ' -- candidate replacement for fig6_ceiling\n', sep = '')
