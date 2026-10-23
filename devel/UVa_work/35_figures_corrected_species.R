## ---------------------------------------------------------------------------
## 35_figures_corrected_species.R
##
## Regenerates figures 4, 5 and 6 on the corrected species grouping
## (26 Pinus sylvestris, 21 Picea abies), plot means, 2000 resamples,
## set.seed(20260824) -- the same basis as Table 5 and Section 3.5.
##
## Supersedes 27_make_figures.R, which drew figure 4 as the corrected state
## only. Figure 4 here has one panel per Table 5 variant:
##   A  foliage corrected, height form 1, pine set 30, spruce set 37
##   B  height form 2, pine set 30 -> 39
##   C  B + boreal temperatures and alphaCx 0.050
##
## Sources
##   validation_out/params_results.csv     script 10, variants A/B/C   -> fig 4
##   validation_out/preflight_results.csv  script 16, observed/obs     -> fig 5, 6
##   validation_out/plotlevel/plotlevel_stats.csv  script 26 intervals -> fig 5, 6
##   results/table6_corrected.csv           script 33, Table 8         -> fig 6b
##
## Output: results/figures/. Does not write to the Overleaf fig/ directory.
##
##   Rscript 35_figures_corrected_species.R
## ---------------------------------------------------------------------------

set.seed(20260824)
IN  <- 'validation_out'
FIG <- 'results/figures'
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
NB  <- 2000

SPP  <- c('Pinus sylvestris', 'Picea abies')
LBL  <- c('Pinus sylvestris' = 'Scots pine', 'Picea abies' = 'Norway spruce')
COLS <- c('Pinus sylvestris' = '#1b6ca8', 'Picea abies' = '#c1410b')
PCH  <- c('Pinus sylvestris' = 21, 'Picea abies' = 24)

agg <- function(x, vars) {
    vars <- intersect(vars, names(x))
    a <- aggregate(x[vars], by = list(Plot_ID = x$Plot_ID), FUN = mean, na.rm = TRUE)
    merge(a, x[!duplicated(x$Plot_ID), c('Plot_ID', 'sp')], by = 'Plot_ID')
}
## Intervals are not bootstrapped here. Model r and ceiling intervals come from
## script 26 (validation_out/plotlevel/plotlevel_stats.csv, observed/obs), the
## single run quoted in the text and tables; the point value computed here
## must equal script 26's or the script stops.
pl26 <- read.csv(file.path(IN, 'plotlevel', 'plotlevel_stats.csv'), stringsAsFactors = FALSE)
pl26 <- pl26[pl26$variant == 'observed/obs', ]
ci26 <- function(s, what, r) {
    z <- pl26[pl26$species == s, ]
    ref <- if (what == 'model') z$r_corrected else z$ceiling_corrected
    if (abs(ref - r) > 1e-9) stop(sprintf('%s %s r %.6f differs from script 26 %.6f', s, what, r, ref))
    if (what == 'model') c(z$r_lo, z$r_hi) else c(z$ceil_lo, z$ceil_hi)
}
lab <- function(txt, cx = 0.95) mtext(txt, side = 3, adj = 0, line = 0.4, font = 2, cex = cx)

pr <- read.csv(file.path(IN, 'params_results.csv'),    stringsAsFactors = FALSE)
pf <- read.csv(file.path(IN, 'preflight_results.csv'), stringsAsFactors = FALSE)

## ===========================================================================
## FIGURE 4  effect of each correction, one panel per Table 5 variant
## ===========================================================================
## Short panel titles; the full variant definitions belong in the caption,
## where they can match Table 5 verbatim without colliding at this width.
VAR <- c(A_control = '(a) A: foliage corrected',
         B_hm2     = '(b) B: + height form 2',
         C_boreal  = '(c) C: + boreal T')

pdf(file.path(FIG, 'fig4_isolation.pdf'), width = 7.5, height = 3.2)
par(mfrow = c(1, 3), mar = c(4.1, 4.1, 2.4, 0.8), mgp = c(2.4, 0.8, 0))
for (v in names(VAR)) {
    a <- agg(pr[pr$variant == v, ], c('sim_gpp', 'obs_gpp'))
    lim <- range(c(a$sim_gpp, a$obs_gpp), na.rm = TRUE)
    plot(a$obs_gpp, a$sim_gpp, xlim = lim, ylim = lim, las = 1,
         pch = PCH[a$sp], bg = adjustcolor(COLS[a$sp], .6), col = COLS[a$sp],
         cex = 1.15, xlab = expression('Measured tree GPP (gC m'^-2*' yr'^-1*')'),
         ylab = expression('Simulated GPP (gC m'^-2*' yr'^-1*')'))
    abline(0, 1, lty = 2, col = 'grey45')
    ln <- character(0)
    for (s in SPP) {
        z <- a[a$sp == s, ]
        r  <- cor(z$sim_gpp, z$obs_gpp, use = 'complete.obs')
        b  <- mean(z$sim_gpp - z$obs_gpp)
        sl <- coef(lm(sim_gpp ~ obs_gpp, z))[2]
        ln <- c(ln, sprintf('%s r = %+.2f, bias %+.0f, slope %.2f', LBL[s], r, b, sl))
    }
    legend('topleft', legend = ln, bty = 'n', cex = 0.72, text.col = COLS[SPP])
    lab(VAR[[v]], 0.85)
}
dev.off()

## ===========================================================================
## FIGURE 5  simulated against measured tree GPP, by species, with ceiling
## ===========================================================================
m <- agg(pf[pf$variant == 'observed/obs', ], c('sim_gpp','obs_gpp','obs_lai'))

pdf(file.path(FIG, 'fig5_validation.pdf'), width = 7.5, height = 3.9)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 0.8), mgp = c(2.5, 0.8, 0))
ann5 <- list()                       # the annotations, persisted below
for (s in SPP) {
    z <- m[m$sp == s, ]
    lim <- range(c(z$sim_gpp, z$obs_gpp), na.rm = TRUE)
    plot(z$obs_gpp, z$sim_gpp, xlim = lim, ylim = lim, las = 1,
         pch = PCH[s], bg = adjustcolor(COLS[s], .6), col = COLS[s], cex = 1.3,
         xlab = expression('Measured tree GPP (gC m'^-2*' yr'^-1*')'),
         ylab = expression('Simulated GPP (gC m'^-2*' yr'^-1*')'))
    abline(0, 1, lty = 2, col = 'grey45')
    r  <- cor(z$sim_gpp, z$obs_gpp, use = 'complete.obs'); ci <- ci26(s, 'model', r)
    ce <- cor(z$obs_gpp, z$obs_lai, use = 'complete.obs'); cc <- ci26(s, 'ceiling', ce)
    b  <- mean(z$sim_gpp - z$obs_gpp)
    sl <- coef(lm(sim_gpp ~ obs_gpp, z))[2]
    legend('topleft', bty = 'n', cex = 0.74, legend = c(
        sprintf('model r = %+.3f [%+.3f, %+.3f]', r, ci[1], ci[2]),
        sprintf('ceiling r = %+.3f [%+.3f, %+.3f]', ce, cc[1], cc[2]),
        sprintf('mean bias %+.0f gC m-2 yr-1', b),
        sprintf('slope of fit %.2f', sl),
        sprintf('%d plots', nrow(z))))
    lab(sprintf('%s. %s', c('(a)','(b)')[match(s, SPP)], LBL[s]))
    ann5[[s]] <- data.frame(species = s, variant = 'observed/obs', n_plots = nrow(z),
        model_r = r, model_r_lo = ci[[1]], model_r_hi = ci[[2]],
        ceiling_r = ce, ceiling_r_lo = cc[[1]], ceiling_r_hi = cc[[2]],
        mean_bias = b, slope_sim_on_obs = sl)
}
write.csv(do.call(rbind, ann5), 'results/fig5_annotations.csv', row.names = FALSE)
dev.off()

## ===========================================================================
## FIGURE 6  the ceiling: measured GPP against measured LAI, and predictors
## ===========================================================================
## Panel (b) draws Table 8: the eight measured attributes on the 47 analysed
## plots (26 pine, 21 spruce), from results/table6_corrected.csv (script 33),
## in the table's row order. The leaf area row is the ceiling, so its interval
## is script 26's, as in panel (a) and Section 3.5; the other intervals are
## script 33's, the only run that computes them.
t8 <- read.csv('results/table6_corrected.csv', stringsAsFactors = FALSE)
T8ORDER <- c('Basal area', 'Leaf area index (max)', 'Soil nitrogen concentration',
             'Soil nitrogen stock', 'Soil carbon', 'Stem density', 'Stand age', 'Elevation')
for (s in SPP) {
    k <- t8$sp == s & t8$column == 'LAImax'
    t8$lo[k] <- NA; t8$hi[k] <- NA
    t8[k, c('lo', 'hi')] <- ci26(s, 'ceiling', t8$r[k])
}
stopifnot(all(t8$n[t8$sp == 'Pinus sylvestris'] == 26),
          all(t8$n[t8$sp == 'Picea abies'] == 21))
write.csv(t8[, c('attribute', 'sp', 'n', 'r', 'lo', 'hi')],
          'results/fig6_values.csv', row.names = FALSE)

pdf(file.path(FIG, 'fig6_ceiling.pdf'), width = 7.5, height = 4.0)
layout(matrix(c(1, 2), 1, 2), widths = c(1.1, 1))
par(mar = c(4.2, 4.2, 2.6, 0.6), mgp = c(2.5, 0.8, 0))
plot(m$obs_lai, m$obs_gpp, las = 1, pch = PCH[m$sp],
     bg = adjustcolor(COLS[m$sp], .6), col = COLS[m$sp], cex = 1.25,
     xlab = 'Measured LAI', ylab = expression('Measured tree GPP (gC m'^-2*' yr'^-1*')'))
ln <- character(0)
for (s in SPP) {
    z <- m[m$sp == s & is.finite(m$obs_lai) & is.finite(m$obs_gpp), ]
    if (nrow(z) > 2) abline(lm(obs_gpp ~ obs_lai, z), col = COLS[s], lwd = 2)
    r <- cor(z$obs_gpp, z$obs_lai, use = 'complete.obs'); ci <- ci26(s, 'ceiling', r)
    ln <- c(ln, sprintf('%s r = %+.3f [%+.3f, %+.3f]', LBL[s], r, ci[1], ci[2]))
}
legend('topleft', legend = ln, bty = 'n', cex = 0.74, text.col = COLS[SPP])
lab('(a) Ceiling', 0.9)

## panel (b): rows bottom-to-top so Table 8's first row is drawn at the top
ROWS <- rev(T8ORDER)
get8 <- function(s, col) t8[[col]][t8$sp == s][match(ROWS, t8$attribute[t8$sp == s])]
rr <- sapply(SPP, function(s) get8(s, 'r'))
lo <- sapply(SPP, function(s) get8(s, 'lo')); hi <- sapply(SPP, function(s) get8(s, 'hi'))
PLAB <- sub('Leaf area index \\(max\\)', 'Leaf area index', ROWS)
PLAB <- sub('Soil nitrogen concentration', 'Soil N concentration', PLAB)
PLAB <- sub('Soil nitrogen stock', 'Soil N stock', PLAB)
par(mar = c(4.2, 8.2, 2.6, 1.0))
bp <- barplot(t(rr), beside = TRUE, horiz = TRUE, las = 1, xlim = c(-1, 1),
              names.arg = PLAB, col = adjustcolor(COLS[SPP], .75), border = NA,
              xlab = 'Correlation with GPP', cex.names = 0.75)
abline(v = 0, col = 'grey40')
arrows(t(lo), bp, t(hi), bp, angle = 90, code = 3, length = 0.02, col = 'grey20', lwd = 0.8)
## No legend: panel (a) carries the species colours, and the panel has no
## free space beside the bars. Plot counts (26, 21) are given in the caption.
lab('(b) Measured predictors', 0.9)
dev.off()

cat('written to ', FIG, ':\n  fig4_isolation.pdf\n  fig5_validation.pdf\n  fig6_ceiling.pdf\n', sep = '')
