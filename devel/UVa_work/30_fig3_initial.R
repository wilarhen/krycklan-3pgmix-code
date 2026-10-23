## ---------------------------------------------------------------------------
## 30_fig3_initial.R   -- run on the server, from ~/datos/3PG
##
## Figure 3: the data cube before the configuration hazards were corrected.
##
##   (a) simulated against measured tree GPP, original configuration
##   (b) simulated GPP for July 2017, showing the spatial field
##   (c) residuals against measured leaf area index
##
## Sources
##   old_outputs_pre_v03/    100 tiles, 240 bands, GPP in t DM/ha/month.
##                           Origin (4, -2): these tiles sit 4 m east and 2 m
##                           south of the structural grid, which is the grid
##                           offset itself, visible in the tile coordinates.
##   validation_pairs.csv    plot-level pairs from the same run
##                           (r = -0.121, R2 = 0.0145)
##
##   Rscript 30_fig3_initial.R
## ---------------------------------------------------------------------------

library(terra)

OLD   <- 'old_outputs_pre_v03'
PAIRS <- '~/Documents/SLU/3PG/devel/UVa_work/validation_out/validation_pairs.csv'
BORD  <- 'krycklan_borders/krycklan_borders.shp'
OUT   <- 'fig'
dir.create(OUT, showWarnings = FALSE)

GDM_TO_GC <- 0.47 * 100        # t DM/ha -> gC/m2
TARGET    <- '2017-07'         # band 199 of 240

hr <- function(x) cat('\n', strrep('-', 66), '\n', x, '\n', sep = '')
lab <- function(s) mtext(s, side = 3, adj = 0, line = 0.5, font = 2)

## ---- plot-level pairs -------------------------------------------------------
hr('PLOT PAIRS, ORIGINAL CONFIGURATION')
d <- read.csv(path.expand(PAIRS))
d <- d[is.finite(d$sim) & is.finite(d$obs_tree), ]

## annual sums per plot-year are already in sim / obs_tree; aggregate to plots
a <- aggregate(cbind(sim, obs_tree, LAImax) ~ Plot_ID, d, mean, na.rm = TRUE)
r_py <- cor(d$sim, d$obs_tree)
r_pl <- cor(a$sim, a$obs_tree)
cat(sprintf('plot-years %d  r %.3f  R2 %.4f\n', nrow(d), r_py, r_py^2))
cat(sprintf('plots      %d  r %.3f  R2 %.4f\n', nrow(a), r_pl, r_pl^2))
a$resid <- a$sim - a$obs_tree

## ---- the July 2017 field ----------------------------------------------------
hr('MOSAIC, JULY 2017')
fs <- list.files(OLD, pattern = '\\.tif$', full.names = TRUE)
cat('tiles: ', length(fs), '\n', sep = '')
b1 <- rast(fs[1])
idx <- which(substr(names(b1), 1, 7) == TARGET)
if (!length(idx)) stop('no band for ', TARGET)
cat('band ', idx, ' = ', names(b1)[idx], '\n', sep = '')
cat('origin of the original tiles: (',
    paste(origin(b1), collapse = ', '), ')  <- the grid offset\n', sep = '')

m <- mosaic(sprc(lapply(fs, function(f) rast(f, lyrs = idx))), fun = 'mean')
m <- m * GDM_TO_GC                                   # gC/m2/month
bord <- vect(BORD)
kry <- disagg(bord)[which.max(expanse(disagg(bord), unit = 'ha'))]
if (!same.crs(kry, m)) kry <- project(kry, crs(m))
m <- mask(crop(m, kry), kry)
print(global(m, c('min', 'mean', 'max'), na.rm = TRUE))

## ---- draw -------------------------------------------------------------------
hr('DRAWING')
pdf(file.path(OUT, 'fig3_initial.pdf'), width = 11, height = 3.9)
op <- par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.4, 1.2), mgp = c(2.5, 0.8, 0))

## (a) simulated against measured ---------------------------------------------
lim <- range(c(a$sim, a$obs_tree), na.rm = TRUE)
lim <- c(0, max(lim) * 1.05)
plot(a$obs_tree, a$sim, xlim = lim, ylim = lim, las = 1, pch = 21,
     bg = adjustcolor('grey30', .5), col = 'grey20', cex = 1.2,
     xlab = expression('Measured tree GPP (gC m'^-2*' yr'^-1*')'),
     ylab = expression('Simulated GPP (gC m'^-2*' yr'^-1*')'))
abline(0, 1, lty = 2, col = 'grey45')
abline(lm(sim ~ obs_tree, a), col = 'firebrick', lwd = 2)
legend('topleft', bty = 'n', cex = 0.85, legend = c(
    sprintf('n = %d plots', nrow(a)),
    sprintf('r = %+.3f', r_pl),
    sprintf('R%s = %.3f', '\u00b2', r_pl^2)))
lab('(a)')

## (b) the field ---------------------------------------------------------------
plot(m, col = hcl.colors(60, 'Batlow'), las = 1,
     plg = list(title = expression('gC m'^-2*' mo'^-1), cex = 0.75),
     mar = c(4.2, 4.4, 2.4, 5.0))
lines(kry, lwd = 1.4)
lab('(b) July 2017')

## (c) residual against measured leaf area -------------------------------------
plot(a$LAImax, a$resid, las = 1, pch = 21, bg = adjustcolor('grey30', .5),
     col = 'grey20', cex = 1.2, xlab = 'Measured LAI',
     ylab = expression('Simulated - measured (gC m'^-2*' yr'^-1*')'))
abline(h = 0, lty = 2, col = 'grey45')
if (sum(is.finite(a$LAImax) & is.finite(a$resid)) > 3)
    abline(lm(resid ~ LAImax, a), col = 'firebrick', lwd = 2)
rr <- cor(a$LAImax, a$resid, use = 'complete.obs')
legend('topright', bty = 'n', cex = 0.85,
       legend = sprintf('r = %+.3f', rr))
lab('(c)')

par(op); dev.off()
cat('\nwritten: ', file.path(OUT, 'fig3_initial.pdf'), '\n', sep = '')

cat('\nNUMBERS FOR THE CAPTION AND TEXT\n')
cat(sprintf('  r on plots      %+.3f   R2 %.4f\n', r_pl, r_pl^2))
cat(sprintf('  r on plot-years %+.3f   R2 %.4f  (the value in the draft)\n',
            r_py, r_py^2))
cat(sprintf('  residual vs measured LAI  r = %+.3f\n', rr))
cat('  tile origin (', paste(origin(b1), collapse = ', '),
    ') against (0, 0) for the corrected cube\n', sep = '')
cat('\nThe correlation is NEGATIVE, not merely weak. Worth saying so:\n')
cat('the original cube ordered the plots slightly against measurement.\n')
