## ---------------------------------------------------------------------------
## 29_fig1_site_v2.R   -- run on the server, from ~/datos/3PG
##
## Fixes over v1, all three visible in the first draft:
##   (a) plain grey rectangle -> real Sweden outline from the maps package,
##       which is installed here.
##   (b) drew kryck_elev.tif, which turned out to be a 1.3 x 1.5 km test tile,
##       so only two plots fell inside it. Now crops the national 2 m DEM to
##       the catchment and aggregates to 10 m, showing all 47 plots.
##   (c) plotted the NMD raster as if continuous, giving no legend and no
##       classes. Now reclassified to the model species using the NMD base
##       layer codes, with a discrete legend.
##   (d) drew Birch_ALS80.tif, first alphabetically. Now PINE_AGE_2025.tif.
##       The age product is species-resolved, so the caption must say which
##       species is shown.
##
##   Rscript 29_fig1_site_v2.R
## ---------------------------------------------------------------------------

library(terra)

BASE <- '.'
OUT  <- 'fig'
dir.create(OUT, showWarnings = FALSE)

F_BORDER <- file.path(BASE, 'krycklan_borders', 'krycklan_borders.shp')
F_DEM    <- file.path(BASE, 'Sweden_dem', 'dem.tif')          # national, 2 m
F_LAND   <- file.path(BASE, 'nationalndcoverdata',
                      'nmd2018bas_ogeneraliserad_v1_0.tif')
F_AGE    <- file.path(BASE, 'Sweden_Forest_Age', 'PINE_AGE_2025.tif')
PLOTS    <- '~/Documents/SLU/3PG/devel/UVa_work/validation_out/itc_vs_measured.csv'
GPPCSV   <- '~/Documents/SLU/3PG/devel/UVa_work/validation_out/preflight_results.csv'

## NMD base layer forest classes, per the class list supplied by SLU.
## 12x are the same classes on wetland and are folded into their 11x
## equivalent here; the distinction is not what this panel is about.
NMD <- data.frame(
  from = c(111,112,113,114,115,116,117,118,
           121,122,123,124,125,126,127,128),
  to   = c(  1,  2,  3,  4,  5,  5,  5,  6,
             1,  2,  3,  4,  5,  5,  5,  6))
LEV <- c('Pine', 'Spruce', 'Mixed conifer', 'Conifer-broadleaf',
         'Broadleaf', 'Temporarily not forest')
COL <- c('#e8c15a', '#2e6f4e', '#7fa86a', '#a8c98a', '#c76b3a', '#d9d9d9')

hr <- function(x) cat('\n', strrep('-', 66), '\n', x, '\n', sep = '')
lab <- function(s) mtext(s, side = 3, adj = 0, line = 0.5, font = 2)

## ---- boundary: keep the Krycklan part only ---------------------------------
v <- disagg(vect(F_BORDER))
kry <- v[which.max(expanse(v, unit = 'ha'))]
cat('Krycklan parcel: ', round(expanse(kry, unit = 'ha')), ' ha\n', sep = '')
E <- ext(kry)

## ---- plots -----------------------------------------------------------------
p <- read.csv(path.expand(PLOTS))
pts <- vect(p[is.finite(p$E) & is.finite(p$N), ], geom = c('E','N'), crs = crs(kry))
g <- read.csv(path.expand(GPPCSV)); g <- g[g$variant == 'observed/obs', ]
ag <- aggregate(obs_gpp ~ Plot_ID, g, mean)
pts$gpp <- ag$obs_gpp[match(pts$Plot_ID, ag$Plot_ID)]
cat('plots: ', nrow(pts), ' | inside the catchment: ',
    sum(!is.na(extract(rasterize(kry, rast(kry, res = 10)), pts)[, 2])), '\n', sep = '')

crop_kry <- function(f, agg = NULL, method = 'near') {
    r <- rast(f)
    k <- if (same.crs(kry, r)) kry else project(kry, crs(r))
    r <- crop(r, ext(k))
    if (!is.null(agg) && agg > 1) r <- aggregate(r, agg, fun = 'mean', na.rm = TRUE)
    r <- mask(r, k)
    if (!same.crs(r, kry)) r <- project(r, crs(kry), method = method)
    r
}

hr('BUILDING PANELS')

pdf(file.path(OUT, 'fig1_site.pdf'), width = 9.5, height = 9)
op <- par(mfrow = c(2, 2), mar = c(3.4, 3.6, 2.4, 1.0), mgp = c(2.2, 0.7, 0))

## (a) location, real outline -------------------------------------------------
cen <- colMeans(crds(centroids(project(kry, 'EPSG:4326'))))
if (requireNamespace('maps', quietly = TRUE)) {
    maps::map('world', regions = c('Sweden','Norway','Finland'),
              xlim = c(4, 32), ylim = c(54.5, 70), fill = TRUE,
              col = c('grey88','grey96','grey96'), border = 'grey55',
              mar = c(3.4, 3.6, 2.4, 1.0))
    box(); axis(1); axis(2, las = 1)
    title(xlab = 'Longitude', ylab = 'Latitude')
} else {
    plot(NA, xlim = c(4, 32), ylim = c(54.5, 70), xlab = 'Longitude',
         ylab = 'Latitude', las = 1)
}
points(cen[1], cen[2], pch = 21, bg = 'firebrick', cex = 2.0, lwd = 1.2)
text(cen[1] + 0.9, cen[2], 'Krycklan', adj = 0, font = 2, cex = 0.95)
lab('(a)')
cat('  (a) done\n')

## (b) elevation over the whole catchment, with all plots ---------------------
e <- crop_kry(F_DEM, agg = 5, method = 'bilinear')   # 2 m -> 10 m
plot(e, col = hcl.colors(60, 'Terrain2'), las = 1,
     plg = list(title = 'm a.s.l.', cex = 0.7), mar = c(3.4, 3.6, 2.4, 4.2))
lines(kry, lwd = 1.5)
sp  <- as.character(pts$sp)
pch <- ifelse(grepl('Pin', sp, ignore.case=TRUE), 21, ifelse(grepl('Spruce|Picea', sp, ignore.case=TRUE), 24, 22))
gp  <- pts$gpp
cx  <- if (any(is.finite(gp)))
    0.7 + 1.5 * (gp - min(gp, na.rm = TRUE)) / diff(range(gp, na.rm = TRUE)) else 1.1
points(pts, pch = pch, bg = 'white', col = 'black', cex = cx, lwd = 1.1)
legend('topleft', inset = c(0.04, 0.02), c('Pinus sylvestris','Picea abies'),
       pch = c(21,24), pt.bg = 'white', bty = 'o', bg = 'white', box.col = 'grey70', cex = 0.7, pt.cex = 1.1)
lab('(b)')
cat('  (b) done, ', nrow(pts), ' plots\n', sep = '')

## (c) species, reclassified with a discrete legend ---------------------------
s <- crop_kry(F_LAND)
s <- classify(s, as.matrix(NMD), others = NA)
levels(s) <- data.frame(value = 1:6, class = LEV)
plot(s, col = COL, las = 1, plg = list(cex = 0.62), mar = c(3.4, 3.6, 2.4, 7.5))
lines(kry, lwd = 1.5)
lab('(c)')
cat('  (c) done\n')

## (d) stand age, pine ---------------------------------------------------------
a <- crop_kry(F_AGE, method = 'bilinear')
plot(a, col = hcl.colors(50, 'YlGnBu', rev = TRUE), las = 1,
     plg = list(title = 'years', cex = 0.7), mar = c(3.4, 3.6, 2.4, 4.2))
lines(kry, lwd = 1.5)
lab('(d)')
cat('  (d) done\n')

par(op); dev.off()
cat('\nwritten: ', file.path(OUT, 'fig1_site.pdf'), '\n', sep = '')

cat('\nCAPTION MUST NOW READ\n')
cat('  (a) location of the Krycklan catchment in northern Sweden\n')
cat('  (b) elevation and the 47 validation plots, symbol shape by dominant\n')
cat('      species and size by measured tree GPP\n')
cat('  (c) dominant species from the national land cover product, NMD 2018\n')
cat('      base layer; wetland variants folded into their forest class\n')
cat('  (d) stand age for Pinus sylvestris from the national age product\n')
cat('\nNote 47, not 57, and note that (d) is species-specific.\n')
