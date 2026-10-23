## ---------------------------------------------------------------------------
## 44_input_variation_cube_plots.R
##
## Spatial variation of the coarse drivers as the cube received them and at the
## 47 analysed plots, for the Methods (Sections 2.4, 2.5) and Section 4.1.
## Re-runs the steps of iteration_by_sectors/polygon_grids_v03.r with the
## functions of devel/R/internals.r. Read-only on the inputs; outputs to
## results/.
##
##   results/input_variation_sectors.csv   per sector of the 10 x 10 grid:
##       reprojected 4 km PTHBV cells in the sector window and distinct 10 m
##       values after disaggregation (January 2017 precipitation)
##   results/input_variation_monthly.csv   per month 2001-2020: range, SD and
##       mean of temperature and precipitation over the reprojected PTHBV cells
##       intersecting the Krycklan parcel and over those holding the plots, and
##       of radiation on the 10 m grid (v03 line 164) within the parcel and at
##       the plot cells
##   results/input_variation_summary.csv   July radiation SD and range/mean,
##       median monthly temperature and precipitation ranges
##
## The Krycklan parcel is part 3 of krycklan_borders.shp, the polygon
## polygon_sector keeps (internals.r:467). The 10 m climate values are
## bilinear combinations of the reprojected cells, so the cell ranges bound
## them. The plot-level simulations do not use any of this: they read one
## monthly climate and radiation series for every plot (16_preflight_corrected_V1.R:46).
##
##   Rscript 44_input_variation_cube_plots.R      (from UVa_work/, T7 mounted)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra); library(sf); library(lubridate); library(stringr); library(magrittr)
})
terraOptions(progress = 0)
sf_use_s2(FALSE)
B   <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
OUT <- 'results'
sink(file.path(OUT, 'input_variation_log.txt'), split = TRUE)
hr <- function(x) cat('\n==', x, '==\n')

ex <- parse('../R/internals.r')
need <- c('polygon_sector', 'crop_raster_fast', 'agg_and_crop..', 'align_to_ref',
          'monthly_clim_roi_from_daily_raster')
for (e in ex) if (is.call(e) && identical(e[[1]], as.name('<-')) &&
                  as.character(e[[2]]) %in% need) eval(e, globalenv())

bord <- vect(file.path(B, 'krycklan_borders', 'krycklan_borders.shp'))
ps   <- suppressWarnings(polygon_sector(bord, 10, 1, 200))
sel  <- ps$kryck_border_sf                                  # part 3
grid <- st_union(ps$grid..)
kb_clim <- suppressWarnings(polygon_sector(bord, 1, 1, 4000))[[1]]
cat('parcel', round(as.numeric(st_area(sel)) / 1e4, 1), 'ha\n')

## plots
xl <- '../validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
co <- as.data.frame(readxl::read_excel(xl, sheet = 'plot coordinates'))
idc <- grep('plot|id|stand', names(co), ignore.case = TRUE, value = TRUE)[1]
co <- data.frame(Plot_ID = suppressWarnings(as.numeric(co[[idc]])),
                 x = suppressWarnings(as.numeric(co$E_Sweref99TM)),
                 y = suppressWarnings(as.numeric(co$N_Sweref99TM)))
ids <- unique(read.csv('validation_out/preflight_results.csv')$Plot_ID)
p47 <- co[co$Plot_ID %in% ids & is.finite(co$x), ]
pts <- vect(p47, geom = c('x', 'y'), crs = crs(bord))
cat('analysed plots:', nrow(p47), '\n')

## ---- 1. climate per sector (January 2017 precipitation) ----------------------
hr('PTHBV per sector')
prf <- file.path(B, 'climate_SW', 'SMHI_pthbv_pr_1961_2025_daily.nc')
taf <- file.path(B, 'climate_SW', 'SMHI_pthbv_tas_1961_2025_daily.nc')
m1  <- monthly_clim_roi_from_daily_raster(prf, kb_clim, fun = 'sum',
                                          t1 = '2017-01-01', t2 = '2017-02-01')
bio <- rast(file.path(B, 'biomass', 'BIO.tif'))
S <- do.call(rbind, lapply(1:100, function(s) tryCatch({
  kb  <- suppressWarnings(polygon_sector(bord, 10, s, 200))[[1]]
  ref <- crop_raster_fast(bio, kb)[[1]]
  a   <- crop_raster_fast(m1, kb)
  b   <- suppressMessages(agg_and_crop..(a, 10, kb))
  v   <- values(align_to_ref(b, ref, method = 'near'))[, 1]
  data.frame(sector = s, coarse_cells = ncell(a), distinct_10m = length(unique(v[is.finite(v)])))
}, error = function(e) NULL)))
S$interpolated <- S$distinct_10m > S$coarse_cells
print(table(coarse_cells = S$coarse_cells, interpolated = S$interpolated))
write.csv(S, file.path(OUT, 'input_variation_sectors.csv'), row.names = FALSE)

## ---- 2. monthly variation 2001-2020 ------------------------------------------
hr('monthly variation 2001-2020')
tas <- monthly_clim_roi_from_daily_raster(taf, kb_clim, fun = 'mean', t1 = '2001-01-01', t2 = '2020-12-31')
prc <- monthly_clim_roi_from_daily_raster(prf, kb_clim, fun = 'sum',  t1 = '2001-01-01', t2 = '2020-12-31')
cells_in <- function(r) unique(cells(r[[1]], vect(st_transform(sel, crs(r))), touches = TRUE)[, 2])
cells_pt <- function(r) unique(cellFromXY(r[[1]], crds(project(pts, crs(r)))))
st <- function(M) cbind(range = apply(M, 2, function(v) diff(range(v, na.rm = TRUE))),
                        sd = apply(M, 2, sd, na.rm = TRUE), mean = colMeans(M, na.rm = TRUE))
T_par <- st(as.matrix(tas[cells_in(tas)])); T_pts <- st(as.matrix(tas[cells_pt(tas)]))
P_par <- st(as.matrix(prc[cells_in(prc)])); P_pts <- st(as.matrix(prc[cells_pt(prc)]))
cat('PTHBV cells: parcel', length(cells_in(tas)), '| holding plots', length(cells_pt(tas)), '\n')

md <- file.path(B, 'MODIS_DSR_Kricklan_solar_radiation')
fs <- list.files(md, pattern = '^DSR_(200[1-9]|201[0-9]|2020)_\\d{2}\\.tif$', full.names = TRUE)
ref  <- crop(bio, vect(grid), snap = 'out')
inp  <- !is.na(values(rasterize(vect(sel), ref))[, 1])
pc   <- cellFromXY(ref, crds(pts))
R <- do.call(rbind, lapply(fs, function(f) {
  v <- values(align_to_ref(rast(f)[[1]], ref))[, 1]           # as v03 line 164
  a <- v[inp]; b <- v[pc]
  data.frame(rad_par_range = diff(range(a, na.rm = TRUE)), rad_par_sd = sd(a, na.rm = TRUE),
             rad_par_mean = mean(a, na.rm = TRUE),
             rad_pts_range = diff(range(b, na.rm = TRUE)), rad_pts_sd = sd(b, na.rm = TRUE),
             rad_pts_mean = mean(b, na.rm = TRUE))
}))
tag <- sub('^DSR_', '', sub('\\.tif$', '', basename(fs)))
M <- data.frame(year = as.integer(substr(tag, 1, 4)), month = as.integer(substr(tag, 6, 7)),
                tas_par_range = T_par[, 'range'], tas_pts_range = T_pts[, 'range'],
                pr_par_range = P_par[, 'range'], pr_par_mean = P_par[, 'mean'],
                pr_pts_range = P_pts[, 'range'], pr_pts_mean = P_pts[, 'mean'], R)
write.csv(M, file.path(OUT, 'input_variation_monthly.csv'), row.names = FALSE)

jul <- M[M$month == 7, ]
Z <- data.frame(
  sectors_replicated = sum(!S$interpolated), sectors_interpolated = sum(S$interpolated),
  tas_range_median_par = median(M$tas_par_range), tas_range_max_par = max(M$tas_par_range),
  tas_range_median_pts = median(M$tas_pts_range),
  pr_range_pct_median_par = 100 * median(M$pr_par_range / M$pr_par_mean),
  pr_range_pct_median_pts = 100 * median(M$pr_pts_range / M$pr_pts_mean),
  july_sd_min_par = min(jul$rad_par_sd), july_sd_max_par = max(jul$rad_par_sd),
  july_range_pct_min_par = 100 * min(jul$rad_par_range / jul$rad_par_mean),
  july_range_pct_max_par = 100 * max(jul$rad_par_range / jul$rad_par_mean),
  july_range_pct_median_par = 100 * median(jul$rad_par_range / jul$rad_par_mean),
  july_sd_min_pts = min(jul$rad_pts_sd), july_sd_max_pts = max(jul$rad_pts_sd),
  july_range_pct_min_pts = 100 * min(jul$rad_pts_range / jul$rad_pts_mean),
  july_range_pct_max_pts = 100 * max(jul$rad_pts_range / jul$rad_pts_mean),
  july_range_pct_median_pts = 100 * median(jul$rad_pts_range / jul$rad_pts_mean))
print(t(round(Z, 3)))
write.csv(Z, file.path(OUT, 'input_variation_summary.csv'), row.names = FALSE)
sink()
