## ---------------------------------------------------------------------------
## 43_tier5_support.R
##
## Recomputes, from the source data on the T7, values that main.tex quotes
## but that no earlier script wrote to a file. Read-only on
## the inputs; outputs to results/. Each section prints what it computes.
##
##   results/v01_object_check.csv     saved v01 sector object (11,120 cells):
##       maximum LAI over all months, initial foliage-to-stem ratio, range of
##       2017 GPP against range of initial stem biomass
##   results/catchment_areas.csv      boundary parcels; NMD class areas within
##       the Krycklan parcel; biomass statistics in the parcel and the box
##   results/relief_vpd_bound.csv     elevation range within the parcel (2 m
##       DEM), temperature difference at 6.5 C/km, saturation vapour pressure
##       spread, effect on the VPD modifier
##   results/radiation_summary.csv    MODIS MCD18A1, 2001-2020, masked to the
##       Krycklan parcel: July within-catchment SD, July range/mean, annual
##       total, November-February share
##   results/initial_cube_cells.csv   initial cube (out_3PG_netcdf/GPP.nc)
##
## The Krycklan parcel is the largest polygon of krycklan_borders.shp.
##
##   Rscript 43_tier5_support.R      (from UVa_work/, T7 mounted)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(terra) })
terraOptions(progress = 0)
B <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
OUT <- 'results'
sink(file.path(OUT, 'tier5_support_log.txt'), split = TRUE)
hr <- function(x) cat('\n==', x, '==\n')

v <- disagg(vect(file.path(B, 'krycklan_borders', 'krycklan_borders.shp')))
ha <- expanse(v, unit = 'ha')
kry <- v[which.max(ha)]

## ---- 1. saved v01 sector object -----------------------------------------------
hr('v01 sector object')
nm <- load(file.path(B, '3PG_10m', 'x731866.017196386_y7132899.71658545.RData'))
x <- get(nm[1]); rm(list = nm)
stand <- function(d, v) { s <- d[d$variable == v, ]
  if (!nrow(s)) return(NULL); tapply(s$value, s$date, sum, na.rm = TRUE) }
rows <- lapply(seq_along(x), function(i) {
  d <- x[[i]]; d0 <- d[d$date == min(d$date), ]
  i0 <- function(v) sum(d0$value[d0$variable == v], na.rm = TRUE)
  lai <- stand(d, 'lai'); g <- d[d$variable == 'gpp' & format(d$date, '%Y') == '2017', ]
  c(init_stem = i0('biom_stem'), init_foliage = i0('biom_foliage'),
    lai_max = if (is.null(lai)) NA else max(lai), gpp_2017 = sum(g$value, na.rm = TRUE))
})
O <- as.data.frame(do.call(rbind, rows)); rm(x); invisible(gc())
O$fs <- O$init_foliage / O$init_stem
q <- quantile(O$init_stem, c(.05, .95), na.rm = TRUE)
in90 <- O$init_stem >= q[1] & O$init_stem <= q[2]
V <- data.frame(
  n_cells = nrow(O),
  lai_max_any_month = max(O$lai_max, na.rm = TRUE),
  foliage_stem_median = median(O$fs, na.rm = TRUE), foliage_stem_mean = mean(O$fs, na.rm = TRUE),
  stem_min = min(O$init_stem), stem_max = max(O$init_stem),
  stem_ratio = max(O$init_stem) / min(O$init_stem),
  gpp2017_min = min(O$gpp_2017), gpp2017_max = max(O$gpp_2017),
  gpp2017_ratio = max(O$gpp_2017) / min(O$gpp_2017),
  stem_ratio_p05_p95 = q[[2]] / q[[1]],
  gpp_ratio_within_p05_p95 = max(O$gpp_2017[in90]) / min(O$gpp_2017[in90]),
  r_stem_gpp = cor(O$init_stem, O$gpp_2017))
print(t(V))
write.csv(V, file.path(OUT, 'v01_object_check.csv'), row.names = FALSE)

## ---- 2. areas, land cover, biomass ----------------------------------------------
hr('areas')
A <- data.frame(item = c(sprintf('parcel %d', seq_along(ha)), 'parcels total'),
                value = c(ha, sum(ha)), unit = 'ha (geodesic)')
nmd <- crop(rast(file.path(B, 'nationalndcoverdata', 'nmd2018bas_ogeneraliserad_v1_0.tif')), kry)
nmd <- mask(nmd, kry)
tab <- freq(nmd); tab$ha <- tab$count * prod(res(nmd)) / 1e4
grp <- function(codes) sum(tab$ha[tab$value %in% codes])
A <- rbind(A, data.frame(
  item = c('NMD pine (111,121)', 'NMD spruce (112,122)', 'NMD mixed (113,114,123,124)',
           'NMD broadleaf (115-117,125-127)', 'NMD forest total',
           'NMD temporarily not forest (118,128)', 'NMD all classes in parcel',
           'NMD pine, non-wetland only (111)', 'NMD spruce, non-wetland only (112)',
           'NMD mixed, non-wetland only (113,114)', 'NMD broadleaf, non-wetland only (115-117)',
           'NMD forest total, non-wetland only', 'NMD temporarily not forest, non-wetland only (118)'),
  value = c(grp(c(111, 121)), grp(c(112, 122)), grp(c(113, 114, 123, 124)),
            grp(c(115:117, 125:127)), grp(c(111:117, 121:127)),
            grp(c(118, 128)), sum(tab$ha),
            grp(111), grp(112), grp(c(113, 114)), grp(115:117), grp(111:117), grp(118)),
  unit = 'ha (NMD 10 m cells)'))
bio <- rast(file.path(B, 'biomass', 'BIO.tif'))
bk <- mask(crop(bio, kry), kry); bk[bk < 0] <- NA
box <- ext(726580, 737150, 7126990, 7138810)             # assembled cube extent
bb <- crop(bio, box); bb[bb < 0] <- NA
st <- function(r) { m <- global(r, c('min', 'max', 'mean'), na.rm = TRUE); unlist(m) }
sk <- st(bk); sb <- st(bb)
A <- rbind(A, data.frame(item = c('biomass parcel min', 'biomass parcel max', 'biomass parcel mean',
                                  'biomass cube box min', 'biomass cube box max', 'biomass cube box mean'),
                         value = c(sk, sb), unit = 't/ha'))
print(A, row.names = FALSE)
write.csv(A, file.path(OUT, 'catchment_areas.csv'), row.names = FALSE)

## ---- 3. relief and the VPD bound ------------------------------------------------------
hr('relief and VPD')
dem <- mask(crop(rast(file.path(B, 'Sweden_dem', 'dem.tif')), kry), kry)
er <- unlist(global(dem, c('min', 'max'), na.rm = TRUE))
dT <- diff(er) / 1000 * 6.5
svp <- function(t) 0.6108 * exp(17.27 * t / (t + 237.3))            # Magnus, kPa
R <- data.frame(elev_min = er[[1]], elev_max = er[[2]], relief = diff(er), dT = dT)
for (t0 in c(10, 15, 20)) R[[sprintf('svp_spread_%d', t0)]] <- svp(t0 + dT / 2) / svp(t0 - dT / 2) - 1
## f_vpd = exp(-CoeffCond * VPD): a fractional VPD change s moves ln f_vpd by s
fv <- 0.800; s <- R$svp_spread_15
R$f_vpd_base <- fv; R$f_vpd_shifted <- exp(log(fv) * (1 + s)); R$gpp_change <- 1 - R$f_vpd_shifted / fv
print(t(R))
write.csv(R, file.path(OUT, 'relief_vpd_bound.csv'), row.names = FALSE)

## ---- 4. radiation ------------------------------------------------------------------------
hr('radiation 2001-2020, Krycklan parcel')
md <- file.path(B, 'MODIS_DSR_Kricklan_solar_radiation')
fs <- list.files(md, pattern = '^DSR_\\d{4}_\\d{2}\\.tif$', full.names = TRUE)
tag <- sub('^DSR_', '', sub('\\.tif$', '', basename(fs)))
yr <- as.integer(substr(tag, 1, 4)); mo <- as.integer(substr(tag, 6, 7))
keep <- yr >= 2001 & yr <= 2020
M <- do.call(rbind, lapply(which(keep), function(i) {
  r <- rast(fs[i])[[1]]; kp <- project(kry, crs(r))          # MODIS grid CRS
  r <- mask(crop(r, kp, snap = 'out'), kp, touches = TRUE)
  g <- unlist(global(r, c('mean', 'sd', 'min', 'max'), na.rm = TRUE))
  data.frame(year = yr[i], month = mo[i], mean = g[[1]], sd = g[[2]], min = g[[3]], max = g[[4]])
}))
dim_m <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
M$total <- M$mean * dim_m[M$month]
jul <- M[M$month == 7, ]
ann <- aggregate(total ~ year, M, sum); ncomp <- aggregate(month ~ year, M, length)
ann <- ann[ncomp$month == 12, ]
clim <- aggregate(total ~ month, M, mean)
Rd <- data.frame(
  july_sd_min = min(jul$sd), july_sd_max = max(jul$sd),
  july_range_over_mean_min = min((jul$max - jul$min) / jul$mean),
  july_range_over_mean_max = max((jul$max - jul$min) / jul$mean),
  july_range_over_mean_median = median((jul$max - jul$min) / jul$mean),
  annual_MJ_mean = mean(ann$total), annual_kWh_mean = mean(ann$total) / 3.6,
  years_complete = nrow(ann),
  novfeb_share = sum(clim$total[clim$month %in% c(11, 12, 1, 2)]) / sum(clim$total),
  n_cells_july = NA)
print(t(Rd))
write.csv(Rd, file.path(OUT, 'radiation_summary.csv'), row.names = FALSE)
write.csv(M, file.path(OUT, 'radiation_monthly.csv'), row.names = FALSE)

## ---- 5. initial cube -----------------------------------------------------------------------
hr('initial cube')
g <- rast(file.path(B, 'out_3PG_netcdf', 'GPP.nc'))
cnt <- function(k) global(!is.na(g[[k]]), 'sum')[[1]]
I <- data.frame(nrow = nrow(g), ncol = ncol(g), ncell = ncell(g), nlyr = nlyr(g),
                nonNA_band1 = cnt(1), nonNA_band199 = cnt(199), nonNA_band240 = cnt(240))
I$share_band199 <- I$nonNA_band199 / I$ncell
print(I, row.names = FALSE)
write.csv(I, file.path(OUT, 'initial_cube_cells.csv'), row.names = FALSE)
sink()
