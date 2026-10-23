## ---------------------------------------------------------------------------
## 48_age_raster_catchment.R
##
## Distribution of the three stand-age rasters over the Krycklan parcel of the
## catchment boundary (the largest polygon of krycklan_borders.shp), for the
## sentence on the age layer's accuracy in Section 2.4 (Forest structure):
## the product description reports it accurate up to about 80 years, with
## older stands increasingly underestimated. Read-only on the inputs; prints
## to the console and writes results/age_raster_catchment.csv.
##
## For each raster: valid cells (finite, >= 0; -1 is no data), maximum, 99th
## percentile, median and the share of cells above 80 years. Values are read
## as delivered, on each raster's own grid, without resampling.
##
##   Rscript 48_age_raster_catchment.R      (from UVa_work/, T7 mounted)
## ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(terra))
T <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
v <- disagg(vect(file.path(T, 'krycklan_borders/krycklan_borders.shp')))
k <- v[which.max(expanse(v, unit = 'ha'))]
out <- do.call(rbind, lapply(c('PINE_AGE_2025.tif', 'SPRUCE_AGE_2025.tif', 'Birch_ALS80.tif'), function(f) {
  r <- rast(file.path(T, 'Sweden_Forest_Age', f)); kk <- project(k, crs(r))
  x <- values(mask(crop(r, kk), kk), mat = FALSE); x <- x[is.finite(x) & x >= 0]
  data.frame(raster = f, cells = length(x), max = max(x), p99 = unname(quantile(x, .99)),
             median = median(x), share_over_80_pct = round(100 * mean(x > 80), 1))
}))
print(out, row.names = FALSE)
dir.create('results', showWarnings = FALSE)
write.csv(out, 'results/age_raster_catchment.csv', row.names = FALSE)
