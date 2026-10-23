## ---------------------------------------------------------------------------
## 56_cube_qc_fail_check.R
##
## Read-only check of the two FAILs in the corrected cube's qc_summary.csv (LAI maximum 14.56,
## height maximum 250.5 m). Rebuilds each cell's inputs on the cube grid from the source rasters
## as polygon_grids_v03_server.r does (NMD species slots via spp_codes_from_land.csv, stand age,
## BIO), then reports per-cell LAI and height maxima against them: the AGB bins of the tall
## mixed cells, GPP in the tall cells, the cells at the validation plots and the composition of
## the LAI exceedance. Writes results/cube_qc_fail_check.csv (quantity, value).
## First run 18 Oct 2026.
##
##   Rscript 56_cube_qc_fail_check.R      (local; reads the checksum-verified copy of the cube on the T7 disk; ~15 min)
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(terra); library(readxl)})
OUT <- if (nzchar(Sys.getenv('OUT'))) Sys.getenv('OUT') else 'results'
B  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
C  <- file.path(B, 'outputs_3PG_corrected_cube')
XL <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
LAI_HI <- 6.75                                   # measured LAImax at the plots, 0.28-6.75
res <- list(); add <- function(q, v) res[[length(res) + 1]] <<- data.frame(quantity = q, value = v)

## ---- inputs on the cube grid ----
g <- rast(file.path(C, 'gpp.tif'))[[1]]
on <- function(r, m = 'near') { r <- crop(r, project(as.polygons(ext(g), crs = crs(g)), crs(r)), snap = 'out')
  if (!compareGeom(r, g, stopOnError = FALSE)) r <- resample(project(r, crs(g), method = m), g, method = m); r }
nmd <- values(on(rast(file.path(B, 'nationalndcoverdata', 'nmd2018bas_ogeneraliserad_v1_0.tif'))), mat = FALSE)
A <- file.path(B, 'Sweden_Forest_Age')
ap <- values(on(rast(file.path(A, 'PINE_AGE_2025.tif'))), mat = FALSE)
as <- values(on(rast(file.path(A, 'SPRUCE_AGE_2025.tif'))), mat = FALSE)
bio <- values(on(rast(file.path(B, 'biomass', 'BIO.tif')), 'bilinear'), mat = FALSE)
lk <- read.csv(file.path(B, 'nationalndcoverdata', 'spp_codes_from_land.csv'))
sl <- sapply(1:3, function(k) { x <- trimws(lk[[paste0('Layer_', k)]])[match(nmd, lk$code)]; x[x %in% c('NA', '')] <- NA; x })
code <- c('Pinus sylvestris' = 'P', 'Picea abies' = 'S', 'Betula' = 'B')
combo <- apply(sl, 1, function(r) { r <- r[!is.na(r)]; if (!length(r)) NA else paste(sort(code[r]), collapse = '') })
cm <- values(rast(file.path(C, 'qc_lai_cellmax.tif')), mat = FALSE)
s <- data.frame(cell = seq_along(nmd), nmd, combo, bio, ap, as)[is.finite(cm), ]
xy <- xyFromCell(g, s$cell); s$x <- xy[, 1]; s$y <- xy[, 2]
add('cells_simulated', nrow(s))

## ---- per-cell LAI and height ----
lv <- as.matrix(rast(file.path(C, 'lai.tif'))[s$cell])
s$lai_max <- apply(lv, 1, max, na.rm = TRUE); s$lai_1 <- lv[, 1]; s$lai_240 <- lv[, 240]
s$lai_bmax <- apply(lv, 1, which.max); rm(lv); invisible(gc())
hv <- as.matrix(rast(file.path(C, 'height.tif'))[s$cell])
s$h_max <- apply(hv, 1, max, na.rm = TRUE); s$h_240 <- hv[, 240]; rm(hv); invisible(gc())
ann <- rast(file.path(C, 'qc_annual_gpp.tif'))                     # gC m-2 yr-1, 2001-2020
av <- as.matrix(ann[s$cell]); s$gpp_mean <- rowMeans(av, na.rm = TRUE); s$gpp_2020 <- av[, nlyr(ann)]; rm(av)

## ---- height ----
add('height_last_median_m', median(s$h_240, na.rm = TRUE))
add('height_last_p99_m', quantile(s$h_240, .99, na.rm = TRUE))
for (t in c(27.16, 40, 50, 60, 100)) add(sprintf('height_cellmax_n_gt_%g', t), sum(s$h_max > t))
add('height_cellmax_max_m', max(s$h_max))
for (k in c('B', 'BPS', 'P', 'PS', 'S')) add(sprintf('height_cellmax_n_gt_40_%s', k), sum(s$h_max > 40 & s$combo %in% k))
for (k in c('P', 'S')) add(sprintf('height_cellmax_max_%s_only_m', k), max(s$h_max[s$combo %in% k]))
t60 <- s$h_max > 60
add('gt60_all_mixed_class114', all(s$nmd[t60] %in% c(114, 124)))
add('gt60_agb_values', paste(sort(unique(s$bio[t60])), collapse = ';'))
m <- s$combo %in% 'BPS'
add('mixed_gt40_agb_values', paste(sort(unique(s$bio[m & s$h_max > 40])), collapse = ';'))
add('mixed_agb189_n', sum(m & s$bio == 189)); add('mixed_agb189_n_gt40', sum(m & s$bio == 189 & s$h_max > 40))
for (r in list(c(181, 188), c(190, 200))) { i <- m & s$bio >= r[1] & s$bio <= r[2]
  add(sprintf('mixed_agb_%d_%d_n', r[1], r[2]), sum(i))
  add(sprintf('mixed_agb_%d_%d_height_cellmax_max_m', r[1], r[2]), max(s$h_max[i]))
  add(sprintf('mixed_agb_%d_%d_noninteger_agb', r[1], r[2]), sum(s$bio[i] != round(s$bio[i]))) }
bv <- values(crop(rast(file.path(B, 'biomass', 'BIO.tif')), ext(g)), mat = FALSE)
for (v in 188:190) add(sprintf('BIO_cells_in_cube_extent_eq_%d', v), sum(bv == v, na.rm = TRUE)); rm(bv)

## ---- GPP in the tall cells (gC m-2 yr-1, mean of 2001-2020 unless stated) ----
add('gpp_mean_all_simulated', mean(s$gpp_mean, na.rm = TRUE))
add('gpp_median_all_simulated', median(s$gpp_mean, na.rm = TRUE))
add('gpp_mean_mixed_all', mean(s$gpp_mean[m], na.rm = TRUE))
i189 <- m & s$bio == 189
add('gpp_mean_mixed_agb189_gt40', mean(s$gpp_mean[i189 & s$h_max > 40]))
add('gpp_median_mixed_agb189_gt40', median(s$gpp_mean[i189 & s$h_max > 40]))
add('gpp_median_mixed_agb189_le40', median(s$gpp_mean[i189 & s$h_max <= 40]))
add('gpp2020_median_mixed_agb189_gt40', median(s$gpp_2020[i189 & s$h_max > 40]))
add('gpp2020_median_mixed_agb189_le40', median(s$gpp_2020[i189 & s$h_max <= 40]))
add('gpp_mean_gt60', mean(s$gpp_mean[t60])); add('gpp_median_gt60', median(s$gpp_mean[t60]))

## ---- validation plots ----
cc <- as.data.frame(suppressMessages(read_excel(XL, sheet = 'plot coordinates')))
names(cc)[1:3] <- c('Plot_ID', 'E', 'N'); cc$E <- suppressWarnings(as.numeric(cc$E)); cc$N <- suppressWarnings(as.numeric(cc$N))
cc <- cc[is.finite(cc$E), ]; cc$cell <- cellFromXY(g, cbind(cc$E, cc$N))
pm <- merge(cc, s, by = 'cell')
add('plots_with_coordinates', nrow(cc)); add('plots_on_simulated_cells', nrow(pm))
add('plot_cells_height_cellmax_max_m', max(pm$h_max))
add('plot_cells_lai_gt_6.75', paste(sort(pm$Plot_ID[pm$lai_max > LAI_HI]), collapse = ';'))
top <- which.max(s$h_max)
dd <- sqrt((cc$E - s$x[top])^2 + (cc$N - s$y[top])^2)
add('tallest_cell_nearest_plot', cc$Plot_ID[which.min(dd)]); add('tallest_cell_nearest_plot_m', min(dd))

## ---- LAI ----
e <- s$lai_max > LAI_HI
add('lai_cellmax_max', max(s$lai_max)); add('lai_n_gt_6.75', sum(e)); add('lai_share_gt_6.75', mean(e))
for (t in c(8, 10, 12, 14)) add(sprintf('lai_n_gt_%g', t), sum(s$lai_max > t))
add('lai_excess_median', median(s$lai_max[e] - LAI_HI)); add('lai_excess_p75', quantile(s$lai_max[e] - LAI_HI, .75))
add('lai_n_gt_6.75_in_jan2001', sum(e & s$lai_1 > LAI_HI)); add('lai_n_gt_6.75_in_dec2020', sum(s$lai_240 > LAI_HI))
add('lai_n_peak_aug2020', sum(e & s$lai_bmax == 236))
for (k in c('S', 'BPS', 'PS', 'B', 'P')) add(sprintf('lai_gt_6.75_n_%s', k), sum(e & s$combo %in% k))
add('lai_gt_6.75_share_with_spruce', mean(grepl('S', s$combo[e])))
add('lai_gt_6.75_share_spruce_only', mean(s$combo[e] %in% 'S'))
add('lai_share_of_spruce_only_cells_gt_6.75', mean(e[s$combo %in% 'S']))
add('lai_gt_6.75_agb_median', median(s$bio[e])); add('lai_le_6.75_agb_median', median(s$bio[!e]))

out <- do.call(rbind, res); out$value <- as.character(out$value)
write.csv(out, file.path(OUT, 'cube_qc_fail_check.csv'), row.names = FALSE)
print(out, right = FALSE)
