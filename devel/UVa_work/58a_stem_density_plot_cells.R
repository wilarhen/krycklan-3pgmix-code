## 58a_stem_density_plot_cells.R
## Stem-density test of the deposited cube (run 21 Oct 2026 on iufor in /tmp/stems_test.qkS5Um; filed unchanged).
## The cube used stems_n 967 (pine), 1162 (spruce), 1000 (birch), species means under the earlier, transposed species
## column; the corrected grouping gives 1094 and 1116. Two arms, everything else identical: arm 1 as the cube, arm 2
## with pine 1094 and spruce 1116 (birch 1000 in both). Inputs are rebuilt by evaluating the pipeline's own lines 1-501
## (polygon_grids_v03_server.r, l.477 removed), so they are the cube's exact per-cell inputs, misindexed site rows included.
##   58a: reproduces the published cube scoring with script 03's 10 m-buffer method and lists the buffer cells
##        (57 plots, 228 cells; 123 plot-years at 41 plots: R2 0.0229, NSE -1.6266, r 0.3760 pine, 0.0994 spruce).
##   58b: one sector per process (40 sectors hold buffer cells, sectors.txt), both arms for the buffer cells; in
##        sector 56 also all 10 934 cells, keeping Dec-2020 height (species mean) and per-cell maximum LAI (species sum).
##   58c: keeps, per cell, the sector whose arm 1 reproduces the cube (177 cells, max |arm1 - cube| 3.7e-7 t DM/ha/month),
##        then scores both arms as 58a: arm 1 reproduces the published values; arm 2 gives R2 0.0210, NSE -1.6557,
##        r 0.3695 pine, 0.0785 spruce. Sector 56 summary for both arms.
## Run on iufor:  RD=$(mktemp -d); Rscript 58a_... $RD; cat $RD/sectors.txt | xargs -P 10 -I{} Rscript 58b_... {} 10 FALSE $RD;
##                Rscript 58c_... $RD          (58a writes plot_buffer_cells_cube_gpp.csv and cube_scoring_pairs.csv)
## Outputs as filed: results/stem_density_test/.
## ---------------------------------------------------------------------------
## Read-only. Reproduce script 03's 10 m-buffer scoring of the corrected cube; record each plot's buffer cells.
suppressPackageStartupMessages({ library(terra); library(sf); library(readxl); library(dplyr); library(tidyr) })
OUT <- commandArgs(TRUE)[1]
cube <- rast('/home/wilar/datos/3PG/outputs_3PG/cube/gpp.tif'); tt <- time(cube); if (any(is.na(tt))) tt <- as.Date(paste0(names(cube), '-01'))
xl <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
rp <- function(sheet) { x <- as.data.frame(suppressMessages(read_excel(xl, sheet = sheet, skip = 2, .name_repair = 'unique')))
  if (sheet == '50 main plots') { i <- match(c('SPM', 'DTS'), names(x)); if (!any(is.na(i))) names(x)[i] <- c('DTS', 'SPM') }; x }
obs <- bind_rows(rp('50 main plots') |> mutate(set = 'main'), rp('7 extra old-stand plots') |> mutate(set = 'oldstand')) |>
  mutate(Plot_ID = as.integer(Plot_ID), Year = as.integer(Year)) |> filter(!is.na(Plot_ID), Year %in% 2016:2018)
co <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |> rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(Plot_ID = as.integer(Plot_ID), E = as.numeric(E), N = as.numeric(N)) |> filter(!is.na(E), !is.na(N))
pts <- st_as_sf(distinct(inner_join(obs, co, by = 'Plot_ID'), Plot_ID, E, N), coords = c('E', 'N'), crs = 3006) |> st_transform(crs(cube))
buf <- vect(st_buffer(pts, 10))
cl <- terra::extract(cube[[1]], buf, cells = TRUE, ID = TRUE); cl$Plot_ID <- pts$Plot_ID[cl$ID]
cl <- cl[, c('Plot_ID', 'cell')]; xy <- xyFromCell(cube, cl$cell); cl$x <- xy[, 1]; cl$y <- xy[, 2]
v <- as.matrix(cube[cl$cell]); keep <- which(format(tt, '%Y') %in% c('2016', '2017', '2018'))
cl_gpp <- cbind(cl, v[, keep]); write.csv(cl_gpp, file.path(OUT, 'plot_buffer_cells_cube_gpp.csv'), row.names = FALSE)
cat('plots with coordinates', nrow(pts), '| buffer cells', nrow(cl), '| distinct cells', length(unique(cl$cell)), '\n')
## script 03: monthly mean over the buffer (na.rm), annual sum of complete years, x 0.47 x 100
mon <- data.frame(Plot_ID = cl$Plot_ID, v[, keep]) |> pivot_longer(-Plot_ID, names_to = 'm', values_to = 'value')
mon$m <- rep(format(tt[keep], '%Y-%m'), times = nrow(cl))
ann <- mon |> group_by(Plot_ID, m) |> summarise(value = mean(value, na.rm = TRUE), .groups = 'drop') |>
  mutate(Year = as.integer(substr(m, 1, 4)), value = ifelse(is.nan(value), NA, value)) |>
  group_by(Plot_ID, Year) |> summarise(n_month = sum(!is.na(value)), sim = sum(value, na.rm = TRUE) * 0.47 * 100, .groups = 'drop') |> filter(n_month == 12)
val <- inner_join(ann, obs |> transmute(Plot_ID, Year, obs = GPP - GPPff, DTS = factor(DTS, 1:2, c('Pine', 'Spruce'))), by = c('Plot_ID', 'Year')) |> filter(!is.na(obs))
e <- val$sim - val$obs
cat(sprintf('cube scoring: n %d plot-years at %d plots | R2 %.4f | NSE %.4f | bias %.1f\n', nrow(val), n_distinct(val$Plot_ID), cor(val$obs, val$sim)^2, 1 - sum(e^2) / sum((val$obs - mean(val$obs))^2), mean(e)))
pm <- val |> group_by(Plot_ID, DTS) |> summarise(sim = mean(sim), obs = mean(obs), .groups = 'drop')
for (s in c('Pine', 'Spruce')) { z <- pm[pm$DTS == s, ]; cat(sprintf('  %s: r %.4f on %d plot means\n', s, cor(z$sim, z$obs), nrow(z))) }
write.csv(val, file.path(OUT, 'cube_scoring_pairs.csv'), row.names = FALSE)
