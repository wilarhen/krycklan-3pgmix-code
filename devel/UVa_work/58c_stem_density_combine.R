## 58c_stem_density_combine.R
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
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
RD <- commandArgs(TRUE)[1]
cub <- read.csv(file.path(RD, 'plot_buffer_cells_cube_gpp.csv'), check.names = FALSE)
mcols <- setdiff(names(cub), c('Plot_ID', 'cell', 'x', 'y'))
cl <- cub |> pivot_longer(all_of(mcols), names_to = 'layer', values_to = 'cube') |> group_by(Plot_ID, cell) |> mutate(month = sprintf('%s', c(sprintf('2016-%02d', 1:12), sprintf('2017-%02d', 1:12), sprintf('2018-%02d', 1:12)))) |> ungroup()
f <- list.files(RD, '^plotcells_s[0-9]+\\.csv$', full.names = TRUE)
sim <- do.call(rbind, lapply(f, read.csv))
a1 <- sim |> filter(arm == 1) |> inner_join(cl |> select(cell, month, cube), by = c('cell', 'month'), relationship = 'many-to-many')
pick <- a1 |> group_by(cell, sector) |> summarise(maxdiff = max(abs(gpp - cube)), .groups = 'drop') |> group_by(cell) |> slice_min(maxdiff, n = 1, with_ties = FALSE) |> ungroup()
cat('buffer cells with cube values:', n_distinct(cl$cell[!is.na(cl$cube)]), '| recovered from the archive rebuild:', nrow(pick), '| max |arm1 - cube| over all cells and months:', signif(max(pick$maxdiff), 3), '\n')
cat('cube cells NA (unsimulated) but rebuilt anyway:', sum(pick$cell %in% unique(cl$cell[is.na(cl$cube)])), '\n')
sel <- sim |> inner_join(pick |> select(cell, sector), by = c('cell', 'sector')) |> distinct(Plot_ID, cell, arm, month, gpp)
xl <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
pairs <- read.csv(file.path(RD, 'cube_scoring_pairs.csv'))
obs <- distinct(pairs, Plot_ID, Year, obs, DTS)
score <- function(m) {
  ann <- m |> group_by(Plot_ID, month) |> summarise(v = mean(value, na.rm = TRUE), .groups = 'drop') |> mutate(Year = as.integer(substr(month, 1, 4)), v = ifelse(is.nan(v), NA, v)) |>
    group_by(Plot_ID, Year) |> summarise(n = sum(!is.na(v)), sim = sum(v, na.rm = TRUE) * 0.47 * 100, .groups = 'drop') |> filter(n == 12)
  val <- inner_join(ann, obs, by = c('Plot_ID', 'Year')); e <- val$sim - val$obs
  pm <- val |> group_by(Plot_ID, DTS) |> summarise(sim = mean(sim), obs = mean(obs), .groups = 'drop')
  data.frame(n = nrow(val), plots = n_distinct(val$Plot_ID), R2 = cor(val$obs, val$sim)^2, NSE = 1 - sum(e^2) / sum((val$obs - mean(val$obs))^2), bias = mean(e),
             r_pine = cor(pm$sim[pm$DTS == 'Pine'], pm$obs[pm$DTS == 'Pine']), n_pine = sum(pm$DTS == 'Pine'),
             r_spruce = cor(pm$sim[pm$DTS == 'Spruce'], pm$obs[pm$DTS == 'Spruce']), n_spruce = sum(pm$DTS == 'Spruce'),
             mean_sim = mean(val$sim)) }
base <- cl |> filter(!is.na(cube)) |> transmute(Plot_ID, cell, month, value = cube)
out <- rbind(cbind(source = 'cube (published)', score(base)),
             cbind(source = 'arm 1: 967 / 1162', score(sel |> filter(arm == 1) |> transmute(Plot_ID, cell, month, value = gpp))),
             cbind(source = 'arm 2: 1094 / 1116', score(sel |> filter(arm == 2) |> transmute(Plot_ID, cell, month, value = gpp))))
print(transform(out, R2 = round(R2, 4), NSE = round(NSE, 4), bias = round(bias, 1), r_pine = round(r_pine, 4), r_spruce = round(r_spruce, 4), mean_sim = round(mean_sim, 1)), row.names = FALSE)
s <- read.csv(file.path(RD, 'sector56_cells.csv'))
for (arm in 1:2) { z <- s[s$arm == arm, ]; cat(sprintf('sector 56 arm %d: cells %d | height Dec 2020 median %.3f m, p99 %.3f m, max %.2f m | LAI max %.3f | cells LAI > 6.75: %d (%.1f %%)\n',
  arm, nrow(z), median(z$h_dec2020), quantile(z$h_dec2020, .99), max(z$h_dec2020), max(z$lai_max), sum(z$lai_max > 6.75), 100 * mean(z$lai_max > 6.75))) }
w <- merge(s[s$arm == 1, c('k','h_dec2020','lai_max')], s[s$arm == 2, c('k','h_dec2020','lai_max')], by = 'k', suffixes = c('_1', '_2'))
cat(sprintf('per-cell change, arm 2 - arm 1: height median %+.3f m (range %+.2f to %+.2f) | LAI max median %+.3f (range %+.3f to %+.3f)\n',
  median(w$h_dec2020_2 - w$h_dec2020_1), min(w$h_dec2020_2 - w$h_dec2020_1), max(w$h_dec2020_2 - w$h_dec2020_1),
  median(w$lai_max_2 - w$lai_max_1), min(w$lai_max_2 - w$lai_max_1), max(w$lai_max_2 - w$lai_max_1)))
