## 58b_stem_density_worker.R
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
## Read-only two-arm stem-density test for one sector. Usage: Rscript worker.R <sector> 10 FALSE <RD>
## Evaluates the pipeline's lines 1-501 (l.477, a write.csv, removed) to rebuild the cube's exact inputs, then
## runs r3PG with stems_n as in the cube (arm 1: 967/1162/1000) and with 1094 pine / 1116 spruce (arm 2).
RD <- commandArgs(TRUE)[4]
PIPE <- '/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r'
src <- readLines(PIPE)[1:501][-477]; eval(parse(text = src))
cl <- read.csv(file.path(RD, 'plot_buffer_cells_cube_gpp.csv'))
g <- specs.[[1]]; k <- cellFromXY(g, cbind(cl$x, cl$y)); cl$k <- as.character(k)
pc <- cl[!is.na(k) & cl$k %in% common_names, ]
st <- list(light_model = 2, transp_model = 2, phys_model = 2, height_model = 2, correct_bias = 0, calculate_d13c = 0)
new_stems <- c('Pinus sylvestris' = 1094, 'Picea abies' = 1116)
run <- function(kk, arm) {
  sp <- filtered_ls[['d_species']][[kk]]
  old <- c('Pinus sylvestris' = 967, 'Picea abies' = 1162, 'Betula' = 1000)
  stopifnot(all(abs(sp$stems_n - old[sp$species]) < 1e-9))
  if (arm == 2) { i <- sp$species %in% names(new_stems); sp$stems_n[i] <- new_stems[sp$species[i]] }
  tryCatch(run_3PG(site = filtered_ls[['d_site']][[kk]], species = sp, climate = filtered_ls[['d_climate']][[kk]], thinning = NULL,
    parameters = filtered_ls[['d_parameters']][[kk]], size_dist = filtered_ls[['d_sizeDist']][[kk]],
    settings = st, check_input = TRUE, df_out = TRUE), error = function(e) NULL) }
stand <- function(o, v, f) { z <- o[o$variable == v, ]; tapply(z$value, format(z$date, '%Y-%m'), f, na.rm = TRUE) }
rows <- list()
for (j in seq_len(nrow(pc))) for (arm in 1:2) {
  o <- run(pc$k[j], arm); if (is.null(o)) next
  gm <- stand(o, 'gpp', sum); gm <- gm[substr(names(gm), 1, 4) %in% c('2016', '2017', '2018')]
  rows[[length(rows) + 1]] <- data.frame(sector = sector_n, Plot_ID = pc$Plot_ID[j], cell = pc$cell[j], arm = arm, month = names(gm), gpp = as.vector(gm)) }
if (length(rows)) write.csv(do.call(rbind, rows), file.path(RD, sprintf('plotcells_s%03d.csv', sector_n)), row.names = FALSE)
cat('sector', sector_n, ': plot-buffer cells in this sector', nrow(pc), '\n')
if (sector_n == 56) {
  t0 <- Sys.time(); res <- list()
  for (kk in common_names) for (arm in 1:2) { o <- run(kk, arm); if (is.null(o)) next
    h <- stand(o, 'height', mean); l <- stand(o, 'lai', sum)
    res[[length(res) + 1]] <- data.frame(k = kk, arm = arm, h_dec2020 = h[['2020-12']], lai_max = max(l)) }
  write.csv(do.call(rbind, res), file.path(RD, 'sector56_cells.csv'), row.names = FALSE)
  cat('sector 56 full two-arm run:', length(common_names), 'cells,', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), 'min\n') }
