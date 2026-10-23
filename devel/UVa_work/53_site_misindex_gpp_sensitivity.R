## ---------------------------------------------------------------------------
## 53_site_misindex_gpp_sensitivity.R
##
## For one sector, reruns every simulated cell twice with the pipeline's settings: with the site
## row the cube used and with the site values at the cell's own position; reports annual GPP
## differences (mean 2001-2020 and 2017). Writes results/site_misindex_sensitivity_sNN.csv.
## First run 17-18 Oct 2026. Sector 56, 10 934 cells, 9.3 min: at most 0.0041 % (2017), 0.0036 % (mean).
##
##   Rscript 53_site_misindex_gpp_sensitivity.R <sector> 10 FALSE   (on iufor)
## ---------------------------------------------------------------------------
## The pipeline's own code, lines 1-501 of polygon_grids_v03_server.r (iufor clone), is evaluated
## first, so every list is built exactly as in the production run. Line 477 (write.csv of d_parameters_adapted) is dropped. The pipeline sets its own
## working directory, so outputs go to the absolute OUT captured here.
OUT  <- normalizePath(if (nzchar(Sys.getenv('OUT'))) Sys.getenv('OUT') else 'results')
PIPE <- '/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r'
src  <- readLines(PIPE)[1:501][-477]
eval(parse(text = src))

## ---- read-only sensitivity: misindexed site (as the cube used it) against each cell's own ----
smp <- common_names
own_row <- function(k) {
  p <- xyFromCell(specs.[[1]], as.integer(k)); v <- terra::extract(sites., p)
  s <- filtered_ls[['d_site']][[k]]
  for (nm in c('latitude', 'altitude', 'soil_class', 'asw_i', 'asw_min', 'asw_max')) s[[nm]] <- v[[nm]][1]
  s }
st <- list(light_model = 2, transp_model = 2, phys_model = 2, height_model = 2, correct_bias = 0, calculate_d13c = 0)
run1 <- function(site, k) tryCatch(run_3PG(site = site, species = filtered_ls[['d_species']][[k]], climate = filtered_ls[['d_climate']][[k]],
  thinning = NULL, parameters = filtered_ls[['d_parameters']][[k]], size_dist = filtered_ls[['d_sizeDist']][[k]],
  settings = st, check_input = TRUE, df_out = TRUE), error = function(e) NULL)
annual <- function(o) { if (is.null(o)) return(NULL); g <- o[o$variable == 'gpp', ]
  m <- tapply(g$value, as.character(g$date), sum, na.rm = TRUE); y <- substr(names(m), 1, 4); tapply(m, y, sum) }
t0 <- Sys.time(); rows <- list()
for (k in smp) {
  a <- filtered_ls[['d_site']][[k]]; b <- own_row(k)
  if (any(is.na(unlist(b[, c('latitude','altitude','soil_class','asw_i','asw_min','asw_max')])))) { rows[[k]] <- data.frame(cell = k, note = 'own position has NA site input'); next }
  ga <- annual(run1(a, k)); gb <- annual(run1(b, k))
  if (is.null(ga) || is.null(gb)) { rows[[k]] <- data.frame(cell = k, note = 'run failed'); next }
  rows[[k]] <- data.frame(cell = k, note = 'ok', dlat = a$latitude - b$latitude, dalt = a$altitude - b$altitude,
    dsoil = a$soil_class != b$soil_class, dasw = a$asw_i - b$asw_i,
    gpp_used_mean = mean(ga), gpp_own_mean = mean(gb), gpp_used_2017 = ga[['2017']], gpp_own_2017 = gb[['2017']])
}
R <- do.call(rbind, lapply(rows, function(x) { for (c in c('dlat','dalt','dsoil','dasw','gpp_used_mean','gpp_own_mean','gpp_used_2017','gpp_own_2017')) if (!c %in% names(x)) x[[c]] <- NA; x }))
cat(sprintf('\nall %d cells of %d simulated in sector %d; %.1f min\n', length(smp), length(common_names), sector_n, as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
print(table(R$note)); ok <- R[R$note == 'ok', ]
cat('site differences used minus own: latitude (deg)\n'); print(summary(ok$dlat)); cat('altitude (m)\n'); print(summary(ok$dalt)); cat('soil class differs:', sum(ok$dsoil), '| ASW differs:', sum(ok$dasw != 0), '\n')
for (w in c('mean', '2017')) { u <- ok[[paste0('gpp_used_', w)]]; o <- ok[[paste0('gpp_own_', w)]]; d <- u - o; p <- 100 * d / o
  cat(sprintf('\nannual GPP (%s), t DM/ha/yr: used median %.4f | own median %.4f\n', if (w == 'mean') 'mean 2001-2020' else '2017', median(u), median(o)))
  cat('difference used - own (t DM/ha/yr):\n'); print(quantile(d, c(0, .05, .25, .5, .75, .95, 1)))
  cat('difference (%):\n'); print(quantile(p, c(0, .05, .25, .5, .75, .95, 1)))
  cat(sprintf('cells with |difference| = 0: %d of %d | max |diff| %.6f t DM/ha/yr (%.4f %%)\n', sum(d == 0), length(d), max(abs(d)), max(abs(p)))) }
write.csv(R, file.path(OUT, sprintf('site_misindex_sensitivity_s%02d.csv', sector_n)), row.names = FALSE)
