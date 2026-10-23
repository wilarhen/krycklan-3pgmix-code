## ---------------------------------------------------------------------------
## 52_site_misindex_check.R
##
## For one sector, rebuilds d_site, d_climate and d_species with the pipeline's own code and
## prints each list's raster geometry, the position three common cell numbers resolve to in each,
## the offsets over all COMMON cells, and latitude, altitude, soil class and climate used against
## the values at each cell's own position. Read-only; prints only.
## First run 17-18 Oct 2026. Sector 37: offsets 0; sector 56: all 10 934 cells displaced, altitude median 4.42 m, max 19.21 m.
##
##   Rscript 52_site_misindex_check.R <sector> 10 FALSE   (on iufor)
## ---------------------------------------------------------------------------
## The pipeline's own code, lines 1-447 of polygon_grids_v03_server.r (iufor clone), is evaluated
## first, so every list is built exactly as in the production run. The pipeline sets its own
## working directory, so outputs go to the absolute OUT captured here.
OUT  <- normalizePath(if (nzchar(Sys.getenv('OUT'))) Sys.getenv('OUT') else 'results')
PIPE <- '/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r'
src  <- readLines(PIPE)[1:447]
eval(parse(text = src))

## ---- read-only: which grid does each list's cell name index? ----
geo <- function(lab, r) { e <- as.vector(ext(r)); cat(sprintf('%-28s ncol %4d nrow %4d cells %6d | ext x %.1f-%.1f y %.1f-%.1f | res %g | origin %s\n',
  lab, ncol(r), nrow(r), ncell(r), e[1], e[2], e[3], e[4], res(r)[1], paste(round(origin(r), 3), collapse = ','))) }
cat('\n== rasters behind each list ==\n')
geo('ref_grid', ref_grid); geo('kryck_elev', kryck_elev)
geo('d_site: sites.', sites.[[1]])
geo('d_climate: clim', rasters[[1]]); geo('d_climate: radiation', monthly_rad_roi_proj[[1]])
geo('d_species: planted', specs.[[1]]); geo('d_species: species', specs.[[2]]); geo('d_species: biom_stem', specs.[[5]])
nS <- names(test_site); nC <- names(test_clim); nP <- names(test_species)
cat(sprintf('\nlist sizes: d_site %d | d_climate %d | d_species %d | COMMON %d\n', length(nS), length(nC), length(nP), length(Reduce(intersect, list(nS, nC, nP)))))
common <- as.integer(Reduce(intersect, list(nS, nC, nP)))
pick <- common[round(quantile(seq_along(common), c(.1, .5, .9)))]
cat('\n== three common cell names and the position each list resolves them to ==\n')
for (k in pick) {
  s <- xyFromCell(sites.[[1]], k); c <- xyFromCell(rasters[[1]], k); p <- xyFromCell(specs.[[1]], k)
  cat(sprintf('cell %6d: d_site (%.0f, %.0f) | d_climate (%.0f, %.0f) | d_species (%.0f, %.0f) | site-species offset dx %+.0f dy %+.0f | climate-species dx %+.0f dy %+.0f\n',
      k, s[1], s[2], c[1], c[2], p[1], p[2], s[1] - p[1], s[2] - p[2], c[1] - p[1], c[2] - p[2]))
}
## all common cells: offsets per list relative to d_species
dS <- xyFromCell(sites.[[1]], common) - xyFromCell(specs.[[1]], common)
dC <- xyFromCell(rasters[[1]], common) - xyFromCell(specs.[[1]], common)
cat('\nsite minus species offsets over all COMMON cells:\n'); print(table(paste0('dx ', dS[, 1], ' dy ', dS[, 2])))
cat('climate minus species offsets over all COMMON cells:\n'); print(table(paste0('dx ', dC[, 1], ' dy ', dC[, 2])))
## consequence: values each simulated cell received against the values at its own position
pos <- xyFromCell(specs.[[1]], common)
alt_used <- values(sites.[['altitude']])[common]; alt_true <- terra::extract(sites.[['altitude']], pos)[[1]]
cat(sprintf('\naltitude used vs at own position: identical %d | differ %d | median |diff| %.2f m | max %.2f m\n',
  sum(alt_used == alt_true, na.rm = TRUE), sum(alt_used != alt_true, na.rm = TRUE), median(abs(alt_used - alt_true), na.rm = TRUE), max(abs(alt_used - alt_true), na.rm = TRUE)))
lat_used <- values(sites.[['latitude']])[common]; lat_true <- terra::extract(sites.[['latitude']], pos)[[1]]
cat(sprintf('latitude used vs own position: differ %d | max |diff| %.5f deg\n', sum(lat_used != lat_true, na.rm = TRUE), max(abs(lat_used - lat_true), na.rm = TRUE)))
sc_used <- values(sites.[['soil_class']])[common]; sc_true <- terra::extract(sites.[['soil_class']], pos)[[1]]
cat(sprintf('soil class used vs own position: differ %d | NA at own position %d\n', sum(sc_used != sc_true, na.rm = TRUE), sum(is.na(sc_true))))
t1 <- rasters[[1]][[1]]; cl_used <- values(t1)[common]; cl_true <- terra::extract(t1, pos)[[1]]
cat(sprintf('climate (%s, first month) used vs own position: differ %d of %d\n', names(rasters)[1], sum(cl_used != cl_true, na.rm = TRUE), length(common)))
