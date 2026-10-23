## ---------------------------------------------------------------------------
## 51_sector_grid_crops.R
##
## For each of the 100 processing sectors, the grid the species, biomass and climate rasters
## use (crop with snap = 'out', as crop_raster_fast) against the grid the site stack is cropped to
## (crop(sites, kryck_border), terra's default snap), and the rows and columns trimmed. A left,
## right or top trim shifts the site table's cell numbers relative to the other tables.
## Writes results/sector_grid_crops.csv.
## First run 17-18 Oct 2026. 80 sectors shift; 20 (bottom row only) do not.
##
##   Rscript 51_sector_grid_crops.R      (on iufor, from UVa_work/)
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(terra); library(sf) })
B <- '/home/wilar/datos/3PG'
OUT <- normalizePath(if (nzchar(Sys.getenv('OUT'))) Sys.getenv('OUT') else 'results')
bio <- rast(file.path(B, 'biomass', 'BIO.tif'))
bsf <- st_as_sf(vect(file.path(B, 'krycklan_borders', 'krycklan_borders.shp')))
kb <- suppressWarnings(st_cast(bsf, 'POLYGON'))[3, ]; grid10 <- st_make_grid(st_buffer(kb, 200), n = 10)
out <- do.call(rbind, lapply(seq_along(grid10), function(s) {
  sq <- project(vect(grid10[s]), crs(bio))
  ref <- crop(bio, sq, snap = 'out')                 # crop_raster_fast(..., snap = 'out'): the pipeline grid
  sit <- crop(rast(ref), sq)                         # crop(sites, kryck_border): terra default snap
  er <- as.vector(ext(ref)); es <- as.vector(ext(sit))
  data.frame(sector = s, ref_ncol = ncol(ref), ref_nrow = nrow(ref), site_ncol = ncol(sit), site_nrow = nrow(sit),
             trim_left = round((es[1] - er[1]) / 10), trim_right = round((er[2] - es[2]) / 10),
             trim_bottom = round((es[3] - er[3]) / 10), trim_top = round((er[4] - es[4]) / 10))
}))
out$indices_shift <- out$trim_left > 0 | out$trim_right > 0 | out$trim_top > 0
cat('sectors by trim pattern (left, right, bottom, top):\n'); print(table(paste(out$trim_left, out$trim_right, out$trim_bottom, out$trim_top)))
cat('sectors whose site-grid cell numbers shift relative to the species/climate grid:', sum(out$indices_shift), '\n')
if (any(out$indices_shift)) print(out[out$indices_shift, ], row.names = FALSE)
write.csv(out, file.path(OUT, 'sector_grid_crops.csv'), row.names = FALSE); cat("\nsector 37:\n"); print(out[out$sector == 37, ], row.names = FALSE)
