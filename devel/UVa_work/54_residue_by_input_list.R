## ---------------------------------------------------------------------------
## 54_residue_by_input_list.R
##
## For one sector, classifies every mapped-forest cell of the parcel that has species inputs but
## no simulated value: own position lacks soil inputs; number absent from d_site (index points
## elsewhere); absent from d_climate; or in every list by name. Run for the 60 sectors holding
## residue and concatenated as results/residue_by_input_list.csv.
## First run 17-18 Oct 2026. All 9 604 residue cells: 7 053 number absent from d_site, 2 551 own position lacks soil inputs.
##
##   Rscript 54_residue_by_input_list.R <sector> 10 FALSE   (on iufor)
## ---------------------------------------------------------------------------
## The pipeline's own code, lines 1-447 of polygon_grids_v03_server.r (iufor clone), is evaluated
## first, so every list is built exactly as in the production run. The pipeline sets its own
## working directory, so outputs go to the absolute OUT captured here.
OUT  <- normalizePath(if (nzchar(Sys.getenv('OUT'))) Sys.getenv('OUT') else 'results')
PIPE <- '/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r'
src  <- readLines(PIPE)[1:447]
eval(parse(text = src))

## ---- read-only: why each residue cell left COMMON, by cell name ----
B <- '/home/wilar/datos/3PG'
g <- rast(file.path(B, 'outputs_3PG', 'cube', 'gpp.tif'))[[1]]
parts <- disagg(vect(file.path(B, 'krycklan_borders', 'krycklan_borders.shp'))); kry <- project(parts[which.max(expanse(parts, unit = 'ha'))], crs(g))
on <- function(r) { r <- crop(r, project(as.polygons(ext(kryck_elev), crs = crs(g)), crs(r)), snap = 'out'); if (!compareGeom(r, kryck_elev, stopOnError = FALSE)) r <- resample(project(r, crs(g), method = 'near'), kryck_elev, method = 'near'); r }
gs <- crop(g, kryck_elev); stopifnot(compareGeom(gs, kryck_elev, stopOnError = FALSE))
inside <- values(!is.na(mask(rast(gs, vals = 1), kry)), mat = FALSE) %in% TRUE
nmdv <- values(on(rast(file.path(B, 'nationalndcoverdata', 'nmd2018bas_ogeneraliserad_v1_0.tif'))), mat = FALSE)
forest <- inside & nmdv %in% c(111:117, 121:127); sim <- !is.na(values(gs, mat = FALSE))
nS <- as.integer(names(test_site)); nC <- as.integer(names(test_clim)); nP <- as.integer(names(test_species))
k <- which(forest & !sim & seq_along(sim) %in% nP)       # species-grid names of unsimulated forest cells that have species inputs
## a cell in this sector's square only (centre in square): drop buffer-overlap cells belonging to neighbours
sq <- project(kryck_border, crs(g)); insq <- seq_len(ncell(gs)) %in% cells(gs, sq)[, 'cell']; k <- k[insq[k]]
own_site_ok <- !is.na(terra::extract(sites., xyFromCell(gs, k))[['soil_class']])   # site values present at the cell's own position
inS <- k %in% nS; inC <- k %in% nC
reason <- ifelse(!inS & !own_site_ok, 'own position lacks soil inputs',
          ifelse(!inS & own_site_ok, 'name absent from d_site (index points elsewhere)',
          ifelse(!inC, 'name absent from d_climate', 'in all lists by name (unexplained)')))
shifted <- ncol(sites.) != ncol(kryck_elev) || ymax(sites.) != ymax(kryck_elev)
cat(sprintf('RESULT sector %d shifted %s cells %d | %s\n', sector_n, shifted, length(k), paste(names(table(reason)), table(reason), sep = ': ', collapse = ' | ')))
write.csv(data.frame(sector = sector_n, cell = k, reason), file.path(OUT, sprintf('residue_by_input_list_s%03d.csv', sector_n)), row.names = FALSE)
