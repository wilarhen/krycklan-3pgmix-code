## ---------------------------------------------------------------------------
## 21_cube_qc.R
##
## Whole-cube sanity checks after 20_assemble_cube.R. This is the cube-wide
## version of what sector 62 passed locally: the point of it is that sector 62
## was one dense interior sector, so LAI and height there were never going to
## expose an open or young stand behaving badly.
##
## Runs the acceptance tests that do not need the plot data. The plot-based
## tests (pine r, spruce r, bias) stay in the extraction script.
##
##   Rscript 21_cube_qc.R
## ---------------------------------------------------------------------------

library(terra)
library(data.table)

base_path <- '/home/wilar/datos/3PG'
out_dir   <- file.path(base_path, 'outputs_3PG')
cube_dir  <- file.path(out_dir, 'cube')
kry_dir   <- file.path(base_path, 'krycklan_borders')

## reference values from the sector 62 validated run and from the measurements
LAI_MEAS_MAX <- 6.75
LAI_S62_MAX  <- 6.6
HGT_S62_MAX  <- 30.5
GDM_TO_GC    <- 0.47 * 100      # tDM/ha -> gC/m2

terra::terraOptions(tempdir = file.path(base_path, 'r_tmp'), progress = 0)
res_tab <- list()
add <- function(test, value, expect, pass)
    res_tab[[length(res_tab) + 1L]] <<-
        data.table(test = test, value = value, expect = expect,
                   result = ifelse(pass, 'PASS', 'FAIL'))

cube <- function(v) {
    f <- file.path(cube_dir, paste0(v, '.tif'))
    if (!file.exists(f)) return(NULL)
    rast(f)
}

hr <- function(x) cat('\n', strrep('-', 75), '\n', x, '\n', sep = '')

## ---------------------------------------------------------------------------
## 1. coverage
## ---------------------------------------------------------------------------
hr('1. COVERAGE')
g <- cube('gpp')
if (is.null(g)) stop('no gpp cube found in ', cube_dir)

cat('grid: ', ncol(g), ' x ', nrow(g), ' cells, res ', res(g)[1],
    ' m, origin (', paste(origin(g), collapse = ', '), '), ',
    nlyr(g), ' layers\n', sep = '')

g1  <- g[[1]]
n_sim <- as.numeric(global(!is.na(g1), 'sum', na.rm = TRUE))
cat('simulated cells in the cube: ', format(n_sim, big.mark = ' '), '\n', sep = '')

shp <- list.files(kry_dir, pattern = '\\.shp$', full.names = TRUE)
if (length(shp)) {
    bord <- vect(shp[1])
    if (!same.crs(bord, g)) bord <- project(bord, crs(g))
    inside <- mask(g1, bord)
    n_in   <- as.numeric(global(!is.na(inside), 'sum', na.rm = TRUE))
    n_poly <- as.numeric(global(!is.na(mask(rast(g1, vals = 1), bord)), 'sum', na.rm = TRUE))
    cat(sprintf('inside the catchment polygon: %s / %s cells (%.2f %%)\n',
                format(n_in, big.mark = ' '), format(n_poly, big.mark = ' '),
                100 * n_in / n_poly))
    add('catchment coverage', sprintf('%.2f %%', 100 * n_in / n_poly),
        'high, gaps only where land cover is non-forest', n_in / n_poly > 0.5)
    ## where are the holes
    holes <- mask(is.na(g1), bord)
    writeRaster(holes, file.path(cube_dir, 'qc_gaps.tif'), overwrite = TRUE)
    cat('gap map written to qc_gaps.tif, inspect it against the land cover map\n')
} else cat('no catchment shapefile found, skipping the polygon check\n')

## ---------------------------------------------------------------------------
## 2. LAI, hazard 2 (branch biomass fed as foliage)
## ---------------------------------------------------------------------------
hr('2. LAI')
lai <- cube('lai')
if (!is.null(lai)) {
    lai_cellmax <- max(lai, na.rm = TRUE)          # per cell over time
    q <- global(lai_cellmax, fun = function(x)
        quantile(x, c(0, 0.01, 0.5, 0.99, 1), na.rm = TRUE))
    print(q)
    lmax <- as.numeric(global(lai_cellmax, 'max', na.rm = TRUE))
    lmin <- as.numeric(global(lai_cellmax, 'min', na.rm = TRUE))
    n_over <- as.numeric(global(lai_cellmax > LAI_MEAS_MAX, 'sum', na.rm = TRUE))
    cat(sprintf('per-cell max LAI: min %.2f, max %.2f | cells above the measured max %.2f: %s\n',
                lmin, lmax, LAI_MEAS_MAX, format(n_over, big.mark = ' ')))
    add('LAI cube max', sprintf('%.2f', lmax),
        sprintf('near %.1f, not 33.6', LAI_S62_MAX), lmax < 10)
    add('LAI low end', sprintf('%.2f', lmin),
        'well below 6.4, open and young stands must appear', lmin < 3)
    writeRaster(lai_cellmax, file.path(cube_dir, 'qc_lai_cellmax.tif'), overwrite = TRUE)
    if (n_over > 0.05 * n_sim)
        cat('NOTE more than 5 % of cells exceed the measured maximum LAI. ',
            'Check the Marklund foliage fractions before trusting the cube.\n', sep = '')
}

## ---------------------------------------------------------------------------
## 3. height, hazard 3 (wrong height equation form)
## ---------------------------------------------------------------------------
hr('3. HEIGHT')
h <- cube('height')
if (!is.null(h)) {
    h_last <- h[[nlyr(h)]]
    hs <- global(h_last, fun = function(x)
        quantile(x, c(0, 0.5, 0.99, 1), na.rm = TRUE))
    print(hs)
    hmax <- as.numeric(global(h_last, 'max', na.rm = TRUE))
    hmed <- as.numeric(global(h_last, fun = function(x) median(x, na.rm = TRUE)))
    add('height max, last month', sprintf('%.1f m', hmax),
        sprintf('near %.1f m, not 1e22', HGT_S62_MAX), hmax < 60)
    add('height median, last month', sprintf('%.1f m', hmed),
        '5 to 30 m', hmed > 5 && hmed < 30)
}

## ---------------------------------------------------------------------------
## 4. GPP level and the seam test
## ---------------------------------------------------------------------------
hr('4. GPP')
yr <- substr(names(g), 1, 4)
if (length(yr) != nlyr(g) || all(is.na(yr))) yr <- substr(names(g), 1, 4)
years <- sort(unique(yr[!is.na(yr)]))
cat('years in the cube: ', paste(range(years), collapse = ' to '),
    ' (', length(years), ')\n', sep = '')

ann <- rast(lapply(years, function(y) sum(g[[which(yr == y)]], na.rm = TRUE) * GDM_TO_GC))
names(ann) <- years
gs <- global(mean(ann), fun = function(x) quantile(x, c(0, 0.5, 1), na.rm = TRUE))
print(gs)
gmed <- as.numeric(global(mean(ann), fun = function(x) median(x, na.rm = TRUE)))
cat(sprintf('median annual GPP over the cube: %.0f gC/m2/yr\n', gmed))
add('annual GPP median', sprintf('%.0f gC/m2/yr', gmed),
    'below the measured level by roughly 440 with published parameters',
    gmed > 100 && gmed < 1500)
writeRaster(ann, file.path(cube_dir, 'qc_annual_gpp.tif'), overwrite = TRUE)

## Sector seams. Neighbouring sectors overlap by 200 m and were simulated from
## identical inputs, so a visible seam in this map means one sector ran with a
## different configuration, which the overlap spread test in script 20 would
## also catch.
png(file.path(cube_dir, 'qc_annual_gpp.png'), width = 1400, height = 1200, res = 130)
plot(mean(ann), main = 'mean annual GPP, gC/m2/yr  (look for sector seams)')
if (exists('bord')) lines(bord)
dev.off()
cat('map written to qc_annual_gpp.png, inspect for a 10 x 10 seam pattern\n')

## ---------------------------------------------------------------------------
## 5. modifiers, is water still non-limiting cube-wide
## ---------------------------------------------------------------------------
hr('5. MODIFIERS')
for (v in c('f_sw','f_tmp','f_frost','f_vpd')) {
    r <- cube(v)
    if (is.null(r)) next
    m <- as.numeric(global(mean(r, na.rm = TRUE), 'mean', na.rm = TRUE))
    cat(sprintf('%-8s cube mean %.3f\n', v, m))
    if (v == 'f_sw')
        add('f_sw cube mean', sprintf('%.3f', m),
            'near 0.996 as in the preflight, water not limiting', m > 0.9)
    if (v == 'f_tmp')
        add('f_tmp cube mean', sprintf('%.3f', m),
            'near 0.265, temperature is the binding constraint', m < 0.5)
}

## ---------------------------------------------------------------------------
## 6. summary
## ---------------------------------------------------------------------------
hr('6. ACCEPTANCE TESTS, CUBE LEVEL')
tab <- rbindlist(res_tab)
print(tab)
fwrite(tab, file.path(cube_dir, 'qc_summary.csv'))
cat('\nsummary written to qc_summary.csv\n')
cat('the remaining tests, pine r near 0.50, spruce near zero, bias near -440,',
    '\nneed the 57 plots and stay in the extraction script\n')
