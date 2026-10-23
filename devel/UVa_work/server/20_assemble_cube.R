## ---------------------------------------------------------------------------
## 20_assemble_cube.R
##
## Assemble the wall-to-wall monthly cube from the per-sector tiles written by
## polygon_grids_v03.r into outputs_3PG/<variable>/.
##
## Three stages, each can be run alone:
##   A. RUN AUDIT      did all 100 sectors finish, and how complete are they
##   B. TILE INVENTORY are the tiles mutually compatible (grid, bands, names)
##   C. ASSEMBLE       one multi-band GeoTIFF per variable in outputs_3PG/cube/
##
## Sectors are built with polygon_sector(..., 200), so neighbours OVERLAP by
## 200 m (20 cells). Overlapping cells were simulated twice from identical
## inputs, so they must agree. Stage C exploits that as a free consistency
## test before the mosaic is trusted.
##
##   Rscript 20_assemble_cube.R            # audit + inventory only, no writing
##   Rscript 20_assemble_cube.R assemble   # everything
## ---------------------------------------------------------------------------

library(terra)
library(data.table)

## ---- configuration ---------------------------------------------------------
base_path <- '/home/wilar/datos/3PG'
work_dir  <- '/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors'
out_dir   <- file.path(base_path, 'outputs_3PG')
cube_dir  <- file.path(out_dir, 'cube')
log_dir   <- file.path(work_dir, 'logs')

N_SECTORS <- 100L
DATE_FROM <- '2001-01'          # the d_climate span in v03
DATE_TO   <- '2020-12'
N_LAYERS  <- 240L

## same split as FIX 7 in the pipeline: fluxes and stocks sum across species,
## state variables and modifiers average. Only affects the overlap check here,
## the tiles are already reduced, but keep the lists in one place.
vars_sum  <- c('gpp','npp','lai','biom_stem','biom_foliage','biom_root',
               'basal_area','stems_n')
vars_mean <- c('height','dbh','f_sw','f_vpd','f_tmp','f_frost','f_nutr','f_age',
               'alpha_c')
keep_vars <- c(vars_sum, vars_mean)

RAST_PAT <- '\\.(tif|tiff|TIF|grd|nc)$'
do_assemble <- 'assemble' %in% commandArgs(trailingOnly = TRUE)

terra::terraOptions(tempdir = file.path(base_path, 'r_tmp'), progress = 0)
dir.create(cube_dir, showWarnings = FALSE, recursive = TRUE)

hr <- function(x) cat('\n', strrep('-', 75), '\n', x, '\n', sep = '')

## ===========================================================================
## A. RUN AUDIT
## ===========================================================================
audit_run <- function() {
    hr('A. RUN AUDIT')

    ## .done may be a directory of markers or a flat list, handle both
    done_path <- file.path(log_dir, '.done')
    done_raw <- if (dir.exists(done_path)) list.files(done_path) else
                if (file.exists(done_path)) readLines(done_path, warn = FALSE) else
                character(0)
    done_n <- suppressWarnings(as.integer(gsub('\\D', '', done_raw)))
    done_n <- sort(unique(done_n[is.finite(done_n)]))
    missing_n <- setdiff(seq_len(N_SECTORS), done_n)

    cat('sectors marked done: ', length(done_n), ' / ', N_SECTORS, '\n', sep = '')
    if (length(missing_n))
        cat('NOT done: ', paste(missing_n, collapse = ', '), '\n', sep = '')

    ## completeness from the logs
    logs <- list.files(log_dir, pattern = '^sector_\\d+\\.log$', full.names = TRUE)
    parse_log <- function(f) {
        ln  <- readLines(f, warn = FALSE)
        hit <- grep('cells simulated:', ln, value = TRUE)
        sec <- suppressWarnings(as.integer(gsub('\\D', '', basename(f))))
        if (!length(hit))
            return(data.table(sector = sec, ok = NA_integer_, total = NA_integer_,
                              fail = NA_integer_, empty = any(grepl('no simulated cells', ln))))
        m <- regmatches(hit[length(hit)],
                        regexec('cells simulated: *([0-9]+) */ *([0-9]+).*failures: *([0-9]+)',
                                hit[length(hit)]))[[1]]
        data.table(sector = sec, ok = as.integer(m[2]), total = as.integer(m[3]),
                   fail = as.integer(m[4]), empty = FALSE)
    }
    tab <- if (length(logs)) rbindlist(lapply(logs, parse_log)) else
           data.table(sector = integer(0), ok = integer(0), total = integer(0),
                      fail = integer(0), empty = logical(0))
    setorder(tab, sector)

    if (nrow(tab)) {
        tab[, pct := 100 * ok / total]
        cat('\ncells simulated over all logged sectors: ',
            format(sum(tab$ok, na.rm = TRUE), big.mark = ' '), ' / ',
            format(sum(tab$total, na.rm = TRUE), big.mark = ' '),
            sprintf(' (%.2f %%)\n', 100 * sum(tab$ok, na.rm = TRUE) /
                                    sum(tab$total, na.rm = TRUE)), sep = '')
        bad <- tab[is.na(pct) | pct < 100 | fail > 0 | empty]
        if (nrow(bad)) { cat('\nsectors not at 100 %:\n'); print(bad) } else
            cat('every logged sector reports 100 % and zero failures\n')
    } else cat('no sector logs found under ', log_dir, '\n', sep = '')

    ## the per-cell failure tables, if any were written
    fails <- list.files(out_dir, pattern = '^failures_.*\\.csv$', full.names = TRUE)
    if (length(fails)) {
        ft <- rbindlist(lapply(fails, fread), fill = TRUE)
        cat('\nfailure csv files: ', length(fails), ', rows: ', nrow(ft), '\n', sep = '')
        print(head(sort(table(ft$msg), decreasing = TRUE), 5))
    } else cat('\nno failures_*.csv written\n')

    fwrite(tab, file.path(cube_dir, 'run_audit.csv'))
    invisible(tab)
}

## ===========================================================================
## B. TILE INVENTORY
## ===========================================================================
tile_files <- function(v)
    list.files(file.path(out_dir, v), pattern = RAST_PAT, full.names = TRUE)

inventory_var <- function(v) {
    fs <- tile_files(v)
    if (!length(fs)) return(data.table(variable = v, tiles = 0L))
    info <- rbindlist(lapply(fs, function(f) {
        r <- rast(f)
        data.table(file = basename(f), nlyr = nlyr(r),
                   resx = res(r)[1], resy = res(r)[2],
                   ox = origin(r)[1], oy = origin(r)[2],
                   crs = crs(r, describe = TRUE)$code,
                   nm_first = names(r)[1], nm_last = names(r)[nlyr(r)])
    }))
    data.table(variable = v,
               tiles      = nrow(info),
               nlyr_min   = min(info$nlyr), nlyr_max = max(info$nlyr),
               res_ok     = all(abs(info$resx - 10) < 1e-6 & abs(info$resy - 10) < 1e-6),
               origin_ok  = all(abs(info$ox) < 1e-6 & abs(info$oy) < 1e-6),
               crs_n      = uniqueN(info$crs),
               names_n    = uniqueN(paste(info$nm_first, info$nm_last)),
               size_gb    = round(sum(file.size(fs)) / 1e9, 2))
}

inventory <- function(vars = keep_vars) {
    hr('B. TILE INVENTORY')
    inv <- rbindlist(lapply(vars, inventory_var), fill = TRUE)
    print(inv)
    if (any(inv$tiles == 0L))
        cat('\nno tiles for: ',
            paste(inv[tiles == 0L, variable], collapse = ', '), '\n', sep = '')
    ok <- inv[tiles > 0L]
    if (nrow(ok)) {
        if (any(!ok$origin_ok))
            cat('\nWARNING origin is not (0,0) for: ',
                paste(ok[origin_ok == FALSE, variable], collapse = ', '),
                '  <- FIX 1 did not propagate\n', sep = '')
        if (any(ok$nlyr_min != N_LAYERS | ok$nlyr_max != N_LAYERS))
            cat('\nWARNING band count is not ', N_LAYERS,
                ' everywhere. Assembly will fall back to layerwise mosaic.\n', sep = '')
        if (any(ok$crs_n > 1L))
            cat('\nWARNING mixed CRS among tiles\n')
    }
    fwrite(inv, file.path(cube_dir, 'tile_inventory.csv'))
    invisible(inv)
}

## ===========================================================================
## C. ASSEMBLE
## ===========================================================================

## dates: take them from the tile band names if they parse, else rebuild the
## known d_climate span. Band descriptions survive GeoTIFF, time() does not.
cube_dates <- function(f) {
    nm <- names(rast(f))
    d  <- suppressWarnings(as.Date(nm, format = '%Y-%m-%d'))
    if (all(is.na(d))) d <- suppressWarnings(as.Date(paste0(nm, '-01'), format = '%Y-%m-%d'))
    if (anyNA(d)) d <- seq(as.Date(paste0(DATE_FROM, '-01')),
                           as.Date(paste0(DATE_TO,   '-01')), by = 'month')
    d
}

## Overlap consistency. Sector buffers are 200 m, so interior cells appear in
## up to four tiles. Identical inputs must give identical output, so the
## spread across tiles is a direct test that no sector ran with a different
## configuration. One layer only, cheap.
overlap_check <- function(v, layer = 7L) {
    fs <- tile_files(v)
    if (length(fs) < 2L) return(invisible(NULL))
    lst <- lapply(fs, function(f) rast(f, lyrs = min(layer, nlyr(rast(f)))))
    sc  <- sprc(lst)
    mn  <- mosaic(sc, fun = 'min')
    mx  <- mosaic(sc, fun = 'max')
    cnt <- mosaic(sprc(lapply(lst, function(r) !is.na(r))), fun = 'sum')
    spread <- mx - mn
    cat(sprintf('  overlap: cells covered by >1 tile = %s, max spread = %.6g\n',
                format(as.numeric(global(cnt > 1, 'sum', na.rm = TRUE)), big.mark = ' '),
                as.numeric(global(spread, 'max', na.rm = TRUE))))
    invisible(list(count = cnt, spread = spread))
}

assemble_vrt <- function(v, fs, dts) {
    vf <- file.path(cube_dir, paste0(v, '.vrt'))
    r  <- terra::vrt(fs, filename = vf, overwrite = TRUE)
    if (nlyr(r) != length(dts)) stop('vrt band count does not match dates')
    names(r) <- format(dts, '%Y-%m')
    r
}

assemble_layerwise <- function(v, fs, dts) {
    tmpd <- file.path(cube_dir, paste0('.tmp_', v))
    dir.create(tmpd, showWarnings = FALSE, recursive = TRUE)
    nm_target <- format(dts, '%Y-%m')
    nms <- lapply(fs, function(f) names(rast(f)))
    parts <- character(0)
    for (i in seq_along(nm_target)) {
        idx <- vapply(nms, function(n) {
            j <- which(substr(n, 1, 7) == nm_target[i])
            if (length(j)) j[1] else NA_integer_
        }, integer(1))
        use <- which(!is.na(idx))
        if (!length(use)) { cat('  no tile carries ', nm_target[i], '\n', sep = ''); next }
        sc <- sprc(lapply(use, function(k) rast(fs[k], lyrs = idx[k])))
        pf <- file.path(tmpd, sprintf('%s_%03d.tif', v, i))
        mosaic(sc, fun = 'mean', filename = pf, overwrite = TRUE,
               gdal = c('COMPRESS=DEFLATE','PREDICTOR=3'))
        parts <- c(parts, pf)
        if (i %% 24L == 0L) { invisible(gc()); cat('  layer', i, 'of', length(nm_target), '\n') }
    }
    r <- rast(parts)
    names(r) <- nm_target[seq_len(nlyr(r))]
    r
}

assemble_var <- function(v) {
    fs <- tile_files(v)
    if (!length(fs)) { cat('skip ', v, ': no tiles\n', sep = ''); return(invisible(NULL)) }
    cat('\n', v, ': ', length(fs), ' tiles\n', sep = '')

    overlap_check(v)

    dts <- cube_dates(fs[1])
    nl  <- vapply(fs, function(f) as.integer(nlyr(rast(f))), integer(1))
    r <- if (length(unique(nl)) == 1L && unique(nl) == length(dts))
             assemble_vrt(v, fs, dts) else assemble_layerwise(v, fs, dts)

    outf <- file.path(cube_dir, paste0(v, '.tif'))
    writeRaster(r, outf, overwrite = TRUE,
                gdal = c('COMPRESS=DEFLATE','PREDICTOR=3','TILED=YES',
                         'BIGTIFF=IF_SAFER'))
    ## GeoTIFF keeps band descriptions but not time(), so the dates are
    ## recovered from the names on read. 21_cube_qc.R does exactly that.
    rr <- rast(outf)
    names(rr) <- format(dts[seq_len(nlyr(rr))], '%Y-%m')
    time(rr)  <- dts[seq_len(nlyr(rr))]
    cat('  wrote ', outf, ' (', round(file.size(outf) / 1e9, 2), ' GB, ',
        nlyr(rr), ' layers, ', ncell(rr), ' cells)\n', sep = '')

    unlink(file.path(cube_dir, paste0('.tmp_', v)), recursive = TRUE)
    invisible(rr)
}

assemble_all <- function(vars = keep_vars) {
    hr('C. ASSEMBLE')
    for (v in vars) { assemble_var(v); invisible(gc()) }
    cat('\ncube written to ', cube_dir, '\n', sep = '')
}

## ===========================================================================
audit_run()
inventory()
if (do_assemble) assemble_all() else
    cat('\nrun again with the argument "assemble" to write the cube\n')
