## ---------------------------------------------------------------------------
## 22_check_srad.R
##
## Diagnostics on the radiation layers already on disk. Answers four questions
## before anything is rewritten:
##
##   Q1  how many DSR files are there, and how many bands each
##   Q2  is the stack 240 layers (2001-2020) or 250-252 (2000 included)
##   Q3  does radiation vary spatially inside the catchment, or is Table 1
##       right that it expresses none
##   Q4  are the values plausible for 64 N, and does the seasonal cycle line
##       up with temperature rather than lagging it
##
## Read-only. Nothing is rewritten here.
## ---------------------------------------------------------------------------

library(terra)

base_path <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
md        <- file.path(base_path, 'MODIS_DSR_Kricklan_solar_radiation')
kry_dir   <- file.path(base_path, 'krycklan_borders')

hr <- function(x) cat('\n', strrep('-', 72), '\n', x, '\n', sep = '')

## ---- Q1  file and band inventory -------------------------------------------
hr('Q1  FILES AND BANDS')
fs <- list.files(md, pattern = 'DSR.*\\.tif$', full.names = TRUE)
cat('DSR files: ', length(fs), '\n', sep = '')

nb <- vapply(fs, function(f) as.integer(nlyr(rast(f))), integer(1))
cat('bands per file: ', paste(sort(unique(nb)), collapse = ', '), '\n', sep = '')
if (any(nb > 1)) {
    cat('\n  Files carry more than one band. The v1 GEE script wrote\n',
        '  DSR_MJ_m2_day AND valid_days into each GeoTIFF. Confirm that the\n',
        '  pipeline reads band 1 only, or band 2 (a day count) is being fed\n',
        '  to 3-PG as MJ/m2/day.\n', sep = '')
    b <- rast(fs[1])
    cat('  band names in ', basename(fs[1]), ': ',
        paste(names(b), collapse = ', '), '\n', sep = '')
    print(global(b, c('mean'), na.rm = TRUE))
}

## ---- Q2  temporal span ------------------------------------------------------
hr('Q2  TEMPORAL SPAN')
tag <- sub('^DSR_', '', sub('\\.tif$', '', basename(fs)))
yr  <- substr(tag, 1, 4)
cat('years present: ', paste(range(yr), collapse = ' to '), '\n', sep = '')
print(table(yr))

n2000 <- sum(yr == '2000')
if (n2000 > 0) {
    cat('\n  ', n2000, ' file(s) from 2000 are present. The climate stack is\n',
        '  240 layers (2001-01 to 2020-12). d_input_from_rast_back joins by\n',
        '  layer POSITION and discards names, so an oversized radiation stack\n',
        '  can offset srad against temperature without any error being raised.\n',
        '  Verify that the radiation stack handed to it is exactly 240 layers\n',
        '  starting 2001-01, or move the 2000 files out of the folder.\n', sep = '')
} else {
    cat('\n  no 2000 files: the stack should align with the 240-layer climate\n',
        '  stack. Confirm the count below.\n', sep = '')
}
cat('expected for 2001-2020: 240 | found: ', sum(yr >= '2001'), '\n', sep = '')

## ---- Q2b  reproduce the stack the pipeline actually builds -----------------
## polygon_grids_v03.r does, in effect:
##   drs_rast <- import_feature_from_folder('rast', md, 'DSR')
##   names(drs_rast) <- dir(md)[grepl('DSR', dir(md))] %>% ... 'X2001.01' ...
##   time(...) <- names(...) %>% ym()
## then hands it to d_input_from_rast_back(from = '2001-01', to = '2020-12'),
## which is documented to join by layer POSITION and discard names. So the
## question is not what the files are called but which layer sits in slot 1.
hr('Q2b  THE STACK AS THE PIPELINE BUILDS IT')
ord   <- dir(md)[grepl('DSR', dir(md))]
ord   <- ord[grepl('\\.tif$', ord)]
nm    <- sub('\\.tif$', '', sub('^DSR_', 'X', ord))
nm    <- sub('_', '.', nm)
dts   <- as.Date(paste0(sub('^X', '', sub('\\.', '-', nm)), '-01'))

cat('folder order is alphabetical: ', identical(ord, sort(ord)), '\n', sep = '')
cat('dates monotonic:              ', !is.unsorted(dts), '\n', sep = '')
cat('layer 1 of the radiation stack: ', format(dts[1], '%Y-%m'), '\n', sep = '')
cat('layer 1 of the climate stack:   2001-01\n')

off <- length(dts[dts < as.Date('2001-01-01')])
if (off > 0) {
    cat('\n  ', off, ' layers precede 2001-01.\n', sep = '')
    cat('  IF d_input_from_rast_back subsets on time(), this is harmless.\n')
    cat('  IF it joins by position, srad is offset by ', off, ' months: the\n', sep = '')
    cat('  cell labelled 2001-01 would carry ', format(dts[1], '%Y-%m'),
        ' radiation.\n', sep = '')
    cat('  Testable: under a ', off, '-month offset the apparent radiation peak\n', sep = '')
    cat('  moves from ', format(dts[which(format(dts, '%m') == '07')[1]], '%m'),
        ' to month ', sprintf('%02d', ((7 - off - 1) %% 12) + 1), '.\n', sep = '')
    cat('  Q4 below reads the FILES, so it shows the true cycle. To test the\n')
    cat('  stack, compare Q4 peak against the f_* modifier seasonality in the\n')
    cat('  saved sector archive, which reflects what the model actually saw.\n')
}

## ---- Q3  spatial variation --------------------------------------------------
hr('Q3  SPATIAL VARIATION WITHIN THE CATCHMENT')
## Table 1 currently reports "none (aggregated)" for solar radiation. The
## pipeline applies align_to_ref with no aggregation step, unlike the climate
## layers, so this should be checked rather than assumed.
shp <- list.files(kry_dir, pattern = '\\.shp$', full.names = TRUE)
pick <- fs[grep('_07\\.tif$', fs)]                    # July, strongest signal
if (!length(pick)) pick <- fs[seq_len(min(3, length(fs)))]
pick <- pick[seq_len(min(3, length(pick)))]

for (f in pick) {
    r <- rast(f)[[1]]
    if (length(shp)) {
        b <- vect(shp[1])
        if (!same.crs(b, r)) b <- project(b, crs(r))
        r <- mask(crop(r, b), b)
    }
    s  <- as.numeric(global(r, 'sd',   na.rm = TRUE))
    mu <- as.numeric(global(r, 'mean', na.rm = TRUE))
    rg <- as.numeric(global(r, 'range', na.rm = TRUE))
    cat(sprintf('%-22s mean %6.2f  sd %6.3f  range %6.2f to %6.2f  (cv %.2f %%)\n',
                basename(f), mu, s, rg[1], rg[2], 100 * s / mu))
}
cat('\n  If sd is not ~0, radiation is a spatially varying driver and\n',
    '  Table 1 needs correcting. That is not merely a fix: it gives a\n',
    '  second varying driver besides structure, which bears on H1 and H5.\n', sep = '')

## ---- Q4  magnitude and seasonality -----------------------------------------
hr('Q4  MAGNITUDE AND SEASONAL CYCLE')
## Expected at 64 N: July roughly 15-22 MJ/m2/day, December near 0-1.
## A cycle that peaks outside June-July indicates a layer offset.
mon <- substr(tag, 6, 7)
mm  <- data.frame(month = sprintf('%02d', 1:12), mean = NA_real_, n = NA_integer_)
for (i in 1:12) {
    sel <- fs[mon == mm$month[i] & yr >= '2001']
    if (!length(sel)) next
    v <- vapply(sel, function(f)
        as.numeric(global(rast(f)[[1]], 'mean', na.rm = TRUE)), numeric(1))
    mm$mean[i] <- mean(v, na.rm = TRUE); mm$n[i] <- length(sel)
}
print(mm)
pk <- mm$month[which.max(mm$mean)]
cat('\npeak month: ', pk, '\n', sep = '')
if (!pk %in% c('06','07')) cat(
    '  Peak is outside June-July. Check for a layer offset between the\n',
    '  radiation stack and the climate stack before trusting the cube.\n', sep = '')

jul <- mm$mean[mm$month == '07']
if (is.finite(jul)) {
    cat(sprintf('July mean: %.2f MJ/m2/day (expect roughly 15-22 at 64 N)\n', jul))
    if (jul < 10)
        cat('  Low. Consistent with the inverted DSR_Quality mask in the v1\n',
            '  GEE script retaining only poor retrievals, or with monthly\n',
            '  means resting on very few days. Check the valid_days band.\n', sep = '')
    if (jul > 26)
        cat('  High. Check the W/m2 to MJ/m2/day conversion and whether a\n',
            '  scale factor was applied twice.\n', sep = '')
}

## ---- valid_days, if it was exported ----------------------------------------
hr('VALID DAY COUNT')
vd <- list.files(md, pattern = 'VALIDDAYS.*\\.tif$', full.names = TRUE)
if (!length(vd) && any(nb > 1)) {
    cat('valid_days is band 2 of the DSR files. Summarising band 2 for July:\n')
    for (f in pick) {
        r <- rast(f)
        if (nlyr(r) < 2) next
        cat(sprintf('  %-22s mean valid days %.1f\n', basename(f),
                    as.numeric(global(r[[2]], 'mean', na.rm = TRUE))))
    }
} else if (length(vd)) {
    for (f in vd[seq_len(min(3, length(vd)))])
        cat(sprintf('  %-24s mean valid days %.1f\n', basename(f),
                    as.numeric(global(rast(f), 'mean', na.rm = TRUE))))
} else {
    cat('no valid_days layer found. Re-export with the v2 script to obtain it;\n',
        'it is the direct measure of how many days each monthly mean rests on.\n', sep = '')
}
cat('\n  A July count well below ~30 means the monthly mean rests on a small\n',
    '  and probably non-random subset of days.\n', sep = '')
