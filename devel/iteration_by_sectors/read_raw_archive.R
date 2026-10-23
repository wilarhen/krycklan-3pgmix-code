## ---------------------------------------------------------------------------
## read_raw_archive.R
##
## The archive written by polygon_grids_v02.r (save_rdata = TRUE) holds the
## complete per-cell model output. Any variable can be mapped from it without
## re-running the model.
##
## One compact data.table per sector, columns:
##   cell      integer raster cell index
##   date      integer YYYYMM
##   species   factor
##   group     factor
##   variable  factor
##   value     numeric
## ---------------------------------------------------------------------------

library(data.table); library(terra)

## qs2 is the current serialization package; qs is deprecated. This reader
## handles whichever format the archive was written in.
read_sector <- function(f) {
  if (grepl('\\.qs2$', f)) {
    if (!requireNamespace('qs2', quietly = TRUE))
      stop('archive is .qs2 but package qs2 is not installed')
    qs2::qs_read(f)
  } else if (grepl('\\.qs$', f)) {
    if (!requireNamespace('qs', quietly = TRUE))
      stop('archive is .qs but package qs is not installed')
    qs::qread(f)
  } else if (grepl('\\.rds$', f, ignore.case = TRUE)) {
    readRDS(f)
  } else {                                  # legacy .RData
    e <- new.env(); load(f, e)
    get(ls(e)[1], envir = e)
  }
}

base_path <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
raw_dir   <- file.path(base_path, '3PG_raw')

## ---- 1. what is in the archive ---------------------------------------------
ARCHIVE_PAT <- '\\.(qs2|qs|rds|RData)$'
files <- list.files(raw_dir, pattern = ARCHIVE_PAT, full.names = TRUE)
cat('sectors archived:', length(files), '\n')
cat('total size:', round(sum(file.size(files))/1e9, 2), 'GB\n')

dt <- read_sector(files[1])
cat('\nvariables available:\n')
print(dt[, .N, by = .(group, variable)][order(group, variable)])
cat('\ndropped as uninformative:',
    paste(attr(dt, 'dropped_variables'), collapse = ', '), '\n')

## ---- 2. build a raster for ANY variable, no re-run needed -------------------
## template must be the same grid the model ran on
make_cube <- function(variable_name, template,
                      agg = c('sum','mean'),
                      files = list.files(raw_dir, ARCHIVE_PAT, full.names = TRUE)) {
  agg <- match.arg(agg)
  f_agg <- if (agg == 'sum') sum else mean

  cells <- list()
  for (f in files) {
    d <- read_sector(f)
    d <- d[variable == variable_name]
    if (!nrow(d)) next
    ## collapse species: sum for fluxes and stocks, mean for state variables
    d <- d[, .(value = f_agg(value, na.rm = TRUE)), by = .(cell, date)]
    cells[[length(cells)+1L]] <- d
    rm(d)
  }
  if (!length(cells)) stop('variable not found: ', variable_name)
  all <- rbindlist(cells); rm(cells); gc()

  dates <- sort(unique(all$date))
  r <- rast(template, nlyrs = length(dates))
  for (i in seq_along(dates)) {
    v <- rep(NA_real_, ncell(template))
    sub <- all[date == dates[i]]
    v[sub$cell] <- sub$value
    r[[i]][] <- v
  }
  time(r) <- as.Date(paste0(substr(dates,1,4), '-', substr(dates,5,6), '-01'))
  names(r) <- format(time(r), '%Y-%m')
  r
}

## ---- 3. examples -----------------------------------------------------------
## template <- rast(file.path(base_path, 'biomass/BIO.tif'))  # the reference grid
##
## npp  <- make_cube('npp',  template, agg = 'sum')
## lai  <- make_cube('lai',  template, agg = 'sum')
## f_sw <- make_cube('f_sw', template, agg = 'mean')
##
## Annual NPP for 2017 in gC/m2/yr:
##   i <- which(format(time(npp), '%Y') == '2017')
##   npp_2017 <- sum(npp[[i]]) * 0.47 * 100

## ---- 4. quick check of what a variable looks like ---------------------------
peek <- function(variable_name, f = files[1]) {
  d <- read_sector(f)[variable == variable_name]
  if (!nrow(d)) { cat('not found\n'); return(invisible(NULL)) }
  cat(variable_name, ': ', nrow(d), ' rows, ',
      uniqueN(d$cell), ' cells, ', uniqueN(d$date), ' dates\n', sep = '')
  print(summary(d$value))
  invisible(d)
}
## peek('lai'); peek('f_sw'); peek('height')
