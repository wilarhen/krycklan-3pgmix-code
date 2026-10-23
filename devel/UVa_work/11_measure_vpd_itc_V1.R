#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 11_measure_vpd_itc_V1.R
##
## Measures the value of two proposed improvements before committing to either.
##
## PART 1  Is tmp_min/tmp_max built correctly?
##         polygon_grids_v01.r takes min() and max() of DAILY MEAN temperature
##         within each month. 3-PG expects the MEAN DAILY MIN and MEAN DAILY
##         MAX. If the source NetCDF is daily mean only, the diurnal range
##         cannot be recovered from it and VPD is being estimated from the
##         month's extreme days instead. f_vpd = 0.80 is currently the binding
##         constraint in the model, so this matters.
##
## PART 2  How much does VPD actually vary across Krycklan?
##         Elevation range is 164-321 m, so ~1 C of lapse and ~5-10% VPD
##         variation is expected. If confirmed, spatial VPD cannot rank plots
##         (observed GPP spans a factor of ten) and is a bias correction only.
##
## PART 3  Do the ITC crowns give usable stem density?
##         Tests ITC counts at the plots against measured Trees. At 1-2 pts/m2
##         suppressed stems are missed, so expect underestimation in dense
##         young stands.
##
##   Rscript 11_measure_vpd_itc_V1.R
## ---------------------------------------------------------------------------

## SPECIES COLUMN CORRECTED: this script reads only '50 main plots',
## whose header transposes SPM and DTS. The identifiers below have been
## exchanged so that species is read from the column holding species.
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(readxl); library(dplyr)
})

bp  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
out <- 'validation_out'; dir.create(out, showWarnings = FALSE)
sink(file.path(out, 'measure_log.txt'), split = TRUE)

svp <- function(T) 0.6108 * exp(17.27 * T / (T + 237.3))   # kPa

## =================== PART 1: what is in the climate NetCDF? ================
cat('\n=========== PART 1: temperature source ===========\n')
f_tas <- file.path(bp, 'climate_SW/SMHI_pthbv_tas_1961_2025_daily.nc')
if (file.exists(f_tas)) {
  nc <- nc_open(f_tas)
  cat('variables:', paste(names(nc$var), collapse = ', '), '\n')
  for (v in names(nc$var)) {
    a <- ncatt_get(nc, v)
    cat(sprintf('  %-12s units=%-12s long_name=%s\n', v,
                ifelse(is.null(a$units), '?', a$units),
                ifelse(is.null(a$long_name), '?', a$long_name)))
  }
  cat('dims:', paste(sapply(nc$dim, function(x) paste0(x$name,'=',x$len)),
                    collapse = ', '), '\n')
  nc_close(nc)
  cat('\nIF THE ONLY VARIABLE IS DAILY MEAN (tas):\n')
  cat('  true daily min/max are NOT recoverable, so 3-PG VPD is being computed\n')
  cat('  from the month extremes. Measured humidity is then the ONLY route to\n')
  cat('  a correct VPD -> see PART 1b.\n')
} else cat('NOT FOUND:', f_tas, '\n')

## ---- 1b: does the SITES station data carry humidity? ----------------------
cat('\n--- SITES Svartberget station files ---\n')
sf <- list.files(file.path(bp, 'climate_SW'), pattern = 'SITES.*\\.csv$',
                 full.names = TRUE)
for (f in sf) {
  h <- tryCatch(names(read.csv(f, nrows = 5)), error = function(e) NULL)
  if (is.null(h)) { cat(basename(f), ': unreadable\n'); next }
  hit <- grep('RH|humid|vpd|VPD|dew|Td|TA|AirTemp|Temp', h, value = TRUE,
              ignore.case = TRUE)
  cat('\n', basename(f), '\n  ncol:', length(h), '\n')
  cat('  candidate columns:', paste(head(hit, 15), collapse = ', '), '\n')
}
cat('\n  RH + air temperature at 30-min resolution would give a TRUE VPD:\n')
cat('      VPD = svp(Tair) * (1 - RH/100), averaged over daylight hours.\n')
cat('  Compare that against what 3-PG currently derives internally.\n')

## =================== PART 2: spatial VPD variation ========================
cat('\n=========== PART 2: does VPD vary across Krycklan? ===========\n')
kb <- tryCatch(vect(list.files(file.path(bp,'krycklan_borders'),
                               pattern='\\.shp$', full.names=TRUE)[1]),
               error = function(e) NULL)
if (!is.null(kb) && file.exists(f_tas)) {
  r <- tryCatch(rast(f_tas), error = function(e) NULL)
  if (!is.null(r)) {
    cat('NetCDF grid:', ncol(r), 'x', nrow(r), '| res', res(r), '\n')
    kbp <- project(kb, crs(r))
    rc  <- crop(r[[1]], ext(kbp) * 1.5)
    cat('grid cells covering Krycklan:', ncell(rc), '\n')
    if (ncell(rc) < 2)
      cat('  -> ONE cell. Climate is spatially uniform BY CONSTRUCTION and\n',
          '     no VPD field can vary within the catchment.\n', sep = '')
  }
}

## elevation-driven bound, computed from the DEM
dm <- rast(file.path(bp, 'Sweden_dem/dem.tif'))
kbd <- project(kb, crs(dm))
el  <- crop(dm[[1]], kbd)
rng <- unlist(global(el, c('min','max'), na.rm = TRUE))
cat('\nDEM elevation range in catchment:', round(rng, 1), 'm\n')
lapse <- 0.0065 * diff(rng)
cat('temperature difference at 6.5 C/km lapse:', round(lapse, 2), 'C\n')
for (Tb in c(10, 15, 20)) {
  v1 <- svp(Tb); v2 <- svp(Tb - lapse)
  cat(sprintf('  at %2d C base: svp varies %.4f -> %.4f kPa = %.1f%% spread\n',
              Tb, v2, v1, 100 * (v1 - v2) / v1))
}
cat('\n  Compare with the observed GPP range: 169-1712 gC/m2/yr (10x).\n')
cat('  A few percent of VPD variation cannot rank these plots.\n')
cat('  f_vpd = exp(-CoeffCond * VPD); at f_vpd = 0.80 a 10% VPD change moves\n')
cat('  f_vpd by ~0.02, i.e. ~2-3% in GPP.\n')

## =================== PART 3: ITC stem density =============================
cat('\n=========== PART 3: ITC crowns vs measured stem density ===========\n')
itc <- list.files(file.path(bp, 'trees'), pattern = '\\.tif$', full.names = TRUE)
cat('ITC tiles found:', length(itc), '\n')
if (length(itc)) {
  t1 <- rast(itc[1])
  cat('\nfirst tile:', basename(itc[1]), '\n')
  print(t1)
  cat('value summary:\n'); print(global(t1[[1]], c('min','max','mean'), na.rm = TRUE))
  fr <- tryCatch(head(freq(t1[[1]]), 10), error = function(e) NULL)
  if (!is.null(fr)) { cat('top values (are these crown IDs or a mask?):\n'); print(fr) }

  ## plots
  o <- suppressMessages(read_excel(xl, sheet='50 main plots', skip=2,
                                   .name_repair='unique')) |> as.data.frame() |>
    mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
           Year = suppressWarnings(as.integer(Year))) |>
    filter(!is.na(Plot_ID), Year == 2016) |>
    transmute(Plot_ID, Trees, BA, Dgv, Hgv,
              sp = ifelse(SPM == 1, 'Pine', 'Spruce'))
  cc <- suppressMessages(read_excel(xl, sheet='plot coordinates')) |>
    rename(Plot_ID=1, E=2, N=3) |>
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |>
    filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))
  d <- inner_join(o, cc, by = 'Plot_ID')

  ## Krycklan plots are ~10 m radius; use a 10 m buffer = 314 m2
  buf_r <- 10
  pv  <- vect(cbind(d$E, d$N), type='points', crs='EPSG:3006')
  bufv <- buffer(pv, buf_r)
  area_ha <- pi * buf_r^2 / 10000

  vr <- vrt(itc)   # mosaic the tiles virtually
  cat('\nvirtual mosaic:', ncol(vr), 'x', nrow(vr), '| res', res(vr), '\n')

  ex <- terra::extract(vr[[1]], bufv, ID = TRUE)
  names(ex)[2] <- 'v'
  cnt <- ex |> group_by(ID) |>
    summarise(n_cell = n(),
              n_valid = sum(!is.na(v)),
              n_unique = n_distinct(v[!is.na(v) & v > 0]),
              .groups = 'drop') |>
    mutate(Plot_ID = d$Plot_ID[ID],
           itc_per_ha = n_unique / area_ha)
  dd <- inner_join(d, cnt, by = 'Plot_ID')

  cat('\nplots with ITC coverage:', sum(dd$n_valid > 0), '/', nrow(dd), '\n')
  ok <- dd[dd$n_valid > 0 & is.finite(dd$Trees), ]
  if (nrow(ok) > 5) {
    cat(sprintf('\n  r(ITC stems/ha, measured Trees) = %+0.3f   n = %d\n',
                cor(ok$itc_per_ha, ok$Trees, use='complete.obs'), nrow(ok)))
    cat(sprintf('  r(ITC stems/ha, BA)             = %+0.3f\n',
                cor(ok$itc_per_ha, ok$BA, use='complete.obs')))
    cat(sprintf('  measured Trees : median %6.0f  range %.0f - %.0f\n',
                median(ok$Trees), min(ok$Trees), max(ok$Trees)))
    cat(sprintf('  ITC stems/ha   : median %6.0f  range %.0f - %.0f\n',
                median(ok$itc_per_ha), min(ok$itc_per_ha), max(ok$itc_per_ha)))
    cat(sprintf('  ratio ITC/measured median: %.2f\n',
                median(ok$itc_per_ha / ok$Trees)))
    cat('\n  A ratio well below 1 confirms suppressed stems are missed.\n')
    cat('  What matters is the CORRELATION: if ITC ranks density correctly it\n')
    cat('  is usable after a calibration factor, even if biased low.\n')
    write.csv(ok, file.path(out, 'itc_vs_measured.csv'), row.names = FALSE)
  } else {
    cat('  too few plots with ITC coverage to correlate.\n')
    cat('  Check whether the tiles cover the catchment and what the values mean.\n')
  }
} else cat('no ITC tiles found under', file.path(bp,'trees'), '\n')

cat('\n=========== DECISION RULES ===========\n')
cat('VPD: if the NetCDF has only daily mean tas, better VPD requires the SITES\n')
cat('     station humidity. Worth doing for BIAS (f_vpd 0.80 -> ~0.95, GPP +15-20%\n')
cat('     against a current bias of -392), NOT for spatial ranking.\n')
cat('ITC: if r(ITC, Trees) > 0.7 it replaces the constant 2000 and gives\n')
cat('     independent height/crown data for the 3-PGmix light sub-model.\n')
cat('     Pine is already at its ceiling (0.641 vs 0.754) so expect little\n')
cat('     correlation gain; the value is internal consistency.\n')
sink()
