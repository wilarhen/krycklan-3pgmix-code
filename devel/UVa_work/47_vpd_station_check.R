## ---------------------------------------------------------------------------
## 47_vpd_station_check.R
##
## Vapour pressure deficit as 3-PG derives it from the monthly temperature
## inputs of the cube against the same formula applied to the diurnal range,
## using the SITES Stortjarn station record inside the catchment (30-min air
## temperature, TA), to test whether the cube's construction inflates the
## derived VPD (Section 2.4, Climate). Read-only on the inputs; prints to the console.
##
## The PTHBV product gives daily mean temperature only, so the cube's monthly
## minimum and maximum are the monthly extremes of DAILY MEAN temperature. 3-PG
## expects the mean daily minimum and maximum. Both pairs are computed here
## from the station record and passed through the r3PG get_vpd formula,
##     VPD = (svp(Tmax) - svp(Tmin)) / 2,  svp(t) = 6.1078 exp(17.269 t / (237.3 + t))
## (mbar). Station days count as complete with >= 46 of 48 half-hours; months
## need >= 25 complete days. The ratio extremes-of-means / diurnal is reported
## for all months and for April-September.
##
##   Rscript 47_vpd_station_check.R      (from UVa_work/, T7 mounted)
## ---------------------------------------------------------------------------

f <- file.path('/run/media/wilar/T7 Shield/slu_data_bases/3PG/climate_SW',
               'SITES_MET_SVB_STO_20160702-20230908_L2_30min.csv')

L <- readLines(f, warn = FALSE); h <- grep('^TIMESTAMP,', L)
x <- read.csv(text = L[h:length(L)], na.strings = c('NaN', 'NA'))
x$t <- as.POSIXct(x$TIMESTAMP, format = '%Y-%m-%d %H:%M', tz = 'Etc/GMT-1')   # file is UTC+1
x$day <- as.Date(format(x$t - 1, '%Y-%m-%d', tz = 'Etc/GMT-1'))               # stamps mark interval ends

d <- aggregate(TA ~ day, x, na.action = na.pass,
               function(v) c(n = sum(is.finite(v)), mean = mean(v), min = min(v), max = max(v)))
d <- data.frame(day = d$day, d$TA); d <- d[d$n >= 46, ]
d$ym <- format(d$day, '%Y-%m'); d$mo <- as.integer(format(d$day, '%m'))

svp <- function(t) 6.10780 * exp(17.2690 * t / (237.30 + t))           # r3PG get_vpd, mbar
vpd <- function(tmin, tmax) (svp(tmax) - svp(tmin)) / 2
m <- do.call(rbind, lapply(split(d, d$ym), function(s) if (nrow(s) >= 25)
  data.frame(ym = s$ym[1], mo = s$mo[1], ext_min = min(s$mean), ext_max = max(s$mean),
             dmin = mean(s$min), dmax = mean(s$max))))
m$vpd_ext <- vpd(m$ext_min, m$ext_max); m$vpd_diurnal <- vpd(m$dmin, m$dmax)
m$ratio <- m$vpd_ext / m$vpd_diurnal
g <- m[m$mo %in% 4:9, ]

cat('months:', nrow(m), sprintf('(%s to %s)', m$ym[1], tail(m$ym, 1)),
    '| April-September months:', nrow(g), '\n')
cat(sprintf(paste0('VPD from monthly extremes of daily means / from mean daily min and max:\n',
  '  all months      median ratio %.2f (range %.2f-%.2f)\n',
  '  April-September median ratio %.2f (range %.2f-%.2f); ratio > 1 in %d of %d months\n'),
  median(m$ratio), min(m$ratio), max(m$ratio),
  median(g$ratio), min(g$ratio), max(g$ratio), sum(g$ratio > 1), nrow(g)))
print(round(aggregate(cbind(vpd_ext, vpd_diurnal, ratio) ~ mo, g, mean), 2))
