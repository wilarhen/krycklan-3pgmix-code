## ---------------------------------------------------------------------------
## 46_radiation_station_check.R
##
## The MODIS-derived radiation used by the model (MCD18A1, eight three-hourly
## estimates per day summed as rectangles; Section 2.4) against the incoming
## shortwave record of the SITES Stortjarn station inside the catchment
## (64.2615 N, 19.7635 E, Apogee PYR-P pyranometer, 30-min means), month by
## month, to test radiation as a source of the negative bias. Read-only on
## the inputs; outputs to results/.
##
##   results/radiation_station_check.csv   one row per month with >= 20
##       complete station days: station and MODIS monthly mean daily
##       radiation (MJ m-2 d-1), ratio MODIS / station
##   results/radiation_station_summary.csv  totals ratio, median ratio and
##       correlation for all months, April-September, and 2016-2018
##
## Station days count as complete with >= 46 of 48 half-hours; negative night
## values are set to zero. MODIS values are the catchment-parcel means written
## by 43_tier5_support.R (results/radiation_monthly.csv).
##
##   Rscript 46_radiation_station_check.R      (from UVa_work/, T7 mounted)
## ---------------------------------------------------------------------------

f <- file.path('/run/media/wilar/T7 Shield/slu_data_bases/3PG/climate_SW',
               'SITES_MET_SVB_STO_20160702-20230908_L2_30min.csv')
OUT <- 'results'

L <- readLines(f, warn = FALSE); h <- grep('^TIMESTAMP,', L)
x <- read.csv(text = L[h:length(L)], na.strings = c('NaN', 'NA'))
x$t <- as.POSIXct(x$TIMESTAMP, format = '%Y-%m-%d %H:%M', tz = 'Etc/GMT-1')   # file is UTC+1
x$day <- as.Date(format(x$t - 1, '%Y-%m-%d', tz = 'Etc/GMT-1'))               # stamps mark interval ends

d <- aggregate(SW_IN ~ day, x, na.action = na.pass,
               function(v) c(n = sum(is.finite(v)), MJ = mean(pmax(v, 0), na.rm = TRUE) * 86400 / 1e6))
d <- data.frame(day = d$day, n = d$SW_IN[, 'n'], MJ = d$SW_IN[, 'MJ'])
d <- d[d$n >= 46, ]
d$year <- as.integer(format(d$day, '%Y')); d$month <- as.integer(format(d$day, '%m'))

s <- aggregate(MJ ~ year + month, d, mean)
s$ndays <- aggregate(MJ ~ year + month, d, length)$MJ
m <- read.csv(file.path(OUT, 'radiation_monthly.csv'))
j <- merge(s, m[, c('year', 'month', 'mean')], by = c('year', 'month'))
names(j)[names(j) == 'MJ'] <- 'station'; names(j)[names(j) == 'mean'] <- 'modis'
j <- j[j$ndays >= 20, ]; j <- j[order(j$year, j$month), ]
j$ratio <- j$modis / j$station
write.csv(j, file.path(OUT, 'radiation_station_check.csv'), row.names = FALSE)

summ <- function(k, lab) data.frame(set = lab, months = sum(k),
    total_ratio = sum(j$modis[k]) / sum(j$station[k]), median_ratio = median(j$ratio[k]),
    r = cor(j$modis[k], j$station[k]))
S <- rbind(summ(rep(TRUE, nrow(j)), 'all months'),
           summ(j$month %in% 4:9, 'April-September'),
           summ(j$year %in% 2016:2018, '2016-2018'),
           summ(j$year %in% 2016:2018 & j$month %in% 4:9, '2016-2018, April-September'))
cat('months compared:', nrow(j), sprintf('(%d-%02d to %d-%02d)', j$year[1], j$month[1],
    tail(j$year, 1), tail(j$month, 1)), '\n')
print(S, row.names = FALSE, digits = 3)
write.csv(S, file.path(OUT, 'radiation_station_summary.csv'), row.names = FALSE)
