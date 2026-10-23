## ---------------------------------------------------------------------------
## 49_saved_sector_full.R
##
## The per-cell correlations of Section 3.2 (Attribution of predictive error)
## and Table 5 (state-to-output row), over ALL cells of the saved sector
## 3PG_10m/x731866.017196386_y7132899.71658545.RData (run of 10 November 2025,
## initial-cube configuration). Same extraction as
## 07_within_object_link_V1.R, which sampled 800 cells (set.seed(1)) and gave
## the 0.901 / -0.901 / -0.331 of within_test.txt; this script takes every cell
## so that Section 3.2 rests on one run. Read-only on the RData.
##
## Writes results/saved_sector_full_cells.csv (one row per cell: initial stem
## and foliage biomass, 2017 GPP sum over months and species, 2017 mean LAI)
## and results/saved_sector_full.csv (the three correlations and n).
## First run 17 Oct 2026: n = 11 120; LAI vs initial foliage +0.8976, LAI vs
## GPP 2017 -0.8968, initial stem biomass vs GPP 2017 -0.3058. Needs about
## 18 GB of memory and five minutes.
##
##   Rscript 49_saved_sector_full.R      (from UVa_work/, T7 mounted)
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr) })
rd <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG/3PG_10m/x731866.017196386_y7132899.71658545.RData'
nm <- load(rd); x <- get(nm[1]); rm(list = nm); gc()
cat('cells in object:', length(x), '\n')
grab <- function(d, v, when = NULL, fun = sum) {
  s <- d[d$variable == v, ]
  if (!nrow(s)) return(NA_real_)
  if (!is.null(when)) s <- s[format(s$date, '%Y') == when, ]
  if (!nrow(s)) return(NA_real_)
  a <- tapply(s$value, s$date, sum, na.rm = TRUE)   # sum across species = stand total
  fun(a)
}
res <- lapply(seq_along(x), function(i) {
  d <- x[[i]]; d0 <- d[d$date == min(d$date), ]
  init <- function(v) { z <- d0$value[d0$variable == v]; if (length(z)) sum(z, na.rm = TRUE) else NA_real_ }
  data.frame(cell = names(x)[i], init_stem = init('biom_stem'), init_foliage = init('biom_foliage'),
             gpp_2017 = grab(d, 'gpp', '2017', sum), lai_2017 = grab(d, 'lai', '2017', mean))
})
r <- bind_rows(res)
write.csv(r, 'results/saved_sector_full_cells.csv', row.names = FALSE)
pairs <- list(c('lai_2017', 'init_foliage'), c('lai_2017', 'gpp_2017'), c('init_stem', 'gpp_2017'))
S <- do.call(rbind, lapply(pairs, function(p) {
  ok <- is.finite(r[[p[1]]]) & is.finite(r[[p[2]]])
  data.frame(x = p[1], y = p[2], n = sum(ok), r = cor(r[[p[1]]][ok], r[[p[2]]][ok]))
}))
print(S, row.names = FALSE, digits = 4)
write.csv(S, 'results/saved_sector_full.csv', row.names = FALSE)
