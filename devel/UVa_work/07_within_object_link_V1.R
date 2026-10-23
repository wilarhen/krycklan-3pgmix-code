#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 07_within_object_link_V1.R
##
## The saved sector object is the run_3PG OUTPUT (inps), keyed by integer cell
## index, holding both the INITIAL STATE and the SIMULATED GPP for every cell.
## That lets us split the two remaining hypotheses without touching a raster:
##
##   GPP correlates with initial biom_stem  -> model is fine; the corruption is
##                                             in timeSeriesRaster() writing
##                                             values into the wrong cells.
##   GPP does NOT correlate                 -> 3-PG discards the initial state;
##                                             the raster writer is innocent.
##
## Also reports the mean physiological modifiers, which name the limiting factor.
##
##   Rscript 07_within_object_link_V1.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr) })

rd <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG/3PG_10m/x731866.017196386_y7132899.71658545.RData'
n_sample <- 800          # cells to sample; raise if fast enough

nm <- load(rd)
x  <- get(nm[1])
cat('cells in object :', length(x), '\n')
cat('index range     :', min(as.integer(names(x))), '-',
    max(as.integer(names(x))), '\n')
cat('gaps (failed)   :', max(as.integer(names(x))) - min(as.integer(names(x))) +
      1 - length(x), '\n\n')

## what is actually in one cell
e1 <- x[[1]]
cat('groups   :', paste(unique(e1$group), collapse = ', '), '\n')
cat('n vars   :', length(unique(e1$variable)), '\n')
cat('variables:', paste(head(sort(unique(e1$variable)), 40), collapse = ', '), '\n\n')

set.seed(1)
idx <- sample(length(x), min(n_sample, length(x)))

## ---- per-cell summary ------------------------------------------------------
grab <- function(d, v, when = NULL, fun = sum) {
  s <- d[d$variable == v, ]
  if (!nrow(s)) return(NA_real_)
  if (!is.null(when)) s <- s[format(s$date, '%Y') == when, ]
  if (!nrow(s)) return(NA_real_)
  ## sum across species = stand total (NOT mean, which is what line 236 does)
  a <- tapply(s$value, s$date, sum, na.rm = TRUE)
  fun(a)
}

res <- lapply(idx, function(i) {
  d <- x[[i]]
  d0 <- d[d$date == min(d$date), ]
  init <- function(v) { z <- d0$value[d0$variable == v]
                        if (length(z)) sum(z, na.rm = TRUE) else NA_real_ }
  data.frame(
    cell        = names(x)[i],
    n_species   = length(unique(d$species)),
    init_stem   = init('biom_stem'),
    init_foliage= init('biom_foliage'),
    init_stems_n= init('stems_n'),
    gpp_2017    = grab(d, 'gpp', '2017', sum),
    npp_2017    = grab(d, 'npp', '2017', sum),
    lai_2017    = grab(d, 'lai', '2017', mean),
    f_sw        = grab(d, 'f_sw',    '2017', mean),
    f_vpd       = grab(d, 'f_vpd',   '2017', mean),
    f_tmp       = grab(d, 'f_tmp',   '2017', mean),
    f_frost     = grab(d, 'f_frost', '2017', mean),
    f_nutr      = grab(d, 'f_nutr',  '2017', mean),
    f_age       = grab(d, 'f_age',   '2017', mean),
    stringsAsFactors = FALSE)
})
r <- bind_rows(res)

cat('=============== THE DECISIVE CORRELATIONS ===============\n')
v <- c('init_stem','init_foliage','init_stems_n','lai_2017','gpp_2017','npp_2017')
v <- v[sapply(v, function(k) sum(is.finite(r[[k]])) > 10)]
print(round(cor(r[, v], use = 'pairwise.complete.obs'), 3))

if (all(c('init_stem','gpp_2017') %in% v)) {
  rr <- cor(r$init_stem, r$gpp_2017, use = 'complete.obs')
  cat(sprintf('\n  cor(initial biom_stem, GPP 2017) = %+0.3f\n', rr))
  if (abs(rr) > 0.5) {
    cat('  >> MODEL IS FINE. The signal exists inside run_3PG output.\n')
    cat('     The corruption is in timeSeriesRaster() / the cell-index\n')
    cat('     mapping in polygon_grids_v01.r lines 206-251.\n')
  } else {
    cat('  >> 3-PG ITSELF loses the initial state. The raster writer is\n')
    cat('     innocent. Investigate initialisation and the 2001 start date.\n')
  }
}

cat('\n=============== LIMITING FACTOR (2017 means) ============\n')
mods <- c('f_sw','f_vpd','f_tmp','f_frost','f_nutr','f_age')
mods <- mods[sapply(mods, function(k) any(is.finite(r[[k]])))]
for (m in mods)
  cat(sprintf('  %-8s mean %.3f   sd %.3f   min %.3f\n', m,
              mean(r[[m]], na.rm = TRUE), sd(r[[m]], na.rm = TRUE),
              min(r[[m]], na.rm = TRUE)))
cat('\n  Lowest mean = the factor throttling growth.\n')
cat('  sd ~ 0 = that factor is spatially uniform and cannot create pattern.\n')

cat('\n=============== SPREAD OF SIMULATED GPP =================\n')
print(summary(r$gpp_2017))
cat('n_species distribution:\n'); print(table(r$n_species))
cat('\nGPP by number of species (tests the mean() vs sum() issue):\n')
print(r |> group_by(n_species) |>
        summarise(n = n(), gpp = mean(gpp_2017, na.rm = TRUE),
                  stem = mean(init_stem, na.rm = TRUE), .groups = 'drop'))

write.csv(r, 'validation_out/within_object.csv', row.names = FALSE)
cat('\nwrote validation_out/within_object.csv\n')
