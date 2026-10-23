## ---------------------------------------------------------------------------
## 33_table6_predictors.R
##
## Rebuilds manuscript Table 6 (tab:predictors): correlation between measured
## tree GPP and measured stand and soil attributes, by dominant species.
##
## The published table had two faults:
##   1. grouping by the column labelled DTS in '50 main plots', which holds
##      soil parent material (header transposition, confirmed by the Zenodo
##      deposit for Peichl et al. 2023);
##   2. correlations on plot-years (48 and 102) while every other result in
##      the manuscript is on plot means.
## This version uses the corrected species column (labelled SPM), one row per
## plot (2016-2018 mean), and the 47 plots of the observed/obs preflight, so
## the LAI row equals the ceiling reported in Section 3.5. The old table is
## reproduced first as a check on the column mapping.
##
## Bootstrap resamples plots, B = 2000, percentile intervals, seed as in
## ceiling_config.R.
##
##   Rscript 33_table6_predictors.R      (from UVa_work/)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(readxl); library(dplyr) })
set.seed(20260824)
NB  <- 2000
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
OUT <- 'results'
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, 'table6_log.txt'), split = TRUE)

## attribute label -> workbook column, in the published row order
attr <- c('Stand age' = 'Age', 'Soil nitrogen stock' = 'SN', 'Basal area' = 'BA',
          'Leaf area index (max)' = 'LAImax', 'Soil nitrogen concentration' = 'SoilN',
          'Stem density' = 'Trees', 'Soil carbon' = 'SoilC', 'Elevation' = 'Elevation')

x <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2)) |>
  filter(!is.na(Plot_ID)) |> mutate(gpp_tree = GPP - GPPff)

## ---- check: the published table, old grouping, plot-years -------------------
cat('Published table reproduced (old grouping = column labelled DTS, plot-years):\n')
old <- sapply(1:2, function(s) {
  d <- x[x$DTS == s, ]
  sapply(attr, function(k) cor(as.numeric(d[[k]]), d$gpp_tree, use = 'complete.obs'))
})
colnames(old) <- c('pine (soil=sediment)', 'spruce (soil=till)')
print(round(old, 3))

## ---- rebuilt: corrected species, plot means, analysed plots -----------------
keep <- read.csv('validation_out/preflight_results.csv', stringsAsFactors = FALSE) |>
  filter(variant == 'observed/obs') |> distinct(Plot_ID, sp)

pm <- x |> group_by(Plot_ID) |>
  summarise(sp_code = first(SPM), gpp_tree = mean(gpp_tree, na.rm = TRUE),
            across(all_of(unname(attr)), ~ mean(as.numeric(.x), na.rm = TRUE)),
            .groups = 'drop') |>
  inner_join(keep, by = 'Plot_ID')

stopifnot(all((pm$sp_code == 1) == (pm$sp == 'Pinus sylvestris')))

boot_r <- function(dd, k) {
  o <- replicate(NB, { b <- dd[sample.int(nrow(dd), replace = TRUE), ]
                       suppressWarnings(cor(b[[k]], b$gpp_tree, use = 'complete.obs')) })
  quantile(o[is.finite(o)], c(.025, .975))
}

res <- list()
for (s in c('Pinus sylvestris', 'Picea abies')) {
  d <- pm[pm$sp == s, ]
  for (lab in names(attr)) {
    k  <- attr[[lab]]; ci <- boot_r(d, k)
    res[[length(res) + 1]] <- data.frame(attribute = lab, column = k, sp = s,
      n = sum(complete.cases(d[[k]], d$gpp_tree)),
      r = cor(d[[k]], d$gpp_tree, use = 'complete.obs'), lo = ci[1], hi = ci[2])
  }
}
R <- do.call(rbind, res); rownames(R) <- NULL
R$old_r <- c(old[, 1], old[, 2])
cat('\nRebuilt (corrected species, plot means, observed/obs plots):\n')
print(transform(R, r = round(r, 3), lo = round(lo, 3), hi = round(hi, 3),
                old_r = round(old_r, 3)), row.names = FALSE)

## sanity: the LAI row is the ceiling of Section 3.5
cat('\nLAI row (should equal the ceiling, 0.465 pine, 0.477 spruce):\n')
print(R[R$column == 'LAImax', c('sp', 'n', 'r')], row.names = FALSE)

write.csv(R, file.path(OUT, 'table6_corrected.csv'), row.names = FALSE)
sink()
