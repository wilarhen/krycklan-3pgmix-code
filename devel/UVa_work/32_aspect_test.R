## ---------------------------------------------------------------------------
## 32_aspect_test.R
##
## Does measured spruce tree GPP vary with aspect?
##
## Hypothesis tested: spring recovery of
## photosynthesis in Norway spruce is governed by local thermal regime, so
## that plots may differ systematically by cold-air pooling and by earlier
## snow-free dates on south-facing slopes. Elevation (r = 0.099) and
## depth-to-water (r = -0.152) were already null for spruce; aspect had not
## been tested.
##
## Aspect is circular, so the primary test uses its cosine (N-S, +1 north)
## and sine (E-W, +1 east), taken from the public deposit for Peichl et al.
## (2023), Zenodo 10.5281/zenodo.7312986 (CC-BY 4.0), file
## validation_data/Peichl2023_zenodo7312986_NEP_deposit.xlsx. The workbook
## carries aspect only as four cardinal classes; the numeric code Aspect_2
## (1 N, 2 E, 3 S, 4 W) is reported as requested but is not an ordered
## variable, so its r has no physical reading.
##
## Plots: the spruce plots of the observed/obs preflight variant (species
## from the corrected column), so the set is the one behind the ceiling
## table. Response: measured tree GPP, GPP - GPPff, averaged 2016-2018.
## Bootstrap resamples plots, B = 2000, percentile intervals, seed as in
## ceiling_config.R. Pine is reported alongside for comparison.
##
##   Rscript 32_aspect_test.R      (from UVa_work/)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(readxl); library(dplyr) })
set.seed(20260824)
NB  <- 2000
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
zen <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Peichl2023_zenodo7312986_NEP_deposit.xlsx'
OUT <- 'validation_out/aspect'
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, 'aspect_log.txt'), split = TRUE)

## ---- data ------------------------------------------------------------------
wb <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2)) |>
  filter(!is.na(Plot_ID)) |>
  group_by(Plot_ID) |>
  summarise(gpp_tree = mean(GPP - GPPff, na.rm = TRUE),
            aspect_card = first(Aspect_1), aspect_code = first(Aspect_2),
            Elevation = first(Elevation), .groups = 'drop')

z <- read_excel(zen, sheet = 'data')
z <- data.frame(Plot_ID = z[['Stand ID']], dep_species = z[['Tree species']],
                ns = z[['N–S']], ew = z[['E–W']])

pf <- read.csv('validation_out/preflight_results.csv', stringsAsFactors = FALSE) |>
  filter(variant == 'observed/obs') |> distinct(Plot_ID, sp)

d <- pf |> inner_join(wb, by = 'Plot_ID') |> inner_join(z, by = 'Plot_ID') |>
  mutate(aspect_deg = (atan2(ew, ns) * 180 / pi) %% 360)

cat('Plots joined:', nrow(d), ' (preflight observed/obs:', nrow(pf), ')\n')
cat('\nCorrected species against the published deposit:\n')
print(table(model = d$sp, deposit = d$dep_species))

## ---- statistics ------------------------------------------------------------
boot <- function(dd, f) {
  o <- replicate(NB, suppressWarnings(f(dd[sample.int(nrow(dd), replace = TRUE), ])))
  quantile(o[is.finite(o)], c(.025, .975), na.rm = TRUE)
}
r_of  <- function(x) function(dd) cor(dd[[x]], dd$gpp_tree)
r2_of <- function(dd) summary(lm(gpp_tree ~ ns + ew, data = dd))$r.squared

res <- list()
for (s in sort(unique(d$sp))) {
  ds <- d[d$sp == s, ]
  for (x in c('ns', 'ew', 'aspect_code', 'Elevation')) {
    ci <- boot(ds, r_of(x))
    ct <- cor.test(ds[[x]], ds$gpp_tree)
    res[[length(res) + 1]] <- data.frame(sp = s, n = nrow(ds), predictor = x,
      stat = 'r', est = unname(ct$estimate), lo = ci[1], hi = ci[2], p = ct$p.value)
  }
  m  <- lm(gpp_tree ~ ns + ew, data = ds); fs <- summary(m)$fstatistic
  ci <- boot(ds, r2_of)
  res[[length(res) + 1]] <- data.frame(sp = s, n = nrow(ds), predictor = 'ns + ew',
    stat = 'R2', est = summary(m)$r.squared, lo = ci[1], hi = ci[2],
    p = pf(fs[1], fs[2], fs[3], lower.tail = FALSE))
  cat('\n----', s, ' n =', nrow(ds), '----\n'); print(round(coef(summary(m)), 3))
  cat('\nMean tree GPP by cardinal class:\n')
  print(as.data.frame(ds |> group_by(aspect_card) |>
          summarise(n = n(), gpp_tree = round(mean(gpp_tree)), .groups = 'drop')))
}
R <- do.call(rbind, res); rownames(R) <- NULL
cat('\n==== Summary (bootstrap 95 % CI, plots resampled, B =', NB, ') ====\n')
print(transform(R, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3),
                p = signif(p, 2)), row.names = FALSE)

write.csv(R, file.path(OUT, 'aspect_test.csv'), row.names = FALSE)
write.csv(d, file.path(OUT, 'aspect_plots.csv'), row.names = FALSE)
sink()
