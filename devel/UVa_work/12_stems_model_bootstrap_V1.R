#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 12_stems_model_bootstrap_V1.R
##
## PART A  A stem-density sub-model to replace the constant stems_n = 2000.
##         Only BIO.tif and the age rasters exist wall-to-wall, so the model
##         may use ONLY biomass and age as predictors. Basal area and measured
##         Trees are plot-only, so using them would be circular for the cube;
##         they appear here as a reference ceiling.
##         Self-thinning gives the log-log form a physical basis: density
##         declines as mean tree size increases.
##         Model choice is by LEAVE-ONE-OUT RMSE, not by in-sample R2, because
##         it will be applied to pixels that were never fitted.
##
## PART B  Bootstrap of the validation correlations, resampling PLOTS not
##         plot-years. Pine rests on 16 plots and the 3 years from one plot are
##         not independent observations, so row-wise resampling would report a
##         falsely narrow interval.
##
##   Rscript 12_stems_model_bootstrap_V1.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl); library(dplyr)
})

xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
out <- 'validation_out'; dir.create(out, showWarnings = FALSE)
cf  <- 0.47
set.seed(42)
sink(file.path(out, 'stems_bootstrap_log.txt'), split = TRUE)

## ============================ PART A =====================================
cat('\n===================== PART A: stem density model =====================\n')

rd <- function(sheet) {
  x <- suppressMessages(read_excel(xl, sheet = sheet, skip = 2,
                                   .name_repair = 'unique')) |> as.data.frame()
  ## CORRECTION: the '50 main plots' sheet header transposes SPM and DTS.
  ## The '7 extra old-stand plots' sheet has them the right way round, and
  ## the marginal counts settle it: position 9 is species (27/23 against the
  ## published 28/22) and position 10 is soil (16/34 against 15/35).
  ## Reading by name is therefore NOT safe for the main sheet - the names are
  ## what is wrong. Corrected here, at read time, for that sheet only.
  if (sheet == '50 main plots') {
    .i <- match(c('SPM', 'DTS'), names(x))
    if (!any(is.na(.i))) names(x)[.i] <- c('DTS', 'SPM')
  }
  x |>
    mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
           Year    = suppressWarnings(as.integer(Year))) |>
    filter(!is.na(Plot_ID), Year == 2016)
}

main <- rd('50 main plots')
old  <- rd('7 extra old-stand plots')
## NOTE: SPM and DTS appear in opposite column ORDER in the two sheets,
## and the '50 main plots' header is the transposed one, so reading by
## NAME does NOT resolve it. rd() corrects that sheet at read time. The main plots carry the primary fit; the old-stand plots are added
## only as a robustness check, so any ambiguity there cannot affect the result.

mk <- function(x, tag) x |>
  transmute(Plot_ID, set = tag, Trees, AGCtree, Age, BA,
            Dgv = if ('Dgv' %in% names(x)) Dgv else NA_real_,
            sp = ifelse(DTS == 1, 'Pine', 'Spruce')) |>
  mutate(agb = AGCtree / cf) |>
  filter(is.finite(Trees), Trees > 0, is.finite(agb), agb > 0,
         is.finite(Age), Age > 0)

d50 <- mk(main, 'main')
d57 <- bind_rows(d50, mk(old, 'oldstand'))
cat('plots usable: main =', nrow(d50), '| main+oldstand =', nrow(d57), '\n')
cat('Trees  : median', round(median(d50$Trees)), ' range',
    round(min(d50$Trees)), '-', round(max(d50$Trees)), '\n')
cat('agb    : median', round(median(d50$agb), 1), ' range',
    round(min(d50$agb), 1), '-', round(max(d50$agb), 1), 'tDM/ha\n')
cat('Age    : median', round(median(d50$Age)), ' range',
    round(min(d50$Age)), '-', round(max(d50$Age)), '\n')

## candidate models. APPLICABLE = uses only agb and Age (both available as
## rasters). REFERENCE = uses BA, which is plot-only and therefore circular.
forms <- list(
  'const (2000)'                = NULL,
  'sp mean'                     = Trees ~ sp,
  'agb'                         = log(Trees) ~ log(agb),
  'age'                         = log(Trees) ~ log(Age),
  'agb + age'                   = log(Trees) ~ log(agb) + log(Age),
  'agb + age + sp'              = log(Trees) ~ log(agb) + log(Age) + sp,
  'agb * age'                   = log(Trees) ~ log(agb) * log(Age),
  'agb + age + sp interact'     = log(Trees) ~ (log(agb) + log(Age)) * sp,
  'REFERENCE: BA + agb'         = log(Trees) ~ log(BA) + log(agb))

loo <- function(fo, dat) {
  ## leave-one-out prediction on the natural scale
  pr <- numeric(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    fit <- try(lm(fo, data = dat[-i, ]), silent = TRUE)
    if (inherits(fit, 'try-error')) { pr[i] <- NA; next }
    p <- try(predict(fit, newdata = dat[i, ]), silent = TRUE)
    if (inherits(p, 'try-error')) { pr[i] <- NA; next }
    pr[i] <- if (identical(all.vars(fo)[1], 'Trees') &&
                 !grepl('log\\(Trees\\)', deparse(fo)[1])) p else exp(p)
  }
  pr
}

cat('\n', sprintf('%-26s %8s %8s %8s %8s', 'model', 'R2', 'LOO_RMSE', 'LOO_MAE', 'LOO_r'), '\n')
res <- list()
for (nm in names(forms)) {
  fo <- forms[[nm]]
  if (is.null(fo)) {
    pr <- rep(2000, nrow(d50)); r2 <- NA
  } else {
    fit <- lm(fo, data = d50); r2 <- summary(fit)$r.squared
    pr  <- loo(fo, d50)
  }
  ok <- is.finite(pr) & is.finite(d50$Trees)
  e  <- pr[ok] - d50$Trees[ok]
  res[[nm]] <- data.frame(model = nm, R2 = r2,
                          RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
                          r = cor(pr[ok], d50$Trees[ok]))
  cat(sprintf('%-26s %8.3f %8.0f %8.0f %8.3f\n', nm,
              ifelse(is.na(r2), NA_real_, r2),
              sqrt(mean(e^2)), mean(abs(e)), cor(pr[ok], d50$Trees[ok])))
}
tabA <- bind_rows(res)

best_app <- tabA |> filter(!grepl('REFERENCE', model)) |>
  slice_min(RMSE, n = 1)
cat('\nbest APPLICABLE model:', best_app$model,
    '| LOO RMSE', round(best_app$RMSE), 'stems/ha\n')
cat('constant 2000 LOO RMSE:', round(tabA$RMSE[tabA$model == 'const (2000)']), '\n')
cat('improvement:', round(100 * (1 - best_app$RMSE /
      tabA$RMSE[tabA$model == 'const (2000)'])), '%\n')

## coefficients of the preferred applicable form
cat('\n--- log(Trees) ~ log(agb) + log(Age) ---\n')
m1 <- lm(log(Trees) ~ log(agb) + log(Age), data = d50)
print(summary(m1)$coefficients, digits = 3)
cat('R2 =', round(summary(m1)$r.squared, 3),
    '| residual SE (log) =', round(summary(m1)$sigma, 3), '\n')
cat('\nSelf-thinning expectation: the log(agb) coefficient should be NEGATIVE\n')
cat('(more biomass in fewer, larger stems) and log(Age) negative too.\n')

cat('\n--- same fit including the 7 old-stand plots ---\n')
m2 <- lm(log(Trees) ~ log(agb) + log(Age), data = d57)
print(summary(m2)$coefficients, digits = 3)
cat('R2 =', round(summary(m2)$r.squared, 3), '\n')

cat('\n--- by species ---\n')
for (s in unique(d50$sp)) {
  ds <- d50[d50$sp == s, ]
  if (nrow(ds) < 8) { cat(s, ': n =', nrow(ds), 'too few\n'); next }
  ms <- lm(log(Trees) ~ log(agb) + log(Age), data = ds)
  cat(sprintf('%-7s n=%2d  R2=%.3f  Trees median=%4.0f\n',
              s, nrow(ds), summary(ms)$r.squared, median(ds$Trees)))
}

write.csv(tabA, file.path(out, 'stems_model_comparison.csv'), row.names = FALSE)

## ============================ PART B =====================================
cat('\n\n===================== PART B: bootstrap by plot =====================\n')
f <- file.path(out, 'params_results.csv')
if (!file.exists(f)) {
  cat('params_results.csv not found - run script 10 first\n')
} else {
  r <- read.csv(f)
  B <- 2000

  boot_r <- function(dat, xcol, ycol) {
    plots <- unique(dat$Plot_ID)
    out <- numeric(B)
    for (b in seq_len(B)) {
      pick <- sample(plots, length(plots), replace = TRUE)
      idx  <- unlist(lapply(pick, function(p) which(dat$Plot_ID == p)))
      s <- dat[idx, ]
      out[b] <- suppressWarnings(cor(s[[xcol]], s[[ycol]],
                                     use = 'complete.obs'))
    }
    out[is.finite(out)]
  }

  cat('\nresampling PLOTS with replacement,', B, 'iterations\n')
  cat(sprintf('\n%-10s %-8s %5s %7s %8s %18s\n',
              'variant','species','plots','r','median','95% CI'))
  rows <- list()
  for (v in unique(r$variant)) {
    for (s in unique(r$sp)) {
      dat <- r[r$variant == v & r$sp == s, ]
      np <- n_distinct(dat$Plot_ID)
      if (np < 5) next
      pt <- cor(dat$sim_gpp, dat$obs_gpp, use = 'complete.obs')
      bb <- boot_r(dat, 'sim_gpp', 'obs_gpp')
      ci <- quantile(bb, c(.025, .975))
      cat(sprintf('%-10s %-8s %5d %7.3f %8.3f   [%+.3f, %+.3f]\n',
                  v, substr(s, 1, 8), np, pt, median(bb), ci[1], ci[2]))
      rows[[length(rows)+1]] <- data.frame(variant = v, sp = s, plots = np,
                                           r = pt, med = median(bb),
                                           lo = ci[1], hi = ci[2])
    }
  }
  tabB <- bind_rows(rows)

  ## the ceiling: how well does measured LAI predict measured GPP?
  cat('\n--- CEILING: r(observed GPP, observed LAI), bootstrapped ---\n')
  base <- r[r$variant == unique(r$variant)[1], ]
  for (s in unique(base$sp)) {
    dat <- base[base$sp == s, ]
    pt <- cor(dat$obs_gpp, dat$obs_lai, use = 'complete.obs')
    bb <- boot_r(dat, 'obs_gpp', 'obs_lai')
    ci <- quantile(bb, c(.025, .975))
    cat(sprintf('%-8s plots=%2d  r=%.3f  95%% CI [%+.3f, %+.3f]\n',
                substr(s, 1, 8), n_distinct(dat$Plot_ID), pt, ci[1], ci[2]))
  }

  cat('\nHOW TO READ THIS\n')
  cat('  Pine: if the lower CI bound stays well above 0, the result holds\n')
  cat('        despite resting on 16 plots.\n')
  cat('  Spruce: if the CI spans 0, spruce has no demonstrable skill at all,\n')
  cat('        which is a stronger statement than a low point estimate.\n')
  cat('  Compare each model r against its species CEILING above: pine was\n')
  cat('        0.641 vs a ceiling of 0.754, i.e. near-optimal.\n')
  write.csv(tabB, file.path(out, 'bootstrap_ci.csv'), row.names = FALSE)
}
sink()
