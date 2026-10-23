## ---------------------------------------------------------------------------
## 26_plotlevel_stats.R
##
## Recomputes every correlation in the validation on the correct unit.
##
## 16_preflight_corrected_V1.R computes r on PLOT-YEARS (line 258, group_by
## variant+sp then cor() on the raw rows) while reporting n as
## n_distinct(Plot_ID). The point estimate and its stated n come from
## different units. The bootstrap at line 262 is NOT affected: it resamples
## unique Plot_ID and carries each plot's years along, which is a correct
## cluster bootstrap, so the published intervals stand.
##
## This script is standalone -- it reads preflight_results.csv and recomputes
## everything, so it can be run without touching the original. Once the
## numbers are confirmed, the fix inside 16_preflight is a single aggregate()
## before the group_by at line 256, plus the same for r_lai (236) and
## r_gpp_apar (238). boot_r stays as it is.
##
##   Rscript 26_plotlevel_stats.R
## ---------------------------------------------------------------------------

set.seed(20260824)
IN  <- 'validation_out/preflight_results.csv'
OUT <- 'validation_out/plotlevel'
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(IN, stringsAsFactors = FALSE)
NB <- 2000

## ---- aggregation to plots ---------------------------------------------------
## Three annual observations from one plot are not independent. Averaging them
## is the same reasoning the existing bootstrap already applies.
agg_plot <- function(x) {
    v <- intersect(c('sim_gpp','obs_gpp','lai','obs_lai','apar','height',
                     'alpha_c','f_sw','f_vpd','f_tmp','f_frost','f_nutr','f_age'),
                   names(x))
    a <- aggregate(x[v], by = list(Plot_ID = x$Plot_ID), FUN = mean, na.rm = TRUE)
    sp <- x[!duplicated(x$Plot_ID), c('Plot_ID','sp')]
    merge(a, sp, by = 'Plot_ID')
}

## ---- cluster bootstrap, unchanged in logic from boot_r ----------------------
## Resamples plots. On plot-level data that is an ordinary plot resample; the
## function is kept in this form so the two units are directly comparable.
boot_ci <- function(dat, f, B = NB) {
    p <- unique(dat$Plot_ID); o <- rep(NA_real_, B)
    for (b in seq_len(B)) {
        idx <- unlist(lapply(sample(p, length(p), TRUE),
                             function(x) which(dat$Plot_ID == x)))
        o[b] <- suppressWarnings(f(dat[idx, ]))
    }
    quantile(o[is.finite(o)], c(.025, .975), na.rm = TRUE)
}

r_of <- function(a, b) function(x) cor(x[[a]], x[[b]], use = 'complete.obs')

## ---- per variant and species ------------------------------------------------
res <- list()
for (v in unique(d$variant)) {
    dv <- d[d$variant == v, ]
    for (s in unique(dv$sp)) {
        x <- dv[dv$sp == s, ]
        if (length(unique(x$Plot_ID)) < 5) next
        a <- agg_plot(x)

        r_py <- cor(x$sim_gpp, x$obs_gpp, use = 'complete.obs')   # as published
        r_pl <- cor(a$sim_gpp, a$obs_gpp, use = 'complete.obs')   # corrected
        ci   <- boot_ci(a, r_of('sim_gpp','obs_gpp'))

        ceil_py <- cor(x$obs_gpp, x$obs_lai, use = 'complete.obs')
        ceil_pl <- cor(a$obs_gpp, a$obs_lai, use = 'complete.obs')
        cci     <- boot_ci(a, r_of('obs_gpp','obs_lai'))

        ## gap, bootstrapped jointly so the dependence is respected
        gci <- boot_ci(a, function(z)
            cor(z$sim_gpp, z$obs_gpp, use = 'complete.obs') -
            cor(z$obs_gpp, z$obs_lai, use = 'complete.obs'))

        res[[length(res) + 1L]] <- data.frame(
            variant = v, species = s,
            n_plotyears = nrow(x), n_plots = nrow(a),
            r_published = r_py, r_corrected = r_pl,
            r_lo = ci[1], r_hi = ci[2],
            ceiling_published = ceil_py, ceiling_corrected = ceil_pl,
            ceil_lo = cci[1], ceil_hi = cci[2],
            gap = r_pl - ceil_pl, gap_lo = gci[1], gap_hi = gci[2],
            bias = mean(a$sim_gpp - a$obs_gpp, na.rm = TRUE),
            RMSE = sqrt(mean((a$sim_gpp - a$obs_gpp)^2, na.rm = TRUE)))
    }
}
R <- do.call(rbind, res)

cat('\n== MODEL SKILL, published unit vs corrected unit ==\n')
print(R[, c('variant','species','n_plotyears','n_plots',
            'r_published','r_corrected','r_lo','r_hi')],
      digits = 3, row.names = FALSE)

cat('\n== CEILING, and the gap between model and ceiling ==\n')
print(R[, c('variant','species','n_plots','ceiling_published','ceiling_corrected',
            'ceil_lo','ceil_hi','gap','gap_lo','gap_hi')],
      digits = 3, row.names = FALSE)

cat('\n== BIAS AND RMSE on plot means ==\n')
print(R[, c('variant','species','n_plots','bias','RMSE')],
      digits = 4, row.names = FALSE)

## ---- whole-sample diagnostics that also used plot-years ---------------------
cat('\n== r_lai and r_gpp_apar, both units (lines 236 and 238) ==\n')
w <- do.call(rbind, lapply(unique(d$variant), function(v) {
    x <- d[d$variant == v, ]; a <- agg_plot(x)
    data.frame(variant = v,
               r_lai_py   = cor(x$lai, x$obs_lai, use = 'complete.obs'),
               r_lai_pl   = cor(a$lai, a$obs_lai, use = 'complete.obs'),
               r_apar_py  = cor(x$sim_gpp, x$apar, use = 'complete.obs'),
               r_apar_pl  = cor(a$sim_gpp, a$apar, use = 'complete.obs'),
               lai_med    = median(a$lai, na.rm = TRUE),
               lai_max    = max(a$lai, na.rm = TRUE))
}))
print(w, digits = 3, row.names = FALSE)

write.csv(R, file.path(OUT, 'plotlevel_stats.csv'), row.names = FALSE)
write.csv(w, file.path(OUT, 'plotlevel_diagnostics.csv'), row.names = FALSE)
cat('\nwritten to ', OUT, '\n', sep = '')

cat('\nNOTE the acceptance test "r(GPP, APAR) > 0.7" should be read off the\n')
cat('plot-year column if the manuscript reports it that way, or restated on\n')
cat('plot means. Pick one unit and use it everywhere.\n')

## ---- variant A of Table 7 (added 2 Oct 2026) ----------------------------------
## Table 7's variant A (foliage corrected, original height form and parameter
## sets) is script 10's A_control and is not in preflight_results.csv. It is
## processed here, after every variant above, so that the random sequence and
## therefore every interval written above are unchanged, and the manuscript
## quotes one bootstrap run for all of Table 7. Written to its own file.
pa <- read.csv('validation_out/params_results.csv', stringsAsFactors = FALSE)
pa <- pa[pa$variant == 'A_control', ]
resA <- list()
for (s in unique(pa$sp)) {
    x <- pa[pa$sp == s, ]; a <- agg_plot(x)
    r_pl <- cor(a$sim_gpp, a$obs_gpp, use = 'complete.obs')
    ci   <- boot_ci(a, r_of('sim_gpp','obs_gpp'))
    resA[[length(resA) + 1L]] <- data.frame(variant = 'A_control', species = s,
        n_plotyears = nrow(x), n_plots = nrow(a), r_corrected = r_pl,
        r_lo = ci[1], r_hi = ci[2],
        bias = mean(a$sim_gpp - a$obs_gpp, na.rm = TRUE),
        RMSE = sqrt(mean((a$sim_gpp - a$obs_gpp)^2, na.rm = TRUE)))
}
RA <- do.call(rbind, resA)
cat('\n== Table 7 variant A (script 10 A_control), plot means ==\n')
print(RA, digits = 3, row.names = FALSE)
write.csv(RA, file.path(OUT, 'plotlevel_stats_variantA.csv'), row.names = FALSE)
