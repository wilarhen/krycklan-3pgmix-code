## ---------------------------------------------------------------------------
## 23_ceiling_multivariate.R
##
## FIX 1 of 3. The ceiling in the manuscript is a single-predictor correlation:
## measured GPP against measured LAI. 3-PG does not predict from LAI alone, so
## a reviewer will ask what the ceiling looks like when every measured stand
## attribute is admissible. If the multivariate ceiling for spruce is much
## higher than 0.268, the central claim does not hold.
##
## Selection is nested inside the cross-validation: for each held-out plot the
## subset is chosen on the remaining plots only. Choosing the subset once on
## all plots and then reporting its in-sample r is the standard way to
## manufacture a ceiling that does not exist, and with n = 14 it would be
## severe.
##
##   Rscript 23_ceiling_multivariate.R
## ---------------------------------------------------------------------------

source('ceiling_config.R')
set.seed(SEED)

d <- to_plot_level(load_plots())

## species actually present, in descending plot count. SPECIES_KEEP may be
## NULL in the config, in which case every label found is analysed.
SPP <- if (is.null(SPECIES_KEEP))
    names(sort(table(d$species), decreasing = TRUE)) else SPECIES_KEEP

cat('plots by species:\n'); print(table(d$species))

pnames <- names(PREDICTORS)

## ---- nested leave-one-out --------------------------------------------------
## returns the out-of-sample predictions of a size-capped exhaustive search
nested_loo <- function(dd, preds, kmax) {
    n <- nrow(dd)
    subsets <- unlist(lapply(1:min(kmax, length(preds)),
                             function(k) combn(preds, k, simplify = FALSE)),
                      recursive = FALSE)
    pred <- rep(NA_real_, n)
    chosen <- character(n)
    for (i in seq_len(n)) {
        tr <- dd[-i, ]; te <- dd[i, , drop = FALSE]
        ## inner LOO over the training set picks the subset
        err <- vapply(subsets, function(s) {
            f <- stats::as.formula(paste('gpp_obs ~', paste(s, collapse = '+')))
            e <- vapply(seq_len(nrow(tr)), function(j) {
                m <- try(stats::lm(f, data = tr[-j, ]), silent = TRUE)
                if (inherits(m, 'try-error')) return(NA_real_)
                tr$gpp_obs[j] - stats::predict(m, tr[j, , drop = FALSE])
            }, numeric(1))
            sqrt(mean(e^2, na.rm = TRUE))
        }, numeric(1))
        best <- subsets[[which.min(err)]]
        chosen[i] <- paste(best, collapse = '+')
        f <- stats::as.formula(paste('gpp_obs ~', paste(best, collapse = '+')))
        m <- stats::lm(f, data = tr)
        pred[i] <- stats::predict(m, te)
    }
    list(pred = pred, chosen = chosen)
}

res <- list()
for (sp in SPP) {
    dd <- d[d$species == sp, ]
    n  <- nrow(dd)
    cat('\n', strrep('-', 68), '\n', sp, '  n = ', n, ' plots\n', sep = '')
    if (n < 8) { cat('  too few plots, skipped\n'); next }

    ## the manuscript's ceiling: univariate, LAI only
    r_uni <- stats::cor(dd$gpp_obs, dd$lai)
    cat(sprintf('  univariate ceiling  r(GPP, LAI)          = %6.3f\n', r_uni))

    ## univariate ceiling, honest out-of-sample version
    p_uni <- vapply(seq_len(n), function(i)
        stats::predict(stats::lm(gpp_obs ~ lai, data = dd[-i, ]),
                       dd[i, , drop = FALSE]), numeric(1))
    r_uni_cv <- stats::cor(dd$gpp_obs, p_uni)
    cat(sprintf('  univariate, LOO out-of-sample            = %6.3f\n', r_uni_cv))

    ## multivariate ceiling with selection nested inside the CV
    kmax <- max(1, min(MAX_SUBSET, floor(n / 5)))
    nl <- nested_loo(dd, pnames, kmax)
    r_mv <- stats::cor(dd$gpp_obs, nl$pred)
    cat(sprintf('  multivariate, nested LOO (kmax = %d)      = %6.3f\n', kmax, r_mv))

    ## the overfitted version, reported only to show the size of the trap
    subsets <- unlist(lapply(1:kmax, function(k) combn(pnames, k, simplify = FALSE)),
                      recursive = FALSE)
    r_in <- max(vapply(subsets, function(s)
        sqrt(summary(stats::lm(
            stats::as.formula(paste('gpp_obs ~', paste(s, collapse = '+'))),
            data = dd))$r.squared), numeric(1)))
    cat(sprintf('  best in-sample fit (NOT a ceiling)       = %6.3f\n', r_in))

    cat('  subsets selected across folds:\n')
    print(sort(table(nl$chosen), decreasing = TRUE))

    ## per-predictor univariate correlations, for the Results table
    uni <- sort(vapply(pnames, function(p)
        stats::cor(dd$gpp_obs, dd[[p]]), numeric(1)), decreasing = TRUE)
    cat('  univariate r by measured attribute:\n')
    print(round(uni, 3))

    res[[sp]] <- data.frame(species = sp, n = n,
                            r_lai_insample = r_uni, r_lai_loo = r_uni_cv,
                            r_multivariate = r_mv, r_best_insample = r_in,
                            kmax = kmax)
}

tab <- do.call(rbind, res)
print(tab)
utils::write.csv(tab, file.path(path.expand(OUT_DIR), 'ceiling_multivariate.csv'),
                 row.names = FALSE)

cat('\n', strrep('=', 68), '\n', sep = '')
cat('READ THIS AS:\n')
cat('  If r_multivariate is close to r_lai_loo, the LAI-only ceiling in the\n')
cat('  manuscript is defensible and should be relabelled a STRUCTURAL ceiling.\n')
cat('  If r_multivariate is much higher, especially for spruce, the ceiling\n')
cat('  claim must be weakened: other measured attributes carry information\n')
cat('  that 3-PG could in principle use.\n')
cat('  r_best_insample is the number you would get by selecting once on all\n')
cat('  plots. It is not a ceiling. It is printed so the gap is visible.\n')
