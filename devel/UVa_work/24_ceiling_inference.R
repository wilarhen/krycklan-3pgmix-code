## ---------------------------------------------------------------------------
## 24_ceiling_inference.R
##
## FIX 2 of 3. The manuscript reports a CI on model skill but not on the
## CEILING, and then compares the two as though they were independent numbers.
## They are not: both are correlations with the same measured GPP, computed on
## the same plots. Three things are needed.
##
##   A  bootstrap CI on each ceiling, resampling plots
##   B  Williams' test for two dependent correlations sharing one variable,
##      which is the correct test of "the model sits at its ceiling"
##   C  disattenuation sensitivity: how far the ceiling moves if measured GPP
##      and LAI carry measurement error
##
## The honest outcome for pine is likely "consistent with, not demonstrated".
## Saying so first is much stronger than having a reviewer say it.
##
##   Rscript 24_ceiling_inference.R
## ---------------------------------------------------------------------------

source('ceiling_config.R')
set.seed(SEED)

d <- to_plot_level(load_plots())

## species actually present, in descending plot count. SPECIES_KEEP may be
## NULL in the config, in which case every label found is analysed.
SPP <- if (is.null(SPECIES_KEEP))
    names(sort(table(d$species), decreasing = TRUE)) else SPECIES_KEEP


## ---- Williams (1959) test for dependent correlations -----------------------
## r_jk : model vs measured GPP        (skill)
## r_jh : LAI   vs measured GPP        (ceiling)
## r_kh : model vs LAI                 (the dependency)
williams <- function(r_jk, r_jh, r_kh, n) {
    R <- 1 - r_jk^2 - r_jh^2 - r_kh^2 + 2 * r_jk * r_jh * r_kh
    num <- (r_jk - r_jh) * sqrt((n - 1) * (1 + r_kh))
    den <- sqrt(2 * ((n - 1) / (n - 3)) * R +
                (((r_jk + r_jh) / 2)^2) * (1 - r_kh)^3)
    t <- num / den
    list(t = t, df = n - 3, p = 2 * stats::pt(-abs(t), n - 3))
}

boot_stat <- function(dd, f, B = N_BOOT) {
    n <- nrow(dd)
    vapply(seq_len(B), function(b) f(dd[sample.int(n, n, replace = TRUE), ]),
           numeric(1))
}
pci <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)

out <- list()
for (sp in SPP) {
    dd <- d[d$species == sp, ]
    n  <- nrow(dd)
    cat('\n', strrep('=', 68), '\n', sp, '  n = ', n, ' plots\n', sep = '')
    if (n < 8) { cat('  too few plots, skipped\n'); next }

    r_sim  <- stats::cor(dd$gpp_obs, dd$gpp_sim)   # model skill
    r_ceil <- stats::cor(dd$gpp_obs, dd$lai)       # ceiling
    r_dep  <- stats::cor(dd$gpp_sim, dd$lai)       # shared structure

    ## ---- A  bootstrap CIs, plots resampled ---------------------------------
    b_sim  <- boot_stat(dd, function(x) stats::cor(x$gpp_obs, x$gpp_sim))
    b_ceil <- boot_stat(dd, function(x) stats::cor(x$gpp_obs, x$lai))
    b_diff <- boot_stat(dd, function(x)
        stats::cor(x$gpp_obs, x$gpp_sim) - stats::cor(x$gpp_obs, x$lai))

    cat(sprintf('  model skill   r = %6.3f   95%% CI [%6.3f, %6.3f]\n',
                r_sim, pci(b_sim)[1], pci(b_sim)[2]))
    cat(sprintf('  ceiling       r = %6.3f   95%% CI [%6.3f, %6.3f]   <- NEW\n',
                r_ceil, pci(b_ceil)[1], pci(b_ceil)[2]))
    cat(sprintf('  difference        %6.3f   95%% CI [%6.3f, %6.3f]\n',
                r_sim - r_ceil, pci(b_diff)[1], pci(b_diff)[2]))
    covers0 <- pci(b_diff)[1] <= 0 && pci(b_diff)[2] >= 0
    cat('  the difference CI ', if (covers0) 'INCLUDES' else 'EXCLUDES',
        ' zero\n', sep = '')

    ## ---- B  Williams test --------------------------------------------------
    w <- williams(r_sim, r_ceil, r_dep, n)
    cat(sprintf("\n  Williams' test (dependent correlations, shared measured GPP)\n"))
    cat(sprintf('    r(model, LAI) = %6.3f   t(%d) = %6.3f   p = %.3f\n',
                r_dep, w$df, w$t, w$p))
    cat('    ', if (w$p > 0.05)
            'model skill is NOT distinguishable from the ceiling'
        else 'model skill DIFFERS from the ceiling', '\n', sep = '')

    ## ---- power: what n would be needed -------------------------------------
    z <- function(r) 0.5 * log((1 + r) / (1 - r))
    dz <- abs(z(r_ceil) - z(r_sim))
    n_need <- if (dz > 1e-8)
        ceiling((stats::qnorm(0.975) + stats::qnorm(0.80))^2 / dz^2 + 3) else NA
    cat(sprintf('    plots needed to resolve a gap this size at 80%% power: %s\n',
                ifelse(is.na(n_need), 'n/a', n_need)))

    ## ---- C  attenuation ----------------------------------------------------
    cat('\n  disattenuated ceiling, if measured GPP and LAI each have reliability:\n')
    for (rel in c(1.00, 0.95, 0.90, 0.85, 0.80)) {
        cat(sprintf('    %.2f  ->  %6.3f\n', rel, min(r_ceil / rel, 0.999)))
    }
    cat('  Measurement error depresses the ceiling. Unless reliability is\n')
    cat('  argued for explicitly, a reviewer can call the ceiling an artefact\n')
    cat('  of measurement precision.\n')

    out[[sp]] <- data.frame(
        species = sp, n = n,
        r_model = r_sim, model_lo = pci(b_sim)[1], model_hi = pci(b_sim)[2],
        r_ceiling = r_ceil, ceil_lo = pci(b_ceil)[1], ceil_hi = pci(b_ceil)[2],
        r_model_lai = r_dep,
        diff = r_sim - r_ceil, diff_lo = pci(b_diff)[1], diff_hi = pci(b_diff)[2],
        williams_t = w$t, williams_p = w$p, n_needed = n_need)
}

tab <- do.call(rbind, out)
print(tab, row.names = FALSE)
utils::write.csv(tab, file.path(path.expand(OUT_DIR), 'ceiling_inference.csv'),
                 row.names = FALSE)
cat('\nwritten to ', file.path(OUT_DIR, 'ceiling_inference.csv'), '\n', sep = '')
