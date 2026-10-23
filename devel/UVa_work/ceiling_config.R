## ---------------------------------------------------------------------------
## ceiling_config.R  --  wired to validation_out
##
## The measured attributes needed for the ceiling are spread over three files,
## so this loader joins them:
##
##   preflight_results.csv   plot-year, per VARIANT: sim_gpp, obs_gpp, obs_lai
##   itc_vs_measured.csv     plot-level field measurements: Trees, BA, Dgv, Hgv
##   validation_pairs.csv    plot-year: Age, LAImax, and terrain indices
##
## ONE thing to set by hand: VARIANT_KEEP. Run the config once and it will
## print the variant labels it found, then set it and rerun.
## ---------------------------------------------------------------------------

## validation_out sits next to these scripts, so a relative path is used and
## the config survives the project being moved or renamed. The fallbacks cover
## being run from the project root or from a sibling directory.
VAL_DIR <- local({
    cand <- c('validation_out',
              'UVa_work/validation_out',
              '../UVa_work/validation_out',
              '~/Documents/SLU/3PG/devel/UVa_work/validation_out')
    hit <- cand[file.exists(file.path(path.expand(cand), 'preflight_results.csv'))]
    if (!length(hit)) stop('validation_out not found. Tried:\n  ',
                           paste(cand, collapse = '\n  '),
                           '\n  set VAL_DIR by hand in ceiling_config.R')
    hit[1]
})
OUT_DIR <- file.path(VAL_DIR, 'ceiling')

F_PREFLIGHT <- file.path(VAL_DIR, 'preflight_results.csv')
F_ITC       <- file.path(VAL_DIR, 'itc_vs_measured.csv')
F_PAIRS     <- file.path(VAL_DIR, 'validation_pairs.csv')

## ---- the one manual setting -------------------------------------------------
## preflight_results.csv holds several configuration variants per plot-year.
## The ceiling must be computed against ONE variant: the corrected one.
VARIANT_KEEP <- 'observed/obs'

## ---- column mapping ---------------------------------------------------------
COL <- c(
  plot    = 'Plot_ID',
  year    = 'Year',
  species = 'sp',
  gpp_obs = 'obs_gpp',
  gpp_sim = 'sim_gpp'
)

## MEASURED stand attributes only. If a mapped raster value gets in here the
## ceiling stops being input-independent and the whole argument collapses.
##   obs_lai  measured LAI, the manuscript's ceiling predictor
##   Trees    measured stem count (itc_per_ha is the ALS-derived one -- excluded)
##   BA       measured basal area
##   Dgv/Hgv  measured mean diameter and height
##   Age      measured stand age
PREDICTORS <- c(
  lai    = 'obs_lai',
  stems  = 'Trees',
  basal  = 'BA',
  dbh    = 'Dgv',
  height = 'Hgv'
# , age  = 'Age'   # drops 6 plots, 3 of 14 pine; script 25 does not use it
)

## Optional second run: terrain and site indices. These are not stand
## structure, so including them tests a BROADER ceiling -- "can any measured
## attribute predict production" rather than "can structure". Running both is
## worthwhile: the broader one is the more generous test, and if spruce stays
## low under it the result is much harder to argue with.
PREDICTORS_EXTENDED <- c(PREDICTORS,
  elev = 'Elevation', slope = 'Slope', spm = 'SPM', dts = 'DTS', sac = 'SAC')
USE_EXTENDED <- FALSE

## ---- analysis settings ------------------------------------------------------
SPECIES_KEEP <- NULL          # NULL = keep every label found, then inspect
MAX_SUBSET   <- 3
N_BOOT       <- 2000
SEED         <- 20260824

## Beer-Lambert extinction coefficients of the parameter sets used in the
## corrected configuration, Forrester et al. 2021 (r3PG database): Picea abies
## set 37, k = 0.6378; Pinus sylvestris set 39, k = 0.3247. Pine was 0.380
## until 2 Oct 2026, which is set 30 (Forrester et al. 2017), the original
## configuration.
## Names must match the labels in the 'sp' column exactly.
K_EXT <- c('Picea abies' = 0.638, 'Pinus sylvestris' = 0.3247)

## ---- loader -----------------------------------------------------------------
rd <- function(p) {
    p <- path.expand(p)
    if (!file.exists(p)) stop('not found: ', p)
    utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
}

load_plots <- function(verbose = TRUE) {
    pf <- rd(F_PREFLIGHT)

    ## --- variant ------------------------------------------------------------
    if ('variant' %in% names(pf)) {
        lv <- sort(unique(pf$variant))
        if (is.null(VARIANT_KEEP)) {
            cat('\npreflight_results.csv variants and their row counts:\n')
            print(table(pf$variant))
            stop('set VARIANT_KEEP in ceiling_config.R to the corrected variant')
        }
        if (!VARIANT_KEEP %in% lv)
            stop('VARIANT_KEEP = "', VARIANT_KEEP, '" not among: ',
                 paste(lv, collapse = ', '))
        pf <- pf[pf$variant == VARIANT_KEEP, ]
    }

    ## --- joins --------------------------------------------------------------
    itc <- rd(F_ITC)
    itc <- itc[, intersect(c('Plot_ID', 'Trees', 'BA', 'Dgv', 'Hgv'), names(itc))]

    prs <- rd(F_PAIRS)
    prs <- prs[, intersect(c('Plot_ID', 'Year', 'Age', 'LAImax', 'Elevation',
                             'Slope', 'SPM', 'DTS', 'SAC'), names(prs))]

    d <- merge(pf, itc, by = 'Plot_ID', all.x = TRUE)
    d <- merge(d, prs, by = c('Plot_ID', 'Year'), all.x = TRUE)

    preds <- if (USE_EXTENDED) PREDICTORS_EXTENDED else PREDICTORS
    want  <- c(COL, preds)
    miss  <- want[!want %in% names(d)]
    if (length(miss)) {
        cat('\ncolumns available after the joins:\n  ',
            paste(sort(names(d)), collapse = ', '), '\n\n', sep = '')
        stop('mapped columns absent: ', paste(miss, collapse = ', '))
    }

    out <- d[, want]; names(out) <- names(want)
    out$species <- as.character(out$species)

    if (verbose) {
        cat('species labels in "', COL['species'], '":\n', sep = '')
        print(table(out$species))
        bad <- setdiff(unique(out$species), names(K_EXT))
        if (length(bad))
            cat('  NOTE no K_EXT entry for: ', paste(bad, collapse = ', '),
                ' -- script 25 will skip these\n', sep = '')
    }
    if (!is.null(SPECIES_KEEP)) out <- out[out$species %in% SPECIES_KEEP, ]
    out
}

## Correlations are computed on plots, not plot-years: the three annual
## observations from a plot are not independent, the same reasoning the
## manuscript already uses for its bootstrap.
to_plot_level <- function(d) {
    num <- setdiff(names(d), c('plot', 'species', 'year'))
    agg <- stats::aggregate(d[num], by = list(plot = d$plot, species = d$species),
                            FUN = function(x) mean(x, na.rm = TRUE))
    keep <- stats::complete.cases(agg)
    if (any(!keep))
        cat('dropped ', sum(!keep), ' plots with incomplete predictors\n', sep = '')
    agg[keep, ]
}

dir.create(path.expand(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
