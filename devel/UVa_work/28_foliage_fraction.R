## ---------------------------------------------------------------------------
## 28_foliage_fraction.R
##
## The pipeline uses CONSTANT foliage fractions of aboveground biomass:
## 0.045 pine, 0.085 spruce, 0.020 birch. The cube shows these are wrong in
## OPPOSITE directions at the plots -- simulated LAI is 1.21x measured for
## spruce and 0.78x for pine -- so a single mis-specified constant cannot
## explain it. Both are poorly chosen for the local diameter range.
##
## This matters beyond QC. The paper's mechanism runs
##     foliage fraction -> LAI -> saturation -> elasticity -> ceiling
## so an inflated spruce foliage fraction inflates the very effect the
## conclusion depends on.
##
## NO ALLOMETRIC COEFFICIENTS ARE HARD-CODED. Reproducing them from memory is
## how the wrong constants arrived in the first place. The fraction is derived
## from measurements instead:
##
##   r3PG computes LAI = biom_foliage * SLA, so the effective SLA of each
##   simulated plot is sim_lai / sim_biom_foliage. The foliage mass needed to
##   reproduce the MEASURED LAI is obs_lai / SLA_eff, and the fraction of
##   aboveground biomass that represents follows directly.
##
## If you want the functional route as a cross-check, transcribe the Marklund
## (1988) coefficients into marklund_coef.csv (species,component,a,b,c) and
## the script will use it. Take them from the original report.
##
##   Rscript 28_foliage_fraction.R
## ---------------------------------------------------------------------------

PLOTS <- '~/Documents/SLU/3PG/devel/UVa_work/validation_out/preflight_results.csv'
VARIANT <- 'observed/obs'
CONST <- c('Pinus sylvestris' = 0.045, 'Picea abies' = 0.085)
K     <- c('Pinus sylvestris' = 0.3247, 'Picea abies' = 0.638)   # sets 39, 37

hr <- function(x) cat('\n', strrep('-', 70), '\n', x, '\n', sep = '')
elas  <- function(k, L) { kl <- k * L; kl * exp(-kl) / (1 - exp(-kl)) }
fapar <- function(k, L) 1 - exp(-k * L)

d <- read.csv(path.expand(PLOTS))
d <- d[d$variant == VARIANT, ]
need <- c('Plot_ID', 'sp', 'lai', 'obs_lai')
if (!all(need %in% names(d))) stop('missing columns: ',
    paste(setdiff(need, names(d)), collapse = ', '))

## plot means, matching the revised Methods
a <- aggregate(d[intersect(c('lai','obs_lai','sim_gpp','obs_gpp'), names(d))],
               by = list(Plot_ID = d$Plot_ID, sp = d$sp), FUN = mean, na.rm = TRUE)

## ---------------------------------------------------------------------------
hr('1. LAI BIAS AT THE PLOTS')
lai <- do.call(rbind, lapply(split(a, a$sp), function(x) data.frame(
    species = x$sp[1], n = nrow(x),
    sim = median(x$lai, na.rm = TRUE),
    obs = median(x$obs_lai, na.rm = TRUE),
    ratio = median(x$lai, na.rm = TRUE) / median(x$obs_lai, na.rm = TRUE))))
lai$fraction_published <- CONST[lai$species]
lai$fraction_implied   <- lai$fraction_published / lai$ratio
print(lai, row.names = FALSE, digits = 3)

## stated from the data rather than assumed (the old grouping gave
## opposite-signed bias; the corrected species do not)
sg <- sign(lai$ratio - 1)
cat(sprintf('\nBias direction: %s.\n', paste(sprintf('%s %s (ratio %.2f)', lai$species,
    ifelse(abs(lai$ratio - 1) < 0.05, 'within 5 %', ifelse(sg > 0, 'high', 'low')), lai$ratio),
    collapse = '; ')))

## ---------------------------------------------------------------------------
hr('2. IMPLIED FOLIAGE FRACTION')
cat('Under LAI = biom_foliage * SLA with SLA unchanged, the foliage fraction\n')
cat('scales inversely with the LAI bias.\n\n')
for (i in seq_len(nrow(lai))) {
    s <- lai$species[i]
    if (!s %in% names(CONST)) next
    cat(sprintf('  %-18s %.4f  ->  %.4f   (ratio %.2f)\n',
                s, CONST[s], CONST[s] / lai$ratio[i], lai$ratio[i]))
}
fi <- function(col, s) lai[[col]][lai$species == s]
cat(sprintf('\nSpruce-to-pine contrast: published %.2f, implied %.2f. Check against\n',
    fi('fraction_published', 'Picea abies') / fi('fraction_published', 'Pinus sylvestris'),
    fi('fraction_implied', 'Picea abies') / fi('fraction_implied', 'Pinus sylvestris')))
cat('Marklund over the local diameter range before adopting it.\n')

## optional functional cross-check, only if the user supplied coefficients
if (file.exists('marklund_coef.csv')) {
    hr('2b. MARKLUND CROSS-CHECK (user-supplied coefficients)')
    mc <- read.csv('marklund_coef.csv')
    itc <- tryCatch(read.csv(path.expand(
        '~/Documents/SLU/3PG/devel/UVa_work/validation_out/itc_vs_measured.csv')),
        error = function(e) NULL)
    if (is.null(itc) || !'Dgv' %in% names(itc)) {
        cat('no measured Dgv available, skipped\n')
    } else {
        g <- function(D, r) exp(r$a + r$b * D / (D + r$c))
        for (s in unique(mc$species)) {
            D <- itc$Dgv[itc$sp == s & is.finite(itc$Dgv)]
            if (!length(D)) next
            comp <- split(mc[mc$species == s, ], mc$component[mc$species == s])
            if (!all(c('foliage','stem','branch') %in% names(comp))) next
            fo <- g(D, comp$foliage); st <- g(D, comp$stem); br <- g(D, comp$branch)
            fr <- fo / (fo + st + br)
            cat(sprintf('  %-18s n=%2d  Dgv %.1f cm  Marklund fraction %.4f (%.4f-%.4f)\n',
                        s, length(D), mean(D), median(fr), min(fr), max(fr)))
        }
    }
} else {
    cat('\n(marklund_coef.csv not present, functional cross-check skipped)\n')
}

## ---------------------------------------------------------------------------
hr('3. CONSEQUENCE FOR THE SATURATION PREDICTION')
cat('The ceiling is empirical and uses MEASURED LAI, so it is unaffected.\n')
cat('The PREDICTED ceiling uses mean LAI and therefore is affected.\n\n')

for (s in names(K)) {
    r <- lai[lai$species == s, ]
    if (!nrow(r)) next
    cat(sprintf('%-18s k = %.3f\n', s, K[s]))
    for (lab in c('sim', 'obs')) {
        L <- r[[lab]]
        cat(sprintf('   LAI from %-3s = %.2f   fAPAR %.3f   elasticity %.3f\n',
                    lab, L, fapar(K[s], L), elas(K[s], L)))
    }
    cat('\n')
}

if (all(names(K) %in% lai$species)) {
    gp <- function(col, s) lai[[col]][lai$species == s]
    rs <- elas(K[1], gp('sim','Pinus sylvestris')) / elas(K[2], gp('sim','Picea abies'))
    ro <- elas(K[1], gp('obs','Pinus sylvestris')) / elas(K[2], gp('obs','Picea abies'))
    cat(sprintf('pine:spruce elasticity ratio, simulated LAI = %.2f\n', rs))
    cat(sprintf('pine:spruce elasticity ratio, MEASURED  LAI = %.2f\n', ro))
    cat('\nReport the measured value. The mechanism is a claim about real\n')
    cat('canopies, and a reviewer recomputing it from measured LAI must get\n')
    cat('the same number you published.\n')
    if (ro < rs)
        cat('\nThe measured ratio is SMALLER: the current draft overstates the\n')
        cat('species contrast and section 25 must be rerun on measured LAI.\n')
}

## results/, not UVa_work/lai_bias.csv, which holds the pre-correction run
write.csv(lai, 'results/lai_bias.csv', row.names = FALSE)
cat('\nwritten: results/lai_bias.csv\n')
