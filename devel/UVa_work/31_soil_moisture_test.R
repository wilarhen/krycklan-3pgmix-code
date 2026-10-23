## ---------------------------------------------------------------------------
## 31_soil_moisture_test.R
##
## SPECIES COLUMN CORRECTED (1 Oct 2026). The first version took species from
## itc_vs_measured.csv, which script 11 wrote on 23 July with the DTS/SPM
## transposition, so its 'Spruce'/'Pine' (33/14) were till/sediment. Species
## now comes from preflight_results.csv, variant observed/obs. The raster is
## read from the T7, so the script runs locally from UVa_work/ rather than
## on the server. Seed added (the first version had none, so its intervals
## cannot be regenerated exactly; its r values reproduce). Output:
## results/section36_corrected.csv.
##
## The Discussion currently names spatially resolved soil water as the one
## proposed refinement the bounding analysis does not exclude, and calls for
## depth-to-water products derived from laser scanning. The SLU wetness map is
## exactly that product and it is already on this machine. Testing it is
## therefore obligatory before the paper claims the test has not been done.
##
## Two outcomes, both useful:
##   near zero for spruce  -> the last open mechanism is closed, and the
##                            Discussion says so rather than calling for it
##   substantial for spruce -> the missing driver is found, which is a larger
##                            result than the one currently reported
##
##   Rscript 31_soil_moisture_test.R
## ---------------------------------------------------------------------------

SM    <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG/SLU_soil_map/SLUMarkfuktighetskarta.tif'
PLOTS <- 'validation_out/itc_vs_measured.csv'    # coordinates only
GPP   <- 'validation_out/preflight_results.csv'
OUT   <- 'results'
dir.create(OUT, showWarnings = FALSE)
set.seed(20260824)

## published values, Section 3.6, grouped by the transposed column
OLD <- data.frame(old_group = c('Spruce', 'Pine'), old_n = c(33, 14),
                  old_r = c(-0.152, 0.551), old_lo = c(-0.410, 0.054),
                  old_hi = c(0.190, 0.864), old_ceiling = c(0.295, 0.685))

hr <- function(x) cat('\n', strrep('-', 66), '\n', x, '\n', sep = '')

r <- terra::rast(SM)
cat('wetness map: res ', paste(terra::res(r), collapse = ' x '),
    '  crs ', terra::crs(r, describe = TRUE)$code, '\n', sep = '')
mm <- terra::minmax(r)
cat('value range ', mm[1], ' to ', mm[2], '  type ', terra::datatype(r),
    '  colour table ', terra::has.colors(r), '\n', sep = '')

p <- read.csv(path.expand(PLOTS))
v <- terra::vect(p[is.finite(p$E) & is.finite(p$N), ],
                 geom = c('E', 'N'), crs = 'EPSG:3006')
if (!terra::same.crs(v, r)) v <- terra::project(v, terra::crs(r))
p$sm <- terra::extract(r, v)[, 2]

g <- read.csv(path.expand(GPP))
g <- g[g$variant == 'observed/obs', ]
p$sp_old <- p$sp                                  # transposed grouping
p$sp <- g$sp[match(p$Plot_ID, g$Plot_ID)]         # corrected species
p <- p[!is.na(p$sp), ]
ag <- aggregate(cbind(obs_gpp, obs_lai) ~ Plot_ID, g, mean)
p$gpp <- ag$obs_gpp[match(p$Plot_ID, ag$Plot_ID)]
p$lai <- ag$obs_lai[match(p$Plot_ID, ag$Plot_ID)]

u <- sort(unique(p$sm[is.finite(p$sm)]))
cat('distinct values at the plots: ', length(u), '  -> ',
    paste(head(u, 12), collapse = ', '), '\n', sep = '')
categorical <- length(u) <= 8

hr('WETNESS AGAINST MEASURED GPP')
res <- list()
for (s in unique(p$sp)) {
    x <- p[p$sp == s & is.finite(p$sm) & is.finite(p$gpp), ]
    if (nrow(x) < 5) next
    if (categorical) {
        cat(s, '  n = ', nrow(x), '\n', sep = '')
        print(round(tapply(x$gpp, x$sm, function(z)
            c(n = length(z), mean = mean(z), sd = stats::sd(z))), 1))
        if (length(unique(x$sm)) > 1) {
            k <- stats::kruskal.test(x$gpp ~ factor(x$sm))
            cat(sprintf('  Kruskal-Wallis p = %.3f\n', k$p.value))
        }
    } else {
        rr <- cor(x$sm, x$gpp)
        b <- replicate(2000, {
            i <- sample.int(nrow(x), nrow(x), TRUE)
            suppressWarnings(cor(x$sm[i], x$gpp[i]))
        })
        ci <- stats::quantile(b[is.finite(b)], c(.025, .975))
        cat(sprintf('%-17s n=%2d  r = %+.3f  95%% CI [%+.3f, %+.3f]\n',
                    s, nrow(x), rr, ci[1], ci[2]))
        res[[s]] <- data.frame(species = s, n = nrow(x), r = rr,
                               lo = ci[[1]], hi = ci[[2]])
    }
}

hr('FOR CONTEXT: THE CEILING, ON THE SAME PLOTS')
for (s in unique(p$sp)) {
    x <- p[p$sp == s & is.finite(p$lai) & is.finite(p$gpp), ]
    if (nrow(x) < 5) next
    cl <- cor(x$lai, x$gpp)
    cat(sprintf('%-17s n=%2d  r(GPP, measured LAI) = %+.3f\n', s, nrow(x), cl))
    if (!is.null(res[[s]])) res[[s]]$ceiling_lai <- cl
}

R <- do.call(rbind, res); rownames(R) <- NULL
R$old_group <- ifelse(R$species == 'Picea abies', 'Spruce', 'Pine')
R <- merge(R, OLD, by = 'old_group')
R$raster <- basename(SM); R$raster_min <- mm[1]; R$raster_max <- mm[2]
write.csv(R, file.path(OUT, 'section36_corrected.csv'), row.names = FALSE)
cat('\nwrote', file.path(OUT, 'section36_corrected.csv'), '\n')

hr('HOW TO READ THIS')
cat('If wetness explains no more of spruce production than leaf area does,\n')
cat('the one refinement the paper does not currently exclude is excluded,\n')
cat('and Section 4.8 should say so instead of calling for the test.\n')
cat('If it explains substantially more, the missing spruce driver is found\n')
cat('and that belongs in the Results, not the Discussion.\n')
cat('\nNote this is the static wetness index. The 1980-2023 soil moisture\n')
cat('series in Soil_moisture_1980_2023_Alex is a separate and stronger test,\n')
cat('since it resolves time as well as space.\n')
