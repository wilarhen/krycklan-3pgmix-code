## ---------------------------------------------------------------------------
## 45_soil_moisture_at_plots.R
##
## Values of the SLU soil moisture map (Agren et al. 2021; 2 m, 0 = dry to
## 100 = wet) at the validation plots, for Section 2.4: "Values at the 57
## validation plots span 1 to 99 with a median of 19". Read-only on the inputs;
## outputs to results/.
##
##   results/soil_moisture_at_plots.csv       one row per plot with coordinates
##       (57): Plot_ID, wetness, in_47 (one of the 47 analysed plots)
##   results/soil_moisture_summary.csv        n, minimum, median, maximum for
##       the 57 plots and for the 47
##
##   Rscript 45_soil_moisture_at_plots.R      (from UVa_work/, T7 mounted)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(terra); library(readxl) })
B   <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
OUT <- 'results'

co <- as.data.frame(read_excel('../validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx',
                               sheet = 'plot coordinates'))
idc <- grep('plot|id|stand', names(co), ignore.case = TRUE, value = TRUE)[1]
co <- data.frame(Plot_ID = suppressWarnings(as.numeric(co[[idc]])),
                 x = suppressWarnings(as.numeric(co$E_Sweref99TM)),
                 y = suppressWarnings(as.numeric(co$N_Sweref99TM)))
co <- co[is.finite(co$x) & is.finite(co$Plot_ID), ]

r <- rast(file.path(B, 'SLU_soil_map', 'SLUMarkfuktighetskarta.tif'))
co$wetness <- terra::extract(r, vect(co, geom = c('x', 'y'), crs = crs(r)), ID = FALSE)[, 1]
ids47 <- unique(read.csv('validation_out/preflight_results.csv')$Plot_ID)
co$in_47 <- co$Plot_ID %in% ids47
write.csv(co[, c('Plot_ID', 'wetness', 'in_47')], file.path(OUT, 'soil_moisture_at_plots.csv'),
          row.names = FALSE)

summ <- function(set, v) data.frame(set = set, n = sum(is.finite(v)), min = min(v, na.rm = TRUE),
                                    median = median(v, na.rm = TRUE), max = max(v, na.rm = TRUE))
S <- rbind(summ('57 validation plots', co$wetness), summ('47 analysed plots', co$wetness[co$in_47]))
print(S, row.names = FALSE)
write.csv(S, file.path(OUT, 'soil_moisture_summary.csv'), row.names = FALSE)
