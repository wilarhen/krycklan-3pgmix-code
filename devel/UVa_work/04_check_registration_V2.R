#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 04_check_registration_V2.R
##
## R2 = 0.015 with a plausible value distribution means the cube is probably
## spatially scrambled rather than physically wrong. Two independent tests:
##
##   TEST A  mirror the PLOT COORDINATES inside the extent. Reflecting the
##           points is equivalent to flipping the raster but costs nothing and
##           cannot crash terra. If a mirrored version scores much better than
##           the original, the cube axes are inverted.
##
##   TEST B  Moran's I. A real landscape is strongly autocorrelated. If the
##           per-cell values were assigned in the wrong order, the histogram
##           still looks fine but Moran's I collapses toward zero.
##
##   Rscript 04_check_registration_V2.R <GPP.nc> <Dataset.xlsx>
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra); library(sf); library(readxl); library(dplyr); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
cube_file <- if (length(args) >= 1) args[1] else '~/3pg_data/out_3PG_netcdf/GPP.nc'
xlsx_path <- if (length(args) >= 2) args[2] else
  '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'

years_val <- 2016:2018; cf <- 0.47
dir.create('validation_out', showWarnings = FALSE)

cube <- rast(path.expand(cube_file))
lt <- time(cube)
i  <- which(as.integer(format(lt, '%Y')) %in% years_val)
e  <- ext(cube)
cat('cube  :', cube_file, '\n')
cat('extent: x', e$xmin, e$xmax, '| y', e$ymin, e$ymax, '\n')

## Subset the layers ONCE. Everything downstream works on this small stack.
sub <- cube[[i]]
key <- data.frame(layer = names(sub), date = lt[i])

## ---- observed -------------------------------------------------------------
obs <- suppressMessages(read_excel(xlsx_path, sheet = '50 main plots',
                                   skip = 2, .name_repair = 'unique')) |>
  as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         Year = as.integer(Year)) |>
  filter(!is.na(Plot_ID), Year %in% years_val) |>
  transmute(Plot_ID, Year, obs = GPP - GPPff)

coords <- suppressMessages(read_excel(xlsx_path, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)),
         E = as.numeric(E), N = as.numeric(N)) |>
  filter(!is.na(E)) |>
  inner_join(distinct(obs, Plot_ID), by = 'Plot_ID')

cat('plots :', nrow(coords), '\n')

## ---- TEST A: mirrored coordinate sets --------------------------------------
mirror <- function(E, N, mx, my) {
  data.frame(E = if (mx) (e$xmin + e$xmax) - E else E,
             N = if (my) (e$ymin + e$ymax) - N else N)
}

score <- function(mx, my, label) {
  m <- mirror(coords$E, coords$N, mx, my)
  v <- vect(cbind(m$E, m$N), type = 'points', crs = crs(cube))
  ex <- terra::extract(sub, v, ID = TRUE)
  s <- ex |> mutate(Plot_ID = coords$Plot_ID[ID]) |> select(-ID) |>
    pivot_longer(-Plot_ID, names_to = 'layer', values_to = 'v') |>
    left_join(key, by = 'layer') |>
    mutate(Year = as.integer(format(date, '%Y'))) |>
    group_by(Plot_ID, Year) |>
    summarise(nm = sum(!is.na(v)), sim = sum(v, na.rm = TRUE) * cf * 100,
              .groups = 'drop') |>
    filter(nm == 12)
  d <- inner_join(s, obs, by = c('Plot_ID','Year'))
  if (nrow(d) < 5)
    return(data.frame(orientation = label, n = nrow(d), r = NA, R2 = NA,
                      bias = NA, RMSE = NA))
  er <- d$sim - d$obs
  data.frame(orientation = label, n = nrow(d),
             r = cor(d$sim, d$obs), R2 = cor(d$sim, d$obs)^2,
             bias = mean(er), RMSE = sqrt(mean(er^2)))
}

cat('\n== TEST A: coordinate orientation ==\n')
resA <- bind_rows(
  score(FALSE, FALSE, 'as written'),
  score(FALSE, TRUE,  'N mirrored (vertical flip)'),
  score(TRUE,  FALSE, 'E mirrored (horizontal flip)'),
  score(TRUE,  TRUE,  'both mirrored (180 rotation)'))
print(resA, digits = 3, row.names = FALSE)

best <- resA$orientation[which.max(replace(resA$R2, is.na(resA$R2), -1))]
cat('\nbest orientation:', best, '\n')
if (best != 'as written')
  cat('>> AXIS INVERSION CONFIRMED. Fix the netcdf write, not the model.\n')

## ---- TEST B: spatial coherence --------------------------------------------
cat('\n== TEST B: spatial autocorrelation ==\n')
j <- which(format(lt, '%Y-%m') == '2017-07')[1]
l <- cube[[j]]
cat('layer:', format(lt[j]), '\n')
print(global(l, c('mean','sd','min','max'), na.rm = TRUE))

## crop to a manageable window so autocor is quick and memory-safe
win <- crop(l, ext(e$xmin, e$xmin + 3000, e$ymin, e$ymin + 3000))
mi  <- tryCatch(autocor(win, method = 'moran', global = TRUE),
                error = function(err) { cat('autocor failed:', conditionMessage(err), '\n'); NA })
cat("Moran's I on a 3x3 km window:", round(mi, 4), '\n')

## same test on the zero mask: mires should form coherent patches
zm  <- tryCatch(autocor(as.numeric(win == 0), method = 'moran', global = TRUE),
                error = function(err) NA)
cat("Moran's I of the zero mask   :", round(zm, 4), '\n')

cat('\n  I > 0.5  -> coherent landscape; any error is a shift/flip (see TEST A)\n')
cat('  I ~ 0    -> per-cell values scrambled. Suspect the ordering of\n')
cat('              common_names (reduce(intersect), line 206) versus the cell\n')
cat('              order assumed by timeSeriesRaster() in polygon_grids_v01.r\n')

## ---- map -------------------------------------------------------------------
png('validation_out/fig_cube_map_2017-07.png', width = 900, height = 950)
plot(l, main = paste('Simulated GPP', format(lt[j]), '(tDM/ha/month)'))
pv <- vect(cbind(coords$E, coords$N), type = 'points', crs = crs(cube))
plot(pv, add = TRUE, pch = 21, bg = 'white', cex = 1.1)
dev.off()
cat('\nwrote validation_out/fig_cube_map_2017-07.png\n')
cat('Look for: catchment outline, stream network, coherent stand patches.\n')
cat('If it looks like television static, TEST B has already told you why.\n')

write.csv(resA, 'validation_out/registration_test.csv', row.names = FALSE)
