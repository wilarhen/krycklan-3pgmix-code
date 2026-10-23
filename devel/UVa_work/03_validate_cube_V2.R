#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## 03_validate_cube_V2.R
##
## Validation of the 3-PG Krycklan GPP cube against the Peichl et al. (2023)
## plot carbon-flux dataset (Global Change Biology 29:1119-1132).
##
## Run:
##   Rscript 03_validate_cube_V2.R
##   Rscript 03_validate_cube_V2.R /path/to/GPP.nc /path/to/Dataset.xlsx
##
## Everything printed is also written to  validation_out/validation_log.txt
## Send that file back. It is self-contained.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra); library(sf); library(readxl)
  library(dplyr); library(tidyr); library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG ----

args <- commandArgs(trailingOnly = TRUE)

cube_file <- if (length(args) >= 1) args[1] else
  '/run/media/wilar/T7 Shield/slu_data_bases/3PG/out_3PG_netcdf/GPP.nc'
xlsx_path <- if (length(args) >= 2) args[2] else
  'Dataset_2016-2018_50_7_plots_compiled.xlsx'

out_dir     <- 'validation_out'
epsg_plots  <- 3006          # SWEREF99TM
years_val   <- 2016:2018
buffer_m    <- 10            # plot footprint > 10 m pixel; 0 disables
carbon_frac <- 0.47          # tDM -> tC ; sensitivity reported below
tha_to_gm2  <- 100           # 1 t/ha = 100 g/m2

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out_dir, 'validation_log.txt')
con <- file(log_file, open = 'wt'); sink(con, split = TRUE)
## NB: stderr is deliberately NOT sunk. sink(type='message') cannot split, so
## it would hide errors from the terminal. Capture both with:
##   Rscript 03_validate_cube_V2.R 2>&1 | tee validation_out/run.txt
on.exit({ sink(); close(con) }, add = TRUE)

rule <- function(x) cat('\n', strrep('=', 72), '\n', x, '\n',
                        strrep('=', 72), '\n', sep = '')

rule('RUN INFO')
cat('date      : ', format(Sys.time()), '\n')
cat('cube      : ', cube_file, '\n')
cat('xlsx      : ', xlsx_path, '\n')
cat('R         : ', R.version.string, '\n')
cat('terra     : ', as.character(packageVersion('terra')), '\n')
cat('buffer_m  : ', buffer_m, ' | carbon_frac: ', carbon_frac, '\n', sep = '')
cat('workdir   : ', getwd(), '\n')

if (!file.exists(cube_file))
  stop('CUBE NOT FOUND: ', cube_file, call. = FALSE)
if (!file.exists(xlsx_path))
  stop('XLSX NOT FOUND: ', xlsx_path,
       '\n  working directory is: ', getwd(),
       '\n  xlsx files here: ',
       paste(list.files(pattern = '\\.xlsx$'), collapse = ', '),
       '\n  pass an absolute path as the 2nd argument.', call. = FALSE)

## ------------------------------------------------------- 1. READ THE CUBE ---

rule('1. CUBE')
cube <- rast(cube_file)
print(cube)

lyr_time <- time(cube)
stopifnot(!any(is.na(lyr_time)), nlyr(cube) == length(lyr_time))
cat('\nCRS code  : ', crs(cube, describe = TRUE)$code, '\n')
cat('time range: ', format(min(lyr_time)), ' to ', format(max(lyr_time)), '\n')

## ------------------------------------------- 2. VALIDATION DATA -> SPATIAL ---

rule('2. VALIDATION DATA')

## Header is on row 3 (rows 1-2 are variable group and units).
## Several columns are duplicated in the sheet (Cdw, Clitter, SOC) so names
## must be repaired, otherwise dplyr fails.
read_plots <- function(sheet) {
  x <- suppressMessages(
    read_excel(xlsx_path, sheet = sheet, skip = 2, .name_repair = 'unique')
  ) |> as.data.frame()
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
  x
}

main <- read_plots('50 main plots')
old  <- read_plots('7 extra old-stand plots')

## NOTE: column ORDER of SPM and DTS differs between the two sheets, and
## the '50 main plots' header is the transposed one, so reading by NAME
## does NOT resolve it. read_plots() corrects that sheet at read time. The primary analysis uses the 50 main plots only
## (GPPff is absent from the old-stand sheet), so this ambiguity does not
## affect the headline result. Cross-tab printed below for inspection.
cat('main sheet cols 9:11 :', paste(names(main)[9:11], collapse = ', '), '\n')
cat('old  sheet cols 9:11 :', paste(names(old)[9:11],  collapse = ', '), '\n')

obs_raw <- bind_rows(main |> mutate(set = 'main'),
                     old  |> mutate(set = 'oldstand')) |>
  mutate(Plot_ID = as.integer(Plot_ID), Year = as.integer(Year)) |>
  filter(!is.na(Plot_ID))

cat('\nrows: main =', nrow(main), '| oldstand =', nrow(old), '\n')
cat('GPPff non-missing by set:\n')
print(obs_raw |> group_by(set) |>
        summarise(n = n(), GPP_ok = sum(!is.na(GPP)),
                  GPPff_ok = sum(!is.na(GPPff)), .groups = 'drop'))

coords <- suppressMessages(read_excel(xlsx_path, sheet = 'plot coordinates')) |>
  rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(Plot_ID = as.integer(Plot_ID), E = as.numeric(E), N = as.numeric(N)) |>
  filter(!is.na(E), !is.na(N))

obs <- obs_raw |> filter(Year %in% years_val) |> inner_join(coords, by = 'Plot_ID')
cat('\nplots with coordinates:', n_distinct(obs$Plot_ID),
    '| plot-years:', nrow(obs), '\n')

miss <- setdiff(unique(obs_raw$Plot_ID), coords$Plot_ID)
if (length(miss)) cat('WARNING plots lacking coordinates:', paste(miss, collapse = ', '), '\n')

pts <- st_as_sf(distinct(obs, Plot_ID, E, N), coords = c('E', 'N'), crs = epsg_plots) |>
  st_transform(crs(cube))

## are all plots inside the cube extent?
e <- ext(cube); xy <- st_coordinates(pts)
inside <- xy[,1] >= e$xmin & xy[,1] <= e$xmax & xy[,2] >= e$ymin & xy[,2] <= e$ymax
cat('plots inside cube extent:', sum(inside), '/', nrow(pts), '\n')
if (any(!inside)) cat('OUTSIDE:', paste(pts$Plot_ID[!inside], collapse = ', '), '\n')

## ---------------------------------------------------- 3. EXTRACT FROM CUBE ---

rule('3. EXTRACTION + PIXEL QUALITY')

geom_pt  <- vect(pts)
geom_buf <- if (buffer_m > 0) vect(st_buffer(pts, buffer_m)) else geom_pt

## pixel-quality check on one growing-season layer
ref_i <- which(format(lyr_time, '%Y-%m') == '2017-07')[1]
ref   <- cube[[ref_i]]
qual <- terra::extract(ref, geom_buf, ID = TRUE) |>
  setNames(c('ID', 'v')) |>
  group_by(ID) |>
  summarise(n_cell = n(), n_na = sum(is.na(v)),
            n_zero = sum(v == 0, na.rm = TRUE), .groups = 'drop') |>
  mutate(Plot_ID = pts$Plot_ID[ID])

cat('reference layer for QC:', format(lyr_time[ref_i]), '\n')
cat('plots whose buffer contains ANY NA  :', sum(qual$n_na   > 0), '\n')
cat('plots whose buffer contains ANY zero:', sum(qual$n_zero > 0), '\n')
cat('plots entirely NA                   :', sum(qual$n_na == qual$n_cell), '\n')

dead <- qual |> filter(n_na == n_cell)
if (nrow(dead)) {
  cat('\nPLOTS ABSENT FROM SIMULATED DOMAIN (investigate: likely dropped by the\n',
      'reduce(intersect) on common_names in polygon_grids_v01.r line 206)\n', sep = '')
  cat('  Plot_ID:', paste(sort(dead$Plot_ID), collapse = ', '), '\n')
  cat('  lost   :', nrow(dead), 'of', nrow(qual),
      sprintf('(%.0f%% of validation plots)\n', 100 * nrow(dead) / nrow(qual)))
}
partial <- qual |> filter(n_na > 0, n_na < n_cell)
if (nrow(partial)) {
  cat('\nPLOTS WITH PARTIAL NA COVERAGE (buffer mean uses fewer cells)\n')
  print(as.data.frame(partial |> select(Plot_ID, n_cell, n_na, n_zero)), row.names = FALSE)
}
if (any(qual$n_zero > 0))
  print(qual |> filter(n_zero > 0) |> select(Plot_ID, n_cell, n_na, n_zero))

extract_annual <- function(g, label) {
  ex <- terra::extract(cube, g, fun = mean, na.rm = TRUE, ID = TRUE)
  key <- data.frame(layer = names(cube), date = lyr_time)
  ex |>
    mutate(Plot_ID = pts$Plot_ID[ID]) |> select(-ID) |>
    pivot_longer(-Plot_ID, names_to = 'layer', values_to = 'value') |>
    left_join(key, by = 'layer') |>
    mutate(Year = as.integer(format(date, '%Y'))) |>
    filter(Year %in% years_val) |>
    group_by(Plot_ID, Year) |>
    summarise(n_month = sum(!is.na(value)),
              sim_tDM = sum(value, na.rm = TRUE), .groups = 'drop') |>
    filter(n_month == 12) |>
    mutate(sim = sim_tDM * carbon_frac * tha_to_gm2, support = label)
}

sim_buf <- extract_annual(geom_buf, sprintf('buffer%dm', buffer_m))
sim_pt  <- extract_annual(geom_pt,  'pixel')
cat('\nannual sims: buffer =', nrow(sim_buf), '| pixel =', nrow(sim_pt), '\n')

## ------------------------------------------------------ 4. MATCH OBSERVED ---

rule('4. MATCHED PAIRS')

## The metadata sheet documents the stand age class as 'SAC', but the DATA
## sheets label the same column 'SDS'. Harmonise to SAC.
if (!'SAC' %in% names(obs) && 'SDS' %in% names(obs)) {
  obs <- dplyr::rename(obs, SAC = SDS)
  cat("\nnote: renamed data-sheet column 'SDS' -> 'SAC' (stand age class)\n")
}

need <- c('GPP','GPPff','NPPtree','Age','LAImax','Elevation','Slope',
          'SPM','DTS','SAC')
if (length(setdiff(need, names(obs))))
  stop('missing expected columns: ', paste(setdiff(need, names(obs)), collapse = ', '),
       '\n  available: ', paste(names(obs), collapse = ', '), call. = FALSE)

obs_sel <- obs |>
  transmute(Plot_ID, Year, set,
            obs_tree = GPP - GPPff,
            obs_tot  = GPP,
            obs_npptree = NPPtree,
            Age, LAImax, Elevation, Slope,
            SPM = factor(SPM, 1:2, c('Sediment', 'Till')),
            DTS = factor(DTS, 1:2, c('Pine', 'Spruce')),
            SAC = factor(SAC, 1:5, c('Initiation','Young','Middle','Mature','Old')))

val_all <- inner_join(sim_buf, obs_sel, by = c('Plot_ID', 'Year'))
val     <- val_all |> filter(!is.na(obs_tree)) |> mutate(obs = obs_tree, resid = sim - obs)
val_tot <- val_all |> filter(!is.na(obs_tot))  |> mutate(obs = obs_tot,  resid = sim - obs)

cat('PRIMARY   tree GPP :', nrow(val),     'pairs,', n_distinct(val$Plot_ID),     'plots\n')
cat('SECONDARY total GPP:', nrow(val_tot), 'pairs,', n_distinct(val_tot$Plot_ID), 'plots\n')
stopifnot(nrow(val) > 0)

## ------------------------------------------------------- 5. ERROR METRICS ---

rule('5. ERROR METRICS')

metrics <- function(obs, sim) {
  ok <- is.finite(obs) & is.finite(sim); obs <- obs[ok]; sim <- sim[ok]
  if (length(obs) < 3) return(data.frame(n = length(obs)))
  e <- sim - obs
  data.frame(n = length(obs),
             obs_mean = mean(obs), sim_mean = mean(sim),
             bias = mean(e), rbias = 100 * mean(e) / mean(obs),
             MAE = mean(abs(e)),
             RMSE = sqrt(mean(e^2)), rRMSE = 100 * sqrt(mean(e^2)) / mean(obs),
             R2 = cor(obs, sim)^2,
             NSE = 1 - sum(e^2) / sum((obs - mean(obs))^2),
             CCC = 2 * cov(obs, sim) /
               (var(obs) + var(sim) + (mean(obs) - mean(sim))^2))
}

cat('\n-- PRIMARY: tree GPP (GPP - GPPff), 50 main plots --\n')
print(round(metrics(val$obs, val$sim), 3))
cat('\n-- SECONDARY: total ecosystem GPP (expect negative bias) --\n')
print(round(metrics(val_tot$obs, val_tot$sim), 3))

cat('\n-- observed vs simulated distributions (primary) --\n')
print(rbind(observed  = summary(val$obs),
            simulated = summary(val$sim)))

cat('\n-- extraction support: buffer vs single pixel --\n')
cmp <- inner_join(sim_buf |> select(Plot_ID, Year, sim_buf = sim),
                  sim_pt  |> select(Plot_ID, Year, sim_pix = sim),
                  by = c('Plot_ID','Year'))
cat('mean |buffer - pixel| =', round(mean(abs(cmp$sim_buf - cmp$sim_pix)), 2),
    'gC/m2/yr | cor =', round(cor(cmp$sim_buf, cmp$sim_pix), 4), '\n')

strat <- function(v, gv, lab)
  v |> mutate(.g = paste0(lab, ':', .data[[gv]])) |> group_by(.g) |>
    group_modify(~metrics(.x$obs, .x$sim)) |> ungroup()

by_strata <- bind_rows(
  val |> mutate(.g = 'ALL') |> group_by(.g) |> group_modify(~metrics(.x$obs, .x$sim)) |> ungroup(),
  strat(val, 'Year', 'Year'), strat(val, 'SPM', 'SPM'),
  strat(val, 'DTS', 'DTS'),   strat(val, 'SAC', 'SAC'))
cat('\n-- STRATIFIED (SPM is the water-limitation test) --\n')
print(as.data.frame(by_strata), digits = 3, row.names = FALSE)

cat('\n-- carbon fraction sensitivity --\n')
print(do.call(rbind, lapply(c(0.45, 0.47, 0.50), function(cf)
  cbind(cf = cf, round(metrics(val$obs, val$sim_tDM * cf * tha_to_gm2), 2)))),
  row.names = FALSE)

cat('\n-- residual correlates (what the model is missing) --\n')
for (v in c('Age','LAImax','obs','Elevation','Slope'))
  cat(sprintf('  cor(resid, %-10s) = %+.3f\n', v,
              cor(val$resid, val[[v]], use = 'complete.obs')))

cat('\n-- how far can the model reach? --\n')
cat('  simulated max :', round(max(val$sim)), 'gC/m2/yr\n')
cat('  observed  max :', round(max(val$obs)), 'gC/m2/yr\n')
cat('  observed plot-years above simulated max:',
    sum(val$obs > max(val$sim)), '/', nrow(val),
    sprintf('(%.0f%%)', 100 * mean(val$obs > max(val$sim))), '\n')

## ---------------------------------------------------------- 6. DIAGNOSTICS ---

rule('6. FIGURES + FILES')

lims <- range(c(val$obs, val$sim), na.rm = TRUE)
p1 <- ggplot(val, aes(obs, sim, colour = DTS, shape = SPM)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2.3, alpha = .85) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = 'Observed tree GPP (gC m-2 yr-1)',
       y = 'Simulated GPP (gC m-2 yr-1)') + theme_bw()

p2 <- ggplot(val, aes(Age, resid, colour = DTS)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_point(size = 2.3, alpha = .85) +
  labs(x = 'Stand age (years)', y = 'Residual (sim - obs)') + theme_bw()

p3 <- ggplot(val, aes(SPM, resid, fill = SPM)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(alpha = .7) +
  facet_wrap(~Year) + labs(y = 'Residual (sim - obs)') +
  theme_bw() + theme(legend.position = 'none')

ggsave(file.path(out_dir, 'fig_1to1.png'),      p1, width = 6,   height = 5.5, dpi = 300)
ggsave(file.path(out_dir, 'fig_resid_age.png'), p2, width = 6.5, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, 'fig_resid_spm.png'), p3, width = 7,   height = 4,   dpi = 300)

write.csv(val,       file.path(out_dir, 'validation_pairs.csv'),  row.names = FALSE)
write.csv(by_strata, file.path(out_dir, 'metrics_by_strata.csv'), row.names = FALSE)
write.csv(qual,      file.path(out_dir, 'pixel_quality.csv'),     row.names = FALSE)

cat('written to ', normalizePath(out_dir), ':\n', sep = '')
cat(paste(' ', list.files(out_dir)), sep = '\n')

rule('DONE')
cat('Send back: validation_out/validation_log.txt\n')
