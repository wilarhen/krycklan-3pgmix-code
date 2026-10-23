## ---------------------------------------------------------------------------
## 55_tall_cell_rerun.R
##
## Reruns the inputs of the corrected cube's tallest cell (cell 665356, 731605 E, 7132515 N,
## NMD class 114: Betula, Pinus sylvestris, Picea abies; stand ages 127, 280, 233 yr) with
## r3PG 0.1.6 and the pipeline's settings, at total aboveground biomass 150, 188, 189, 190 and
## 250 t DM/ha, and keeps the per-species output that the cube does not store (the cube holds
## stand sums and unweighted species means only). Writes results/tall_cell_rerun.csv.
## First run 18 Oct 2026.
##
##   Rscript 55_tall_cell_rerun.R      (local; needs the T7)
##
## APPROXIMATE, not a rebuild of the cell: climate is the saved v01 sector's (one 4 km PTHBV
## cell, 2001-2020), and the site row is fixed (latitude 64.24, altitude 250 m, soil class 2,
## ASW 120 in 0-150), not the cell's own. It gives a stand height of 67.8 m in Nov 2020 where
## the cube has 250.5 m, so it shows the mechanism, not the cube's values.
## Species inputs follow polygon_grids_v03_server.r: AGB split 1/n among the species present,
## foliage 0.045 / 0.085 / 0.020 of each share, root 0.25/0.76, 0.25/0.75, 0.25/0.75,
## stems 967 / 1162 / 1000 per ha, fertility 0.5, planted = 2025-07 minus stand age,
## parameter sets Betula 33, Picea abies 37, Pinus sylvestris 39, height_model 2.
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(r3PG); library(dplyr); library(tidyr)})
OUT <- 'results'
rd  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG/3PG_10m/x731866.017196386_y7132899.71658545.RData'

nm <- load(rd); xx <- get(nm[1])
cl <- xx[[1]] |> filter(group == 'climate') |> select(date, variable, value) |>
  distinct(date, variable, .keep_all = TRUE) |>
  pivot_wider(names_from = variable, values_from = value) |> arrange(date)
rm(xx); invisible(gc())
clim <- cl |> transmute(year = as.integer(format(date, '%Y')), month = as.integer(format(date, '%m')),
                        tmp_min, tmp_max, tmp_ave, prcp, srad = solar_rad, frost_days) |> as.data.frame()

sets <- c('Betula 33', 'Picea abies 37', 'Pinus sylvestris 39')
par <- get_parameters('parameters', sp_names = sets) |> rename_with(~ sub('\\s*\\d+$', '', .x))
szd <- get_parameters('sizeDist', sp_names = sets) |> rename_with(~ sub('\\s*\\d+$', '', .x))
ff <- c('Pinus sylvestris' = .045, 'Picea abies' = .085, 'Betula' = .020)
rf <- c('Pinus sylvestris' = .25/.76, 'Picea abies' = .25/.75, 'Betula' = .25/.75)
st <- c('Pinus sylvestris' = 967, 'Picea abies' = 1162, 'Betula' = 1000)
## planting date as compute_plantation_raster: reference 2025-07, fractional years to months
pl <- function(age) { a <- 2025.5 - age; y <- floor(a); m <- round((a - y) * 12); if (m == 0) m <- 1; sprintf('%04d-%02d', y, m) }
ages <- c('Betula' = 127, 'Pinus sylvestris' = 280, 'Picea abies' = 233)
sp <- c('Betula', 'Pinus sylvestris', 'Picea abies')     # slot order of class 114

run <- function(agb) {
  sh <- agb / length(sp)
  spec <- data.frame(species = sp, planted = sapply(ages[sp], pl), fertility = .5, stems_n = st[sp],
                     biom_stem = sh * (1 - ff[sp]), biom_root = sh * rf[sp], biom_foliage = sh * ff[sp],
                     stringsAsFactors = FALSE)
  site <- data.frame(latitude = 64.24, altitude = 250, soil_class = 2, asw_i = 120, asw_min = 0,
                     asw_max = 150, from = '2001-01', to = '2020-12')
  run_3PG(site, spec, clim, NULL, prepare_parameters(par, sp_names = sp), prepare_sizeDist(szd, sp_names = sp),
          settings = list(light_model = 2, transp_model = 2, phys_model = 2, height_model = 2,
                          correct_bias = 0, calculate_d13c = 0),
          check_input = TRUE, df_out = TRUE)
}

keep <- c('stems_n', 'biom_stem', 'dbh', 'basal_area', 'height', 'lai')
res <- lapply(c(150, 188, 189, 190, 250), function(agb) {
  o <- run(agb)
  o[o$variable %in% keep & format(o$date, '%Y-%m') %in% c('2001-01', '2020-11'), c('date', 'species', 'variable', 'value')] |>
    pivot_wider(names_from = variable, values_from = value) |>
    mutate(agb = agb, date = format(date, '%Y-%m'), .before = 1)
}) |> bind_rows() |> as.data.frame()
## stand value as the cube forms it: unweighted mean of height over the species present
stand <- res |> group_by(agb, date) |> summarise(species = 'stand (unweighted mean)', height = mean(height),
                                                n_species = n(), .groups = 'drop')
out <- bind_rows(res, stand) |> arrange(agb, date, species)
write.csv(out, file.path(OUT, 'tall_cell_rerun.csv'), row.names = FALSE)
print(out, digits = 4)
