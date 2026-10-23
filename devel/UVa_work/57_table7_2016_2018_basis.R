## ---------------------------------------------------------------------------
## 57_table7_2016_2018_basis.R
##
## Reruns the plot-level variants of Table 7 (script 10's A_control, B_hm2, C_boreal; script 16's
## observed/obs/a070 as D) with the same code and inputs, keeping every year of GPP, APAR, height
## and realised alpha_c. Scripts 10 and 16 store APAR, height and the modifiers for 2017 only
## (pk's default yr = '2017'), so Table 7's r(GPP, APAR) paired 2016-2018 mean GPP with 2017 APAR.
## This script reports r(GPP, APAR) and median height on three bases (Table 7 as first printed,
## 2016-2018 GPP with 2017 APAR; both 2017; both 2016-2018 plot means) and checks that it reproduces the stored
## sim_gpp and 2017 apar of params_results.csv and results/preflight_results.csv.
## Writes results/table7_2016_2018_basis.csv. First run 21 Oct 2026.
##
##   Rscript 57_table7_2016_2018_basis.R      (from devel/UVa_work; needs the T7)
## ---------------------------------------------------------------------------
## Read-only reproduction of script 10's plot runs, keeping every year of APAR, GPP and height.
suppressPackageStartupMessages({ library(terra); library(readxl); library(dplyr); library(tidyr); library(r3PG) })
bp  <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
rd  <- file.path(bp, '3PG_10m/x731866.017196386_y7132899.71658545.RData')
xl  <- '/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx'
cf  <- 0.47
fol_marklund <- c('Pinus sylvestris' = 0.045, 'Picea abies' = 0.085)
root_frac    <- c('Pinus sylvestris' = 0.25/0.76, 'Picea abies' = 0.25/0.75)
nm <- load(rd); xx <- get(nm[1])
cl <- xx[[1]] |> filter(group == 'climate') |> select(date, variable, value) |> distinct(date, variable, .keep_all = TRUE) |>
  pivot_wider(names_from = variable, values_from = value) |> arrange(date)
rm(xx); invisible(gc())
d_climate <- cl |> transmute(year = as.integer(format(date, '%Y')), month = as.integer(format(date, '%m')),
                             tmp_min, tmp_max, tmp_ave, prcp, srad = solar_rad, frost_days) |> as.data.frame()
obs_all <- suppressMessages(read_excel(xl, sheet = '50 main plots', skip = 2, .name_repair = 'unique')) |> as.data.frame() |>
  mutate(Plot_ID = suppressWarnings(as.integer(Plot_ID)), Year = suppressWarnings(as.integer(Year))) |>
  filter(!is.na(Plot_ID), Year %in% 2016:2018)
d <- obs_all |> filter(Year == 2016) |> transmute(Plot_ID, Age, Trees, AGCtree, LAImax, SPM, sp = ifelse(SPM == 1, 'Pinus sylvestris', 'Picea abies'))
cc <- suppressMessages(read_excel(xl, sheet = 'plot coordinates')) |> rename(Plot_ID = 1, E = 2, N = 3) |>
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) |> filter(!is.na(E)) |> mutate(Plot_ID = as.integer(Plot_ID))
d <- inner_join(d, cc, by = 'Plot_ID')
pv <- vect(cbind(d$E, d$N), type = 'points', crs = 'EPSG:3006')
d$latitude <- crds(project(pv, 'EPSG:4326'))[, 2]
dm <- rast(file.path(bp, 'Sweden_dem/dem.tif'))
d$altitude <- terra::extract(dm[[1]], project(pv, crs(dm)), ID = FALSE)[, 1]
bor_spruce <- c(Tmin = -3, Topt = 14, Tmax = 30, alphaCx = 0.050); bor_pine <- c(Tmin = -5, Topt = 13, Tmax = 30, alphaCx = 0.050)
variants <- list(
  A_control = list(hm = 1, pine = 'Pinus sylvestris 30', spruce = 'Picea abies 37', ov = NULL),
  B_hm2     = list(hm = 2, pine = 'Pinus sylvestris 39', spruce = 'Picea abies 37', ov = NULL),
  C_boreal  = list(hm = 2, pine = 'Pinus sylvestris 39', spruce = 'Picea abies 37', ov = list('Pinus sylvestris' = bor_pine, 'Picea abies' = bor_spruce)),
  D_a070    = list(hm = 2, pine = 'Pinus sylvestris 39', spruce = 'Picea abies 37', ov = list('Picea abies' = c(alphaCx = 0.070))))
build_par <- function(sp, pset, ov) { p <- get_parameters('parameters', sp_names = pset); names(p)[2] <- sp
  if (!is.null(ov) && sp %in% names(ov)) { o <- ov[[sp]]; for (k in names(o)) if (k %in% p$parameter) p[p$parameter == k, sp] <- o[[k]] }; p }
run_one <- function(i, v) {
  sp <- d$sp[i]; agb <- d$AGCtree[i] / cf; age <- d$Age[i]
  if (!is.finite(agb) || agb <= 0 || !is.finite(age) || age <= 0) return(NULL)
  fol <- agb * fol_marklund[[sp]]
  sprow <- data.frame(species = sp, planted = sprintf('%04d-01', 2016 - round(age)), fertility = 0.5, stems_n = d$Trees[i],
                      biom_stem = agb * (1 - fol_marklund[[sp]]), biom_root = agb * root_frac[[sp]], biom_foliage = fol, stringsAsFactors = FALSE)
  site <- data.frame(latitude = d$latitude[i], altitude = d$altitude[i], soil_class = 2, asw_i = 150, asw_min = 0, asw_max = 200, from = '2016-01', to = '2018-12')
  pset <- if (sp == 'Pinus sylvestris') v$pine else v$spruce
  par <- build_par(sp, pset, v$ov); sz <- get_parameters('sizeDist', sp_names = pset); names(sz)[2] <- sp
  tryCatch(run_3PG(site = site, species = sprow, climate = d_climate, thinning = NULL,
                   parameters = prepare_parameters(par, sp_names = sp), size_dist = prepare_sizeDist(sz, sp_names = sp),
                   settings = list(light_model = 2, transp_model = 2, phys_model = 2, height_model = v$hm, correct_bias = 0, calculate_d13c = 0),
                   check_input = TRUE, df_out = TRUE), error = function(e) NULL)
}
rows <- list()
for (vn in names(variants)) for (i in seq_len(nrow(d))) {
  o <- run_one(i, variants[[vn]]); if (is.null(o)) next
  o$Year <- as.integer(format(o$date, '%Y'))
  g <- tapply(o$value[o$variable == 'gpp'], o$Year[o$variable == 'gpp'], sum) * cf * 100
  a <- tapply(o$value[o$variable == 'apar'], o$Year[o$variable == 'apar'], sum)
  h <- tapply(o$value[o$variable == 'height'], o$Year[o$variable == 'height'], mean); ac <- tapply(o$value[o$variable == 'alpha_c'], o$Year[o$variable == 'alpha_c'], mean)
  rows[[length(rows) + 1]] <- data.frame(variant = vn, Plot_ID = d$Plot_ID[i], sp = d$sp[i], Year = as.integer(names(g)), sim_gpp = as.vector(g), apar = as.vector(a[names(g)]), height = as.vector(h[names(g)]), alpha_c = as.vector(ac[names(g)]))
}
r <- do.call(rbind, rows)
cat('apar months per plot-variant:', unique(table(paste(r$variant, r$Plot_ID)) ), 'years each;  years present:', sort(unique(r$Year)), '\n')
## reproduction check against the stored files (2017 values)
p10 <- read.csv('validation_out/params_results.csv')
pf  <- read.csv('results/preflight_results.csv')
for (vn in c('A_control','B_hm2','C_boreal')) { m <- merge(r[r$variant == vn, ], p10[p10$variant == vn, ], by = c('Plot_ID','Year'))
  cat(sprintf('check %-9s n %d  max|sim diff| %.2e  max|apar2017 diff| %.2e\n', vn, nrow(m), max(abs(m$sim_gpp.x - m$sim_gpp.y)), max(abs(m$apar.x[m$Year == 2017] - m$apar.y[m$Year == 2017])))) }
m <- merge(r[r$variant == 'D_a070', ], pf[pf$variant == 'observed/obs/a070', ], by = c('Plot_ID','Year'))
cat(sprintf('check D_a070    n %d  max|sim diff| %.2e  max|apar2017 diff| %.2e\n', nrow(m), max(abs(m$sim_gpp.x - m$sim_gpp.y)), max(abs(m$apar.x[m$Year == 2017] - m$apar.y[m$Year == 2017]))))
## three bases, per variant and species, plus pooled and the variant A means quoted in the text
out <- NULL
add <- function(vn, grp, x) {
  pm <- aggregate(cbind(sim_gpp, apar, height, alpha_c) ~ Plot_ID, x, mean); a17 <- x[x$Year == 2017, ]
  mixed <- merge(aggregate(sim_gpp ~ Plot_ID, x, mean), a17[, c('Plot_ID', 'apar')], by = 'Plot_ID')
  data.frame(variant = vn, group = grp, plots = nrow(pm),
             r_gpp_apar_as_printed = cor(mixed$sim_gpp, mixed$apar), r_gpp_apar_both_2017 = cor(a17$sim_gpp, a17$apar),
             r_gpp_apar_both_2016_2018 = cor(pm$sim_gpp, pm$apar),
             median_h_2017 = median(a17$height), median_h_2016_2018 = median(pm$height),
             mean_apar_2017 = mean(a17$apar), mean_apar_2016_2018 = mean(pm$apar),
             mean_alpha_c_2017 = mean(a17$alpha_c), mean_alpha_c_2016_2018 = mean(pm$alpha_c),
             mean_sim_gpp_2016_2018 = mean(pm$sim_gpp)) }
for (vn in names(variants)) {
  for (s in c('Pinus sylvestris', 'Picea abies')) { x <- r[r$variant == vn & r$sp == s, ]; if (nrow(x)) out <- rbind(out, add(vn, s, x)) }
  x <- r[r$variant == vn, ]
  pm <- aggregate(cbind(sim_gpp, apar) ~ Plot_ID + sp, x, mean); a17 <- x[x$Year == 2017, ]
  mixed <- merge(aggregate(sim_gpp ~ Plot_ID, x, mean), a17[, c('Plot_ID', 'apar')], by = 'Plot_ID')
  out <- rbind(out, data.frame(variant = vn, group = 'pooled', plots = nrow(pm),
    r_gpp_apar_as_printed = cor(mixed$sim_gpp, mixed$apar), r_gpp_apar_both_2017 = cor(a17$sim_gpp, a17$apar),
    r_gpp_apar_both_2016_2018 = cor(pm$sim_gpp, pm$apar), median_h_2017 = NA, median_h_2016_2018 = NA,
    mean_apar_2017 = NA, mean_apar_2016_2018 = NA, mean_alpha_c_2017 = NA, mean_alpha_c_2016_2018 = NA, mean_sim_gpp_2016_2018 = NA)) }
write.csv(out, 'results/table7_2016_2018_basis.csv', row.names = FALSE)
print(out[, c('variant','group','plots','r_gpp_apar_as_printed','r_gpp_apar_both_2017','r_gpp_apar_both_2016_2018','median_h_2016_2018')], digits = 4, row.names = FALSE)
