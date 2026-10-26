## go to line 332
library('ncdf4')
library('terra')
library('wilaRtools')
require('R.utils')
require('lidR')
require('sf')
require('r3PG')
require('dplyr')
require('ggplot2')
require('tidyverse')
require('tidyr')
require('stringr')
require('purrr')
require('geodata')
library('lubridate')
require('patchwork')
require('parallel')
require('units')
require('data.table')
## require('terra')
wd('~','3PG/devel/R', setwd = TRUE)
sourceDirectory('.', modifiedOnly = TRUE, verbose = TRUE) #
wd(,'3PG/devel')
## temporary files
base_path <- '/media/wihe/LaCie/projects/3PG'
my_temp <- file.path(base_path, 'r_tmp')
dir.create(my_temp, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = my_temp)

wd. <- '/media/wihe/LaCie/projects/3PG/outputs_3PG/'
fl <- dir(wd.)[grepl('.tif', dir(wd.))]
fl. <- file.path(wd., fl)
vrtl <- vrt(fl., filename='file.vrt', overwrite = TRUE)
n. <- sample(1:length(fl), 1)
plot(vrtl[[n.]], main = names(vrtl[[n.]]))
## tmp <- terra::rast(file.path(wd., fl[n.]))

library('zoo')
n  <- nlyr(vrtl)
end_m <- as.yearmon("2020-12")               # last layer month
t_m   <- seq(end_m - (n - 1)/12, end_m, by = 1/12)
terra::time(vrtl) <- as.Date(t_m, frac = 1)

writeCDF(
  vrtl, "/media/wihe/LaCie/projects/3PG/out_3PG_netcdf/GPP.nc",
  overwrite = TRUE,
  varname   = "GPP",
  longname  = "Gross Primary Productivity",
  unit      = "tons dry matter ha-1 month-1",
  comments  = "These estimates are preliminary and must be validated against independent observations before use."
)



kryck_area_km2 <- st_area(as.polygons(ext(vrtl)) |> st_as_sf()) |> set_units("m^2") |> set_units("km^2")
## krycklan polygon
## start.time <- Sys.time()
sector_n  <- as.numeric(readLines('sector_n'))
kry_folder <- file.path(base_path, 'krycklan_borders')
kryck_bord <- import_feature_from_folder('vect', kry_folder, '.shp')
kryck_border_sf <- st_cast(st_as_sf(kryck_bord), 'POLYGON')[3,2]

plot(kryck_border_sf)


## kryck_bord. <- polygon_sector(kryck_bord,10,65,200)

r <- vrtl
# --- Quantiles across time ---
q_rasters <- quantile(r, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)
# Convert to data frame
df_q <- as.data.frame(q_rasters, xy = TRUE)
# Reshape long
df_long <- df_q %>%
  pivot_longer(
    cols = starts_with("q"),
    names_to = "quantile",
    values_to = "value"
  )
# Clean quantile names
df_long$quantile <- gsub("q", "Q", df_long$quantile)
# --- Compute statistics ---
r_mean <- mean(r, na.rm = TRUE)
q_rasters <- quantile(r, probs = c(0.05, 0.95), na.rm = TRUE)
r_range <- q_rasters[[2]] - q_rasters[[1]]   # p95 - p05
# --- Convert to data frames ---
df_mean <- as.data.frame(r_mean, xy = TRUE) %>%
  rename(value = names(r_mean)) %>%
  mutate(stat = "Mean")
df_range <- as.data.frame(r_range, xy = TRUE) %>%
  rename(value = names(r_range)) %>%
  mutate(stat = "Range (P95 - P05)")
df_plot <- bind_rows(df_mean, df_range)
## # --- Convert to data frames ---
## df_mean <- as.data.frame(r_mean, xy = TRUE) %>%
##   rename(value = names(r_mean)) %>%
##   mutate(stat = "Mean")
## df_range <- as.data.frame(r_range, xy = TRUE) %>%
##   rename(value = names(r_range)) %>%
##   mutate(stat = "Range (P95 - P05)")
## df_plot <- bind_rows(df_mean, df_range)
# --- Convert polygon to coords for plotting ---
poly_sf <- st_as_sf(kryck_border_sf)     # SpatVector → sf
poly_df <- st_coordinates(poly_sf) %>%
  as.data.frame() %>%
  rename(x = X, y = Y)
# Build grouping var to avoid connecting rings wrongly
poly_df$group <- interaction(poly_df$L1, poly_df$L2, drop = TRUE)

# --- Plot ---
gpp_plot <- ggplot() +
  # Raster quantiles
  geom_raster(data = df_plot, aes(x = x, y = y, fill = value)) +
  # Polygon outline
  geom_path(data = poly_df, aes(x = x, y = y),#, group = group),
            color = "grey30", linewidth = 0.8) +
    scale_fill_viridis_c(direction = 1,
                         option = "turbo",
                         begin = 0.07, end = 1) +
                         ## begin = 0.0, end = 0.7) +
  ## scale_fill_viridis_c(option = "plasma") +
  coord_equal(expand = FALSE) +
  ## scale_y_reverse() +
  facet_wrap(~stat, ncol = 2, nrow = 1) +
  theme_minimal() +
    labs(fill = expression("t ha"^-1*" mo"^-1*""),
         x = 'Longitude',
         y = 'Latitude',
        title = "GPP (2001–2020)"
         ) +
    theme(legend.position = 'top')
gpp_plot

ggsave(gpp_plot,
       file = '/home/wihe/Documents/slu_locals/3PG/devel/mosaic_Kryckland_3PG_GPP_output.pdf')



























kryck_bord. <- polygon_sector(kryck_bord,10,sector_n,200)
kryck_border_for_coarse_climate <- polygon_sector(kryck_bord, 1, 1, 4000)[[1]]
kryck_border <- kryck_bord.[[1]]
## DEM raster
dem_folder <- file.path(base_path, 'Sweden_dem')
dem_dt <- import_feature_from_folder('rast', dem_folder, 'dem.tif')
kryck_dem <- crop_raster_fast(dem_dt, kryck_border)
## Scale factor
## fact <- res_target / res(kryck_dem)[1] 
kryck_elev <- agg_and_crop.(kryck_dem, 10, kryck_border)
## CLIMATE DATA
nc_folder <- file.path(base_path,'climate_SW')
daily_r <- import_feature_from_folder('rast',nc_folder, '.nc')
nc_paths <- import_feature_from_folder(NULL,nc_folder, '.nc')
monthly_clim <-  Map(function(x,y)
    monthly_clim_roi_from_daily_raster(nc_paths[grepl(y,nc_paths)],
                                       ## kryck_border,
                                       kryck_border_for_coarse_climate,
                                       ## vect(kryck_border_for_coarse_climate),
                                        fun = x,
                                        t2 = '2020-12-31'),
    c('sum','min','mean','max','frost_days'), c('_pr_', rep('_tas_', 4)))
names(monthly_clim) <- c('prcp','tmp_min','tmp_ave','tmp_max','frost_days')
monthly_clim. <- Map(function(x)
    crop_raster_fast(x, kryck_border), monthly_clim)
monthly_clim.. <- Map(function(x)
                      agg_and_crop..(x, 10, kryck_border), monthly_clim.)
monthly_clim_roi_proj.  <- Map(function(x)
    align_to_ref(x, kryck_elev, method = 'near'), monthly_clim..)
## Solar Radiation
md <- file.path(base_path, 'MODIS_DSR_Kricklan_solar_radiation')
drs_rast <- import_feature_from_folder('rast', md, 'DSR')
names(drs_rast)  <- drs <- dir(md)[grepl('DSR', dir(md))] %>%
    str_replace('DSR_','X') %>% str_remove('\\.tif$') %>% str_replace('_','.')
monthly_rad_roi_proj  <- align_to_ref(drs_rast, kryck_elev)
time(monthly_rad_roi_proj) <- names(monthly_rad_roi_proj) %>% str_remove('X') %>% str_replace('\\.','-') %>% ym() 
rasters <- c(monthly_clim_roi_proj., monthly_rad_roi_proj)
blueprint <- 'd_climate'
colnames. <- c(names(monthly_clim_roi_proj.),'srad')
test_clim <- d_input_from_rast_back(rasters = rasters,
                                      colnames = colnames.,
                                 blueprint = blueprint,
                                      from = '2001-01',
                               to = '2020-12')#,
## Texture map
text_folder <- file.path(base_path,'ESDB-Raster-Library-1k-GeoTIFF-20240507/TEXT')
text_dt <- import_feature_from_folder('rast', text_folder, '.tif')
kryck_text. <- crop_raster_fast(text_dt, kryck_border)
kryck_text <- agg_and_crop.(kryck_text., 10, kryck_border) 
kryck_text_proj  <- align_to_ref(kryck_text, kryck_elev)
tex_matrix <- matrix(c(
  0, 0,   # no info → others
  1, 1,   # coarse → sandy
  2, 2,   # medium → sandy loam
  3, 3,   # medium fine → clay loam
  4, 4,   # fine → clay
  5, 4,   # very fine → clay
  6, 0,   # no texture → others
  7, 0,   # no texture → others
  8, 0    # no texture → others
), ncol = 2, byrow = TRUE)
kryck_text_proj. <- classify(kryck_text_proj[[1]], tex_matrix)
## soil water
asw_folder <- file.path(base_path,'ESDB-Raster-Library-1k-GeoTIFF-20240507/AWC_TOP')
asw_dt <- import_feature_from_folder('rast', asw_folder, '.tif')
kryck_asw. <- crop_raster_fast(asw_dt, kryck_border)
kryck_asw <- agg_and_crop.(kryck_asw., 10, kryck_border) 
kryck_asw_proj..  <- align_to_ref(kryck_asw, kryck_elev)
asw_matrix <- matrix(c(
  1, 75,
  2, 120,
  3, 165,
  4, 210,
  5, NA
), ncol = 2, byrow = TRUE)
kryck_asw_proj. <- classify(kryck_asw_proj..[[1]], asw_matrix)
## raster of latitudes
kryck_lat <- latitude_raster(kryck_elev)
## age ranges
from_age <- make_date_raster(kryck_elev, '2001-01')
to_age <- make_date_raster(kryck_elev, '2020-12')
## data list of sites
sites <- c(kryck_lat, kryck_elev, kryck_text_proj., kryck_asw_proj.,
  kryck_asw_proj. - kryck_asw_proj., kryck_asw_proj. + 30, from_age, to_age)
sites. <- crop(sites, kryck_border)
sites_nm <- c('latitude','altitude',
              'soil_class', 'asw_i',
              'asw_min', 'asw_max',
              'from','to')
names(sites.) <- sites_nm
    test_site  <- d_input_from_rast_back(terra::as.list(sites.),
                                         names(sites.), 'd_site')
## Forest species
land_folder <- file.path(base_path,'nationalndcoverdata')
land_dt <- import_feature_from_folder('rast', land_folder, '.tif$')
kryck_land. <- crop_raster_fast(land_dt, kryck_border)
kryck_land_proj  <- align_to_ref(kryck_land., kryck_elev, method = 'mode')
lookup <- read.csv(file.path(land_folder,'spp_codes_from_land.csv'), sep = ',') %>%
select(code, starts_with('Layer'))
r_kryck_spp <- reclass_char_stack(kryck_land_proj, lookup)
## Age map
age_folder <- file.path(base_path, 'Sweden_Forest_Age')
age_dt <- import_feature_from_folder('rast', age_folder, '(PINE|SPRUCE).*\\.tif$')
age_dt. <- import_feature_from_folder('rast', age_folder, '(Birch).*\\.tif$')
kryck_age  <- crop_raster_with_polygon(age_dt, kryck_border, mask = FALSE, buffer_width = 0)
kryck_age.  <- crop_raster_with_polygon(age_dt., kryck_border, mask = FALSE, buffer_width = 0)
kryck_age_proj <- align_to_ref(c(kryck_age, kryck_age.), kryck_elev)
# Create a raster of plantation dates as "YYYY-mm"
## next code chunk is very slow
ref_date <- ymd("2025-07-01")  # lubridate date
kryck_planted  <- compute_plantation_raster(kryck_age_proj, ref_date)
## Fertilization raster
kryck_fert <- kryck_elev
kryck_fert[] <- 0.5
## Number of stems
lookup_stem <- read.csv(file.path(land_folder,'spp_codes_from_land.csv'), sep = ',') %>%
select(code, starts_with('Planted'))
kryck_stems <- reclass_stack(kryck_land_proj, lookup_stem)
names(kryck_stems)  <- names(r_kryck_spp)
## Biomass
biom_folder <- file.path(base_path, 'biomass')
biom_dt <- import_feature_from_folder('rast', biom_folder, '.*\\.tif$')
kryck_biom <- crop_raster_fast(biom_dt, kryck_border)
kryck_biom_proj <- align_to_ref(kryck_biom, kryck_elev)
lookup_props_ <- read.csv(file.path(biom_folder,'biom_stocks_from_Hans_Peterson.csv'),
                         sep = ',') %>%
    select(Tree_parts,ends_with('prop')) %>% slice(-4)
lookup_props <- lookup_props_ %>% pivot_longer(cols = -1) %>%
    pivot_wider(names_from = Tree_parts, values_from = value)  %>% rename(code = name)
prop <- Map(function(x) replace_species_with_proportions(
                            lookup,
                            lookup_props[lookup_props$code == x,]),
            lookup_props$code)
biom_perc <- Map(function(x)
    reclass_stack(kryck_land_proj, x), prop)
kryck_biom_comps. <- Map(function(x)
           kryck_biom_proj * x, as.list(biom_perc))
kryck_biom_comps <- lapply(kryck_biom_comps., function(x){names(x) <- names(r_kryck_spp)
    return(x)})
## d_species inputs
specs. <- c(list(rep(kryck_planted[[1L]],3),r_kryck_spp, c(rep(kryck_fert,3)), kryck_stems), kryck_biom_comps)
nm_specs.  <- c('planted','species','fertility','stems_n', 'biom_stem','biom_foliage', 'biom_root')
names(specs.) <- nm_specs.
test_species  <- d_input_from_rast_back(specs., nm_specs., 'd_species')
## Parameters input
d_parameters_adapted <- get_parameters(mode = 'parameters', sp_names = c('Betula', 'Picea abies', 'Pinus sylvestris 30')) %>% rename_with(~ str_remove(.x, "\\s*\\d+$"))
d_parameters_adapted[,'Betula']  <- d_parameters_adapted[,'Pinus sylvestris'] 
test_parameters <- Map(function(x)
prepare_parameters(parameters = d_parameters_adapted, sp_names = x$'species'), test_species)
## Size-distribution input
d_sizeDist_adapted <- get_parameters(mode = 'sizeDist', sp_names = c('Betula', 'Picea abies', 'Pinus sylvestris 30')) %>% rename_with(~ str_remove(.x, "\\s*\\d+$"))
test_sizeDist <- Map(function(x)
    prepare_sizeDist(size_dist = d_sizeDist_adapted, sp_names = x$'species'), test_species)
## list of data inputs
all_d <- list(test_site, test_species, test_clim, test_parameters, test_sizeDist) 
names(all_d)  <- paste0('d_', c('site', 'species', 'climate', 'parameters', 'sizeDist'))
## common and not-common elements in the data lists
common_names <- all_d %>% map(names) %>% reduce(intersect)
not_common_names  <- all_d %>% map(names) %>% reduce(union) %>% setdiff(common_names)
filtered_ls  <- all_d %>% map(~ .x[common_names])
marg.  <- list(function(t,w,x,y,z)
    tryCatch(run_3PG(site = t,
                     species = w,
                     climate = x,
                     thinning = NULL,
                     parameters = y,
                     size_dist = z,
                     settings  = list(
                         light_model = 2,
                         transp_model = 2,
                         phys_model = 2, 
                         height_model = 1,
                         correct_bias = 0,
                         calculate_d13c = 0),
                     check_input = TRUE,
                     df_out = TRUE),
             error = function(e) NULL),
    t = filtered_ls[['d_site']],
    w = filtered_ls[['d_species']],
    x = filtered_ls[['d_climate']],
    y = filtered_ls[['d_parameters']],
    z = filtered_ls[['d_sizeDist']])## ,
inps <- do.call('Map', marg.)
## select GPP outputs & summarise mixed species
pg_lst <- inps %>% map(function(x) filter(x, variable %in% "gpp")) %>%
    map(~ .x %>% group_by(date) %>%
            summarise(species = if (n_distinct(species) == 1) first(species) else "mixed",
value = mean(value, na.rm = TRUE), .groups = 'drop'))
# Step 1: extract first species name per pixel
species_vec <- map_chr(pg_lst, ~ .x$species[1])
# Step 2: make a r aster with the same dimensions as a reference raster
ref_raster <- kryck_elev
species_rast <- rast(ref_raster)
# Step 3: convert to factor and assign numeric codes
species_fac <- as.factor(species_vec)
values(species_rast) <- as.integer(species_fac)
# Step 4: add categorical levels (species names as labels)
levels(species_rast) <- data.frame(
  ID = seq_along(levels(species_fac)),
  species = levels(species_fac)
)
## Create and print time-series taster 
tsr <- timeSeriesRaster(kryck_stems, pg_lst, names(pg_lst))
names(tsr) <- time(tsr)
end.time <- Sys.time()
save_timeseries_raster(tsr, file.path(base_path, 'outputs_3PG'))
time.taken <- end.time - start.time
time.taken
## global_stats(tmp)
## tmp <- rast("/media/wihe/LaCie/projects/3PG/outputs_3PG/3PG_735034_7126998.tif")
## plot(tmp[[239]])


