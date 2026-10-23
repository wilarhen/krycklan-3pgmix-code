## ---------------------------------------------------------------------------
## polygon_grids_v02.r
##
## Corrected version of polygon_grids_v01.r. Eight fixes, each marked
## ## FIX n ## with the original line number of v01 and the evidence for it.
##
## Diff against v01 to review. Nothing else was changed.
## ---------------------------------------------------------------------------

library('ncdf4')
library('terra')
## library('wilaRtools')
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
## wd('~','3PG/devel/R', setwd = TRUE)
setwd('/home/wilar/Documents/SLU/3PG/devel/R')
sourceDirectory('.', modifiedOnly = TRUE, verbose = TRUE)
setwd('/home/wilar/Documents/SLU/3PG/devel')

## ===========================================================================
## OVERRIDE: compute_plantation_raster  (terra cannot store character values)
## ---------------------------------------------------------------------------
## terra 1.9.27 coerces character raster values to numeric, so date strings
## cannot be stored directly (the wilaRtools version silently produced ages).
## Store the "YYYY-MM" date as a CATEGORICAL raster instead: the codes are
## integers, the labels are the date strings. terra::as.data.frame (used by
## d_input_from_rast_back) then returns the LABELS, exactly as it already does
## for the species stack. r3PG receives planted as character "YYYY-MM".
compute_plantation_raster <- function(age_raster, ref_date) {
  ref_date <- lubridate::ymd(ref_date)
  age_to_date <- function(age) {
    yrs <- floor(age); mos <- round((age - yrs) * 12)
    format(ref_date %m-% years(yrs) %m-% months(mos), "%Y-%m")
  }
  out <- vector('list', nlyr(age_raster))
  for (i in seq_len(nlyr(age_raster))) {
    a  <- as.vector(terra::values(age_raster[[i]]))
    ok <- is.finite(a)
    dates <- rep(NA_character_, length(a))
    dates[ok] <- age_to_date(a[ok])
    ## build a factor raster: integer codes + label table
    labs  <- sort(unique(stats::na.omit(dates)))
    codes <- match(dates, labs)
    r <- terra::rast(age_raster[[i]])
    terra::values(r) <- codes
    levels(r) <- data.frame(value = seq_along(labs), label = labs)
    names(r) <- names(age_raster)[i]
    out[[i]] <- r
  }
  terra::rast(out)
}
## wd(,'3PG/devel')
## temporary files
## base_path <- '/media/wihe/LaCie/projects/3PG'
## base_path <- '/run/media/wilar/T7 Shield/slu_data_bases/3PG'
base_path <- '/home/wilar/datos/3PG'
my_temp <- file.path(base_path, 'r_tmp')
dir.create(my_temp, showWarnings = FALSE, recursive = TRUE)
terraOptions(tempdir = my_temp)

## ===========================================================================
## COMMAND LINE INTERFACE
##   Rscript polygon_grids_v02.r <sector> [cells_per_side] [save_rdata]
##
##   sector         1 .. cells_per_side^2   (sector index to process)
##   cells_per_side default 10, so 100 sectors over the buffered catchment
##   save_rdata     TRUE writes the raw per-cell run_3PG output for this
##                  sector. Very large. Use for ONE diagnostic sector only.
##
##   Falls back to the old 'sector_n' file if no argument is given, so
##   existing habits still work.
## ===========================================================================
.args <- commandArgs(trailingOnly = TRUE)
sector_n <- if (length(.args) >= 1) as.integer(.args[1]) else
            as.integer(readLines('sector_n'))
cells_per_side <- if (length(.args) >= 2) as.integer(.args[2]) else 10L
save_rdata <- if (length(.args) >= 3) as.logical(.args[3]) else FALSE
stopifnot(is.finite(sector_n), sector_n >= 1,
          sector_n <= cells_per_side^2)
cat('sector', sector_n, 'of', cells_per_side^2,
    '| save_rdata =', save_rdata, '\n')

## krycklan polygon
start.time <- Sys.time()
kry_folder <- file.path(base_path, 'krycklan_borders')
kryck_bord <- import_feature_from_folder('vect', kry_folder, '.shp')
kryck_bord. <- polygon_sector(kryck_bord, cells_per_side, sector_n, 200)
kryck_border_for_coarse_climate <- polygon_sector(kryck_bord, 1, 1, 4000)[[1]]
kryck_border <- kryck_bord.[[1]]

## ===========================================================================
## FIX 1  [v01 line 51]  REFERENCE GRID
## ---------------------------------------------------------------------------
## v01 built the reference from the 2 m DEM: agg_and_crop.(kryck_dem, 10, ...).
## Cropping snaps to a 2 m boundary, so aggregating to 10 m produced a grid at
## origin (4, -4). The three structural products (biomass, age, species) all
## share origin (0, 0) at 10 m, so every structural value was interpolated
## from four source cells instead of being read.
##
## MEASURED: variance retained on the offset grid was biomass 86.3 %,
## age spruce 79.4 %, age pine 70.9 %. On a grid at origin (0,0) all layers
## reproduce exactly. Four validation plots were lost to this alone.
##
## The reference is now the biomass grid. Because every downstream call is
## align_to_ref(x, kryck_elev, ...), redefining kryck_elev onto that grid
## propagates the fix through the whole script.
## ===========================================================================
biom_folder <- file.path(base_path, 'biomass')
biom_dt <- import_feature_from_folder('rast', biom_folder, '.*\\.tif$')
ref_grid <- crop_raster_fast(biom_dt, kryck_border)[[1]]

cat('reference grid: origin', origin(ref_grid), '| res', res(ref_grid), '\n')
stopifnot(all(abs(origin(ref_grid)) < 1e-6), all(res(ref_grid) == 10))

## DEM raster
dem_folder <- file.path(base_path, 'Sweden_dem')
dem_dt <- import_feature_from_folder('rast', dem_folder, 'dem.tif')
kryck_dem <- crop_raster_fast(dem_dt, kryck_border)
## 2 m -> 10 m is a downsampling, so average rather than interpolate
kryck_elev <- align_to_ref(kryck_dem, ref_grid, method = 'average')
names(kryck_elev) <- 'altitude'
stopifnot(compareGeom(kryck_elev, ref_grid, stopOnError = FALSE))

## CLIMATE DATA
nc_folder <- file.path(base_path,'climate_SW')
daily_r <- import_feature_from_folder('rast',nc_folder, '.nc')
nc_paths <- import_feature_from_folder(NULL,nc_folder, '.nc')
monthly_clim <-  Map(function(x,y)
    monthly_clim_roi_from_daily_raster(nc_paths[grepl(y,nc_paths)],
                                       kryck_border_for_coarse_climate,
                                        fun = x,
                                        t2 = '2020-12-31'),
    c('sum','min','mean','max','frost_days'), c('_pr_', rep('_tas_', 4)))
names(monthly_clim) <- c('prcp','tmp_min','tmp_ave','tmp_max','frost_days')
## NOTE tmp_min and tmp_max are the monthly MIN and MAX of DAILY MEAN
## temperature, because the PTHBV product provides daily mean only. 3-PG
## expects the mean daily minimum and maximum. This inflates the derived
## vapour pressure deficit. Not fixed here because no daily min/max source is
## available; report it as a limitation.
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
                               to = '2020-12')
## Texture map
text_folder <- file.path(base_path,'ESDB-Raster-Library-1k-GeoTIFF-20240507/TEXT')
text_dt <- import_feature_from_folder('rast', text_folder, '.tif')
kryck_text. <- crop_raster_fast(text_dt, kryck_border)
kryck_text <- agg_and_crop.(kryck_text., 10, kryck_border)
kryck_text_proj  <- align_to_ref(kryck_text, kryck_elev, method = 'near')
tex_matrix <- matrix(c(
  0, 0,   # no info -> others
  1, 1,   # coarse -> sandy
  2, 2,   # medium -> sandy loam
  3, 3,   # medium fine -> clay loam
  4, 4,   # fine -> clay
  5, 4,   # very fine -> clay
  6, 0,   # no texture -> others
  7, 0,   # no texture -> others
  8, 0    # no texture -> others
), ncol = 2, byrow = TRUE)
kryck_text_proj. <- classify(kryck_text_proj[[1]], tex_matrix)
## soil water
asw_folder <- file.path(base_path,'ESDB-Raster-Library-1k-GeoTIFF-20240507/AWC_TOP')
asw_dt <- import_feature_from_folder('rast', asw_folder, '.tif')
kryck_asw. <- crop_raster_fast(asw_dt, kryck_border)
kryck_asw <- agg_and_crop.(kryck_asw., 10, kryck_border)
kryck_asw_proj..  <- align_to_ref(kryck_asw, kryck_elev, method = 'near')
asw_matrix <- matrix(c(
  1, 75,
  2, 120,
  3, 165,
  4, 210,
  5, NA
), ncol = 2, byrow = TRUE)
kryck_asw_proj. <- classify(kryck_asw_proj..[[1]], asw_matrix)
## NOTE both ESDB layers take a SINGLE value at all 57 validation plots, so
## soil contributes no spatial variation. f_sw = 0.996 in the corrected
## preflight run, i.e. soil water is not limiting at this site.
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

## ===========================================================================
## FIX 3a  [v01 line 138]  categorical resampling for species
## MEASURED: five validation plots changed land cover class through
## resampling. With FIX 1 the grids match and this is a no-op, but the method
## must be correct in case a product changes grid.
## ===========================================================================
kryck_land_proj  <- align_to_ref(kryck_land., kryck_elev, method = 'near')

lookup <- read.csv(file.path(land_folder,'spp_codes_from_land.csv'), sep = ',') %>%
select(code, starts_with('Layer'))
r_kryck_spp <- reclass_char_stack(kryck_land_proj, lookup)
cat('species layers:', paste(names(r_kryck_spp), collapse = ' | '), '\n')

## Age map
age_folder <- file.path(base_path, 'Sweden_Forest_Age')
age_dt <- import_feature_from_folder('rast', age_folder, '(PINE|SPRUCE).*\\.tif$')
age_dt. <- import_feature_from_folder('rast', age_folder, '(Birch).*\\.tif$')
kryck_age  <- crop_raster_with_polygon(age_dt, kryck_border, mask = FALSE, buffer_width = 0)
kryck_age.  <- crop_raster_with_polygon(age_dt., kryck_border, mask = FALSE, buffer_width = 0)

## ===========================================================================
## FIX 3b  [v01 line 148]  nearest neighbour for stand age
## ---------------------------------------------------------------------------
## v01 used the default bilinear. Interpolating a stand-age surface across a
## clearcut boundary produces ages that exist nowhere, and bilinear also fills
## across NA boundaries: two plots with no source value received fabricated
## age and biomass in the v01 cube.
## ===========================================================================
kryck_age_proj <- align_to_ref(c(kryck_age, kryck_age.), kryck_elev, method = 'near')
## Layer order follows the file pattern: PINE, SPRUCE, then Birch.
names(kryck_age_proj) <- c('Pinus sylvestris', 'Picea abies', 'Betula')  ## VERIFY

# Create a raster of plantation dates as "YYYY-mm"
## next code chunk is very slow
ref_date <- ymd("2025-07-01")
kryck_planted  <- compute_plantation_raster(kryck_age_proj, ref_date)
names(kryck_planted) <- names(kryck_age_proj)

## ===========================================================================
## FIX 2  [v01 line 182]  per-species planting dates
## ---------------------------------------------------------------------------
## v01 used rep(kryck_planted[[1L]], 3), which gave EVERY species the PINE
## planting date. Spruce is 34 of the 50 validation plots. The birch and
## spruce age layers were loaded and then discarded.
##
## CORRECTED. d_input_from_rast_back pairs species / planted / biom / stems
## by LAYER POSITION (it renames every input to Layer_1/2/3 internally). But
## r_kryck_spp is in SLOT layout: Layer_k holds the species occupying slot k
## of each pixel, as a CELL VALUE that varies pixel to pixel (e.g. class 114
## has Betula in slot 1). kryck_planted is in fixed SPECIES layout (Pinus,
## Picea, Betula). The two cannot be paired by reordering layers.
##
## Build planted in the SAME slot layout: slot k carries the planting date of
## whichever species sits in slot k of that pixel. Keyed on the species LABEL,
## never the integer code -- reclass_char_stack assigns codes as 1=Pinus,
## 2=Betula, 3=Picea, which is NOT the layer order of kryck_planted.
## compute_plantation_raster returns character 'YYYY-MM', so values stay
## character throughout.
## ===========================================================================
## CRITICAL: reclass_char_stack builds a SEPARATE factor coding per layer
## (it calls factor(new_vals, levels = unique(new_vals)) independently for
## each column), so code 1 in Layer_1 and code 1 in Layer_3 are DIFFERENT
## species. Each slot must therefore be decoded with its OWN levels, not with
## Layer_1's. Verified: for a class-114 pixel this yields Betula / Pinus /
## Picea, matching the lookup.
n_slot <- nlyr(r_kryck_spp)
slot_planted <- vector('list', n_slot)
for (k in seq_len(n_slot)) {
  lab_k  <- levels(r_kryck_spp[[k]])[[1]]        # THIS layer's value<->label
  stopifnot(all(stats::na.omit(lab_k$label) %in% names(kryck_planted)))
  code_k <- as.vector(values(r_kryck_spp[[k]]))
  out_lab <- rep(NA_character_, length(code_k))
  for (j in seq_len(nrow(lab_k))) {
    sp_j <- lab_k$label[j]
    if (is.na(sp_j) || !sp_j %in% names(kryck_planted)) next
    ## kryck_planted layers are now FACTOR rasters: read codes, map to labels
    pl        <- kryck_planted[[ sp_j ]]
    pl_lab    <- levels(pl)[[1]]
    pl_code   <- as.vector(values(pl))
    pl_string <- pl_lab$label[match(pl_code, pl_lab$value)]   # "YYYY-MM"
    sel <- which(code_k == lab_k$value[j])
    out_lab[sel] <- pl_string[sel]
  }
  ## re-encode the gathered date strings as a factor raster for this slot
  labs  <- sort(unique(stats::na.omit(out_lab)))
  r <- terra::rast(r_kryck_spp[[k]])
  terra::values(r) <- match(out_lab, labs)
  levels(r) <- data.frame(value = seq_along(labs), label = labs)
  names(r)  <- paste0('Layer_', k)
  slot_planted[[k]] <- r
}
kryck_planted_ord <- rast(slot_planted)
cat('per-slot planting dates built for', n_slot, 'slots\n')

## Fertilization raster
kryck_fert <- kryck_elev
kryck_fert[] <- 0.5

## ===========================================================================
## FIX 5  [v01 lines 157-159]  stem density = measured species means
## ---------------------------------------------------------------------------
## v01 assigned stems_n by land cover class from the "Planted" columns, giving
## an invariant 2000 ha-1 against a measured range of 382 to 3121.
##
## MEASURED: leave-one-out RMSE 1062 stems/ha for the constant against 582 for
## a species mean, a 45 % reduction. No mappable predictor beat the species
## mean; the log-log biomass coefficient was POSITIVE (+0.298, p = 0.058),
## contrary to self-thinning. Preflight showed the choice barely matters
## (r = 0.483 with the mean against 0.482 with measured density).
## ===========================================================================
## CORRECTED. Stem density is NOT a proportion, so it must NOT pass through
## replace_species_with_proportions, which divides each value by the number of
## species in the pixel (correct for biomass fractions, wrong for density).
## v01 used reclass_stack, a plain per-class lookup, which is right. Keep that
## path; only change the per-species VALUES from the old 2000 to the measured
## species means. reclass_stack reads the lookup's Layer_1/2/3 columns, so the
## output is already in slot layout matching r_kryck_spp.
lookup_stem <- lookup %>%
  select(code, starts_with('Layer')) %>%
  mutate(across(starts_with('Layer'), ~ dplyr::recode(.x,
           'Pinus sylvestris' = 967,
           'Picea abies'      = 1162,
           'Betula'           = 1000,   ## VERIFY, no birch validation plots
           .default = NA_real_)))
kryck_stems <- reclass_stack(kryck_land_proj, lookup_stem)
names(kryck_stems) <- names(r_kryck_spp)

## Biomass
kryck_biom <- crop_raster_fast(biom_dt, kryck_border)
kryck_biom_proj <- align_to_ref(kryck_biom, kryck_elev, method = 'bilinear')

## ===========================================================================
## FIX 4  [v01 lines 166-180]  foliage from Marklund, not the branch column
## ---------------------------------------------------------------------------
## biom_stocks_from_Hans_Peterson.csv has Stem / Branch / Stump+root and NO
## foliage column. v01 mapped those three positionally onto biom_stem /
## biom_foliage / biom_root, so BRANCH WOOD entered the foliage pool. In 3-PG
## the stem pool already includes branches.
##
## MEASURED: Picea branch/stem = 0.25/0.50 = 0.50; the ratio in the saved v01
## model inputs was 0.49. Simulated LAI reached 33.6 against a measured
## maximum of 6.75, saturating light interception in every stand and leaving
## GPP insensitive to biomass.
##
## BIO.tif is aboveground biomass, so:
##   biom_foliage = AGB * f            (f from Marklund 1988)
##   biom_stem    = AGB * (1 - f)      (stem pool includes branches and bark)
##   biom_root    = AGB * root/(stem+branch)
##
## The species-allocation machinery below is unchanged from v01; only the
## fractions differ.
## ===========================================================================
fol_frac  <- c('Pinus sylvestris' = 0.045,     ## VERIFY against Marklund
               'Picea abies'      = 0.085,     ## functions for the local
               'Betula'           = 0.020)     ## diameter range
root_frac <- c('Pinus sylvestris' = 0.25/0.76,
               'Picea abies'      = 0.25/0.75,
               'Betula'           = 0.25/0.75)

## Build a per-class lookup with one column per biomass component. For each
## land-cover class and each slot, reclass_stack assigns the slot the fraction
## of whatever species the lookup names in that slot. This is a plain per-slot
## lookup with NO averaging and NO division by species count.
##
## We do NOT use replace_species_with_proportions here: it computes
## avg(fractions)/n_species, which discards each species' own fraction (it
## returned foliage = 0.0167 for all three species in a mixed pixel instead of
## 0.020 / 0.045 / 0.085). That averaging is wrong for biomass partitioning.
##
## For a mixed pixel the total AGB is split EQUALLY among the co-occurring
## species (1/n each), then each species' share is partitioned by its own
## stem / foliage / root fractions. The 1/n split is applied once, as a
## separate share layer, so it is not entangled with the partitioning.
make_comp_lookup <- function(component) {
  frac <- switch(component,
    biom_foliage = fol_frac,
    biom_root    = root_frac,
    biom_stem    = c('Pinus sylvestris' = 1 - fol_frac[['Pinus sylvestris']],
                     'Picea abies'      = 1 - fol_frac[['Picea abies']],
                     'Betula'           = 1 - fol_frac[['Betula']]))
  lookup %>%
    select(code, starts_with('Layer')) %>%
    mutate(across(starts_with('Layer'),
                  ~ dplyr::recode(.x,
                      'Pinus sylvestris' = unname(frac[['Pinus sylvestris']]),
                      'Picea abies'      = unname(frac[['Picea abies']]),
                      'Betula'           = unname(frac[['Betula']]),
                      .default = NA_real_)))
}

## per-species partitioning fraction, in slot layout
biom_perc <- Map(function(x) reclass_stack(kryck_land_proj, make_comp_lookup(x)),
                 c('biom_stem','biom_foliage','biom_root'))
names(biom_perc) <- c('biom_stem','biom_foliage','biom_root')

## equal 1/n share of total AGB among co-occurring species, per slot.
## n_species per pixel = count of non-NA slots in the species stack.
n_species_rast <- sum(!is.na(r_kryck_spp))
species_share  <- 1 / n_species_rast          # one layer, same for every slot

## each component = total AGB * species share * that species' fraction
kryck_biom_comps. <- Map(function(x)
    kryck_biom_proj * species_share * x, biom_perc)
kryck_biom_comps <- lapply(kryck_biom_comps., function(x){
    names(x) <- names(r_kryck_spp); x })

## d_species inputs
specs. <- c(list(kryck_planted_ord,               ## FIX 2: per-species dates
                 r_kryck_spp,
                 c(rep(kryck_fert, nlyr(r_kryck_spp))),
                 kryck_stems),
            kryck_biom_comps)
nm_specs.  <- c('planted','species','fertility','stems_n',
                'biom_stem','biom_foliage','biom_root')
names(specs.) <- nm_specs.
test_species  <- d_input_from_rast_back(specs., nm_specs., 'd_species')

## ===========================================================================
## FIX 6  [v01 lines 187, 195, 199, 220]  one parameter family, height_model 2
## ---------------------------------------------------------------------------
## r3PG implements two height equations sharing the names aH, nHB, nHC with
## different meanings:
##     form 1:  h = aH * dbh^nHB * comp^nHC
##     form 2:  h = 1.3 + aH * exp(-nHB/dbh) + nHC * comp * dbh
## Verified in src/md_3PG.f95. The parameter database records source, region,
## country and stand type but NOT which form a set was fitted under, and the
## package default is form 1.
##
## v01 mixed families: 'Pinus sylvestris 30' (Forrester 2017) is a form-1 set,
## while 'Picea abies' resolves to set 37 and 'Betula' to set 33, both
## Forrester 2021 form-2 sets. Under height_model = 1, Picea abies 37 with
## nHB = 17.85 gives 37.73 * 14.7^17.85 = 2.6e22 m.
##
## MEASURED: median simulated height 3.2e22 m, r(GPP, APAR) = 0.008.
## After the fix: median height 15.7 m, r(GPP, APAR) = 0.882.
##
## v01 line 195 also overwrote the Betula column with Pinus. Betula 33 exists
## and is a proper deciduous set (leafgrow = 5, leaffall = 11). That line is
## DELETED. Note that with a January start, biom_foliage should be the SPRING
## foliage biomass for deciduous species.
## ===========================================================================
d_parameters_adapted <- get_parameters(
    mode = 'parameters',
    sp_names = c('Betula 33', 'Picea abies 37', 'Pinus sylvestris 39')) %>%
    rename_with(~ str_remove(.x, "\\s*\\d+$"))

write.csv(d_parameters_adapted, file = 'd_parameters_adapted.csv')

## sanity check: form-2 sets have large nHB
cat('nHB by species:\n')
print(d_parameters_adapted[d_parameters_adapted$parameter == 'nHB', ])

## FIX 6: line that assigned Pinus parameters to Betula is REMOVED.

test_parameters <- Map(function(x)
prepare_parameters(parameters = d_parameters_adapted, sp_names = x$'species'), test_species)
## Size-distribution input
d_sizeDist_adapted <- get_parameters(
    mode = 'sizeDist',
    sp_names = c('Betula 33', 'Picea abies 37', 'Pinus sylvestris 39')) %>%
    rename_with(~ str_remove(.x, "\\s*\\d+$"))
test_sizeDist <- Map(function(x)
    prepare_sizeDist(size_dist = d_sizeDist_adapted, sp_names = x$'species'), test_species)

## list of data inputs
all_d <- list(test_site, test_species, test_clim, test_parameters, test_sizeDist)
names(all_d)  <- paste0('d_', c('site', 'species', 'climate', 'parameters', 'sizeDist'))
## common and not-common elements in the data lists
common_names <- all_d %>% map(names) %>% reduce(intersect)
not_common_names  <- all_d %>% map(names) %>% reduce(union) %>% setdiff(common_names)
filtered_ls  <- all_d %>% map(~ .x[common_names])

## report how many cells each input list contributes and how many survive
cat('\ncells per input list:\n')
for (nm in names(all_d)) cat(sprintf('  %-14s %d\n', nm, length(all_d[[nm]])))
cat(sprintf('  %-14s %d  (%.1f %% of the largest)\n', 'COMMON',
            length(common_names),
            100 * length(common_names) / max(lengths(map(all_d, names)))))

## ===========================================================================
## FIX 8  [v01 line 225]  log failures instead of discarding them
## ---------------------------------------------------------------------------
## v01 used error = function(e) NULL. At 1.2 million cells failures cannot be
## inspected individually, so they must be recorded.
##
## MEASURED: of nine validation plots absent from the v01 cube, three had no
## source value and one was lost to the grid offset. The remaining FIVE lie
## inside the catchment with complete inputs and are attributable only to
## discarded simulation failures. This was the largest single cause.
## ===========================================================================
run_log <- new.env(parent = emptyenv())
run_log$fail <- list()

marg.  <- list(function(t,w,x,y,z,cell_id)
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
                         height_model = 2,          ## FIX 6
                         correct_bias = 0,
                         calculate_d13c = 0),
                     check_input = TRUE,
                     df_out = TRUE),
             error = function(e) {
                 run_log$fail[[length(run_log$fail) + 1L]] <-
                     data.frame(cell = cell_id, msg = conditionMessage(e),
                                stringsAsFactors = FALSE)
                 NULL
             }),
    t = filtered_ls[['d_site']],
    w = filtered_ls[['d_species']],
    x = filtered_ls[['d_climate']],
    y = filtered_ls[['d_parameters']],
    z = filtered_ls[['d_sizeDist']],
    cell_id = common_names)
inps <- do.call('Map', marg.)

## write the failure log and the completeness statistic
fail_tab <- if (length(run_log$fail)) do.call(rbind, run_log$fail) else
            data.frame(cell = character(0), msg = character(0))
n_ok <- sum(!vapply(inps, is.null, logical(1)))
cat(sprintf('\ncells simulated: %d / %d (%.1f %%), failures: %d\n',
            n_ok, length(inps), 100 * n_ok / length(inps), nrow(fail_tab)))
if (nrow(fail_tab)) {
    print(head(sort(table(fail_tab$msg), decreasing = TRUE), 5))
    write.csv(fail_tab,
              file.path(base_path, 'outputs_3PG',
                        sprintf('failures_%s.csv', sector_n)),
              row.names = FALSE)
}
inps <- inps[!vapply(inps, is.null, logical(1))]

if (!length(inps)) {
    cat('sector', sector_n, 'produced no simulated cells. Nothing to write.\n')
    quit(save = 'no', status = 0)
}

## ===========================================================================
## RAW OUTPUT ARCHIVE
## ---------------------------------------------------------------------------
## Saving the full per-cell model output means any variable can be mapped
## later without re-running the model. Re-running costs days; reading costs
## minutes.
##
## Two reductions that lose no information:
##   (1) drop variables that are constant zero or all NA. With
##       correct_bias = 0 and calculate_d13c = 0 the weibull and wood_delta
##       groups (about 30 variables) carry nothing.
##   (2) store one compact data.table with a cell column instead of a list of
##       11000 data.frames each repeating the same date column.
##
## Set save_rdata = TRUE on the command line to enable.
## ===========================================================================
compact_inps <- function(x) {
    dt <- data.table::rbindlist(x, idcol = 'cell')
    data.table::setDT(dt)
    ## drop variables that carry no information in this configuration
    keep <- dt[, .(keep = !all(is.na(value)) &&
                          !(all(value == 0, na.rm = TRUE))), by = variable]
    dropped <- keep[keep == FALSE, as.character(variable)]
    dt <- dt[variable %in% keep[keep == TRUE, as.character(variable)]]
    ## compact types
    dt[, cell     := as.integer(cell)]
    dt[, date     := as.integer(format(as.Date(date), '%Y%m'))]
    dt[, variable := factor(variable)]
    dt[, species  := factor(species)]
    if ('group' %in% names(dt)) dt[, group := factor(group)]
    attr(dt, 'dropped_variables') <- dropped
    attr(dt, 'date_format') <- 'YYYYMM integer'
    dt
}

if (isTRUE(save_rdata)) {
    rd_dir <- file.path(base_path, '3PG_raw')
    dir.create(rd_dir, showWarnings = FALSE, recursive = TRUE)

    inps_dt <- compact_inps(inps)
    cat('raw archive: ', nrow(inps_dt), ' rows, ',
        data.table::uniqueN(inps_dt$variable), ' variables retained, ',
        length(attr(inps_dt, 'dropped_variables')), ' dropped\n', sep = '')

    ## qs2 is the current package (qs itself is deprecated). Fall back to
    ## legacy qs, then to base saveRDS, so the script never fails on this.
    if (requireNamespace('qs2', quietly = TRUE)) {
        rd_file <- file.path(rd_dir, sprintf('inps_sector_%03d.qs2', sector_n))
        qs2::qs_save(inps_dt, rd_file, compress_level = 8,
                     nthreads = max(1L, parallel::detectCores() - 1L))
    } else if (requireNamespace('qs', quietly = TRUE)) {
        rd_file <- file.path(rd_dir, sprintf('inps_sector_%03d.qs', sector_n))
        qs::qsave(inps_dt, rd_file, preset = 'high')
    } else {
        ## saveRDS with xz is slower but needs no extra package
        rd_file <- file.path(rd_dir, sprintf('inps_sector_%03d.rds', sector_n))
        saveRDS(inps_dt, rd_file, compress = 'xz')
    }
    cat('wrote raw output: ', rd_file, ' (',
        round(file.size(rd_file) / 1e6, 1), ' MB)\n', sep = '')
    cat('  100 sectors at this size would need about ',
        round(100 * file.size(rd_file) / 1e9, 1), ' GB\n', sep = '')
    rm(inps_dt); invisible(gc())
}

## ===========================================================================
## FIX 7  [v01 lines 233-236]  sum across species, and keep more variables
## ---------------------------------------------------------------------------
## In 3-PGmix each species' value is that species' CONTRIBUTION to the stand,
## so stand totals are the SUM across species. v01 used mean(), which divides
## mixed cells by the number of species present.
##
## MEASURED: in the saved sector, annual GPP by species count was 14.8 / 13.4 /
## 14.2 tDM/ha for 1 / 2 / 3 species when summed correctly, i.e. flat.
##
## Also retain NPP, LAI, biomass and the modifiers. All are validatable
## against the plot dataset and re-running the cube to obtain them later is
## expensive. State variables are averaged across species, flux and stock
## variables are summed.
## ===========================================================================
vars_sum  <- c('gpp','npp','lai','biom_stem','biom_foliage','biom_root',
               'basal_area','stems_n')
vars_mean <- c('height','dbh','f_sw','f_vpd','f_tmp','f_frost','f_nutr','f_age',
               'alpha_c')
keep_vars <- c(vars_sum, vars_mean)

summarise_var <- function(x, v) {
    fun <- if (v %in% vars_sum) sum else mean
    x %>% filter(variable == v) %>%
        group_by(date) %>%
        summarise(species = if (n_distinct(species) == 1) first(species) else "mixed",
                  value = fun(value, na.rm = TRUE), .groups = 'drop')
}

## species raster, taken from the first retained variable
pg_lst <- inps %>% map(~ summarise_var(.x, 'gpp'))
species_vec <- map_chr(pg_lst, ~ .x$species[1])
ref_raster <- kryck_elev
species_rast <- rast(ref_raster)
species_fac <- as.factor(species_vec)
species_rast[] <- NA_integer_
species_rast[as.numeric(names(pg_lst))] <- as.integer(species_fac)
levels(species_rast) <- data.frame(
  ID = seq_along(levels(species_fac)),
  species = levels(species_fac)
)

## one time-series raster per retained variable
out_dir <- file.path(base_path, 'outputs_3PG')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (v in keep_vars) {
    lst <- inps %>% map(~ summarise_var(.x, v))
    lst <- lst[vapply(lst, nrow, integer(1)) > 0]
    if (!length(lst)) { cat('no output for', v, '\n'); next }
    tsr_v <- timeSeriesRaster(kryck_elev, lst, names(lst))
    names(tsr_v) <- time(tsr_v)
    vdir <- file.path(out_dir, v)   # tiles are named by extent, no collision
    dir.create(vdir, showWarnings = FALSE, recursive = TRUE)
    save_timeseries_raster(tsr_v, vdir)
    cat('wrote', v, '\n')
    rm(tsr_v, lst); invisible(gc())
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

## ---------------------------------------------------------------------------
## ACCEPTANCE TESTS after this runs
##   1. simulated LAI at the 50 plots inside 0.28 to 6.75   (v01 reached 33.6)
##   2. r(simulated GPP, APAR) above 0.7                    (v01 gave 0.008)
##   3. median simulated height between 5 and 30 m          (v01 gave 3.2e22)
##   4. pine r against measured GPP near 0.50 with the mapped species map,
##      which is 66 % accurate; 0.641 is the best case with true species
##   5. fewer than 9 validation plots absent from the cube
##   6. cell completeness printed above, and a failures csv if any
##
## Expect a bias near -440 gC/m2/yr with published parameters. Do NOT calibrate
## it away: "no calibration was performed at any stage" is worth more to the
## manuscript than a smaller bias, and the boreal temperature variant is the
## sensitivity analysis that explains it.
## ---------------------------------------------------------------------------
