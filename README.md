# Krycklan 3-PGmix: analysis and data-cube code

**DOI.** 10.5281/zenodo.23614212 (reserved; it resolves on publication of the record).

Code for Lara Henao, W., Wallerman, J., Gessler, A., Martínez-García, E., Luković, M., Forsman, M., Peichl, M., Valbuena, R., "Where high-resolution inputs stop producing high-resolution predictions: diagnosing information loss in a gridded process-based forest productivity model", in review at Environmental Modelling & Software. The corrected data cube it builds is deposited separately (doi:10.5281/zenodo.23472784, reserved).

Licence: the code is licensed under the GNU General Public License version 3 or later (GPL-3.0-or-later), with the version 3 text in `LICENSE`.

Paths below are relative to the folder holding this README and keep the layout of the authors' repository. Scripts in `devel/UVa_work/` are run from that folder (`Rscript 03_validate_cube_V2.R`); the cube pipeline is run from `devel/iteration_by_sectors/`.

For publication, comments and five printed messages were edited to remove references to the authors' internal records; the code is otherwise unchanged from the version that produced the results.

## What can be reproduced

- **The corrected cube is reproducible** from this code and the inputs listed below: `devel/iteration_by_sectors/server/polygon_grids_v03_server.r` (run per sector by `run_sectors.sh`), then `devel/UVa_work/server/20_assemble_cube.R` and `21_cube_qc.R`. These are the copies that ran on the processing server.
- **The initial cube is the 100 per-sector GeoTIFF tiles written by `devel/polygon_grids_v01.r`** on 21–22 November 2025 (09:23 to 00:14). The file here is byte-identical to the version committed on 21 November 2025, the day the tiles were written. The tiles are not included in this record; the authors hold three copies, verified against a SHA-256 list. `out_3PG_netcdf/GPP.nc` is a one-off NetCDF conversion of those tiles, made on 15 December 2025 to share the cube, not a separate product; its values equal the mosaic of the tiles cell for cell. It was written by `devel/mosaic_grids.r`, whose lines 31 to 52 build a virtual mosaic of the tiles and write it with `terra::writeCDF`; the rest of that file is plotting and, from line 164, a copy of `polygon_grids_v01.r` lines 36 to 248. Scripts 03, 04, 06, 13, 14 and 43 and `devel/UVa_diagnostics.R` read `GPP.nc`; script 30 reads the tiles.

## Data the scripts read (not included)

- Plot data: plot attributes from Martínez-García et al. (2022) (doi:10.5281/zenodo.5730154); annual ecosystem GPP from Peichl et al. (2026) (doi:10.5281/zenodo.17581090); the tree and forest-floor GPP split in the cube record (doi:10.5281/zenodo.23472784). The scripts read these from the authors' compiled workbook, `devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx`, which is not included; script 32 also reads the Peichl et al. (2023) deposit (doi:10.5281/zenodo.7312986).
- Model inputs, under one base folder (`bp`, `base_path`, `B` or `BASE` in the scripts): aboveground biomass and stand age from the National Forest Attribute Maps (Sweden), species from the NMD base layer, elevation from the Lantmäteriet national DEM, temperature and precipitation from SMHI PTHBV, solar radiation from MODIS MCD18A1 (exported with `devel/modis_dsr_export.js` in Google Earth Engine), soil texture and available water capacity from the European Soil Database, and the SLU soil moisture map (scripts 31 and 38). Vintages are given in the paper (Section 2.4) and in the cube record's README.
- `devel/UVa_work/marklund_coef.csv` (included) is the authors' transcription of the Marklund (1988) coefficients read by script 28.

## Run order

From the authors' pipeline notes, with the scripts written after them placed in their phase.

0. **Build the cube** (processing server). `cd devel/iteration_by_sectors; seq 1 100 | xargs -P 6 -I {} bash -c 'SCRIPT=server/polygon_grids_v03_server.r CELLS_PER_SIDE=10 ./run_sectors.sh {} {} TRUE'` (the runner's default script name, `polygon_grids_v02.r`, is not the one used). Then `Rscript devel/UVa_work/server/20_assemble_cube.R assemble` and `Rscript devel/UVa_work/server/21_cube_qc.R`. `read_raw_archive.R` reads the per-sector archives written with the third argument `TRUE`. The pipeline loads the helper functions in `devel/R/internals.r` with `R.utils::sourceDirectory()` over `devel/R/` and reads `devel/spp_codes_from_land.csv`.
1. **Validate the cube.** `03_validate_cube_V2.R <cube> <workbook>` (pass both arguments: the default workbook name is wrong). Its `validation_pairs.csv` feeds phase 7.
2. **Eliminate explanations.** `04_check_registration_V2.R <cube> <workbook>`, `05_check_inputs_V1.R <base_path> <workbook>`, `06_input_output_link_V1.R`, `07_within_object_link_V1.R`.
3. **Locate the faults at the plots.** `08_pointmode_57plots_V1.R`, `09_light_model_V1.R`, `10_params_heightmodel_V1.R`.
4. **Bound the proposed refinements.** `11_measure_vpd_itc_V1.R`, `12_stems_model_bootstrap_V1.R` (after 10: Part B reads `validation_out/params_results.csv`).
5. **Grid alignment and missing plots.** `13_alignment_audit_V1.R`, `14_missing_plots_V1.R`, `15_unexplained_plots_V1.R`.
6. **Plot-level runs with the corrections.** `16_preflight_corrected_V1.R` (feeds phase 7); `36_alphacx_zero_bias.R`; `57_table7_2016_2018_basis.R`.
7. **Ceilings.** Set `VARIANT_KEEP` in `ceiling_config.R` (run it once to list the variants), then `bash run_ceiling_analysis.sh` (scripts 23, 24, 25).
8. **Plot-mean statistics and tables.** `26_plotlevel_stats.R`, then `33_table6_predictors.R`, `34_variant_tables.R`, `37_elasticity_interval.R`, `38_section4_supporting.R`, `39_alphacx_slope_from_sweep.R`, `40_dispersion_vs_elasticity.R`, `41_tier2_support.R`, `43_tier5_support.R`.
9. **Figures.** Figure 1: `29_fig1_site_v2.R` (`server/29_fig1_site_v2.R` is the copy run on the server, which holds the national DEM). Figure 2 has its own source and PDF in the manuscript repository (`fig/fig2_workflow.tex`), not in this record. Figure 3: `30_fig3_initial.R`, then `42_fig3_annotations.R`. Figures 4 to 6: `35_figures_corrected_species.R` (it replaces the authors' earlier `27_make_figures.R`, not included).

**Standalone checks**, cited for individual statements in the paper: `22_check_srad.R`, `28_foliage_fraction.R`, `31_soil_moisture_test.R`, `32_aspect_test.R`, `44` to `50`, `56_cube_qc_fail_check.R`, `devel/UVa_diagnostics.R`, `devel/UVa_work/patches/which_column_is_species.R` (the species-column check). On the processing server, with a sector number as first argument: `51_sector_grid_crops.R`, `52_site_misindex_check.R`, `53_site_misindex_gpp_sensitivity.R`, `54_residue_by_input_list.R`, and the stem-density test `58a_stem_density_plot_cells.R`, `58b_stem_density_worker.R <sector> 10 FALSE <dir>`, `58c_stem_density_combine.R`. `55_tall_cell_rerun.R` reruns the tall cells of one saved sector. `devel/get_segmented_trees.r` is the individual-tree crown segmentation behind the airborne-laser stem density (script 11).

## Software

| Group | Where it ran | Versions |
|---|---|---|
| Corrected-cube pipeline, assembly and QC; scripts 29 (server copy), 30, 51 to 54 and 58 | Processing server | R 4.5.2, terra 1.8-86, sf 1.0.23, r3PG 0.1.6 (terra logged; R and sf are the server's installed versions, not logged) |
| Initial cube: per-sector tiles (`polygon_grids_v01.r`) | Remote workstation under Windows Subsystem for Linux, accessed over SSH | terra, R and r3PG not recorded |
| Initial cube: NetCDF conversion of the tiles (`mosaic_grids.r`) | Laptop | terra 1.8-80 (from the metadata of `GPP.nc`) |
| Cube scoring (03) | Server and local machine | as above, or the local versions |
| Plot-level runs, statistics, tables, figures 1 (local copy), 4 to 6, and the other checks | Local machine | R 4.6.1, terra 1.9-27, sf 1.1.1, r3PG 0.1.6 |
| `modis_dsr_export.js` | Google Earth Engine | not recorded |
| `get_segmented_trees.r` (crown segmentation) | Local machine | lasR, with lasRtools, lidR, terra, sf and R.utils; versions not recorded |

R packages the code loads, none of them included:

- On CRAN: data.table, dplyr, geodata, ggplot2, lidR, lubridate, magrittr, maps, ncdf4, patchwork, purrr, qs2, r3PG, readxl, rlang, R.utils, scales, sf, stringr, terra, tidyr, tidyverse, units and zoo, with parallel, stats and utils from base R.
- lasR, not on CRAN but publicly available through r-universe (r-lidar.r-universe.dev): the crown segmentation in `get_segmented_trees.r` is built from its functions.
- lasRtools, the authors' helper package, not publicly available: `get_segmented_trees.r` uses it only to locate and group the input tiles, which sets each sector's tiles and buffer, and to plot.
- wilaRtools, the authors' package, not publicly available: loaded only by the initial-cube script, `polygon_grids_v01.r`, and by `mosaic_grids.r`, which repeats that script's code.

The corrected cube, its assembly and quality control, and every result behind Figures 4 to 6 and Table 7 of the paper need only base R and packages on CRAN.

The pipeline writes its per-sector archives with qs2. The package qs is used only as a fallback for legacy archives and was removed from CRAN in January 2026, so a reader needs qs2.

## Hard-coded paths a reader must change

Every absolute path in the code, by file and line. Most point to the authors' input folder (`/run/media/wilar/T7 Shield/slu_data_bases/3PG`, or `/home/wilar/datos/3PG` on the processing server) or to the plot workbook; set them to the local copies. Scripts 52, 53, 54 and 58b point `PIPE` at `devel/iteration_by_sectors/polygon_grids_v03_server.r`, the location of the pipeline on the processing server; here it is in `server/`.

| File | Line | Path |
|---|---|---|
| `devel/UVa_diagnostics.R` | 6 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/out_3PG_netcdf/GPP.nc` |
| `devel/UVa_diagnostics.R` | 14 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/out_3PG_netcdf/` |
| `devel/UVa_work/03_validate_cube_V2.R` | 26 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/out_3PG_netcdf/GPP.nc` |
| `devel/UVa_work/04_check_registration_V2.R` | 25 | `~/3pg_data/out_3PG_netcdf/GPP.nc` |
| `devel/UVa_work/04_check_registration_V2.R` | 27 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/05_check_inputs_V1.R` | 28 | `~/3pg_data` |
| `devel/UVa_work/05_check_inputs_V1.R` | 30 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/06_input_output_link_V1.R` | 24 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/06_input_output_link_V1.R` | 26 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/07_within_object_link_V1.R` | 22 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/3PG_10m/x731866.017196386_y7132899.71658545.RData` |
| `devel/UVa_work/08_pointmode_57plots_V1.R` | 41 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/08_pointmode_57plots_V1.R` | 43 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/09_light_model_V1.R` | 27 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/09_light_model_V1.R` | 29 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/10_params_heightmodel_V1.R` | 40 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/10_params_heightmodel_V1.R` | 42 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/11_measure_vpd_itc_V1.R` | 35 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/11_measure_vpd_itc_V1.R` | 36 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/12_stems_model_bootstrap_V1.R` | 27 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/13_alignment_audit_V1.R` | 36 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/13_alignment_audit_V1.R` | 38 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/14_missing_plots_V1.R` | 32 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/14_missing_plots_V1.R` | 34 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/15_unexplained_plots_V1.R` | 27 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/15_unexplained_plots_V1.R` | 28 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/16_preflight_corrected_V1.R` | 45 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/16_preflight_corrected_V1.R` | 47 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/22_check_srad.R` | 19 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/28_foliage_fraction.R` | 31 | `~/Documents/SLU/3PG/devel/UVa_work/validation_out/preflight_results.csv` |
| `devel/UVa_work/28_foliage_fraction.R` | 89 | `~/Documents/SLU/3PG/devel/UVa_work/validation_out/itc_vs_measured.csv` |
| `devel/UVa_work/29_fig1_site_v2.R` | 28 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/30_fig3_initial.R` | 24 | `~/Documents/SLU/3PG/devel/UVa_work/validation_out/validation_pairs.csv` |
| `devel/UVa_work/31_soil_moisture_test.R` | 28 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/SLU_soil_map/SLUMarkfuktighetskarta.tif` |
| `devel/UVa_work/32_aspect_test.R` | 33 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/32_aspect_test.R` | 34 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Peichl2023_zenodo7312986_NEP_deposit.xlsx` |
| `devel/UVa_work/33_table6_predictors.R` | 27 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/36_alphacx_zero_bias.R` | 45 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/36_alphacx_zero_bias.R` | 47 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/38_section4_supporting.R` | 88 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/SLU_soil_map/SLUMarkfuktighetskarta.tif` |
| `devel/UVa_work/38_section4_supporting.R` | 93 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/41_tier2_support.R` | 99 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/43_tier5_support.R` | 28 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/44_input_variation_cube_plots.R` | 35 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/45_soil_moisture_at_plots.R` | 18 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/46_radiation_station_check.R` | 24 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/climate_SW` |
| `devel/UVa_work/47_vpd_station_check.R` | 22 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/climate_SW` |
| `devel/UVa_work/48_age_raster_catchment.R` | 18 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/49_saved_sector_full.R` | 22 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/3PG_10m/x731866.017196386_y7132899.71658545.RData` |
| `devel/UVa_work/50_forest_coverage_crosscheck.R` | 28 | `/home/wilar/datos/3PG` |
| `devel/UVa_work/51_sector_grid_crops.R` | 14 | `/home/wilar/datos/3PG` |
| `devel/UVa_work/52_site_misindex_check.R` | 16 | `/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r` |
| `devel/UVa_work/53_site_misindex_gpp_sensitivity.R` | 15 | `/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r` |
| `devel/UVa_work/54_residue_by_input_list.R` | 16 | `/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r` |
| `devel/UVa_work/54_residue_by_input_list.R` | 21 | `/home/wilar/datos/3PG` |
| `devel/UVa_work/55_tall_cell_rerun.R` | 24 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG/3PG_10m/x731866.017196386_y7132899.71658545.RData` |
| `devel/UVa_work/56_cube_qc_fail_check.R` | 16 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/56_cube_qc_fail_check.R` | 18 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/57_table7_2016_2018_basis.R` | 17 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/57_table7_2016_2018_basis.R` | 19 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/58a_stem_density_plot_cells.R` | 21 | `/home/wilar/datos/3PG/outputs_3PG/cube/gpp.tif` |
| `devel/UVa_work/58a_stem_density_plot_cells.R` | 22 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/58b_stem_density_worker.R` | 22 | `/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors/polygon_grids_v03_server.r` |
| `devel/UVa_work/58c_stem_density_combine.R` | 30 | `/home/wilar/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/ceiling_config.R` | 22 | `~/Documents/SLU/3PG/devel/UVa_work/validation_out` |
| `devel/UVa_work/patches/which_column_is_species.R` | 22 | `~/Documents/SLU/3PG/devel/validation_data/Dataset_2016-2018_50+7_plots_compiled.xlsx` |
| `devel/UVa_work/patches/which_column_is_species.R` | 23 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/UVa_work/patches/which_column_is_species.R` | 25 | `~/Documents/SLU/3PG/devel/validation_data/nmd_krycklan_window.tif` |
| `devel/UVa_work/server/20_assemble_cube.R` | 25 | `/home/wilar/datos/3PG` |
| `devel/UVa_work/server/20_assemble_cube.R` | 26 | `/home/wilar/Documents/SLU/3PG/devel/iteration_by_sectors` |
| `devel/UVa_work/server/21_cube_qc.R` | 18 | `/home/wilar/datos/3PG` |
| `devel/UVa_work/server/29_fig1_site_v2.R` | 31 | `~/Documents/SLU/3PG/devel/UVa_work/validation_out/itc_vs_measured.csv` |
| `devel/UVa_work/server/29_fig1_site_v2.R` | 32 | `~/Documents/SLU/3PG/devel/UVa_work/validation_out/preflight_results.csv` |
| `devel/get_segmented_trees.r` | 36 | `/media/wihe/LaCie/projects/wihe_lasR/2_pc_source` |
| `devel/get_segmented_trees.r` | 41 | `/media/wihe/LaCie/0-projects/wihe_lasR/ofile/ITC` |
| `devel/get_segmented_trees.r` | 65 | `/media/wihe/LaCie/0-projects/wihe_lasR/ofile/ITC` |
| `devel/get_segmented_trees.r` | 66 | `/media/wihe/LaCie/0-projects/wihe_lasR/ofile/ITC` |
| `devel/get_segmented_trees.r` | 104 | `/media/wihe/LaCie/2_pc_source` |
| `devel/get_segmented_trees.r` | 105 | `/media/wihe/LaCie/0-projects/dem.tif` |
| `devel/get_segmented_trees.r` | 126 | `/media/wihe/LaCie/2_pc_source` |
| `devel/get_segmented_trees.r` | 127 | `/media/wihe/LaCie/0-projects/dem.tif` |
| `devel/get_segmented_trees.r` | 128 | `/media/wihe/LaCie/2_pc_source/ofile/itc_*.tif` |
| `devel/get_segmented_trees.r` | 152 | `/media/wihe/LaCie/Rutor_231215/Rutor_231215.shp` |
| `devel/get_segmented_trees.r` | 184 | `/media/wihe/LaCie/Rutor_231215/Rutor_231215.shp` |
| `devel/get_segmented_trees.r` | 203 | `/media/wihe/LaCie/2_pc_source` |
| `devel/get_segmented_trees.r` | 253 | `/media/wihe/LaCie/2_pc_source/ofile/itc_*.tif` |
| `devel/get_segmented_trees.r` | 254 | `/media/wihe/LaCie/2_pc_source/ofile` |
| `devel/get_segmented_trees.r` | 260 | `/media/wihe/LaCie/2_pc_source/` |
| `devel/get_segmented_trees.r` | 263 | `/media/wihe/LaCie/2_pc_source` |
| `devel/get_segmented_trees.r` | 264 | `/media/wihe/LaCie/2_pc_source/ofile` |
| `devel/get_segmented_trees.r` | 322 | `/media/wihe/LaCie/Rutor_231215/Rutor_231215.shp` |
| `devel/get_segmented_trees.r` | 417 | `/media/wihe/LaCie/2_pc_source/ofile/segmentation.vrt` |
| `devel/iteration_by_sectors/read_raw_archive.R` | 38 | `/run/media/wilar/T7 Shield/slu_data_bases/3PG` |
| `devel/iteration_by_sectors/server/polygon_grids_v03_server.r` | 30 | `/home/wilar/Documents/SLU/3PG/devel/R` |
| `devel/iteration_by_sectors/server/polygon_grids_v03_server.r` | 32 | `/home/wilar/Documents/SLU/3PG/devel` |
| `devel/iteration_by_sectors/server/polygon_grids_v03_server.r` | 70 | `/home/wilar/datos/3PG` |
| `devel/mosaic_grids.r` | 26 | `/media/wihe/LaCie/projects/3PG` |
| `devel/mosaic_grids.r` | 31 | `/media/wihe/LaCie/projects/3PG/outputs_3PG/` |
| `devel/mosaic_grids.r` | 46 | `/media/wihe/LaCie/projects/3PG/out_3PG_netcdf/GPP.nc` |
| `devel/mosaic_grids.r` | 136 | `/home/wihe/Documents/slu_locals/3PG/devel/mosaic_Kryckland_3PG_GPP_output.pdf` |
| `devel/polygon_grids_v01.r` | 26 | `/media/wihe/LaCie/projects/3PG` |

## Not included

Correspondence, meeting notes and the authors' working records; the plot workbook and other third-party data; published papers; the manuscript; literature-search tooling; outputs (`results/`, `validation_out/`, figures); presentation material; superseded or unused scripts (earlier pipeline versions, `27_make_figures.R`, the laptop copies of the server scripts); and run by-products.
