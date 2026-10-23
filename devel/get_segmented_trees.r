## ## install`reticulate`, only one time
## if (!("lidR" %in% installed.packages()[,1])) {
##  print("Installing package `lidR`...")
##  install.packages("lidR")
## } else {
##     print("Package `lidR` already installed")
## require('lidR')}
## source('~/helpers.r')
## require('rgl') # to render pointclouds
## require('sf')
## require('jsonlite') # reads json data
## require('geojson') 
## require('ecochange')
## require('terra')
## require('lidR')
## require('lasR')
## install.packages('lidR')
require('lasR')
require('lasRtools')
require('parallel')
require('R.utils')
require('lidR')
require('terra')
require('sf')
## require('terra')
locate_directory('~','lasRtools_dev/lasRtools/R', setwd = TRUE)
sourceDirectory('.', modifiedOnly = TRUE, verbose = TRUE) #
locate_directory('~','lasRtools_dev/devel')
dir()
load('sysdata.rda')



## names(from_tif) <- NULL
length_pos <- 7
source <- "/media/wihe/LaCie/projects/wihe_lasR/2_pc_source"
sector <- names(record_list[length_pos])
sector <- '184'
sector <- c('184','185')
sector <- '224'
fit_source <- paste0('/media/wihe/LaCie/0-projects/wihe_lasR/ofile/ITC',sector,'SW')
from_tif <- replace_tif_to_laz(fit_source) 
move_files_to_matching_folders(dir(source)[grepl('la[zx]$', dir(source))], source)
move_files_to_parent(from_tif, source)
ctg <- readLAScatalog(source)
plot(st_geometry(buffered_lattice_400_1000), col = NA, border = 'grey', main = paste0('sector: ',sector))
plot(ctg, add = TRUE, col = NA, border = 'blue')


non_buffered_lattice <- create_buffered_sf_lattice(rutor_table, 20, 0)
## plot(st_geometry(non_buffered_lattice), col = NA, border = 'grey', main = paste0('sector: ',sector))
## plot(ctg, add = TRUE, col = NA, border = 'blue')
non_buffered_sector <- subset(non_buffered_lattice, id %in% sector)
buffered_sector <- subset(buffered_lattice_400_1000, id %in% sector)
st_crs(ctg) <- st_crs(non_buffered_sector) 
int <- st_intersection(non_buffered_sector, st_as_sf(ctg))
int_buf <- st_intersection(buffered_sector, st_as_sf(ctg))
plot(st_geometry(buffered_sector), col = NA, border = 'black', main = paste0('sector: ',sector))
plot(st_geometry(int_buf), add = TRUE, col = NA, border = 'red')
plot(st_geometry(int), add = TRUE, col = NA, border = 'black')
plot(st_geometry(non_buffered_sector), add = TRUE, col = NA, border = 'blue')



sec88 <- replace_tif_to_laz(paste0('/media/wihe/LaCie/0-projects/wihe_lasR/ofile/ITC',88,'SW'))
sec89 <- replace_tif_to_laz(paste0('/media/wihe/LaCie/0-projects/wihe_lasR/ofile/ITC',89,'SW'))
secs <- c(sec88, sec89)
secs <- sub('.laz','', secs)

laz_in_processed_sectors <- unlist(all_files[c('88','89')], use.names = FALSE)
new_base <- sub("^[^_]*_", "", laz_in_processed_sectors)
new_base <- sub("\\.la[zx]$", "", new_base)


laz_non_buffered <- subset(rutor_table, Las_Namn%in%new_base)





indx <- sample(1:465, 1)
tmp <- rast(file.path(fit_source,dir(fit_source)[indx]))


plot(tmp)
ans_roi  <- terra::sel(tmp)
ans_spat <- Map(function(x)
    crop_roi(x, ans_roi), tmp)
ans_spat_col  <- assign_random_colors(ans_spat[[1]])
plot(ans_spat_col)
## plot(ctg, add = TRUE, mapview = FALSE)



length(record_list)


## buffered_lattice_400_1000 <- create_buffered_sf_lattice(rutor_table, 20, buffer = 1000)




length_pos <- 5
source <- "/media/wihe/LaCie/2_pc_source"
dem <- '/media/wihe/LaCie/0-projects/dem.tif'
ofile <- file.path(source,'ofile')
sector <- names(record_list[length_pos])
ofdir <- file.path(ofile, paste0('sector_', sector))
dir.create(ofdir, showWarnings = FALSE)
of <- paste0(ofdir, '/itc_*.tif')

move_files_to_matching_folders(dir(source)[grepl('la[zx]$', dir(source))], source)
move_files_to_parent(all_files[[length_pos]], source)
ctg <- readLAScatalog(source)
plot(st_geometry(buffered_lattice_400_1000), col = NA, border = 'grey', main = paste0('sector: ',sector))
plot(ctg, add = TRUE, col = NA, border = 'blue')

## Indexing
pipeline = write_lax(embedded = FALSE)
start.time <- Sys.time()
ans = exec(pipeline, on = source, ncores = round(parallel::detectCores()/3))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

source <- "/media/wihe/LaCie/2_pc_source"
dem <- '/media/wihe/LaCie/0-projects/dem.tif'
of <- "/media/wihe/LaCie/2_pc_source/ofile/itc_*.tif"
del <- load_raster(dem)
norm <- transform_with(del, "-")
pix <- 1
chm <- lasR::rasterize(pix, 'max')
chm2 <- pit_fill(chm)
## chm2 <- chm
seed <- local_maximum_raster(chm2, ws = pix * 3, min_height = 2)
## seed = local_maximum(ws = pix * 3, min_height = 2)
tree <- region_growing(chm2, seed, ofile = of)
pipeline <- del + norm + chm + chm2 + seed + tree
start.time <- Sys.time()
print(start.time)
ans <- exec(pipeline, on = source, ncores = round(parallel::detectCores()/5))
## ans <- exec(pipeline, on = source, ncores = 5)
## ans <- exec(pipeline, on = ctg)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
record_list[[sec]] <- TRUE
save(all_files, buffered_lattice_400_1000, rutor_table, record_list, file = 'sysdata.rda')


## rutor_table <- sf::st_read('/home/wihe/Documents/lasRtools_dev/devel/Rutor_231215/Rutor_231215.shp')
rutor_table <- sf::st_read('/media/wihe/LaCie/Rutor_231215/Rutor_231215.shp')
## crd <- st_coordinates(rutor_table)
## tail(crd)
lat0 <- create_sf_lattice(rutor_table, 5)
lat1 <- create_buffered_sf_lattice(rutor_table, 5, buffer = 2000)
## st_area(lat0[1L,2L])
## st_area(lat1[1L,2L])
## plot(st_geometry(lat1), col = NA, border = 'red', main = "Buffered Lattice")
## plot(st_geometry(rutor_table), add = TRUE)
## plot(st_geometry(lat0), add = TRUE)

## f_sector  <- function(pol_nr){
## pol. <- lat1[lat1$'id' == pol_nr,]
## section <- st_crop(rutor_table, pol.)
## ## plot(section)
## ## dirs <- paste0(section$'Las_Namn',c('.laz', '.lax'))
## dirs_laz <- paste0(section$'Las_Namn',c('.laz'))
## dirs_lax <- paste0(section$'Las_Namn',c('.lax'))
## return(c(dirs_laz,dirs_lax))}
## dirs_ <- Map(function(x)
##             f_sector(x),1:25)
## names(dirs_) <- 1:25
## filtered_laz_dirs <- Filter(function(x)length(x) > 2, dirs_) 

## plot(st_geometry(lat1), col = NA, border = 'red', main = "Buffered Lattice")
## pol_nr <- 1
## pol. <- lat1[lat1$'id' == pol_nr,]
## section. <- st_crop(rutor_table, pol.)
## plot(st_geometry(section.), add = TRUE)

## save(filtered_laz_dirs, file = 'filtered_laz_dirs.RData')

rutor_table <- sf::st_read('/media/wihe/LaCie/Rutor_231215/Rutor_231215.shp')

lat1 <- create_buffered_sf_lattice(rutor_table, 100, buffer = 100)
plot(st_geometry(lat1), col = NA, border = 'grey')


get_laz_in_buffered_sf <- function(rutor_table, sector = 1, n = 10, buffer = 100){
lat1 <- create_buffered_sf_lattice(rutor_table, n = n, buffer = buffer)
    ## pol_nr <- sector
pol. <- subset(lat1, id == sector)
## pol. <- lat1[lat1$'id' == pol_nr,]
section <- st_crop(rutor_table, pol.)
## plot(section)
## dirs <- paste0(section$'Las_Namn',c('.laz', '.lax'))
dirs_laz <- paste0(section$'Las_Namn',c('.laz'))
dirs_lax <- paste0(section$'Las_Namn',c('.lax'))
dirs_both <- c(dirs_laz, dirs_lax)
return(dirs_both)}

source <- "/media/wihe/LaCie/2_pc_source"
buffered_lattice_400_1000 <- create_buffered_sf_lattice(rutor_table, 20, buffer = 1000)
lazfiles <- get_laz_in_buffered_sf(rutor_table, sector = 11, n = 10, buffer = 1000)

all_laz_files <- Map(function(x)
                     get_laz_in_buffered_sf(rutor_table, sector = x, n = 20, buffer = 1000), 1:400)
names(all_laz_files)  <- 1:400
all_files <- Filter(function(x)
                    length(x)>2, all_laz_files)
all_files <- all_files[order(sapply(all_files,length), decreasing = TRUE)]
## str(all_files)

getwd()
save(all_files, buffered_lattice_400_1000, rutor_table, record_list, file = 'sysdata.rda')



getwd()


sec <- sample(1:length(all_files), 1)
names(all_files[sec])
move_files_to_matching_folders(dir(source)[grepl('la[zx]$', dir(source))], source)
move_files_to_parent(all_files[[sec]], source)
ctg <- readLAScatalog(source)
plot(st_geometry(buffered_lattice_400_1000), col = NA, border = 'grey')
plot(ctg, add = TRUE)


record_list <- as.list(rep(FALSE, length(all_files)))
names(record_list)  <- names(all_files)


## move_files_to_matching_folders(lazfiles, source)
## move_files_to_parent(lazfiles, source)

## writeLines(c(dirs_laz,dirs_lax), '/media/wihe/LaCie/2-pc-source/file.txt')

find_zero_byte_files <- function(path) {
  if (!dir.exists(path)) stop("Path does not exist or is not a directory.")
  # List all files recursively
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  # Get file sizes
  file_info <- file.info(files)
  # Filter files with 0 bytes
  zero_byte_files <- rownames(file_info[file_info$size == 0, , drop = FALSE])
  zero_byte_files  <- basename(zero_byte_files)
  return(zero_byte_files)
}

of <- "/media/wihe/LaCie/2_pc_source/ofile/itc_*.tif"
of. <- "/media/wihe/LaCie/2_pc_source/ofile"
zero_byte_files <- find_zero_byte_files(of.)

not_sub <- sub('itc_', '', zero_byte_files)
add_laz_ext <- sub('tif', 'laz', not_sub)

dir("/media/wihe/LaCie/2_pc_source/")

#++++++++++++++++++++++++++++++++++++++++++++++++++++
source <- "/media/wihe/LaCie/2_pc_source"
ofiles <- "/media/wihe/LaCie/2_pc_source/ofile"
find_laz_lax_files <- function(path, recursive = FALSE) {
  if (!dir.exists(path)) stop("Path does not exist or is not a directory.")
  files <- list.files(path, recursive = recursive, full.names = TRUE)
  laz_lax_files <- files[grepl("\\.(laz|lax)$", files, ignore.case = TRUE)]
  laz_lax_files <- basename(laz_lax_files)
  return(laz_lax_files)
}
source_laz_files <- find_laz_lax_files(source)
generate_laz_lax_names <- function(path, recursive = FALSE) {
  if (!dir.exists(path)) stop("Path does not exist or is not a directory.")
  # List .tif files starting with "itc_"
  tif_files <- list.files(path, recursive = recursive, full.names = TRUE, pattern = "^itc_.*\\.tif$", ignore.case = TRUE)
 # Generate .laz and .lax counterparts
  tif_files <- sub('itc_', '', tif_files)
  laz_files <- sub("\\.tif$", ".laz", tif_files, ignore.case = TRUE)
  lax_files <- sub("\\.tif$", ".lax", tif_files, ignore.case = TRUE)
  laz_files.  <- c(laz_files, lax_files)
  laz_files. <- basename(laz_files.) 
  return(laz_files.)
}
ofiles_with_laz_lax_names  <- generate_laz_lax_names(ofiles)
## writeLines(ofiles_with_laz_lax_names, file.path(source, 'ofiles_with_laz_lax_names.txt'))
##
## files_to_process <- setdiff(source_laz_files, ofiles_with_laz_lax_names)
## writeLines(files_to_process, file.path(source, 'files_to_process.txt'))

move_files_to_matching_folders <- function(file_names, parent_dir) {
  if (!dir.exists(parent_dir)) stop("Parent directory does not exist.")
  for (file_name in file_names) {
    # Extract the prefix before the first underscore
    prefix <- sub("(_.*)$", "", file_name)
    target_dir <- file.path(parent_dir, prefix)
    src_path <- file.path(parent_dir, file_name)
    dest_path <- file.path(target_dir, file_name)
    # Check that source file exists
    if (!file.exists(src_path)) {
      warning(paste("Source file not found:", src_path))
      next
    }
    # Check that target directory exists
    if (!dir.exists(target_dir)) {
      warning(paste("Target directory does not exist:", target_dir))
      next
    }
    # Move the file
    success <- file.rename(src_path, dest_path)
    if (!success) {
      warning(paste("Failed to move file:", file_name))
    }
  }
}
## move_files_to_matching_folders(ofiles_with_laz_lax_names, source)

## setwd(source)
## smp <- dir()[grepl('.la[zx]$', dir())][sample(1:2000, 5)]
## move_files_to_matching_folders(smp, source)

rutor_table <- sf::st_read('/media/wihe/LaCie/Rutor_231215/Rutor_231215.shp')
ctg <- readLAScatalog(source)
## plot(ctg)

section.

laz_in_source  <- dir(source)[grepl('la[zx]$', dir(source))]

length(laz_in_source)

move_files_to_matching_folders(laz_in_source, source)


move_files_to_parent <- function(file_names, parent_dir) {
## parent_dir <- source
## file_names  <- dirs_both
  if (!dir.exists(parent_dir)) stop("Parent directory does not exist.")
  # Recursively search for the files in subdirectories
  all_files <- list.files(path = parent_dir, recursive = TRUE, full.names = TRUE)
  # Match files by their basename
  matched_files <- all_files[basename(all_files) %in% file_names]
 if (length(matched_files) == 0) {
    message("No matching files found.")
    return(invisible(NULL))
  }
  # Move each file to the parent directory
  for (src in matched_files) {
    dest <- file.path(parent_dir, basename(src))
    success <- file.rename(src, dest)
    if (!success) {
      warning(paste("Failed to move file:", src))
    }
  }
}

move_files_to_parent(dirs_both, source)



warnings()

source_name  <- gsub('.la[zx]','',laz_in_source)
rutor_table_sub <- subset(section., Las_Namn%in%source_name) 

# Function to buffer polygons
buffer_polygons <- function(sf_object, buffer_dist = 10) {
  st_buffer(sf_object, dist = buffer_dist)
}
rutor_table_sub_buffered <- buffer_polygons(rutor_table_sub, buffer_dist = 1)
## plot(st_geometry(rutor_table_sub_buffered))
intersections <- st_intersection(section., rutor_table_sub_buffered)
length(unique(intersections$'Las_Namn'))
textsect <- subset(section., Las_Namn%in%unique(intersections$'Las_Namn'))
plot(st_geometry(rutor_table_sub_buffered), col = NA, border = 'grey')
plot(st_geometry(rutor_table_sub), col = NA, border = 'red', add = TRUE)
plot(st_geometry(textsect), col = NA, border = 'black', add = TRUE)


plot(st_geometry(rutor_table_sub), col = NA, border = 'blue', add = TRUE)

joined_sf  <- st_join(rutor_table_sub_buffered, section.)


plot(st_geometry(joined_sf))

st_crs(ctg) <- st_crs(rutor_table) 


st_crs(rutor_table)
st_crs(ctg)

ctg_intersection <- st_intersection(rutor_table, st_as_sf(ctg), sparce = FALSE)


plot(st_geometry(ctg_intersection))

ctg_non_intersecting <- st_filter(rutor_table, st_as_sf(ctg), .predicate = st_disjoint)

section.. <- st_crop(ctg_non_intersecting, pol.)

plot(st_as_sf(section..[,4L]))


plot(st_geometry(st_as_sf(ctg)))
plot(st_geometry(st_as_sf(ctg_intersection)), col = NA, border = 'red', add = TRUE)


plot(st_geometry(lat1), col = NA, border = 'red', main = "Buffered Lattice")




plot(ctg)

## NEXT CODE CREATES VIRTUAL RASTERS
vrt_output <- "/media/wihe/LaCie/2_pc_source/ofile/segmentation.vrt"
# Check if gdalbuildvrt is available
## if (Sys.which("gdalbuildvrt") == "") {
##   stop("gdalbuildvrt not found. Please ensure GDAL is installed and accessible in your PATH.")
## }


# Construct the command
vrt_command <- sprintf("gdalbuildvrt %s %s", shQuote(vrt_output), of)

# Run the system command to create the VRT
system(vrt_command)


tmp <- terra::rast(vrt_output)


plot(tmp)

