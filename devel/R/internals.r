one_agg_cell_that_matches_polygon <- function(r,
                                             pol,
                                             target_res = 10,
                                             fun = mean,
                                             snap_if_misaligned = TRUE,
                                             tol = 1e-6) {
  stopifnot(requireNamespace("terra", quietly = TRUE))
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  if (!inherits(pol, "SpatVector")) pol <- terra::vect(pol)

  # CRS
  if (!terra::same.crs(r, pol)) pol <- terra::project(pol, r)

  # Factor check (must be integer multiple)
  res_in <- terra::res(r)
  if (abs(res_in[1] - res_in[2]) > tol) stop("Raster must have square cells for this helper.")
  fact <- target_res / res_in[1]
  if (abs(fact - round(fact)) > tol) stop("target_res must be an integer multiple of raster resolution.")
  fact <- as.integer(round(fact))

  # Use raster grid origin as reference (this preserves global alignment)
  e <- terra::ext(r)
  xmin0 <- e$xmin
  ymax0 <- e$ymax

  # Choose the cell by polygon centroid (robust when polygon slightly off-grid)
  cxy <- terra::crds(terra::centroids(pol))[1, ]
  x <- cxy[1]; y <- cxy[2]

  # 10 m grid indices (1-based)
  col10 <- floor((x - xmin0) / target_res + 1e-12) + 1L
  row10 <- floor((ymax0 - y) / target_res + 1e-12) + 1L

  # Exact 10 m cell extent
  x0 <- xmin0 + (col10 - 1L) * target_res
  x1 <- x0 + target_res
  y1 <- ymax0 - (row10 - 1L) * target_res
  y0 <- y1 - target_res
  e10 <- terra::ext(x0, x1, y0, y1)

  # Check whether polygon already matches this cell exactly (by extent)
  ep <- terra::ext(pol)
  aligned_exact <- all(abs(c(ep$xmin - x0, ep$xmax - x1, ep$ymin - y0, ep$ymax - y1)) < tol)
  if (!aligned_exact && !isTRUE(snap_if_misaligned)) {
    off <- c(xmin = ep$xmin - x0, xmax = ep$xmax - x1, ymin = ep$ymin - y0, ymax = ep$ymax - y1)
    stop("Polygon does not align to the aggregated grid. Offsets (m): ",
         paste(names(off), signif(off, 6), collapse = ", "))
  }

  # Read only the tiny window from the big raster (for 2m DEM and 10m target => 5x5 cells)
  r_win <- terra::crop(r, e10, snap = "near")

  # Aggregate the small window to 1x1 at 10 m
  r_10m_1cell <- terra::aggregate(r_win, fact = fact, fun = fun, na.rm = TRUE)

  # Build the exact matching 10 m polygon (fast)
  cell_poly <- terra::vect(list(matrix(
    c(x0, y0,
      x0, y1,
      x1, y1,
      x1, y0,
      x0, y0),
    ncol = 2, byrow = TRUE
  )), type = "polygons", crs = terra::crs(r))

  list(
    raster_cell_10m = r_10m_1cell,          # 1x1 SpatRaster at 10 m
    cell_polygon_10m = cell_poly,           # polygon that matches that cell exactly
    aligned_exact = aligned_exact,          # TRUE if your input polygon already matched
    window_2m = r_win,                      # the 5x5 (or fact x fact) window read from disk
    cell_extent = e10,
    rowcol_10m = c(row = row10, col = col10)
  )
}


make_10m_centered_squares <- function(df,
                                      x_col = "E_Sweref99TM",
                                      y_col = "N_Sweref99TM",
                                      crs = 3006,              # SWEREF99 TM (EPSG:3006). Change if needed.
                                      size = 10,               # side length in meters
                                      output = c("sf", "terra")) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(x_col, y_col) %in% names(df)))

  output <- match.arg(output)

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Please install 'sf' (install.packages('sf')).")
  }
  if (output == "terra" && !requireNamespace("terra", quietly = TRUE)) {
    stop("Please install 'terra' (install.packages('terra')).")
  }

  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  half <- size / 2

  # Build one polygon per row (square centered at (x, y))
  polys_sfc <- sf::st_sfc(lapply(seq_along(x), function(i) {
    if (is.na(x[i]) || is.na(y[i])) return(sf::st_polygon())
    coords <- matrix(
      c(x[i] - half, y[i] - half,
        x[i] - half, y[i] + half,
        x[i] + half, y[i] + half,
        x[i] + half, y[i] - half,
        x[i] - half, y[i] - half),
      ncol = 2, byrow = TRUE
    )
    sf::st_polygon(list(coords))
  }), crs = crs)

  out <- sf::st_sf(df, geometry = polys_sfc)

  if (output == "sf") return(out)

  # Convert to terra SpatVector if requested
  terra::vect(out)
}


snap_points_to_raster_cells <- function(df,
                                        r,
                                        x_col = "E_Sweref99TM",
                                        y_col = "N_Sweref99TM",
                                        keep_outside = TRUE,
                                        output = c("sf", "terra")) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c(x_col, y_col) %in% names(df)))

  output <- match.arg(output)

  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Please install 'terra' (install.packages('terra')).")
  }
  if (output == "sf" && !requireNamespace("sf", quietly = TRUE)) {
    stop("Please install 'sf' (install.packages('sf')).")
  }

  # Load raster if a path was provided
  if (is.character(r)) r <- terra::rast(r)
  stopifnot(inherits(r, "SpatRaster"))

  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])

  # Track points outside raster extent
  e <- terra::ext(r)
  inside <- !is.na(x) & !is.na(y) &
    (x >= e$xmin) & (x <= e$xmax) &
    (y >= e$ymin) & (y <= e$ymax)

  # For outside points, either snap to nearest edge cell (keep_outside=TRUE)
  # or return NA geometry (keep_outside=FALSE)
  x_use <- x
  y_use <- y

  if (keep_outside) {
    # clamp to raster extent (use a tiny epsilon so xmax/ymax don't fall on boundary edge-cases)
    epsx <- terra::res(r)[1] * 1e-9
    epsy <- terra::res(r)[2] * 1e-9
    x_use <- pmin(pmax(x_use, e$xmin + epsx), e$xmax - epsx)
    y_use <- pmin(pmax(y_use, e$ymin + epsy), e$ymax - epsy)
  }

  # Cell index for each point (containing cell; for outside points: nearest edge cell if clamped)
  cells <- rep(NA_integer_, length(x_use))
  ok <- !is.na(x_use) & !is.na(y_use) & (keep_outside | inside)
  if (any(ok)) {
    cells[ok] <- terra::cellFromXY(r, cbind(x_use[ok], y_use[ok]))
  }

  # Row/col per cell
  rc <- matrix(NA_integer_, nrow = length(cells), ncol = 2)
  has_cell <- !is.na(cells)
  if (any(has_cell)) {
    rc[has_cell, ] <- terra::rowColFromCell(r, cells[has_cell])
  }

  # Build one polygon per input point (duplicates allowed, order preserved)
  xr <- terra::res(r)[1]
  yr <- terra::res(r)[2]

  xmin_r <- e$xmin
  ymax_r <- e$ymax

  poly_mats <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    if (is.na(cells[i])) {
      poly_mats[[i]] <- NULL
      next
    }
    row <- rc[i, 1]
    col <- rc[i, 2]

    xleft  <- xmin_r + (col - 1) * xr
    xright <- xleft + xr
    ytop   <- ymax_r - (row - 1) * yr
    ybot   <- ytop - yr

    poly_mats[[i]] <- matrix(
      c(xleft,  ybot,
        xleft,  ytop,
        xright, ytop,
        xright, ybot,
        xleft,  ybot),
      ncol = 2, byrow = TRUE
    )
  }

  # Create SpatVector polygons (one per point), then attach attributes
  crs_r <- terra::crs(r)

  # terra::vect expects a list of rings (each polygon can be list of matrices)
  geom_list <- lapply(poly_mats, function(m) if (is.null(m)) NULL else list(m))
  v <- terra::vect(geom_list, type = "polygons", crs = crs_r)

  # Add attributes (preserve input rows)
  out_df <- df
  out_df$.cell <- cells
  out_df$.row  <- rc[, 1]
  out_df$.col  <- rc[, 2]
  out_df$.inside_extent <- inside

  # Add snapped cell center (useful for QA)
  center_xy <- matrix(NA_real_, nrow = length(cells), ncol = 2)
  if (any(has_cell)) center_xy[has_cell, ] <- terra::xyFromCell(r, cells[has_cell])
  out_df$.cell_center_x <- center_xy[, 1]
  out_df$.cell_center_y <- center_xy[, 2]

  terra::values(v) <- out_df

  if (output == "terra") return(v)
  sf::st_as_sf(v)
}



## align_to_vector <- function(r, v, res, method = "bilinear", origin = NULL, buffer = 0) {
##   tmpl <- template_from_vector(v, res = res, origin = origin, buffer = buffer)
##   # One step: reproject + resample onto the template grid
##   ## terra::warp(r, tmpl, method = method)
##   warp(r, tmpl, method = method)
## }


template_from_vector <- function(v, res, crs = NULL, origin = NULL, buffer = 0) {
  stopifnot(inherits(v, "SpatVector"))
  e <- ext(v)
  if (buffer > 0) e <- ext(xmin(e) - buffer, xmax(e) + buffer, ymin(e) - buffer, ymax(e) + buffer)
  crs_use <- if (!is.null(crs)) crs else terra::crs(v)
  tmpl <- rast(e, resolution = res, crs = crs_use)
  if (!is.null(origin)) {
    origin(tmpl) <- origin  # origin is c(x0, y0)
  }
  tmpl
}



utm_points_to_pixel_polygons_per_point <- function(r, x, y = NULL, crs = NULL) {
  stopifnot(inherits(r, "SpatRaster"))
  if (is.data.frame(x)) {
    stopifnot(all(c("x", "y") %in% names(x)))
    xy <- as.matrix(x[, c("x", "y")])
  } else {
    stopifnot(is.numeric(x))
    if (is.null(y)) stop("Provide y if x is not a data.frame.")
    stopifnot(is.numeric(y), length(x) == length(y))
    xy <- cbind(x, y)
    colnames(xy) <- c("x", "y")
  }
  if (!is.null(crs) && !is.na(terra::crs(r)) && crs != terra::crs(r)) {
    v  <- terra::vect(xy, geom = c("x", "y"), crs = crs)
    v2 <- terra::project(v, terra::crs(r))
    xy <- terra::crds(v2)
  }
  ids <- terra::cellFromXY(r, xy)  # same length as points; may contain NA
  e <- terra::ext(r); rx <- terra::res(r)[1]; ry <- terra::res(r)[2]
  rows <- terra::rowFromCell(r, ids)
  cols <- terra::colFromCell(r, ids)
  xmin <- terra::xmin(e) + (cols - 1) * rx
  xmax <- xmin + rx
  ymax <- terra::ymax(e) - (rows - 1) * ry
  ymin <- ymax - ry
  # build polygons; keep NA for points outside raster
  polys <- vector("list", length(ids))
  ok <- !is.na(ids)
  polys[ok] <- lapply(which(ok), function(i) {
    sf::st_polygon(list(matrix(
      c(xmin[i], ymin[i],
        xmin[i], ymax[i],
        xmax[i], ymax[i],
        xmax[i], ymin[i],
        xmin[i], ymin[i]),
      ncol = 2, byrow = TRUE
    )))
  })
  sfc <- sf::st_sfc(polys, crs = sf::st_crs(terra::crs(r)))
  # return an sf with one row per point
  sf::st_sf(point_id = seq_along(ids), cell = ids, geometry = sfc)
}



## utm_points_to_pixel_multipolygon_fast <- function(r, x, y = NULL, crs = NULL, unique_cells = TRUE) {
##   stopifnot(inherits(r, "SpatRaster"))
##   # Accept df with x/y or numeric vectors
##   if (is.data.frame(x)) {
##     stopifnot(all(c("x", "y") %in% names(x)))
##     xy <- as.matrix(x[, c("x", "y")])
##   } else {
##     stopifnot(is.numeric(x))
##     if (is.null(y)) stop("Provide y if x is not a data.frame.")
##     stopifnot(is.numeric(y), length(x) == length(y))
##     xy <- cbind(x, y)
##     colnames(xy) <- c("x", "y")
##   }
##   # Reproject points if needed
##   if (!is.null(crs) && !is.na(terra::crs(r)) && crs != terra::crs(r)) {
##     v  <- terra::vect(xy, geom = c("x", "y"), crs = crs)
##     v2 <- terra::project(v, terra::crs(r))
##     xy <- terra::crds(v2)
##   }
##   ids <- terra::cellFromXY(r, xy)
##   ids <- ids[!is.na(ids)]
##   if (unique_cells) ids <- unique(ids)
##   if (!length(ids)) {
##     return(sf::st_sfc(sf::st_multipolygon(), crs = sf::st_crs(terra::crs(r))))
##   }
##   e <- terra::ext(r)
##   rx <- terra::res(r)[1]
##   ry <- terra::res(r)[2]
##   ncols <- terra::ncol(r)
##   nrows <- terra::nrow(r)
##   cols <- terra::colFromCell(r, ids)
##   rows <- terra::rowFromCell(r, ids)
##   # For north-up rasters: y decreases as row increases
##   xmin <- terra::xmin(e) + (cols - 1) * rx
##   xmax <- xmin + rx
##   ymax <- terra::ymax(e) - (rows - 1) * ry
##   ymin <- ymax - ry
##   polys <- lapply(seq_along(ids), function(i) {
##     sf::st_polygon(list(matrix(
##       c(xmin[i], ymin[i],
##         xmin[i], ymax[i],
##         xmax[i], ymax[i],
##         xmax[i], ymin[i],
##         xmin[i], ymin[i]),
##       ncol = 2, byrow = TRUE
##     )))
##   })
##   sfc <- sf::st_sfc(polys, crs = sf::st_crs(terra::crs(r)))
##   sf::st_cast(sf::st_union(sfc), "MULTIPOLYGON")
## }


## utm_points_to_pixel_multipolygon <- function(r, x, y = NULL, crs = NULL, unique_cells = TRUE) {
##   stopifnot(inherits(r, "SpatRaster"))
##   # Accept (x,y) vectors or a data.frame with x/y
##   if (is.data.frame(x)) {
##     stopifnot(all(c("x", "y") %in% names(x)))
##     xy <- as.matrix(x[, c("x", "y")])
##   } else {
##     stopifnot(is.numeric(x))
##     if (is.null(y)) stop("Provide y if x is not a data.frame.")
##     stopifnot(is.numeric(y), length(x) == length(y))
##     xy <- cbind(x, y)
##     colnames(xy) <- c("x", "y")
##   }
##   # If points CRS differs from raster CRS, project via SpatVector, then extract coords
##   if (!is.null(crs) && !is.na(terra::crs(r)) && crs != terra::crs(r)) {
##     v  <- terra::vect(xy, geom = c("x", "y"), crs = crs)
##     v2 <- terra::project(v, terra::crs(r))
##     xy <- terra::crds(v2)  # matrix with x/y in raster CRS
##   }
##   # IMPORTANT: pass matrix/data.frame, not SpatVector
##   ids <- terra::cellFromXY(r, xy)
##   ids <- ids[!is.na(ids)]
##   if (unique_cells) ids <- unique(ids)
##   if (!length(ids)) {
##     return(sf::st_sfc(sf::st_multipolygon(), crs = sf::st_crs(terra::crs(r))))
##   }
##   # Use a 1-layer template to polygonize selected cells
##   tmpl <- r[[1]]
##   terra::values(tmpl) <- seq_len(terra::ncell(tmpl))
##   sel <- terra::as.polygons(tmpl, values = TRUE)[ids]
##   sf::st_cast(sf::st_union(sf::st_geometry(sf::st_as_sf(sel))), "MULTIPOLYGON")
## }



## utm_points_to_pixel_multipolygon <- function(r, x, y = NULL, crs = NULL, unique_cells = TRUE) {
##   stopifnot(inherits(r, "SpatRaster"))
##   # Accept (x,y) vectors or a data.frame with x/y
##   if (is.data.frame(x)) {
##     stopifnot(all(c("x","y") %in% names(x)))
##     df <- x[, c("x","y")]
##   } else {
##     stopifnot(is.numeric(x))
##     if (is.null(y)) stop("Provide y if x is not a data.frame.")
##     stopifnot(is.numeric(y), length(x) == length(y))
##     df <- data.frame(x = x, y = y)
##   }
##   v <- terra::vect(df, geom = c("x","y"), crs = if (!is.null(crs)) crs else terra::crs(r))
##   if (!is.null(crs) && !is.na(terra::crs(r)) && terra::crs(v) != terra::crs(r)) {
##     v <- terra::project(v, terra::crs(r))
##   }
##   ids <- terra::cellFromXY(r, v)
##   ids <- ids[!is.na(ids)]
##   if (unique_cells) ids <- unique(ids)
##   if (!length(ids)) {
##     return(sf::st_sfc(sf::st_multipolygon(), crs = sf::st_crs(terra::crs(r))))
##   }
##   # Use a single-layer template with cell IDs
##   tmpl <- r[[1]]
##   terra::values(tmpl) <- seq_len(terra::ncell(tmpl))
##   sel <- terra::as.polygons(tmpl, values = TRUE)[ids]
##   sf::st_cast(sf::st_union(sf::st_geometry(sf::st_as_sf(sel))), "MULTIPOLYGON")
## }

## global_stats <- function(r) {
##     r = tmp
##   stopifnot(inherits(r, "SpatRaster"))
##   n <- terra::ncell(r)
##   stats <- terra::global(r, c("mean","sd"), na.rm=TRUE)
##   gmean <- mean(stats$mean, na.rm=TRUE)
##   gsd <- sqrt((sum((n-1)*stats$sd^2, na.rm=TRUE) + sum(n*(stats$mean-gmean)^2, na.rm=TRUE)) / (n*nrow(stats)-1))
##   mm <- terra::minmax(r)
##   gmin <- min(mm[1,], na.rm=TRUE)
##   gmax <- max(mm[2,], na.rm=TRUE)
##   qs <- terra::quantile(r, c(0.05,0.95), na.rm=TRUE)
##   q5 <- terra::global(qs[[1]], "max", na.rm=TRUE)[1,1]
##   q95 <- terra::global(qs[[2]], "max", na.rm=TRUE)[1,1]
##   list(global_mean=gmean, global_sd=gsd, global_min=gmin, global_max=gmax, q05=q5, q95=q95)
## }


save_timeseries_raster <- function(r, dir_path) {
  # Create directory if missing
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  # Extract xmin and ymin from raster extent
  ext <- ext(r)
  xmin_val <- ext$xmin
  ymin_val <- ext$ymin
  # Build file name: 3PG_<xmin>_<ymin>.tif
  file_name <- sprintf("3PG_%s_%s.tif", xmin_val, ymin_val)
  # Full output path
  out_path <- file.path(dir_path, file_name)
  # Write multilayer raster into a single GeoTIFF
  writeRaster(
    r,
    filename = out_path,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "TILED=YES")
  )
  return(out_path)
}

polygon_sector <- function(kryck_bord,
                              cells_per_side = 10,
                              sector = 62,
                              buffer_length = 200){
kryck_border_sf. <- st_as_sf(kryck_bord)
kryck_border_sf <- st_cast(kryck_border_sf., 'POLYGON')[3,]
kryck_border_buff  <- st_buffer(kryck_border_sf, buffer_length)
grid.. <- st_make_grid(kryck_border_buff, n = cells_per_side)
grid_target <- grid..[sector,]
## grid_target_area_km2 <- set_units(st_area(st_as_sf(grid_target)), km^2)
## st_area(st_as_sf(grid_target))
kryck_border <- vect(grid_target)##<------ This is the ROI
poly <- list(kryck_border = kryck_border, grid.. = grid.., kryck_border_sf = kryck_border_sf)
class(poly)  <- append(class(poly), 'polygon_sector')
return(poly)}

plot.polygon_sector <- function(x, y = "grey", ...) {
  # 1. Get geometry/extents
  g_grid   <- x[["grid.."]]            # sf or sfc
  g_border <- x[["kryck_border_sf"]][,1]
  v_terra  <- x[["kryck_border"]]      # SpatVector (terra)
  # 2. Common xlim/ylim from the grid (largest)
  bb <- sf::st_bbox(g_grid)
  xlim <- c(bb["xmin"], bb["xmax"])
  ylim <- c(bb["ymin"], bb["ymax"])
  # 3. First layer: grid
  plot(g_grid,
       border = "grey",
       xlim = xlim,
       ylim = ylim,
       ...)
  # 4. Second layer: sf border
  plot(g_border, add = TRUE)
  # 5. Third layer: terra SpatVector
  terra::plot(v_terra, col = y, add = TRUE)
}


## plot.polygon_sector <- function(x,y = 'grey', ...){
## plot(x[['grid..']], border = 'grey', ...)
## plot(x[['kryck_border_sf']][1,], add = TRUE)
## plot(x[['kryck_border']], col = y, add = TRUE)}

## plot.polygon_sector <- function(x, fill_border = "grey", ...) {
##   # 1) Plot grid (sfc or sf)
##   sf::plot(
##     x[["grid.."]],
##     border = "grey",
##     col    = NA,
##     ...
##   )
##   # 2) Add the sf border polygon (first feature only, as you did)
##   sf::plot(
##     x[["kryck_border_sf"]][1, ],
##     add   = TRUE,
##     border = "black",
##     col    = NA
##   )
##   # 3) Convert SpatVector -> sf and add it
##   kryck_border_sf2 <- sf::st_as_sf(x[["kryck_border"]])
##   sf::plot(
##     kryck_border_sf2,
##     add = TRUE,
##     col = fill_border,
##     border = NA
##   )
## }



create_multi_variable_netcdf<-function(inps,base_raster,outfile,time_origin=NULL,cores=4){
  ## library(ncdf4)
  ## library(terra)
  ## library(parallel)
  first_df<-inps[[1]]
  all_dates<-sort(unique(first_df$date))
  n_time<-length(all_dates)
  all_species<-sort(unique(unlist(lapply(inps,function(x)x$species))))
  n_spp<-length(all_species)
  all_vars<-sort(unique(unlist(lapply(inps,function(x)x$variable))))
  n_var<-length(all_vars)
  n_cell<-length(inps)
  if(!is.null(names(inps))){
    cell_ids<-as.integer(names(inps))
    if(any(is.na(cell_ids)))cell_ids<-seq_len(n_cell)
  }else{
    cell_ids<-seq_len(n_cell)
  }
  coords<-xyFromCell(base_raster,cell_ids)
  x_coord<-coords[,1]
  y_coord<-coords[,2]
  if(is.null(time_origin))time_origin<-min(all_dates)
  time_origin<-as.Date(time_origin)
  time_num<-as.numeric(all_dates-time_origin)
  time_dim<-ncdim_def("time",paste0("days since ",format(time_origin,"%Y-%m-%d")),time_num)
  cell_dim<-ncdim_def("cell","index",seq_len(n_cell))
  spp_dim<-ncdim_def("species","index",seq_len(n_spp))
  var_dim<-ncdim_def("variable","index",seq_len(n_var))
  value_var<-ncvar_def("value","unknown",list(time_dim,cell_dim,spp_dim,var_dim),missval=-9999,prec="float",shuffle=TRUE,deflate=4)
  x_var<-ncvar_def("x","m",list(cell_dim),prec="double")
  y_var<-ncvar_def("y","m",list(cell_dim),prec="double")
  id_var<-ncvar_def("cell_id","index",list(cell_dim),prec="int")
  nc<-nc_create(outfile,list(value_var,x_var,y_var,id_var))
  ncvar_put(nc,x_var,x_coord)
  ncvar_put(nc,y_var,y_coord)
  ncvar_put(nc,id_var,cell_ids)
  date_index<-match(all_dates,all_dates)
  species_index<-match(all_species,all_species)
  var_index<-match(all_vars,all_vars)
  pb<-txtProgressBar(min=0,max=n_cell,style=3)
  worker<-function(i){
    df<-inps[[i]]
    ti<-match(df$date,all_dates)
    si<-match(df$species,all_species)
    vi<-match(df$variable,all_vars)
    out<-array(-9999,dim=c(n_time,n_spp,n_var))
    out[cbind(ti,si,vi)]<-df$value
    list(i=i,out=out)
  }
  results<-mclapply(seq_along(inps),worker,mc.cores=cores)
  for(res in results){
    ncvar_put(nc,value_var,res$out,start=c(1,res$i,1,1),count=c(n_time,1,n_spp,n_var))
    setTxtProgressBar(pb,res$i)
  }
  close(pb)
  ncatt_put(nc,0,"species_names",paste(all_species,collapse=","))
  ncatt_put(nc,0,"variable_names",paste(all_vars,collapse=","))
  ncatt_put(nc,0,"crs","EPSG:3006 SWEREF99 TM")
  nc_close(nc)
  invisible(outfile)
}

## create_multi_variable_netcdf<-function(inps,base_raster,outfile,time_origin=NULL){
##   ## library(dplyr)
##   ## library(tidyr)
##   ## library(purrr)
##   ## library(ncdf4)
##   ## library(terra)
##   ## library(tibble)
##   first_df<-inps[[1]]
##   all_dates<-sort(unique(first_df$date))
##   n_time<-length(all_dates)
##   all_species<-sort(unique(unlist(lapply(inps,function(x)unique(x$species)))))
##   n_spp<-length(all_species)
##   all_vars<-sort(unique(unlist(lapply(inps,function(x)unique(x$variable)))))
##   n_var<-length(all_vars)
##   n_cell<-length(inps)
##   if(!is.null(names(inps))){
##     cell_ids<-as.integer(names(inps))
##     if(any(is.na(cell_ids)))cell_ids<-seq_len(n_cell)
##   }else{
##     cell_ids<-seq_len(n_cell)
##   }
##   arr<-array(NA_real_,dim=c(n_time,n_cell,n_spp,n_var),
##              dimnames=list(time=as.character(all_dates),
##                            cell=as.character(cell_ids),
##                            species=all_species,
##                            variable=all_vars))
##   for(i in seq_along(inps)){
##     df<-inps[[i]]%>%
##       mutate(species=factor(species,levels=all_species),
##              variable=factor(variable,levels=all_vars))%>%
##       select(date,species,variable,value)%>%
##       right_join(tibble(date=all_dates),by="date")%>%
##       complete(date,species,variable)%>%
##       arrange(date)
##     cube<-df%>%
##       pivot_wider(names_from=variable,values_from=value)%>%
##       pivot_wider(names_from=species,values_from=all_vars)
##     mat<-array(as.numeric(unlist(cube[,-1])),dim=c(n_time,n_spp,n_var))
##     arr[,i,,]<-mat
##   }
##   coords<-terra::xyFromCell(base_raster,cell_ids)
##   x_coord<-coords[,1]
##   y_coord<-coords[,2]
##   if(is.null(time_origin))time_origin<-min(all_dates)
##   time_origin<-as.Date(time_origin)
##   time_num<-as.numeric(all_dates-time_origin)
##   time_dim<-ncdim_def("time",paste0("days since ",format(time_origin,"%Y-%m-%d")),time_num)
##   cell_dim<-ncdim_def("cell","index",seq_len(n_cell))
##   spp_dim<-ncdim_def("species","index",seq_len(n_spp))
##   var_dim<-ncdim_def("variable","index",seq_len(n_var))
##   value_var<-ncvar_def("value","unknown",list(time_dim,cell_dim,spp_dim,var_dim),missval=-9999,prec="float")
##   x_var<-ncvar_def("x","m",list(cell_dim),prec="double")
##   y_var<-ncvar_def("y","m",list(cell_dim),prec="double")
##   id_var<-ncvar_def("cell_id","index",list(cell_dim),prec="int")
##   nc<-nc_create(outfile,list(value_var,x_var,y_var,id_var))
##   ncvar_put(nc,value_var,arr)
##   ncvar_put(nc,x_var,x_coord)
##   ncvar_put(nc,y_var,y_coord)
##   ncvar_put(nc,id_var,cell_ids)
##   ncatt_put(nc,0,"species_names",paste(all_species,collapse=","))
##   ncatt_put(nc,0,"variable_names",paste(all_vars,collapse=","))
##   ncatt_put(nc,0,"crs","EPSG:3006 SWEREF99 TM")
##   nc_close(nc)
##   invisible(outfile)
## }

latitude_raster <- function(r) {
  # Check input
  stopifnot(inherits(r, "SpatRaster"))
  # Get cell centers in source CRS
  xy <- terra::xyFromCell(r, 1:terra::ncell(r))
  # Project those centers to WGS84 (lon/lat)
  ll <- terra::project(xy, from = terra::crs(r), to = "EPSG:4326")
  # Create a new raster with same geometry
  lat_r <- terra::rast(r)
  terra::values(lat_r) <- ll[, 2]  # latitude (degrees)
  names(lat_r) <- "latitude_deg"
  return(lat_r)
}

## water_availability_mm <- function(wa_dt, depth = "0_289") {
##   # Validate input
##   if (!inherits(wa_dt, "SpatRaster")) {
##     stop("Input must be a terra SpatRaster object.")
##   }
##   # Extract NetCDF file path from source
##   src <- sources(wa_dt)$source[1]
##   ncfile <- sub("NETCDF:\"", "", src)
##   ncfile <- sub("\":.*$", "", ncfile)
##   # Read layer depths from NetCDF metadata
##   nc <- nc_open(ncfile)
##   d1 <- ncatt_get(nc, "swvl1", "GRIB_depthBelowLandLayer")$value
##   d2 <- ncatt_get(nc, "swvl2", "GRIB_depthBelowLandLayer")$value
##   d3 <- ncatt_get(nc, "swvl3", "GRIB_depthBelowLandLayer")$value
##   d4 <- ncatt_get(nc, "swvl4", "GRIB_depthBelowLandLayer")$value
##   nc_close(nc)
##   depths <- c(d1, d2, d3, d4)
##   thickness <- diff(c(depths, 289)) / 100  # convert cm → m
##   names(thickness) <- c("swvl1", "swvl2", "swvl3", "swvl4")
##   # Layer selection
##   i1 <- grep("^swvl1", names(wa_dt))
##   i2 <- grep("^swvl2", names(wa_dt))
##   i3 <- grep("^swvl3", names(wa_dt))
##   i4 <- grep("^swvl4", names(wa_dt))
##   sw1 <- wa_dt[[i1]] * (thickness["swvl1"] * 1000)
##   sw2 <- wa_dt[[i2]] * (thickness["swvl2"] * 1000)
##   sw3 <- wa_dt[[i3]] * (thickness["swvl3"] * 1000)
##   sw4 <- wa_dt[[i4]] * (thickness["swvl4"] * 1000)
##   # Depth selection
##   depth <- match.arg(depth, c("0_7", "0_28", "0_100", "0_289"))
##   wa_total <- switch(
##     depth,
##     "0_7"   = sw1,
##     "0_28"  = sw1 + sw2,
##     "0_100" = sw1 + sw2 + sw3,
##     "0_289" = sw1 + sw2 + sw3 + sw4
##   )
##   # Fix time
##   if (!all(is.na(time(wa_dt)))) {
##     time(wa_total) <- time(wa_dt)[i1]
##   }
##   nm_depth <- paste0("WA_", depth, "_mm_")
##   if (!all(is.na(time(wa_total)))) {
##     names(wa_total) <- paste0(nm_depth, format(time(wa_total), "%Y-%m-%d"))
##   } else {
##     names(wa_total) <- paste0(nm_depth, seq_len(nlyr(wa_total)))
##   }
##   return(wa_total)
## }

## agg_and_crop <- function(rast, fact, pol){
## rast <- c(rast)
##     kryck_elev. <- aggregate(rast, fact = fact, fun = mean)  # mean for DEM
##     kryck_elev <- crop(kryck_elev., pol)
## return(kryck_elev)
## }

agg_and_crop. <- function(r, target_res, pol = NULL, mask = FALSE) {
  # Ensure inputs are terra objects
  if (!inherits(r, "SpatRaster")) r <- rast(r)
  if (!is.null(pol) && !inherits(pol, "SpatVector")) pol <- vect(pol)
  # Input resolution (assumes square cells)
  res_in <- res(r)[1]
  # Compute ratio between target and input resolution
  ratio <- target_res / res_in
  # Dispatch automatically
  if (ratio > 1) {
    message(sprintf("Aggregating raster from %.2f to %.2f (factor = %.2f)", res_in, target_res, ratio))
    r_out <- aggregate(r, fact = ratio, fun = mean, na.rm = TRUE)
  } else if (ratio < 1) {
    fact_d <- 1 / ratio
    message(sprintf("Disaggregating raster from %.2f to %.2f (factor = %.2f)", res_in, target_res, fact_d))
    r_out <- disagg(r, fact = fact_d, method = "bilinear")
  } else {
    message("No resampling applied (resolutions are equal).")
    r_out <- r
  }
  # Crop and optionally mask
  if (!is.null(pol)) {
    r_out <- crop(r_out, pol)
    if (isTRUE(mask)) r_out <- mask(r_out, pol)
  }
  return(r_out)
}

agg_and_crop.. <- function(r, target_res, pol = NULL, mask = FALSE, snap = "out") {
  # Ensure raster is a SpatRaster
  if (!inherits(r, "SpatRaster")) r <- rast(r)
  # Handle polygon or raster input for cropping
  if (!is.null(pol)) {
    if (inherits(pol, "SpatVector")) {
      # Project polygon to raster CRS if needed (cheap)
      if (!identical(crs(pol), crs(r))) pol <- project(pol, crs(r))
    } else if (!inherits(pol, "SpatRaster")) {
      stop("`pol` must be a SpatVector (polygon) or SpatRaster (for extent).")
    }
  }
  # Original resolution
  res_in <- res(r)[1]
  ratio  <- target_res / res_in
  # Initialize output raster
  r_out <- r
  # Determine operation type
  if (ratio > 1) {
    # Target resolution is coarser → aggregate
    fact <- max(1L, round(ratio))
    message(sprintf("Aggregating: %.2f → %.2f (fact = %d)", res_in, target_res, fact))
    r_out <- aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
  } else if (ratio < 1) {
    # Target resolution is finer → disaggregate
    fact <- max(1L, round(1 / ratio))
    message(sprintf("Disaggregating: %.2f → %.2f (fact = %d)", res_in, target_res, fact))
    r_out <- disagg(r, fact = fact, method = "bilinear")
  } else {
    message("No resampling (same resolution).")
  }
  # Crop and optionally mask
  if (!is.null(pol)) {
    r_out <- crop(r_out, pol, snap = snap)
    if (isTRUE(mask) && inherits(pol, "SpatVector")) {
      r_out <- mask(r_out, pol)
    }
  }
  return(r_out)
}

import_feature_from_folder <- function(fun, folder, pattern){
    ## nc_folder <- '/media/wihe/LaCie/projects/wihe_3PG/climate_SW'
    nc_files <- dir(folder)[grepl(pattern, dir(folder))]
    nc_paths <- file.path(folder, nc_files)
    if(is.null(fun))
        return(nc_paths)
    fun.  <- match.fun(fun)
    imp <- fun.(nc_paths)
    return(imp)
}

map_coarse_to_fine <- function(coarse, fine) {
    coarse <- terra::project(coarse, terra::crs(fine), threads = FALSE)
  ## coarse <- align_to_ref(coarse, fine, method = "near")
  values(coarse) <- seq_len(ncell(coarse))
  ids <- terra::extract(coarse, terra::xyFromCell(fine, 1:terra::ncell(fine)))[, 2]
  fine_out <- fine
  values(fine_out) <- ids
  fine_out
}

add_title_first <- function(lst, title) {
# helper that adds title only to the first plot in the list
  lst2 <- lst
  lst2[[1]] <- lst2[[1]] +
    ggtitle(title) +
    theme(
      plot.title = element_text(
        hjust = 0.5,       # center
        ## face  = "bold",    # bold text
        size  = 12,        # bigger font
        margin = margin(b = 20) # add space below title
      )
    )
  Reduce(`/`, lst2)
}


warp_to_ref <- function(r, ref, method = "bilinear") {
  stopifnot(inherits(r, "SpatRaster"), inherits(ref, "SpatRaster"))

  out <- if (!terra::same.crs(r, ref)) {
    terra::project(r, ref, method = method)   # uses ref grid => origin aligned
  } else {
    terra::resample(r, ref, method = method)  # uses ref grid => origin aligned
  }

  # Force exact geometry (extent/rows/cols) to match ref
  out <- terra::crop(out, ref)

  out
}


## warp_to_ref <- function(r, ref, method = "bilinear") {
##   stopifnot(inherits(r, "SpatRaster"), inherits(ref, "SpatRaster"))

##   # If CRS differs, do everything in one shot on ref grid
##   if (!terra::same.crs(r, ref)) {
##     r <- terra::warp(r, ref, method = method)
##   } else {
##     # Same CRS: just resample to ref grid
##     r <- terra::resample(r, ref, method = method)
##   }

##   r
## }


align_to_ref <- function(r, ref, method = 'bilinear') {
# Function to align one raster to the reference
  # Reproject if CRS differs
  if (!terra::compareGeom(r, ref, stopOnError = FALSE)) {
    ## r <- terra::project(r, terra::crs(ref))
    r <- terra::project(r, terra::crs(ref), threads = FALSE)
  }
  # Resample to match resolution, extent, origin
  r <- terra::resample(r, ref, method = method)  # or "near" for categorical
  return(r)
}

biomass_components <- function(total_biomass_rast, species_rast, prop_table) {
  if (nlyr(species_rast) > 1) {
    species_rast <- species_rast[[1]]  # use first species layer
  }
  # Get categories (mapping of codes -> species names)
  cat_table <- cats(species_rast)[[1]]
  colnames(cat_table)[1] <- "code"
  colnames(cat_table)[2] <- "Tree_parts"
  # Drop "All" row from proportions table
  prop_table <- subset(prop_table, Tree_parts != "All")
  # Merge category table with lookup
  merged <- merge(cat_table, prop_table, by = "Tree_parts", all.x = TRUE)
  # Build reclass matrices
  stem_mat   <- merged[, c("code", "Stem_prop")]
  branch_mat <- merged[, c("code", "Branch_prop")]
  stump_mat  <- merged[, c("code", "Stump_and_root_prop")]
  # Apply classify + multiply by total biomass
  stem_r   <- classify(species_rast, stem_mat)   * total_biomass_rast
  branch_r <- classify(species_rast, branch_mat) * total_biomass_rast
  stump_r  <- classify(species_rast, stump_mat)  * total_biomass_rast
  names(stem_r)   <- "Stem_biomass"
  names(branch_r) <- "Branch_biomass"
  names(stump_r)  <- "Root_stump_biomass"
  return(c(stem_r, branch_r, stump_r))
}

central_latitude_raster <- function(rast) {
  # Get coordinates of pixel centers
  coords <- crds(rast)          # matrix of x (lon) and y (lat)
  lat_values <- coords[,2]      # extract latitude
  # Create new raster with same geometry
  lat_rast <- rast
  values(lat_rast) <- lat_values
  return(lat_rast)
}

clim_to_dt_long <- function(clim, nm = 'value'){
terra::as.data.frame(clim, xy = TRUE, cells = TRUE) %>%
  pivot_longer(
    cols = starts_with("X20"),      # or use: names(s)
    names_to = "date",
    values_to = nm) %>%
    separate(date, c('Xyear', 'month')) %>%
    mutate(year = as.numeric(str_remove(Xyear, 'X')),
           month = as.numeric(month)) %>%
    ## select(cell, year, month, nm) %>%
    select(cell, year, month, (nm)) %>%
    arrange(cell, year, month) 
}

compute_plantation_raster <- function(age_raster, ref_date) {
  ref_date <- ymd(ref_date)
  # Function to compute plantation date for a numeric vector
  calc_dates <- function(age_vec) {
    sapply(age_vec, function(age) {
      yrs <- floor(age)
      months <- round((age - yrs) * 12)
      plant_date <- ref_date %m-% years(yrs) %m-% months(months)
      format(plant_date, "%Y-%m")
    })
  }
  if (nlyr(age_raster) == 1) {
    # Single-layer raster
    char_vals <- calc_dates(values(age_raster))
    plantation_rast <- rast(age_raster)
    values(plantation_rast) <- char_vals
  } else {
    # Multi-layer raster: process each layer separately
    plantation_rast <- age_raster
    for (i in 1:nlyr(age_raster)) {
      char_vals <- calc_dates(values(age_raster[[i]]))
      values(plantation_rast[[i]]) <- char_vals
    }
  }
  return(plantation_rast)
}

crop_raster_with_polygon <- function(raster_path, polygon_path, mask = FALSE, output_path = NULL, buffer_width = NULL) {
  # Load raster and polygon
  r <- if (is.character(raster_path)) rast(raster_path) else raster_path
  p <- if (is.character(polygon_path)) vect(polygon_path) else polygon_path
  # Check if polygon is in geographic coordinates (degrees)
  ## if (terra::is.lonlat(p)) {
  ##   # Get UTM zone from centroid
  ##   centroid <- centroids(p)[1, ]
  ##   lon <- crds(centroid)[1]
  ##   zone <- floor((lon + 180) / 6) + 1
  ##   utm_crs <- paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=m +no_defs")
    ## utm_crs  <- crs(r)
  ##   # Reproject to UTM
    ## p_proj <- project(p, utm_crs)
    p_proj <- project(p, crs(r))
  ## } else {
  ##   p_proj <- p
  ## }
  # Buffer in meters
p_buffered <- p_proj
    if(!is.null(buffer_width))
        p_buffered <- buffer(p_proj, width = buffer_width)
  # Reproject back to raster CRS
  p_final <- project(p_buffered, terra::crs(r))
  # Crop and mask
  r_cropped <- crop(r, p_final)
  r_masked  <- if (mask) mask(r_cropped, p_final) else r_cropped
  # Save output if requested
  if (!is.null(output_path)) {
    writeRaster(r_masked, output_path, overwrite = TRUE)
  }
  return(r_masked)
}


crop_raster_fast <- function(r, p,
                             mask = TRUE,
                             snap = "out",
                             buffer_width = NULL,
                             return_in_polygon_crs = TRUE) {
  if (!inherits(r, "SpatRaster")) r <- terra::rast(r)
  if (!inherits(p, "SpatVector")) p <- terra::vect(p)

  # 1) Project polygon -> raster CRS
  p_r <- if (!terra::same.crs(p, r)) terra::project(p, terra::crs(r)) else p

  # 2) Optional buffer
  if (!is.null(buffer_width) && buffer_width != 0) {
    if (terra::is.lonlat(p_r)) {
      cen <- terra::crds(terra::centroids(p, inside = TRUE))[1, ]
      aeqd <- sprintf(
        "+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
        cen[2], cen[1]
      )
      p_m <- terra::project(p, aeqd)
      p_m <- terra::buffer(p_m, width = buffer_width)
      p_r <- terra::project(p_m, terra::crs(r))
    } else {
      p_r <- terra::buffer(p_r, width = buffer_width)
    }
  }

  # 3) Crop and mask
  r_chunk <- terra::crop(r, p_r, snap = snap)
  if (isTRUE(mask)) r_chunk <- terra::mask(r_chunk, p_r)

  # 4) Reproject back to polygon CRS if requested
  if (isTRUE(return_in_polygon_crs) && !terra::same.crs(r_chunk, p)) {
    method <- if (any(terra::is.factor(r_chunk))) "near" else "bilinear"
    r_chunk <- terra::project(r_chunk, terra::crs(p), method = method)
  }

  r_chunk
}

## crop_raster_fast <- function(r, p,
##                              mask = TRUE,
##                              snap = "out",
##                              buffer_width = NULL,      # meters; negative = inward buffer
##                              return_in_polygon_crs = TRUE) {
##   # Accept paths or terra objects
##   if (!inherits(r, "SpatRaster")) r <- rast(r)
##   if (!inherits(p, "SpatVector")) p <- vect(p)
##   # 1) Project polygon -> raster CRS (cheap)
##   p_r <- if (!identical(crs(p), crs(r))) project(p, crs(r)) else p
##   # 2) Optional buffer (in meters)
##   if (!is.null(buffer_width) && buffer_width != 0) {
##     if (is.lonlat(p_r)) {
##       # Raster CRS uses degrees: buffer in a temporary local metric CRS, then return
##       cen <- crds(centroids(p, inside = TRUE)[1, ])
##       aeqd <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
##                       cen[2], cen[1])
##       p_m  <- project(p, aeqd)
##       p_m  <- buffer(p_m, width = buffer_width)
##       p_r  <- project(p_m, crs(r))
##     } else {
##       # Projected CRS (typically meters): buffer directly
##       p_r <- buffer(p_r, width = buffer_width)
##     }
##   }
##   # 3) Crop (and mask) in raster CRS
##   r_chunk <- crop(r, p_r, snap = snap)
##   if (isTRUE(mask)) r_chunk <- mask(r_chunk, p_r)
##   # 4) Optionally reproject the small chunk back to polygon CRS
##   if (isTRUE(return_in_polygon_crs) && !identical(crs(r_chunk), crs(p))) {
##     method <- if (any(is.factor(r_chunk))) "near" else "bilinear"
##     r_chunk <- project(r_chunk, crs(p), method = method)
##   }
##   r_chunk
## }

d_input_from_rast_back  <- function(rasters, colnames, blueprint, ...){
    dots <- rlang::list2(...)
    rasters <- Map(function(x){names(x)  <- paste0('Layer_',1:nlyr(x))
                          return(x)}, rasters)
    names(rasters) <- colnames
    joinBy  <- c('cell','layer')
    bprint <- sub("^d_", "prepare_", blueprint)
    bprint_fun <- match.fun(bprint)
    time  <- rasters %>%
        imap_lgl(~ !is.null(time(.x)) && all(!is.na(time(.x)))) %>% all()
    if(time){
        joinBy  <- c('time', 'cell')
        colnames  <- c('year','month', colnames)
    }
    var_from_rs <-  rasters %>%
    map(function(x){
        terra::as.data.frame(x, xy = FALSE, cells = TRUE, wide = FALSE, time = time) %>%
            { if (time) dplyr::select(., -layer) else . }
        }) %>%
    purrr::imap(~ .x %>% rename(!!.y := 'values')) %>%
        purrr::reduce(function(x, y)
            data.table::merge.data.table(
    data.table::as.data.table(x),
    data.table::as.data.table(y),
    by = joinBy,
    all = FALSE
  )) %>%
  {if ("time" %in% names(.)) {
      dplyr::mutate(., 
             year = lubridate::year(time),
             month = lubridate::month(time)) %>%
          dplyr::select(-time)
    ## } else {.}} #%>%
   } else {
       dplyr::distinct(.)}} %>%
    as.data.table() %>%
{
  # 1) ensure numeric columns are <dbl>
    num_cols  <- names(.)[vapply(., is.numeric, logical(1)) | vapply(., inherits, logical(1), what="integer")]
    if (length(num_cols)) {
      .[, (num_cols) := lapply(.SD, as.double), .SDcols = num_cols]
    }
    # 2) normalize fake NA strings in character columns (preserves class)
    char_cols <- names(.)[vapply(., is.character, logical(1))]
    if (length(char_cols)) {
      .[, (char_cols) := lapply(.SD, function(x) {
        x[x %in% c("NA", "<NA>", "NaN", "")] <- NA_character_
        x
      }), .SDcols = char_cols]
    }
    # 3) remove 'layer' column if present
    if ("layer" %in% names(.)) .[, layer := NULL]
    # 4) DROP any row with NA in ANY column  (this fixes your example)
    . <- .[stats::complete.cases(.), ]
# 5) reorder columns using blueprint data frame
  blueprint_df <- get(blueprint, envir = parent.frame())  # get the DF by name
  desired_order <- intersect(names(blueprint_df), names(.))  # keep only matching columns
  other_cols <- setdiff(names(.), desired_order)            # preserve extra cols at end
  data.table::setcolorder(., c(desired_order, other_cols))  # reorder in-place
    .
} %>%
split(by = "cell", keep.by = TRUE) %>%
    imap(~ .x[, !"cell"])  %>%
map(function(x) {
      rlang::exec(bprint_fun, x, !!!dots)
})    
    var_from_rs
    }
    
d_input_from_rast  <- function(rasters, colnames, blueprint, ...){
    rasters <- Map(function(x){names(x)  <- paste0('Layer_',1:nlyr(x))
                          return(x)}, rasters)
    names(rasters) <- colnames
    joinBy  <- c('cell','layer')
    ## joinBy  <- 'cell'
    bprint <- sub("^d_", "prepare_", blueprint)
    ## return(blueprint%in%c('d_climate', 'd_species'))
    ## time  <- rasters %>% imap_lgl(~ !is.null(time(.x)) && all(!is.null(time(.x)))) %>% all()
    time  <- rasters %>%
        imap_lgl(~ !is.null(time(.x)) && all(!is.na(time(.x)))) %>% all()
    if(time){
        joinBy  <- c('time', 'cell')
        colnames  <- c('year','month', colnames)
    }
    var_from_rs <-  rasters %>%
    map(function(x){
        terra::as.data.frame(x, xy = FALSE, cells = TRUE, wide = FALSE, time = time) %>%
            { if (time) dplyr::select(., -layer) else . }
        ## dplyr::select(-layer)
        ## %>%
        ## dplyr::select(-layer)
        }) %>%
    purrr::imap(~ .x %>% rename(!!.y := 'values')) %>%
## reduce(inner_join, by = joinBy)%>%
purrr::reduce(function(x, y) data.table::merge.data.table(
    data.table::as.data.table(x),
    data.table::as.data.table(y),
    by = joinBy,
    all = FALSE
  )) %>%
  {if ("time" %in% names(.)) {
      dplyr::mutate(., 
             year = lubridate::year(time),
             month = lubridate::month(time)) %>%
          dplyr::select(-time)
    ## } else {.}} #%>%
   } else {dplyr::distinct(.)}} %>%
    dplyr::group_split(cell, .keep = TRUE) %>%
    rlang::set_names(map_chr(., ~ unique(as.character(.x$cell)))) %>%
    purrr::map( ~ dplyr::select(.x, -cell) %>%
mutate(across(where(is.character), ~ {
    nums <- suppressWarnings(as.numeric(.x))
    if (any(!is.na(nums))) nums else .x
  }))
        )
bp_template_dt <- data.table::as.data.table(get(blueprint)[1L, ])
bprint_fun <- match.fun(bprint)
var_from_rs %>%
  imap(function(rast_df, nm) {
    n_pixels <- nrow(rast_df)
    # Replicate template efficiently
    bp_expanded <- bp_template_dt[rep(1L, n_pixels)]
    # Assign multiple columns at once
    bp_expanded[, (colnames) := rast_df[, colnames, drop = FALSE]]
    # Drop incomplete rows efficiently
    bp_expanded <- bp_expanded[stats::complete.cases(bp_expanded)]
    bprint_fun(bp_expanded, ...)
  })
}

get_cell_df <- function(df_plot, cell_id) {
  # Ensure cell_id is numeric
  cell_id_num <- as.numeric(cell_id)
  # Defensive check
  if (cell_id_num < 1 || cell_id_num > nrow(df_plot)) {
    stop("cell_id is outside the range of df_plot rows")
  }
  # Extract coordinates
  coords <- df_plot[cell_id_num, c("x", "y")]
  # Return as data frame
  data.frame(x = coords$x, y = coords$y)
}

make_date_raster <- function(ref_rast, date_str) {
  # Create empty raster with same geometry as reference
  r <- rast(ref_rast)
  # Fill all pixels with the same date string
  values(r) <- date_str
  return(r)
}

make_parameters_from_species <- function(species.){
i_parameters_lit %>%
    filter(species %in% species.) %>%
 filter(!(species == "Pinus sylvestris") | str_detect(country, "Sweden")) %>%
    dplyr::select(species, type, country, any_of(d_parameters$parameter)) %>% 
   mutate(species = species) %>%
    dplyr::select(-type, -country) %>%
     pivot_longer(
         cols = -c(species),
                 names_to = "parameter",
                 values_to = "value")%>%
    pivot_wider(
        names_from = species,
        values_from = value
    ) %>%  prepare_parameters(sp_names = species.)
}


## monthly_clim_roi_from_daily_raster <- function(daily_raster, roi, ..., fun = 'mean', t1 = '2000-01-01', t2 = '2014-12-31'){
monthly_clim_roi_from_daily_raster. <- function(daily_raster, roi, fun = 'mean', t1 = '2000-01-01', t2 = '2014-12-31'){
    r_masked  <- daily_raster
    if(is.character(r_masked))
        r_masked  <- rast(r_masked)
    if(is.character(roi))
        roi <- vect(roi)
r_time <- time(r_masked)
r_masked  <- subset(r_masked, as.Date(time(r_masked)) >= as.Date(t1) &
                              as.Date(time(r_masked)) < as.Date(t2))
    if(!is.null(roi)){
## crs. <- crs(roi)
crs. <- terra::crs(roi)
p <- terra::project(roi, r_masked)
## r_masked <- crop_raster_with_polygon(r_masked, p, ...)
r_masked  <- crop(r_masked, p)
    r_masked <- terra::project(r_masked, crs.)
    }
r_time <- time(r_masked)
group_index <- format(r_time, '%Y-%m')
    if(fun == 'frost_days'){
        fun <- function(x)length(x[x<(0)])
monthly_max <- tapp(r_masked, index = group_index, fun = fun)
        }else{
monthly_max <- tapp(r_masked, index = group_index, fun = fun, na.rm = TRUE)}
time(monthly_max) <- names(monthly_max) %>% str_remove('X') %>% str_replace('\\.','-') %>% ym() #%>% format('%Y-%m')
    ## time(monthly_max)  <- as.Date(time(group_index))
return(monthly_max)
}


## monthly_clim_roi_from_daily_raster <- function(daily_raster, roi, ..., fun = 'mean', t1 = '2000-01-01', t2 = '2014-12-31'){
monthly_clim_roi_from_daily_raster <- function(daily_raster, roi, fun = 'mean', t1 = '2000-01-01', t2 = '2014-12-31'){
    r_masked  <- daily_raster
    if(is.character(r_masked))
        r_masked  <- rast(r_masked)
    if(is.character(roi))
        roi <- vect(roi)
r_time <- time(r_masked)
r_masked  <- subset(r_masked, as.Date(time(r_masked)) >= as.Date(t1) &
                              as.Date(time(r_masked)) < as.Date(t2))
    if(!is.null(roi)){
## crs. <- crs(roi)
crs. <- terra::crs(roi)
p <- terra::project(roi, r_masked)
## r_masked <- crop_raster_with_polygon(r_masked, p, ...)
r_masked  <- crop(r_masked, p)
    r_masked <- terra::project(r_masked, crs.)
    }
r_time <- time(r_masked)
group_index <- format(r_time, '%Y-%m')
    if(fun == 'frost_days'){
        fun <- function(x)length(x[x<(0)])
monthly_max <- tapp(r_masked, index = group_index, fun = fun)
        }else{
monthly_max <- tapp(r_masked, index = group_index, fun = fun, na.rm = TRUE)}
time(monthly_max) <- names(monthly_max) %>% str_remove('X') %>% str_replace('\\.','-') %>% ym() #%>% format('%Y-%m')
    ## time(monthly_max)  <- as.Date(time(group_index))
return(monthly_max)
}

normalize_raster <- function(kryck_lat){
    app(kryck_lat, fun = function(x) {
  (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
})
}

plot_rotate_data <- function(r,
                             label = "",
                             option = "H",
                             limits = NULL,          # pass c(min, max) to keep same scale
                             pad_x = 0.4,            # right margin for label
                             label_dx = 0.05,
                             label_size = 6,
                             label_col = "black") {
  # vectorize + rotate
  g <- as.polygons(r, values = TRUE, na.rm = TRUE, , dissolve = FALSE) |>
    st_as_sf() |>
    rotate_data()
  # bounding box and label point
  bb <- st_bbox(g)
  pt <- st_as_sf(
    data.frame(
      x = bb["xmax"] + label_dx,
      y = (bb["ymin"] + bb["ymax"]) / 2
    ),
    coords = c("x","y"), crs = st_crs(g)
  )
  # first data column to fill
  val_col <- names(g)[1]
  vals <- g[[val_col]]
  # base plot
  p <- ggplot(g) +
    geom_sf(aes(fill = .data[[val_col]]), color = NA, show.legend = FALSE)
  # choose scale depending on data type
  if (is.numeric(vals)) {
    p <- p + scale_fill_viridis_c(
      option = option,
      begin = 0, end = 1,
      limits = limits,
      oob = scales::squish
    )
  } else {
    p <- p + scale_fill_viridis_d(option = option)
  }
  # add label and formatting
  p <- p +
    geom_sf_text(data = pt, aes(label = label),
                 color = label_col, size = label_size, hjust = 0) +
    theme_void() +
    coord_sf(
      xlim = c(bb["xmin"], bb["xmax"] + pad_x),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    )
  return(p)
}

reclass_char_stack <- function(r_codes, lookup) {
  # Ensure lookup is a data frame
  lookup <- as.data.frame(lookup)
  # Check that 'lookup' has a "code" column or at least 2 columns
  stopifnot("code" %in% names(lookup) || ncol(lookup) >= 2)
  # Ensure first column is named "code"
  names(lookup)[1] <- "code"
  # Internal function to reclassify one column from lookup
  reclass_layer <- function(colname) {
    # Extract raster values
    vals <- values(r_codes)
    # Match lookup values by "code" and extract new values for the column
    new_vals <- lookup[[colname]][match(vals, lookup$code)]
    if (all(!is.na(suppressWarnings(as.numeric(new_vals))))) {
      # Case 1: Column is numeric → assign directly as raster values
      new_vals <- as.numeric(new_vals)
      r_new <- r_codes
      values(r_new) <- new_vals
      names(r_new) <- colname
    } else {
      # Case 2: Column is categorical → convert to factor, assign labels
      fac <- factor(new_vals, levels = unique(new_vals))
      r_new <- r_codes
      values(r_new) <- as.integer(fac)
      # Build category table for terra::SpatRaster
      cats <- data.frame(
        value = seq_along(levels(fac)),
        label = levels(fac)
      )
      levels(r_new) <- cats
      names(r_new) <- colname
    }
    r_new
  }
  # Apply reclassification to all lookup columns except "code"
  out <- lapply(setdiff(names(lookup), "code"), reclass_layer)
  # Return as a multi-layer SpatRaster (stack)
  rast(out)
}

reclass_stack <- function(r_codes, lookup) {
  lookup <- as.data.frame(lookup)
  stopifnot("code" %in% names(lookup) || ncol(lookup) >= 2)
  # Ensure first column is named "code"
  names(lookup)[1] <- "code"
  # Force numeric when possible
  lookup[] <- lapply(lookup, function(x) {
    if (all(!is.na(as.numeric(as.character(x))))) {
      as.numeric(as.character(x))
    } else {
      as.character(x)
    }
  })
  reclass_layer <- function(colname) {
    vals <- values(r_codes)
    new_vals <- lookup[[colname]][match(vals, lookup$code)]
    r_new <- r_codes
    values(r_new) <- as.numeric(new_vals)
    names(r_new) <- colname
    # ensure raster is numeric, not categorical
    ## levels(r_new) <- NULL
    r_new
  }
  out <- lapply(setdiff(names(lookup), "code"), reclass_layer)
  rast(out)
}

replace_and_rename_column <- function(df_replace, df_main,
                                      keys = c("year", "month"),
                                      tmp_map = c("tmp_ave" = "tmp"),
                                      keep_cols = NULL) {
  # Rename columns in df_replace according to tmp_map if applicable
  to_rename <- tmp_map[names(tmp_map) %in% names(df_replace)]
  if (length(to_rename) > 0) {
    df_replace <- df_replace %>%
      rename(!!!to_rename)
  }
  # Determine columns to replace (common to both but not keys)
  common_cols <- intersect(setdiff(names(df_replace), keys), names(df_main))
  # Join data
  if (is.null(keys)) {
    # Rename only overlapping columns to avoid conflicts (not keep_cols)
    rename_these <- intersect(names(df_replace), names(df_main))
    df_replace_renamed <- df_replace
    names(df_replace_renamed)[names(df_replace_renamed) %in% rename_these] <-
      paste0(names(df_replace_renamed)[names(df_replace_renamed) %in% rename_these], ".new")
    df_joined <- bind_cols(df_main, df_replace_renamed)
  } else {
    df_joined <- left_join(df_main, df_replace, by = keys, suffix = c("", ".new"))
  }
  # Replace columns from df_replace into df_main
  for (col in common_cols) {
    new_col <- paste0(col, ".new")
    if (new_col %in% names(df_joined)) {
      df_joined[[col]] <- coalesce(df_joined[[new_col]], df_joined[[col]])
    }
  }
  # Columns to keep in final output
  final_cols <- names(df_main)
  if (!is.null(keep_cols)) {
    # Always keep the original name (not .new)
    final_cols <- union(final_cols, keep_cols)
  }
  # Drop all .new columns and return clean result
  df_final <- df_joined %>%
    select(any_of(final_cols))
  return(df_final)
}

replace_species_with_proportions <- function(table1, table2) {
  # build lookup vector from table2
  lookup <- unlist(table2[ , -1], use.names = TRUE)
  names(lookup) <- gsub("_", " ", names(table2)[-1])  # match species names with spaces
  table1 %>%
    # replace species names with their proportions
    mutate(across(starts_with("Layer"),
                  ~ ifelse(is.na(.x), NA, lookup[.x]))) %>%
    # compute scaled proportions
    rowwise() %>%
    mutate(
      avg_value = mean(c_across(starts_with("Layer")), na.rm = TRUE),
      n_species = sum(!is.na(c_across(starts_with("Layer")))),
      across(starts_with("Layer"),
             ~ ifelse(!is.na(.x), avg_value / n_species, NA))
    ) %>%
    ungroup() %>%
    select(-avg_value, -n_species)
}

rotate_data <- function(data, x_add = 0, y_add = 0) {
## Functions to tilt sf
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  rotate_matrix <- function(x){ 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geometry = .$geometry * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}
rotate_data_geom <- function(data, x_add = 0, y_add = 0) {
  shear_matrix <- function(){ matrix(c(2, 1.2, 0, 1), 2, 2) }
  rotate_matrix <- function(x) { 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  data %>% 
    dplyr::mutate(
      geom = .$geom * shear_matrix() * rotate_matrix(pi/20) + c(x_add, y_add)
    )
}

timeSeriesRaster <- function(template, ts_list, pixel_ids = NULL,
                             date_col = "date", value_col = "value") {
  n_cells <- ncell(template)
  if (is.null(pixel_ids)) {
    pixel_ids <- seq_along(ts_list)
  }
  if (any(as.numeric(pixel_ids) > n_cells)) {
    stop("pixel_ids contain values larger than the number of raster cells.")
  }
  # ---- 1. Build master date vector from all pixels
  all_dates <- sort(unique(do.call(c, lapply(ts_list, function(df) as.Date(df[[date_col]])))))
  n_time <- length(all_dates)
  # ---- 2. Initialize raster value matrix: rows=pixels, cols=time steps
  ts_mat <- matrix(NA_real_, nrow = n_cells, ncol = n_time,
                   dimnames = list(NULL, as.character(all_dates)))
  # ---- 3. Fill matrix with pixel data
  for (i in seq_along(ts_list)) {
    df <- ts_list[[i]]
    pix <- as.numeric(pixel_ids[i])  # convert to numeric raster index
    df[[date_col]] <- as.Date(df[[date_col]])
    idx <- match(df[[date_col]], all_dates)       # index in master dates
    valid <- !is.na(idx)
    if (any(valid)) {
      ts_mat[pix, idx[valid]] <- df[[value_col]][valid]
    }
  }
  # ---- 4. Build SpatRaster stack
  r_stack <- rast(template, nlyrs = n_time)
  for (t in seq_len(n_time)) {
    values(r_stack[[t]]) <- ts_mat[, t]
  }
  time(r_stack) <- all_dates
  return(r_stack)
}

