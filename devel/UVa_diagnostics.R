library(terra); library(ncdf4)




nc <- nc_open('/run/media/wilar/T7 Shield/slu_data_bases/3PG/out_3PG_netcdf/GPP.nc')
print(nc)
nc_close(nc)

class(nc)
print(nc)
nc_close(nc)

f <- normalizePath(list.files('/run/media/wilar/T7 Shield/slu_data_bases/3PG/out_3PG_netcdf/',
                pattern = '\\.(tif|nc)$', full.names = TRUE))
length(f); head(basename(f), 3)


r <- rast(f[1])
r
crs(r, describe = TRUE)
head(names(r), 5); head(time(r), 5)
global(r[[100]], c('mean','min','max'), na.rm = TRUE)






r <- rast(shQuote(f[1]))
r                      # dims, resolution, extent, nlyr
crs(r, describe = TRUE)
head(names(r), 5); head(time(r), 5); class(time(r))
global(r[[100]], c('mean','min','max'), na.rm = TRUE)

# if NetCDF:
# nc <- nc_open(f[1]); print(nc); nc_close(nc)

## Run this before the full script
r <- rast(f[1])
i <- which(format(time(r), '%Y') %in% c('2016','2017','2018'))

## decisive: annual sums, compare against ~16 tDM/ha/yr
ann <- tapp(r[[i]], format(time(r)[i], '%Y'), sum, na.rm = TRUE)
global(ann, c('mean','min','max'), na.rm = TRUE)

## are failed pixels NA or silently 0?
l <- r[[100]]
c(ncell = ncell(l),
  valid = unlist(global(!is.na(l), 'sum', na.rm = TRUE)),
  zeros = unlist(global(l == 0, 'sum', na.rm = TRUE)))

## seasonal cycle, checks the spring question
mm <- tapp(r[[i]], format(time(r)[i], '%m'), mean, na.rm = TRUE)
round(unlist(global(mm, 'mean', na.rm = TRUE)), 3)
