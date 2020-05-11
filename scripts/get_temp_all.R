# THIS SCRIPT TAKES HOURS TO RUN, so don't re-run the whole thing if you don't have to! the OISST datasets in particular are slow.

# get all temperature datasets from ERDDAP for all regions:
# https://coastwatch.pfeg.noaa.gov/erddap/index.html
library(here)
source(here("functions","get_rerddap.R"))

# select bounding boxes for all regions 
neus_latrange <- c(35, 45)
neus_lonrange <- c(-78, -66) 
wc_latrange <- c(30, 50)
wc_lonrange <- c(-126, -117)
ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

# program call for each dataset; get from the ERDDAP pages and update end years if more data is added
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdHadISST.html
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst2Agg_LonPM180.html

hadisst <- "erdHadISST"
neus_hadisst_years <- c(1967, 2018)
wc_hadisst_years <- c(1976, 2018) # don't bother downloading early years for regions with no trawl data 
ebs_hadisst_years <- c(1981, 2018)

neus_hadisst_time <- c("1967-01-16", "2018-12-16")

hadisst_fields <- "sst"
oisst <- "ncdcOisst2Agg_LonPM180"
oisst_years <- c(1982, 2018)
oisst_fields <- "sst"

neus_hadisst_grid <- griddap(hadisst, time=neus_hadisst_time, latitude = neus_latrange, longitude = neus_lonrange, fields=hadisst_fields)
# 
# testgrid <- griddap(hadisst, time=c("1967-01-16","1967-12-16"), latitude = neus_latrange, longitude = neus_lonrange, fields=hadisst_fields)
# testfile <- testgrid$summary$filename
# testnc <- nc_open(testfile)
# plot(testnc)
# sst <- ncvar_get(testnc, "sst")
# 
# # get filename
# neus_hadisst_nc_file <- neus_hadisst_grid$summary$filename
# 
# # import .nc file 
# testnc <- nc_open(testfile)
# 
# # get raster from netcdf by extracting dimensions and variables 
# 
# # get dimensions and variables
# lon <- ncvar_get(testnc, "longitude")
# lat <- ncvar_get(testnc, "latitude")
# time <- ncvar_get(testnc, "time")
# tunits <- ncatt_get(testnc,"time","units")
# sstarray <- ncvar_get(testnc, "sst")
# sstname <- ncatt_get(testnc, "sst", "long_name")
# sstunits <- ncatt_get(testnc,"sst","units")
# sst_fillvalue <- ncatt_get(testnc,"sst","_FillValue")
# 
# # get global attributes
# title <- ncatt_get(testnc,0,"title")
# institution <- ncatt_get(testnc,0,"institution")
# datasource <- ncatt_get(testnc,0,"source")
# references <- ncatt_get(testnc,0,"references")
# history <- ncatt_get(testnc,0,"history")
# Conventions <- ncatt_get(testnc,0,"Conventions")
# 
# plot(sstarray[,,1])
# 
# testraster <- raster::stack(testfile)
# 
# # save raster stack 
# 
# # let's try a tidy version of that with tidync 
# library(tidync)
# library(stars)

# 

# neus_hadisst <- get_rerddap(datasetID = hadisst, latrange = neus_latrange, lonrange = neus_lonrange, startyear=neus_hadisst_years[1], endyear=neus_hadisst_years[2], fields=hadisst_fields, verbose=TRUE)
# 
# neus_oisst <- get_rerddap(datasetID = oisst, latrange = neus_latrange, lonrange = neus_lonrange, startyear=oisst_years[1], endyear=oisst_years[2], fields=oisst_fields, verbose=TRUE)
# 
# wc_hadisst <- get_rerddap(datasetID = hadisst, latrange = wc_latrange, lonrange = wc_lonrange, startyear=wc_hadisst_years[1], endyear=wc_hadisst_years[2], fields=hadisst_fields, verbose=TRUE)
# wc_oisst <- get_rerddap(datasetID = oisst, latrange = wc_latrange, lonrange = wc_lonrange, startyear=oisst_years[1], endyear=oisst_years[2], fields=oisst_fields, verbose=TRUE)
# 
# ebs_hadisst <- get_rerddap(datasetID = hadisst, latrange = ebs_latrange, lonrange = ebs_lonrange, startyear=ebs_hadisst_years[1], endyear=ebs_hadisst_years[2], fields=hadisst_fields, verbose=TRUE)
# ebs_oisst <- get_rerddap(datasetID = oisst, latrange = ebs_latrange, lonrange = ebs_lonrange, startyear=oisst_years[1], endyear=oisst_years[2], fields=oisst_fields, verbose=TRUE)
# 
# # these are dataframes by default, but we'd like to keep them as netcdf files, to facilitate cropping/resampling
# 
# # get file names:
# neus_hadisst_file <- neus_hadisst$summary$filename
# 
# # if there is an error that says "Error in R_nc4_open: NetCDF: Unknown file format"
# # find the directory C:\Users\alexafh\AppData\Local\Cache/R/rerddap/ and delete its contents 
# 
# dfs <- list(neus_oisst, neus_hadisst, wc_hadisst, wc_oisst,ebs_hadisst, ebs_oisst)
# names <- c("neus_oisst", "neus_hadisst", "wc_hadisst", "wc_oisst","ebs_hadisst", "ebs_oisst")
# names(dfs) <- names
# 
# for(i in 1:length(dfs)){
#   saveRDS(dfs[i], here("processed-data",paste0(names[i],"_raw.rds")))
# }
# rm(list=ls())
