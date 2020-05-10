# this script takes large temperature dataframes and crops them to a shapefile generated for each region
# the "masks" are created by getting a bathymetric map within the bounding box for each region, cropping it to a depth cutoff, and then ensuring it falls within the US EEZ (which eliminates lakes, Canadian waters, etc.)

# later on I added the COBE data because there are major discrepancies between OISST and HadISST in the years when they overlap; HadISST is systematically warmer, at least in the Northeast, where the comparison matters most (because we have the longest time-series). COBE is not a NOAA product and as of now there is no API to download it. It can be manually accessed here:
# https://psl.noaa.gov/data/gridded/data.cobe.html

# workflow:

# select bounding boxes, bathymetry cutoffs, and resolution of bathymetric maps
# make bathymetric masks
# import raw temperature datasets and confirm that they are all reported as points in the center of spatial cells (for 1-degree resolution datasets like COBE and Hadley, this will be in half-degree values; for 0.25-degree resolution like OISST, points should end in 0.125, 0.375, etc.) 
# crop temperature datasets to the extent of each region's bathymetric mask; note that this effectively retains temperature readings from spatial cells with centers that fell within the mask

##############
## load packages and functions
##############

library(here)
library(tidyverse)
library(raster)
library(sp)
library(sf)
library(oceanmap)
library(data.table)
library(lubridate)
source(here("functions","sfc_as_cols.R"))
here <- here::here # fix conflict with lubridate 
select <- dplyr::select

##############
## make bathymetric mask for each region 
##############

# define bounding boxes for masks 
neus_latrange <- c(35, 45)
neus_lonrange <- c(-78, -66) 
wc_latrange <- c(30, 50)
wc_lonrange <- c(-126, -117)
ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

# choose how far out into the ocean you want temperature data
ebs.depth.cutoff <- 300
neus.depth.cutoff <- 300
wc.depth.cutoff <- 600 # WC shelf is very steep so I increased this from 300m in 100m increments until the bathymetric mask did not have big gaps along the coast 

# get masks for each region 
wc.bathy <- get.bathy(lon = wc_lonrange, lat = wc_latrange, visualize = F, res = 4) 
neus.bathy <- get.bathy(lon = neus_lonrange, lat = neus_latrange, visualize = F, res = 4) 
ebs.bathy <- get.bathy(lon = ebs_lonrange, lat = ebs_latrange, visualize = F, res = 4) 

# get CRS for future reference
bathy.crs <- wc.bathy %>% # works for all regions 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>% 
  st_crs() 

# get shapefile of the US EEZ, reproject to match bathymetry 
eezs <- st_read(here("raw-data/World_EEZ_v10_20180221","eez_v10.shp")) # download from http://www.marineregions.org/downloads.php and move to raw-data folder
useez <- eezs %>% 
  dplyr::filter(Sovereign1 == "United States") %>% 
  st_transform(crs=bathy.crs) 

# get bathymetric masks; slow 
wc.bathy.mask <- wc.bathy %>% 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>% # retains CRS 
  dplyr::filter(layer <= wc.depth.cutoff) %>% # get rid of points over X m deep
  st_intersection(st_union(useez)) %>% # keep only points within the EEZ; crop out lakes, Canada 
  st_union() # merge polygons into one 
# plot(wc.bathy.mask)
# note that this still has gaps along the West Coast

ebs.bathy.mask <- ebs.bathy %>% 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>%
  dplyr::filter(layer <= ebs.depth.cutoff) %>% 
  st_intersection(st_union(useez)) %>%
  st_union()
# plot(ebs.bathy.mask) 

neus.bathy.mask <- neus.bathy %>% 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>%
  dplyr::filter(layer <= neus.depth.cutoff) %>% 
  st_intersection(st_union(useez)) %>%
  st_union()
# plot(neus.bathy.mask)

##############
## rasterize temperature data and resample coarser rasters
##############

proj4 <- CRS("+proj=longlat +ellps=WGS84 +no_defs")

# function to convert SST into a raster brick 
rasterFromSST <- function(sstdf){
  out <- sstdf %>% 
    mutate(time = as.Date(time)) %>% 
    rename("x"=lon,"y"=lat) %>% 
    select(x, y, everything()) %>% # put x before y for rasterFromXYZ
    pivot_wider(names_from="time",values_from="sst") %>%
    rasterFromXYZ(crs=proj4)
  return(out)
}

# get raw datasets 

# should have columns x, y, time, sst
neus.hadisst.df <- read_rds(here("processed-data","neus_hadisst_raw.rds"))$neus_hadisst %>% 
  filter(!is.na(sst)) 
wc.hadisst.df <- read_rds(here("processed-data","wc_hadisst_raw.rds")) %>% 
  filter(!is.na(sst)) 
ebs.hadisst.df <- read_rds(here("processed-data","ebs_hadisst_raw.rds")) %>% 
  filter(!is.na(sst)) 

neus.oisst.df <- read_rds(here("processed-data","neus_oisst_raw.rds"))$neus_oisst %>% 
  filter(!is.na(sst)) %>%
  select(-altitude)
wc.oisst.df <- read_rds(here("processed-data","wc_oisst_raw.rds"))$wc_oisst %>% 
  filter(!is.na(sst)) %>%
  select(-altitude)
ebs.oisst.df <- read_rds(here("processed-data","ebs_oisst_raw.rds"))$ebs_oisst %>% 
  filter(!is.na(sst)) %>%
  select(-altitude)

# convert to raster brick 
neus.hadisst.raster <- rasterFromSST(neus.hadisst.df)
wc.hadisst.raster <- rasterFromSST(wc.hadisst.df)
ebs.hadisst.raster <- rasterFromSST(ebs.hadisst.df)

neus.oisst.raster <- rasterFromSST(neus.oisst.df)
wc.oisst.raster <- rasterFromSST(wc.oisst.df)
ebs.oisst.raster <- rasterFromSST(ebs.oisst.df)

# resample the HadISST data to the same resolution as the OISST data 

##############
## crop all to bathymetric masks 
##############

# hadisst
wc_hadisst <- read_rds(here("processed-data","wc_hadisst_raw.rds")) %>% 
  filter(!is.na(sst)) %>% # do early to make the object smaller
  st_as_sf(coords=c("lon","lat"), crs = bathy.crs) %>% # make spatial object with same CRS as mask 
  st_join(wc.bathy.mask, left=FALSE) %>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) 

neus_hadisst <- read_rds(here("processed-data","neus_hadisst_raw.rds"))$neus_hadisst %>% 
  filter(!is.na(sst)) %>% 
  st_as_sf(coords=c("lon","lat"), crs = bathy.crs) %>% 
  st_intersection(neus.bathy.mask) %>% 
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) 

ebs_hadisst <- read_rds(here("processed-data","ebs_hadisst_raw.rds"))$ebs_hadisst %>% 
  filter(!is.na(sst)) %>% 
  st_as_sf(coords=c("lon","lat"), crs = bathy.crs) %>%  
  st_join(ebs.bathy.mask, left=FALSE) %>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  )

# oisst
wc_oisst <- read_rds(here("processed-data","wc_oisst_raw.rds"))$wc_oisst %>% 
  filter(!is.na(sst)) %>% 
  st_as_sf(coords=c("lon","lat"), crs = bathy.crs) %>%
  st_join(wc.bathy.mask, left=FALSE) %>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) 

neus_oisst <- read_rds(here("processed-data","neus_oisst_raw.rds"))$neus_oisst %>% 
  filter(!is.na(sst)) %>% 
  st_as_sf(coords=c("lon","lat"), crs = bathy.crs) %>%
  st_join(neus.bathy.mask, left=FALSE) %>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  )

ebs_oisst <- read_rds(here("processed-data","ebs_oisst_raw.rds"))$ebs_oisst %>% 
  filter(!is.na(sst)) %>% 
  st_as_sf(coords=c("lon","lat"), crs = bathy.crs) %>% 
  st_join(ebs.bathy.mask, left=FALSE) %>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  )

dfs <- list(neus_oisst, neus_hadisst, wc_oisst, wc_hadisst, ebs_oisst, ebs_hadisst)
names <- c("neus_oisst", "neus_hadisst", "wc_oisst", "wc_hadisst", "ebs_oisst", "ebs_hadisst")
names(dfs) <- names

for(i in 1:length(dfs)){
  saveRDS(dfs[i], here("processed-data",paste0(names[i],".rds")))
}
rm(list=ls())
