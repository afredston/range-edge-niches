# this script takes large temperature dataframes and crops them to a shapefile generated for each region
# the "masks" are created by getting a bathymetric map within the bounding box for each region, cropping it to a depth cutoff, and then ensuring it falls within the US EEZ (which eliminates lakes, Canadian waters, etc.)

library(here)
library(tidyverse)
library(raster)
library(sf)
library(oceanmap)
library(data.table)
source(here("functions","sfc_as_cols.R"))

neus_latrange <- c(35, 45)
neus_lonrange <- c(-78, -66) 
wc_latrange <- c(30, 50)
wc_lonrange <- c(-126, -117)
ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

# choose how far out into the ocean you want temperature data
ebs.depth.cutoff <- 300
neus.depth.cutoff <- 300
wc.depth.cutoff <- 400

# get masks for each region 
wc.bathy <- get.bathy(lon = wc_lonrange, lat = wc_latrange, visualize = F, res = 15) 
neus.bathy <- get.bathy(lon = neus_lonrange, lat = neus_latrange, visualize = F, res = 15) 
ebs.bathy <- get.bathy(lon = ebs_lonrange, lat = ebs_latrange, visualize = F, res = 15) 

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
  st_intersection(st_union(useez)) # keep only points within the EEZ; crop out lakes, Canada 

ebs.bathy.mask <- ebs.bathy %>% 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>%
  dplyr::filter(layer <= ebs.depth.cutoff) %>% 
  st_intersection(st_union(useez))

neus.bathy.mask <- neus.bathy %>% 
  as("SpatialPolygonsDataFrame") %>% 
  st_as_sf() %>%
  dplyr::filter(layer <= neus.depth.cutoff) %>% 
  st_intersection(st_union(useez))

# filter all temperature dataframes by the mask, and correct the time column 
# can be very slow 

# hadisst
wc_hadisst <- read_rds(here("processed-data","wc_hadisst_raw.rds"))$wc_hadisst %>% 
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
