# ULTRA SLOW - LITERALLY TAKES DAYS - BEWARE 

library(rerddap)
library(rnaturalearth)
library(sf)
library(here)
library(tidyverse)
library(sp)
library(raster)

WGS84 <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
usoutline <- rnaturalearth::ne_states("united states of america", returnclass = "sf") %>% 
  # st_union() %>% #just to combine all polygons into a large one
  st_sf()
eezs <- st_read(here("raw-data/World_EEZ_v10_20180221","eez_v10.shp")) # http://www.marineregions.org/downloads.php
useez <- eezs %>% 
  dplyr::filter(Sovereign1 == "United States") %>% 
  st_transform(crs=WGS84) 

# set bounding boxes
neus_latrange <- c(34, 46)
neus_lonrange <- c(-78, -66) 
wc_latrange <- c(30, 50)
wc_lonrange <- c(-126, -116)
ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

ebs.depth.cutoff <- 300
neus.depth.cutoff <- 300
wc.depth.cutoff <- 400

# get bathymetry
neus_bathy_raw <- rerddap::griddap("etopo180", longitude=neus_lonrange, latitude = neus_latrange, fields = "altitude")$data
wc_bathy_raw <- rerddap::griddap("etopo180", longitude=wc_lonrange, latitude = wc_latrange, fields = "altitude")$data
ebs_bathy_raw <- rerddap::griddap("etopo180", longitude=ebs_lonrange, latitude = ebs_latrange, fields = "altitude")$data

# generate bathymetry sf objects that are cropped to X m and to the EEZ (really slow) 
neus_bathy <- SpatialPointsDataFrame(coords=neus_bathy_raw[,c("longitude","latitude")], data=neus_bathy_raw, proj4string=crs(WGS84)) %>%
  st_as_sf() %>%
  filter(altitude <= 0) %>%
  rename(depth = altitude) %>%
  mutate(depth=abs(depth)) %>%
  filter(depth <= neus.depth.cutoff) %>%
  st_intersection(st_union(useez)) # keep only points within the EEZ 

wc_bathy <- SpatialPointsDataFrame(coords=wc_bathy_raw[,c("longitude","latitude")], data=wc_bathy_raw, proj4string=crs(WGS84)) %>%
  st_as_sf() %>%
  filter(altitude <= 0) %>%
  rename(depth = altitude) %>%
  mutate(depth=abs(depth)) %>%
  filter(depth <= wc.depth.cutoff) %>%
  st_intersection(st_union(useez)) 

ebs_bathy <- SpatialPointsDataFrame(coords=ebs_bathy_raw[,c("longitude","latitude")], data=ebs_bathy_raw, proj4string=crs(WGS84)) %>%
  st_as_sf() %>%
  filter(altitude <= 0) %>%
  rename(depth = altitude) %>%
  mutate(depth=abs(depth)) %>%
  filter(depth <= ebs.depth.cutoff) %>%
  st_intersection(st_union(useez)) 

saveRDS(neus_bathy, here("processed-data","neus_bathy_300m.rds"))
saveRDS(ebs_bathy, here("processed-data","ebs_bathy_300m.rds"))
saveRDS(wc_bathy, here("processed-data","wc_bathy_400m.rds"))