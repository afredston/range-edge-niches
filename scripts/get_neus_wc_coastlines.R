# this script uses maps of the US to generate a smoothed coastline for NEUS, WC
# in each case, a smoother function was applied incrementally until unwanted coastal features (bays, estuaries, etc.) vanished from the smoothed coastline

# https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html

library(here)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(raster)
library(smoothr)

usamap <- rnaturalearth::ne_countries(scale = "small", country = "united states of america", returnclass = "sf")[1] %>% 
  st_cast("MULTILINESTRING") # get basic map of the country 

## Northeast:

# set bounding boxes
neus.xmin=-78
neus.xmax=-66
neus.ymin=35
neus.ymax=45

neus.bbox1 <- st_set_crs(st_as_sf(as(raster::extent(neus.xmin, neus.xmax, neus.ymin, neus.ymax), "SpatialPolygons")), st_crs(usamap))
neus.bbox2 <- st_set_crs(st_as_sf(as(raster::extent(-78, -74, 42, 45), "SpatialPolygons")), st_crs(usamap)) # smaller bounding box to get rid of extra lines on the map 

neusmap <- usamap %>% 
  st_intersection(neus.bbox1) %>%  
  st_difference(neus.bbox2) # gets rid of extra non coastal line 

neus.smoothgeom <- neusmap %>% 
  smoothr::smooth(method="ksmooth", smoothness=8) %>% # smoother was applied incrementally more until the Chesapeake went away 
  as("Spatial") %>% 
  geom()

neus.geomdists <- pointDistance(neus.smoothgeom[-nrow(neus.smoothgeom), c("x", "y")], neus.smoothgeom[-1, c("x", "y")], lonlat=TRUE)
neus.coastdistdat <- data.frame(neus.smoothgeom[, c('x','y')], seglength=c(0, neus.geomdists))
neus.coastdistdat$lengthfromhere <- rev(cumsum(rev(neus.coastdistdat[,"seglength"])))
# first row should match st_length(smoothmap)

write_rds(neus.coastdistdat, here("processed-data","neus_coastdistdat.rds"))

## West Coast:

wc.ymin <- 30
wc.ymax <- 50
wc.xmin <- -126
wc.xmax <- -117

wc.bbox1 <- st_set_crs(st_as_sf(as(raster::extent(wc.xmin, wc.xmax, wc.ymin, wc.ymax), "SpatialPolygons")), st_crs(usamap))
wc.bbox2 <- st_set_crs(st_as_sf(as(raster::extent(-123.5, wc.xmax, 45, wc.ymax), "SpatialPolygons")), st_crs(usamap))
wc.bbox3 <- st_set_crs(st_as_sf(as(raster::extent(-124.5, wc.xmax, 48, wc.ymax), "SpatialPolygons")), st_crs(usamap))

wcmap <- usamap %>% 
  st_intersection(wc.bbox1) %>%
  st_difference(wc.bbox2) %>% # crop out Puget Sound
  st_difference(wc.bbox3)

wc.smoothgeom <- wcmap %>% 
  smoothr::smooth(method="ksmooth", smoothness=1)  %>% 
  as("Spatial") %>% 
  geom()

wc.geomdists <- pointDistance(wc.smoothgeom[-nrow(wc.smoothgeom), c("x", "y")], wc.smoothgeom[-1, c("x", "y")], lonlat=TRUE)
wc.coastdistdat <- data.frame(wc.smoothgeom[, c('x','y')], seglength=c(0, wc.geomdists))
wc.coastdistdat$lengthfromhere <- cumsum(rev(wc.coastdistdat[,"seglength"]))
saveRDS(wc.coastdistdat, here("processed-data","wc_coastdistdat.rds"))

rm(list=ls())