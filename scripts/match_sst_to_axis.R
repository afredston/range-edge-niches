# this script matches temperature data to the coastal distance (NEUS, WC) or rotated axis (EBS) metric output from VAST. it takes about 30 minutes on my machine.
library(here)
library(tidyverse)

#####
# starting with EBS because it's much faster
#####
ebs.survey.earliest <- 7
ebs.sst <- read_rds(here("processed-data","ebs_sst_corrected.rds"))
ebs.vast.axes <- read_rds(here("processed-data","ebs_coords_conversion.rds")) %>%
  rename(vastLon=Lon, vastLat=Lat)

# isolate coordinates 
ebs.xy <- ebs.sst %>%
  dplyr::select(x, y) %>%
  distinct()

# function to get all VAST coordinates based on minimizing Euclidean distance between lat/lon of temperature and VAST grid
get_axes <- function(lon, lat, axesdf) {
  tmp <- axesdf %>% 
    mutate(abs.diff.x2 = abs(vastLon-lon)^2,
           abs.diff.y2 = abs(vastLat-lat)^2,
           abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
           )) %>% 
    filter(abs.diff.xy == min(abs.diff.xy)) %>% 
    dplyr::select(vastLon, vastLat, E_km, N_km, line_km) 
  return(tmp)
}

ebs.sst.axes <- NULL
for(i in 1:nrow(ebs.xy)){
  tmp <- ebs.xy[i,]
  out <- get_axes(lon=tmp$x, lat=tmp$y, axesdf=ebs.vast.axes)
  out <- cbind(tmp$x, tmp$y, out)
  ebs.sst.axes <- rbind(ebs.sst.axes, out)
}
ebs.sst.axes <- ebs.sst.axes %>% rename("x"=`tmp$x`, "y"=`tmp$y`)

# join VAST axes to temperature data 
ebs.sst.rotated <- ebs.sst %>%
  left_join(ebs.sst.axes, by=c("x","y")) %>%
  mutate(year_measured = ifelse(month < ebs.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rename(lat=y, lon=x) 

saveRDS(ebs.sst.rotated, file=here("processed-data","ebs_sst_linedist.rds"))


#####
# add coastal distance measurement to temperature data for NEUS, WC
#####
neus.survey.earliest <- 3
neus.vast.coords <- readRDS(here("processed-data","neus_coords_conversion.rds"))%>%
  rename("vastLon"=Lon,
         "vastLat"=Lat)

wc.survey.earliest <- 5
wc.vast.coords <- readRDS(here("processed-data","wc_coords_conversion.rds"))%>%
  rename("vastLon"=Lon,
         "vastLat"=Lat)

# NOTE that this matches the points of the survey to points on the coast using simple distance minimization. This means that a point that is close to a coastline to the north and slightly further offshore to the west will be assigned a coastal distance further north than if I had just measured the latitude of the observation and matched it to the coast. However, I think this is more accurate because it assigns observations to the closest part of the shelf. 

get_length <- function(lon, lat, distdf) {
  tmp <- distdf %>% 
    mutate(abs.diff.x2 = abs(vastLon-lon)^2,
           abs.diff.y2 = abs(vastLat-lat)^2,
           abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
           )) %>% 
    filter(abs.diff.xy == min(abs.diff.xy)) %>% 
    dplyr::select(coast_km) %>% 
    pull()
  return(tmp)
}

# load raw SST data, assign each point a coastal distance. slow!
wc.sst.coastdist <- readRDS(here("processed-data","wc_sst_corrected.rds")) %>%
  mutate(year_measured = ifelse(month < wc.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rowwise() %>%
  mutate(coast_km = (get_length(lon=x, lat=y, distdf = wc.vast.coords))) %>% 
  ungroup()

neus.sst.coastdist <- readRDS(here("processed-data","neus_sst_corrected.rds")) %>%
  mutate(year_measured = ifelse(month < neus.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rowwise() %>%
  mutate(coast_km = (get_length(lon=x, lat=y, distdf = neus.vast.coords))) %>% 
  ungroup() 

saveRDS(neus.sst.coastdist, here("processed-data","neus_sst_coastdist.rds"))
saveRDS(wc.sst.coastdist, here("processed-data","wc_sst_coastdist.rds"))

# rm(list=ls())
