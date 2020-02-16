# this script matches temperature data to the coastal distance (NEUS, WC) or rotated axis (EBS) metric output from VAST. it is slow! 
library(here)
library(tidyverse)

#####
# starting with EBS because it's much faster
#####
ebs.survey.earliest <- 7
ebs.oisst <- read_rds(here("processed-data","ebs_oisst.rds"))$ebs_oisst
ebs.hadisst <- read_rds(here("processed-data","ebs_hadisst.rds"))$ebs_hadisst
ebs.vast.axes <- read_rds(here("processed-data","ebs_coords_conversion.rds")) %>%
  rename(vastLon=Lon, vastLat=Lat)

# isolate coordinates of each temperature dataset 
ebs.oisst.xy <- ebs.oisst %>%
  dplyr::select(x, y) %>%
  distinct()

ebs.hadisst.xy <- ebs.hadisst %>%
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
    dplyr::select(vastLon, vastLat, E_km, N_km, NE_km, NW_km) 
  return(tmp)
}

ebs.oisst.axes <- NULL
for(i in 1:nrow(ebs.oisst.xy)){
  tmp <- ebs.oisst[i,]
  out <- get_axes(lon=tmp$x, lat=tmp$y, axesdf=ebs.vast.axes)
  out <- cbind(tmp$x, tmp$y, out)
  ebs.oisst.axes <- rbind(ebs.oisst.axes, out)
}
ebs.oisst.axes <- ebs.oisst.axes %>% rename("x"=`tmp$x`, "y"=`tmp$y`)

ebs.hadisst.axes <- NULL
for(i in 1:nrow(ebs.hadisst.xy)){
  tmp <- ebs.hadisst[i,]
  out <- get_axes(lon=tmp$x, lat=tmp$y, axesdf=ebs.vast.axes)
  out <- cbind(tmp$x, tmp$y, out)
  ebs.hadisst.axes <- rbind(ebs.hadisst.axes, out)
}
ebs.hadisst.axes <- ebs.hadisst.axes %>% rename("x"=`tmp$x`, "y"=`tmp$y`)

# join VAST axes to temperature datasets 
ebs.oisst.rotated <- ebs.oisst %>%
  left_join(ebs.oisst.axes, by=c("x","y")) %>%
  mutate(year_measured = ifelse(month < ebs.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rename(lat=y, lon=x)%>%
  mutate(dataset="oisst") %>%
  select(-altitude)

ebs.hadisst.rotated <- ebs.hadisst %>%
  left_join(ebs.hadisst.axes, by=c("x","y")) %>%
  mutate(year_measured = ifelse(month < ebs.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rename(lat=y, lon=x)%>%
  mutate(dataset="hadisst")


# create joint dataframe using hadISST to populate years where there's no OISST data 
ebs.sst.rotated <- ebs.hadisst.rotated %>%
  filter(time < min(ebs.oisst.rotated$time)) %>% 
  full_join(ebs.oisst.rotated) # replace with rbind so it runs faster? 

saveRDS(ebs.sst.rotated, file=here("processed-data","ebs_sst_rotated.rds"))


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

# load clean but raw OISST data, assign each point a coastal distance, calculate annual summary stats for SST
# SLOW!
wc.oisst.coastdist <- readRDS(here("processed-data","wc_oisst.rds"))$wc_oisst %>%
  mutate(year_measured = ifelse(month < wc.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rowwise() %>%
  mutate(coast_km = (get_length(lon=x, lat=y, distdf = wc.vast.coords))) %>% 
  ungroup()  %>%
  mutate(dataset="oisst") %>%
  select(-altitude)

neus.oisst.coastdist <- readRDS(here("processed-data","neus_oisst.rds"))$neus_oisst %>%
  mutate(year_measured = ifelse(month < neus.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rowwise() %>%
  mutate(coast_km = (get_length(lon=x, lat=y, distdf = neus.vast.coords))) %>% 
  ungroup()  %>%
  mutate(dataset="oisst")%>%
  select(-altitude)

# do the same for hadISST
neus.hadisst.coastdist <- readRDS(here("processed-data","neus_hadisst.rds"))$neus_hadisst %>%
  mutate(year_measured = ifelse(month < neus.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rowwise() %>%
  mutate(coast_km = (get_length(lon=x, lat=y, distdf = neus.vast.coords))) %>% 
  ungroup() %>%
  mutate(dataset="hadisst")

wc.hadisst.coastdist <- readRDS(here("processed-data","wc_hadisst.rds"))$wc_hadisst %>%
  mutate(year_measured = ifelse(month < wc.survey.earliest, year-1, year), 
         year_match = year_measured + 1) %>% 
  filter(year_match > min(year_match)) %>% 
  dplyr::select(-year) %>%
  rowwise() %>%
  mutate(coast_km = (get_length(lon=x, lat=y, distdf = wc.vast.coords))) %>% 
  ungroup() %>%
mutate(dataset="hadisst")

# create joint dataframe using hadISST to populate years where there's no OISST data 
neus.sst.coastdist <- neus.hadisst.coastdist %>%
  filter(time < min(neus.oisst.coastdist$time)) %>% 
  full_join(neus.oisst.coastdist)

wc.sst.coastdist <- wc.hadisst.coastdist %>%
  filter(time < min(wc.oisst.coastdist$time)) %>% 
  full_join(wc.oisst.coastdist)

saveRDS(neus.sst.coastdist, here("processed-data","neus_sst_coastdist.rds"))
saveRDS(wc.sst.coastdist, here("processed-data","wc_sst_coastdist.rds"))

rm(list=ls())
