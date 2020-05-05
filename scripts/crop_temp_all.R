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
## rasterize temperature data
##############

proj4 <- CRS("+proj=longlat +ellps=WGS84 +no_defs")

# get raw datasets and convert to raster brick 
# note that spatialpixelsdataframe assumes points are centers of cells
neus.hadisst.df <- read_rds(here("processed-data","neus_hadisst_raw.rds"))$neus_hadisst %>% 
  filter(!is.na(sst)) 

testdf <- neus.hadisst.df %>% 
  filter(time==max(time)) %>%
  mutate(time=as.Date(time))

testspdf <- SpatialPixelsDataFrame(points=testdf[,c('lon','lat')], data=testdf %>% select(-lat, -lon), proj4string=proj4)

testraster <- raster(testspdf)
teststack <- stack(testspdf)

spg <- testdf
coordinates(spg) <- ~lon+lat
gridded(spg) <- TRUE
rasterdf <- raster(spg)
rasterdf

jc <- neus.hadisst.df %>% 
  mutate(time = as.Date(time)) %>% 
  rename("x"=lon,"y"=lat) %>% 
  select(x, y, everything()) %>% # put x before y for rasterFromXYZ
  pivot_wider(names_from="time",values_from="sst") %>%
  rasterFromXYZ(crs=proj4)


# this works to produce a single raster 

stack_time_slices <- function(sstdf){
  slices <- unique(sstdf$time)
  slice1 <- slices[1]
  dat <- sstdf %>% filter(time==slice1) %>% select(lat, lon, sst)
  pixels <- SpatialPixelsDataFrame(points=dat[,c('lon','lat')], data=dat %>% select(-lat, -lon), proj4string=proj4)
  out <- stack(pixels)
  names(out) <- paste0(slice1)
  for(i in tail(slices, -1)){
    dat <- sstdf %>% filter(time==i) %>% select(lat, lon, sst)
    pixels <- SpatialPixelsDataFrame(points=dat[,c('lon','lat')], data=dat %>% select(-lat, -lon), proj4string=proj4)
    slice <- stack(pixels)
    names(slice) <- paste0(i)
    out <- stack(out, slice)
  }
  return(out)
}

# testdf <- neus.hadisst.df %>% filter(time==max(time)) %>% select(-time)
# testspdf <- SpatialPixelsDataFrame(points=testdf[,c('lon','lat')], data=testdf %>% select(-lat, -lon), proj4string=proj4)
# testraster <- raster(testspdf)
# plot(testraster)
# teststack <- stack(testspdf)
# names(teststack) <- "timeee"

# define bounding boxes for masks 
neus_latrange <- c(35, 45)
neus_lonrange <- c(-78, -66) 
wc_latrange <- c(30, 50)
wc_lonrange <- c(-126, -117)
ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

# compare hadisst and oisst without bathymetric masks (since they should have been downloaded for the same bounding box)
neus_oisst_all <- read_rds(here("processed-data","neus_oisst_raw.rds"))$neus_oisst %>% 
  filter(!is.na(sst)) %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) %>% 
  dplyr::select(-altitude) 

neus_oisst_summary <- neus_oisst_all %>%
  group_by(lat, lon, year, month) %>%
  mutate(sst.month.mean = mean(sst)) %>%
  ungroup() %>%
  select(lat, lon, year, month, sst.month.mean) %>%
  distinct() %>%
  group_by(lat, lon, year) %>%
  mutate(cell.min = min(sst.month.mean),
         cell.mean = mean(sst.month.mean),
         cell.max = max(sst.month.mean)) %>%
  ungroup() %>%
  select(lat, lon, year, cell.mean, cell.max) %>%
  distinct()%>%
  group_by(year) %>%
  mutate(oisst_only_mean = mean(cell.mean),
         oisst_only_max = mean(cell.max)
  ) %>%
  ungroup() %>%
  select(year, oisst_only_mean, oisst_only_max) %>%
  distinct() 

neus_hadisst_all <- read_rds(here("processed-data","neus_hadisst_raw.rds"))$neus_hadisst %>% 
  filter(!is.na(sst)) %>% 
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) 

neus_hadisst_summary <- neus_hadisst_all %>%
  group_by(lat, lon, year) %>%
  mutate(cell.min = min(sst),
         cell.mean = mean(sst),
         cell.max = max(sst)) %>%  
  ungroup() %>%
  select(lat, lon, year, cell.mean, cell.max) %>%
  distinct()%>%
  group_by(year) %>%
  mutate(hadisst_only_mean = mean(cell.mean),
         hadisst_only_max = mean(cell.max)) %>%
  ungroup() %>%
  select(year, hadisst_only_mean, hadisst_only_max) %>%
  distinct() 


neustest <- neus_hadisst_summary %>%
  left_join(neus_oisst_summary, by="year") %>%
  ggplot() +
  geom_line(aes(x=year, y=oisst_only_max), color="magenta") +
  geom_line(aes(x=year, y=hadisst_only_max), color="darkorange") +
  scale_x_continuous(limits = c(1967, 2018), breaks=seq(1968, 2017, 4), expand = c(0, 0)) +
  labs(x="Year", y="Temperature (°C)") +
  theme_bw() +
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.position="none",
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title=element_text(family="sans",size=12,color="black"),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  NULL

# choose how far out into the ocean you want temperature data
ebs.depth.cutoff <- 300
neus.depth.cutoff <- 300
wc.depth.cutoff <- 600 # WC shelf is very steep so I increased this from 300m in 100m increments until the bathymetric mask did not have big gaps along the coast 

# get masks for each region 
wc.bathy <- get.bathy(lon = wc_lonrange, lat = wc_latrange, visualize = F, res = 4) 
neus.bathy <- get.bathy(lon = neus_lonrange, lat = neus_latrange, visualize = F, res = 4) 
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

# filter all temperature dataframes by the mask, and correct the time column 
# can be very slow 

# import the COBE data manually 
# note that it's ginormous because I haven't cropped it at all yet 
cobe.raw <- raster::stack(here("raw-data","sst.mon.mean.nc"))
cobe.tidy <- raster::as.data.frame(cobe.raw, xy=TRUE)  %>% 
  pivot_longer(cols=c(-x, -y), names_to="date", values_to="sst") %>%
  filter(!is.na(sst)) 

cobe.sf <- cobe.tidy %>% 
  mutate(x= ifelse(x >= 180, x-360, x)) %>% # change lon to -180/+180 from +360
  st_as_sf(coords=c("x","y"), crs = bathy.crs) 

# cobe
neus_cobe <- cobe.sf %>%
  st_intersection(neus.bathy.mask)%>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(time = gsub("X","",date)) %>%
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) %>%
  dplyr::select(-date) %>% 
  filter(year>=1967)

wc_cobe <- cobe.sf %>%
  st_intersection(wc.bathy.mask)%>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(time = gsub("X","",date)) %>%
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) %>%
  dplyr::select(-date) %>% 
  filter(year>=1967)

ebs_cobe <- cobe.sf %>%
  st_intersection(ebs.bathy.mask)%>% 
  dplyr::select(-layer) %>%
  sfc_as_cols() %>% 
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  distinct() %>% 
  mutate(time = gsub("X","",date)) %>%
  mutate(
    time = lubridate::ymd(str_sub(time, 1, 10)),
    year = lubridate::year(time),
    month = lubridate::month(time)
  ) %>%
  select(-date) %>% 
  filter(year>=1967)

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
