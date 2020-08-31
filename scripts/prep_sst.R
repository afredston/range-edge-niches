# this script will get gridded netcdf SST data from NOAA's ERDDAP server, resample the lower-resolution datasets, crop them to the extent of the study regions, and write SST out as dataframes

library(rerddap)
library(raster)
library(sf)
library(tidyverse)
library(here)
library(tabularaster)
library(lubridate)
library(purrr)
library(rnaturalearth)

# explicitly assign functions found in multiple packages
map <- purrr::map
here <- here::here
select <- dplyr::select

source(here("functions","sfc_as_cols.R"))

generate_supplementary_plots = FALSE # toggle this off if you don't want the script to print out a bunch of extra plots exploring anomalies and climatologies of both SST datasets

#####################
### get ERDDAP data
#####################

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
hadisst_fields <- "sst"

neus_hadisst_time <- c("1967-01-16", "2018-12-16")
wc_hadisst_time <- c("1976-01-16","2018-12-16")
ebs_hadisst_time <- c("1981-01-16","2018-12-16") # leaving this in for completeness but we actually don't need HadISST for this region since the data starts in 1989
min_ebs_year_needed <- 1988

oisst <- "ncdcOisst2Agg_LonPM180"
# split up into decades for easier downloading
oisst_time1 <- c("1982-01-01","1989-12-31") # same time interval for all regions
oisst_time2 <- c("1990-01-01","1999-12-31") 
oisst_time3 <- c("2000-01-01","2009-12-31") 
oisst_time4 <- c("2010-01-01","2018-12-31") 

oisst_fields <- "sst"

# get gridded datasets from ERDDAP
neus_hadisst_grid <- griddap(hadisst, time=neus_hadisst_time, latitude = neus_latrange, longitude = neus_lonrange, fields=hadisst_fields)
wc_hadisst_grid <- griddap(hadisst, time=wc_hadisst_time, latitude = wc_latrange, longitude = wc_lonrange, fields=hadisst_fields)
ebs_hadisst_grid <- griddap(hadisst, time=ebs_hadisst_time, latitude = ebs_latrange, longitude = ebs_lonrange, fields=hadisst_fields)

neus_oisst_grid1 <- griddap(oisst, time=oisst_time1, latitude = neus_latrange, longitude = neus_lonrange, fields=oisst_fields)
neus_oisst_grid2 <- griddap(oisst, time=oisst_time2, latitude = neus_latrange, longitude = neus_lonrange, fields=oisst_fields)
neus_oisst_grid3 <- griddap(oisst, time=oisst_time3, latitude = neus_latrange, longitude = neus_lonrange, fields=oisst_fields)
neus_oisst_grid4 <- griddap(oisst, time=oisst_time4, latitude = neus_latrange, longitude = neus_lonrange, fields=oisst_fields)

wc_oisst_grid1 <- griddap(oisst, time=oisst_time1, latitude = wc_latrange, longitude = wc_lonrange, fields=oisst_fields)
wc_oisst_grid2 <- griddap(oisst, time=oisst_time2, latitude = wc_latrange, longitude = wc_lonrange, fields=oisst_fields)
wc_oisst_grid3 <- griddap(oisst, time=oisst_time3, latitude = wc_latrange, longitude = wc_lonrange, fields=oisst_fields)
wc_oisst_grid4 <- griddap(oisst, time=oisst_time4, latitude = wc_latrange, longitude = wc_lonrange, fields=oisst_fields)

ebs_oisst_grid1 <- griddap(oisst, time=oisst_time1, latitude = ebs_latrange, longitude = ebs_lonrange, fields=oisst_fields)
ebs_oisst_grid2 <- griddap(oisst, time=oisst_time2, latitude = ebs_latrange, longitude = ebs_lonrange, fields=oisst_fields)
ebs_oisst_grid3 <- griddap(oisst, time=oisst_time3, latitude = ebs_latrange, longitude = ebs_lonrange, fields=oisst_fields)
ebs_oisst_grid4 <- griddap(oisst, time=oisst_time4, latitude = ebs_latrange, longitude = ebs_lonrange, fields=oisst_fields)

#####################
### resample HadISST to same spatial resolution as OISST
#####################

# find file paths to .nc files
neus_hadisst_nc_file <- neus_hadisst_grid$summary$filename
wc_hadisst_nc_file <- wc_hadisst_grid$summary$filename
ebs_hadisst_nc_file <- ebs_hadisst_grid$summary$filename

neus_oisst_nc_file1 <- neus_oisst_grid1$summary$filename
wc_oisst_nc_file1 <- wc_oisst_grid1$summary$filename
ebs_oisst_nc_file1 <- ebs_oisst_grid1$summary$filename

neus_oisst_nc_file2 <- neus_oisst_grid2$summary$filename
wc_oisst_nc_file2 <- wc_oisst_grid2$summary$filename
ebs_oisst_nc_file2 <- ebs_oisst_grid2$summary$filename

neus_oisst_nc_file3 <- neus_oisst_grid3$summary$filename
wc_oisst_nc_file3 <- wc_oisst_grid3$summary$filename
ebs_oisst_nc_file3 <- ebs_oisst_grid3$summary$filename

neus_oisst_nc_file4 <- neus_oisst_grid4$summary$filename
wc_oisst_nc_file4 <- wc_oisst_grid4$summary$filename
ebs_oisst_nc_file4 <- ebs_oisst_grid4$summary$filename

# read .nc files in as raster bricks

neus_hadisst_brick <- brick(neus_hadisst_nc_file)
wc_hadisst_brick <- brick(wc_hadisst_nc_file)
ebs_hadisst_brick <- brick(ebs_hadisst_nc_file)

neus_oisst_brick1 <- brick(neus_oisst_nc_file1)
neus_oisst_brick2 <- brick(neus_oisst_nc_file2)
neus_oisst_brick3 <- brick(neus_oisst_nc_file3)
neus_oisst_brick4 <- brick(neus_oisst_nc_file4)

wc_oisst_brick1 <- brick(wc_oisst_nc_file1)
wc_oisst_brick2 <- brick(wc_oisst_nc_file2)
wc_oisst_brick3 <- brick(wc_oisst_nc_file3)
wc_oisst_brick4 <- brick(wc_oisst_nc_file4)

ebs_oisst_brick1 <- brick(ebs_oisst_nc_file1)
ebs_oisst_brick2 <- brick(ebs_oisst_nc_file2)
ebs_oisst_brick3 <- brick(ebs_oisst_nc_file3)
ebs_oisst_brick4 <- brick(ebs_oisst_nc_file4)

# resample hadISST to resolution of OISST
neus_hadisst_resample <- resample(neus_hadisst_brick, neus_oisst_brick1, method="ngb") # nearest neighbor method does no interpolation and instead just pastes the value of the nearest point from the coarser raster 
wc_hadisst_resample <- resample(wc_hadisst_brick, wc_oisst_brick1, method="ngb")
ebs_hadisst_resample <- resample(ebs_hadisst_brick, ebs_oisst_brick1, method="ngb")

#####################
### create bathymetric masks for cropping SST data
#####################

# these take hours to generate so don't re-create them unless settings have changed or it's the first time! 
neus.bathy.file <- here("processed-data","neus_bathy_mask.shp")
wc.bathy.file <- here("processed-data","wc_bathy_mask.shp")
ebs.bathy.file <- here("processed-data","ebs_bathy_mask.shp")

# choose how far out into the ocean you want temperature data
ebs.depth.cutoff <- -300
neus.depth.cutoff <- -300
wc.depth.cutoff <- -600 # WC shelf is very steep so I increased this from 300m in 100m increments until the bathymetric mask did not have big gaps along the coast 

# get a bunch of shapefiles to use to crop data 
# US EEZ 
eezs <- st_read(here("raw-data/World_EEZ_v10_20180221","eez_v10.shp")) # download from http://www.marineregions.org/downloads.php and move to raw-data folder
useez <- eezs %>% 
  dplyr::filter(Sovereign1 == "United States") # need to get the new bathy masks in the same CRS

# coarse country outlines  
countries <- ne_countries(type = 'countries', scale = 'small')

# get bathymetric masks; slow 
if(!file.exists(wc.bathy.file)) {
  
  # get masks for each region; same bounding boxes as SST data above. using 4-minute resolution (default)
 # wc.bathy <- get.bathy(lon = wc_lonrange, lat = wc_latrange, visualize = F, res = 4) 
  bathy <- "etopo180" # https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.html 
    
  wc.bathy.grid <- griddap(bathy, latitude=wc_latrange, longitude=wc_lonrange) # this includes land topography in addition to ocean bathymetry 
  
  # get file name
  wc.topo.file <- wc.bathy.grid$summary$filename
  
  # read in as a raster 
  wc.bathy.raster <- raster(wc.topo.file)
  
  # filter out values on land or below depth cutoff 
  wc.bathy.crop <- wc.bathy.raster
  wc.bathy.crop[wc.bathy.crop > 0] <- NA
  wc.bathy.crop[wc.bathy.crop < wc.depth.cutoff] <- NA
  
  # make a bounding box for Puget Sound which is hard to crop out (raster of random values)
  wc.bbox <- st_crop(useez, y=c(xmin=-124, xmax=-121, ymin=46, ymax=49.9)) # crop(useez, extent(-123,-121,46, 49))
  
  # crop to extent of EEZ 
  wc.bathy.mask <- mask(wc.bathy.crop, useez) %>% # KEEP points in the US EEZ
    mask(mask=countries, inverse=TRUE) %>% # DON'T KEEP points in the coarse country map (bays/estuaries/etc)
    mask(mask=wc.bbox, inverse=TRUE) %>% # DON'T KEEP points in the extra bounding box for Puget Sound
    as(., "SpatialPolygonsDataFrame") %>% # make into a polygon (necessary intermediate step)
    st_as_sf() # make into sf object
  
  #plot(wc.bathy.mask)
    
  st_write(wc.bathy.mask, wc.bathy.file) }else {
    wc.bathy.mask <- st_read(wc.bathy.file)
  }

if(!file.exists(ebs.bathy.file)) {
  
  bathy <- "etopo180" # https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.html 
  
  ebs.bathy.grid <- griddap(bathy, latitude=ebs_latrange, longitude=ebs_lonrange) 
  
  ebs.topo.file <- ebs.bathy.grid$summary$filename
  
  ebs.bathy.raster <- raster(ebs.topo.file)
  
  ebs.bathy.crop <- ebs.bathy.raster
  ebs.bathy.crop[ebs.bathy.crop > 0] <- NA
  ebs.bathy.crop[ebs.bathy.crop < ebs.depth.cutoff] <- NA
  
  ebs.bathy.mask <- mask(ebs.bathy.crop, useez) %>%  
    mask(mask=countries, inverse=TRUE) %>%  
    as(., "SpatialPolygonsDataFrame") %>%  
    st_as_sf() 
  
  # plot(ebs.bathy.mask) 
  st_write(ebs.bathy.mask, ebs.bathy.file) }else {
    ebs.bathy.mask <- st_read(ebs.bathy.file)
  }

if(!file.exists(neus.bathy.file)) {
  
  bathy <- "etopo180" # https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.html 
  
  neus.bathy.grid <- griddap(bathy, latitude=neus_latrange, longitude=neus_lonrange) 
  
  neus.topo.file <- neus.bathy.grid$summary$filename
  
  neus.bathy.raster <- raster(neus.topo.file)
  
  neus.bathy.crop <- neus.bathy.raster
  neus.bathy.crop[neus.bathy.crop > 0] <- NA
  neus.bathy.crop[neus.bathy.crop < wc.depth.cutoff] <- NA
  
  neus.bathy.mask <- mask(neus.bathy.crop, useez) %>%  
    mask(mask=countries, inverse=TRUE) %>%  
    as(., "SpatialPolygonsDataFrame") %>%  
    st_as_sf() 
  
  st_write(neus.bathy.mask, neus.bathy.file) } else {
    neus.bathy.mask <- st_read(neus.bathy.file)
  }

#####################
### crop SST datasets to extent of masks 
#####################

neus_hadisst_crop <- mask(neus_hadisst_resample, as_Spatial(neus.bathy.mask)) 
neus_oisst_crop1 <- mask(neus_oisst_brick1, as_Spatial(neus.bathy.mask))
neus_oisst_crop2 <- mask(neus_oisst_brick2, as_Spatial(neus.bathy.mask))
neus_oisst_crop3 <- mask(neus_oisst_brick3, as_Spatial(neus.bathy.mask))
neus_oisst_crop4 <- mask(neus_oisst_brick4, as_Spatial(neus.bathy.mask))

wc_hadisst_crop <- mask(wc_hadisst_resample, as_Spatial(wc.bathy.mask))
# %>%
#   setZ(z=wc_hadisst_times, name="time")
wc_oisst_crop1 <- mask(wc_oisst_brick1, as_Spatial(wc.bathy.mask))
wc_oisst_crop2 <- mask(wc_oisst_brick2, as_Spatial(wc.bathy.mask))
wc_oisst_crop3 <- mask(wc_oisst_brick3, as_Spatial(wc.bathy.mask))
wc_oisst_crop4 <- mask(wc_oisst_brick4, as_Spatial(wc.bathy.mask))

ebs_hadisst_crop <- mask(ebs_hadisst_resample, as_Spatial(ebs.bathy.mask))
# %>%
#   setZ(z=ebs_hadisst_times, name="time")
ebs_oisst_crop1 <- mask(ebs_oisst_brick1, as_Spatial(ebs.bathy.mask))
ebs_oisst_crop2 <- mask(ebs_oisst_brick2, as_Spatial(ebs.bathy.mask))
ebs_oisst_crop3 <- mask(ebs_oisst_brick3, as_Spatial(ebs.bathy.mask))
ebs_oisst_crop4 <- mask(ebs_oisst_brick4, as_Spatial(ebs.bathy.mask))

# HadISST does have dates in the @z dimension, but for some reason they get dropped in resample() (maybe because they're in character format?). let's fix that now, by setting the @z dimension of the final cropped object to the dates from the original file: 
neus_hadisst_crop <- setZ(x=neus_hadisst_crop, z=as_datetime(unlist(neus_hadisst_brick@z)))
wc_hadisst_crop <- setZ(x=wc_hadisst_crop, z=as_datetime(unlist(wc_hadisst_brick@z)))
ebs_hadisst_crop <- setZ(x=ebs_hadisst_crop, z=as_datetime(unlist(ebs_hadisst_brick@z)))

#####################
### convert rasters to tidy dataframes 
#####################

# convert to dataframe, fix column names, and make date into date format 
neus_hadisst_df <- tabularaster::as_tibble(neus_hadisst_crop, cell=FALSE, dim=TRUE, values=TRUE, xy=TRUE) %>% 
  filter(!is.na(cellvalue))%>%
  rename("sst" = cellvalue,
         "time" = dimindex)%>%
  mutate(time = as_date(time),
         year = year(time),
         month=month(time))
wc_hadisst_df <- tabularaster::as_tibble(wc_hadisst_crop, cell=FALSE, dim=TRUE, values=TRUE, xy=TRUE) %>% 
  filter(!is.na(cellvalue))%>%
  rename("sst" = cellvalue,
         "time" = dimindex)%>%
  mutate(time = as_date(time),
         year = year(time),
         month=month(time))
ebs_hadisst_df <- tabularaster::as_tibble(ebs_hadisst_crop, cell=FALSE, dim=TRUE, values=TRUE, xy=TRUE) %>% 
  filter(!is.na(cellvalue),
         cellvalue > -999)%>% # filter out weird cells with a value of -1000
  rename("sst" = cellvalue,
         "time" = dimindex) %>%
  mutate(time = as_date(time),
         year = year(time),
         month=month(time))

# fewer issues with OISST dates, can just convert them to date format here 

# convert all the OISST bricks to data frames 

# make function that takes cropped OISST raster brick and converts it to a dataframe using tabularaster, with some tidying 
oisst_to_df <- function(oisst){
  out <- tabularaster::as_tibble(oisst, cell=FALSE, dim=TRUE, values=TRUE, xy=TRUE) %>% 
    filter(!is.na(cellvalue)) %>%
    rename("sst" = cellvalue,
           "time" = dimindex) %>%
    mutate(time = as_datetime(as.integer(time)),
           time = as_date(time),
           year = year(time),
           month = month(time)) 
  return(out)
}

# get a list of all the rasters 
oisst_crop_list <- c(neus_oisst_crop1,neus_oisst_crop2,neus_oisst_crop3,neus_oisst_crop4,wc_oisst_crop1,wc_oisst_crop2,wc_oisst_crop3,wc_oisst_crop4,ebs_oisst_crop1,ebs_oisst_crop2,ebs_oisst_crop3,ebs_oisst_crop4)

# apply oisst_to_df to all rasters
oisst_df_list <- map(oisst_crop_list, oisst_to_df)

# give them sensible names 
names(oisst_df_list) <- c("neus_oisst_df1","neus_oisst_df2","neus_oisst_df3","neus_oisst_df4","wc_oisst_df1","wc_oisst_df2","wc_oisst_df3","wc_oisst_df4","ebs_oisst_df1","ebs_oisst_df2","ebs_oisst_df3","ebs_oisst_df4")

# unlist dfs into environment 
list2env(oisst_df_list, .GlobalEnv)

# make summary OISST dfs
neus_oisst_df <- bind_rows(neus_oisst_df1, neus_oisst_df2, neus_oisst_df3, neus_oisst_df4)
wc_oisst_df <- bind_rows(wc_oisst_df1, wc_oisst_df2, wc_oisst_df3, wc_oisst_df4)
ebs_oisst_df <- bind_rows(ebs_oisst_df1, ebs_oisst_df2, ebs_oisst_df3, ebs_oisst_df4)

#####################
### confirm grids are identical
#####################

# check that both datasets for a region contain the exact same cells 
# as shown in the plots here, there are minor differences between the two grids in grid cells adjacent to land. I'm just dropping these points for now, before calcluating any climatology
# don't need to worry about this for EBS because we aren't actually using hadISST

neus_hadisst_coords <- neus_hadisst_df %>%
  select(x, y) %>%
  distinct() %>%
  mutate(coords = paste0(x, ",", y))

neus_oisst_coords <- neus_oisst_df %>%
  select(x, y) %>%
  distinct() %>%
  mutate(coords = paste0(x, ",", y))

setdiff(neus_hadisst_coords$coords, neus_oisst_coords$coords)# check for differences
setdiff(neus_oisst_coords$coords, neus_hadisst_coords$coords)

# plot neus differences
usoutline <- rnaturalearth::ne_states("united states of america", returnclass = "sf") %>% 
  st_sf()

ggplot() + 
  geom_sf(data=usoutline, color="#999999") +
  geom_point(data=neus_hadisst_coords %>% filter(coords %in% setdiff(neus_hadisst_coords$coords, neus_oisst_coords$coords)), aes(x=x, y=y), color="pink") +
  geom_point(data=neus_oisst_coords %>% filter(coords %in% setdiff(neus_oisst_coords$coords, neus_hadisst_coords$coords)), aes(x=x, y=y), color="purple") +
  scale_x_continuous(limits=neus_lonrange) +
  scale_y_continuous(limits=neus_latrange)

wc_hadisst_coords <- wc_hadisst_df %>%
  select(x, y) %>%
  distinct() %>%
  mutate(coords = paste0(x, ",", y))

wc_oisst_coords <- wc_oisst_df %>%
  select(x, y) %>%
  distinct() %>%
  mutate(coords = paste0(x, ",", y))

setdiff(wc_hadisst_coords$coords, wc_oisst_coords$coords)
setdiff(wc_oisst_coords$coords, wc_hadisst_coords$coords)

ggplot() + 
  geom_sf(data=usoutline, color="#999999") +
  geom_point(data=wc_hadisst_coords %>% filter(coords %in% setdiff(wc_hadisst_coords$coords, wc_oisst_coords$coords)), aes(x=x, y=y), color="pink") +
  geom_point(data=wc_oisst_coords %>% filter(coords %in% setdiff(wc_oisst_coords$coords, wc_hadisst_coords$coords)), aes(x=x, y=y), color="purple") +
  scale_x_continuous(limits=wc_lonrange) +
  scale_y_continuous(limits=wc_latrange)

#####################
### calculate climatologies 
#####################

# in order to join these datasets, we need to do a mean bias correction, and also aggregate OISST up to monthly resolution (from daily) 

neus_hadisst_df_clim <- neus_hadisst_df %>%
  mutate(coords = paste0(x, ",", y)) %>%
  filter(!coords %in% setdiff(neus_hadisst_coords$coords, neus_oisst_coords$coords)) %>% # get rid of grid cells that are missing from OISST
  select(-coords) %>% 
  group_by(x, y, month) %>%
  mutate(sst_month_clim = mean(sst[time >= min(neus_oisst_df$time)])) %>% # get climatology by month (average conditions across all instances of that month in that cell) from the same time span used for the OISST climatology
  ungroup() %>%
  rename("date_join"=time) %>%
  mutate(sst_month_anom = sst-sst_month_clim, # get anomaly by month (how much hotter/colder that month is relative to the average )
         dataset = "HadISST",
         sst_month_clim = NA) %>% # going to use the OISST climatology 
  select(-sst)

wc_hadisst_df_clim <- wc_hadisst_df %>%
  mutate(coords = paste0(x, ",", y)) %>%
  filter(!coords %in% setdiff(wc_hadisst_coords$coords, wc_oisst_coords$coords)) %>%
  select(-coords) %>%
  group_by(x, y, month) %>%
  mutate(sst_month_clim = mean(sst[time >= min(wc_oisst_df$time)])) %>% 
  ungroup() %>%
  rename("date_join"=time) %>%
  mutate(sst_month_anom = sst-sst_month_clim,  
         dataset = "HadISST",
         sst_month_clim = NA) %>%  
  select(-sst)

# ebs_hadisst_df_clim <- ebs_hadisst_df %>%
#   mutate(coords = paste0(x, ",", y)) %>%
#   filter(!coords %in% setdiff(ebs_hadisst_coords$coords, ebs_oisst_coords$coords)) %>%
#   select(-coords) %>%
#   group_by(x, y, month) %>%
#   mutate(sst_month_clim = mean(sst)) %>% 
#   ungroup() %>%
#   rename("date_join"=time) %>%
#   mutate(sst_month_anom = sst-sst_month_clim,  
#          dataset = "HadISST",
#          sst_month_clim = NA) %>%  
#   select(-sst)

neus_oisst_df_clim <- neus_oisst_df %>% 
  mutate(coords = paste0(x, ",", y)) %>%
  filter(!coords %in% setdiff(neus_oisst_coords$coords, neus_hadisst_coords$coords)) %>% # get rid of grid cells missing from hadISST
  select(-coords) %>%
  group_by(x, y, year, month)  %>%
  summarise(sst_month = median(sst)) %>% # convert into monthly medians for comparability with HadISST before doing anything else; using median instead of mean because it is less sensitive to outliers
  group_by(x, y, month) %>%
  mutate(sst_month_clim = mean(sst_month)) %>%
  ungroup() %>%
  mutate(sst_month_anom = sst_month-sst_month_clim,
         date_join = as_date(paste0(year, "-", month, "-16")),
         dataset="OISST") %>% # prepare for joining to hadISST 
  select(-sst_month)

wc_oisst_df_clim <- wc_oisst_df %>% 
  mutate(coords = paste0(x, ",", y)) %>%
  filter(!coords %in% setdiff(wc_oisst_coords$coords, wc_hadisst_coords$coords)) %>%
  select(-coords) %>%
  group_by(x, y, year, month)  %>%
  summarise(sst_month = median(sst)) %>% 
  group_by(x, y, month) %>%
  mutate(sst_month_clim = mean(sst_month)) %>%
  ungroup() %>%
  mutate(sst_month_anom = sst_month-sst_month_clim,
         date_join = as_date(paste0(year, "-", month, "-16")),
         dataset="OISST") %>%  
  select(-sst_month)

#####################
### merge datasets and write out 
#####################

neus_df_clim <- neus_hadisst_df_clim %>%
  filter(date_join < min(neus_oisst_df_clim$date_join)) %>% # keep only early dates 
  bind_rows(neus_oisst_df_clim) %>%
  group_by(x, y, month) %>%
  arrange(year) %>%
  tidyr::fill(sst_month_clim, .direction="up") %>% # populate NA climatologies in earlier years (hadISST) with OISST climatologies (CELL-SPECIFIC)
  ungroup() %>%
  mutate(sst = sst_month_clim + sst_month_anom) 

wc_df_clim <- wc_hadisst_df_clim %>%
  filter(date_join < min(wc_oisst_df_clim$date_join)) %>% # keep only early dates 
  bind_rows(wc_oisst_df_clim) %>%
  group_by(x, y, month) %>%
  arrange(year) %>%
  tidyr::fill(sst_month_clim, .direction="up") %>% 
  ungroup() %>%
  mutate(sst = sst_month_clim + sst_month_anom) 

# EBS is treated differently because we don't need hadISST years 
# some of this script is redundant (don't actually need to calculate climatology + anomaly here because we aren't combining datasets) but keeping it so columns/dfs are comparable to other regions
ebs_df_clim <- ebs_oisst_df %>%
  filter(year >= min_ebs_year_needed) %>% # DROP EARLY YEARS OF OISST 
  group_by(x, y, year, month)  %>%
  summarise(sst_month = median(sst)) %>% 
  group_by(x, y, month) %>%
  mutate(sst_month_clim = mean(sst_month)) %>%
  ungroup() %>%
  mutate(sst_month_anom = sst_month-sst_month_clim,
         date_join = as_date(paste0(year, "-", month, "-16")),
         dataset="OISST") %>%  
  select(-sst_month) %>%
  mutate(sst = sst_month_clim + sst_month_anom) 

# check there are no NAs
neus_df_clim  %>% filter(is.na(sst))
wc_df_clim  %>% filter(is.na(sst))
ebs_df_clim  %>% filter(is.na(sst))

# save all dfs
saveRDS(neus_df_clim, here("processed-data","neus_sst_corrected.rds"))
saveRDS(wc_df_clim, here("processed-data","wc_sst_corrected.rds"))
saveRDS(ebs_df_clim, here("processed-data","ebs_sst_corrected.rds"))

#####################
### generate supplementary plots
#####################

if(generate_supplementary_plots==TRUE) {
  library(RColorBrewer)
  
  # monthly climatology of each dataset over time, for 1982 onwards; need this to make a plot of climatologies and anomalies over time (since above, we filter out HadISST years after 1982)
  neus.reg.clim.hadisst <- neus_hadisst_df %>%
    filter(year>=1982) %>%
    group_by(x, y, month) %>%
    summarise(sst_month_clim = mean(sst)) %>% # get climatology by month (average conditions across all instances of that month in that cell) 
    ungroup() %>%
    group_by(month) %>%
    summarise(reg_month_clim = mean(sst_month_clim)) %>%
    mutate(dataset="HadISST")
  
  wc.reg.clim.hadisst <- wc_hadisst_df %>%
    filter(year>=1982) %>%
    group_by(x, y, month) %>%
    summarise(sst_month_clim = mean(sst)) %>% 
    ungroup() %>%
    group_by(month) %>%
    summarise(reg_month_clim = mean(sst_month_clim)) %>%
    mutate(dataset="HadISST")
  
  # plot monthly climatologies over time from both datasets 
  neus.clim.time.gg <- neus_oisst_df_clim %>%
    group_by(month) %>%
    summarise(reg_month_clim = mean(sst_month_clim)) %>%
    mutate(dataset="OISST") %>%
    bind_rows(neus.reg.clim.hadisst) %>%
    ggplot() +
    geom_line(aes(x=month, y=reg_month_clim, group=dataset, color=dataset), size=0.8) +
    theme_bw() + 
    scale_color_manual(values=c("#41C03F","mediumblue")) + 
    scale_x_continuous(breaks=seq(1, 12, 1)) +
    scale_y_continuous(breaks=seq(6, 22, 2)) +
    labs(x="Month", y="Climatological Mean SST (°C)", title="Northeast") +
    theme(legend.title=element_blank(),
          legend.position=c(0.2, 0.8))
  neus.clim.time.gg
  
  wc.clim.time.gg <- wc_oisst_df_clim %>%
    group_by(month) %>%
    summarise(reg_month_clim = mean(sst_month_clim)) %>%
    mutate(dataset="OISST") %>%
    bind_rows(wc.reg.clim.hadisst) %>%
    ggplot() +
    geom_line(aes(x=month, y=reg_month_clim, group=dataset, color=dataset), size=0.8) +
    theme_bw() + 
    scale_color_manual(values=c("#41C03F","mediumblue")) + 
    scale_x_continuous(breaks=seq(1, 12, 1)) +
    scale_y_continuous(breaks=seq(11, 17, 1)) +
    labs(x="Month", y="Climatological Mean SST (°C)", title="West Coast") +
    theme(legend.title=element_blank(),
          legend.position=c(0.2, 0.8))
  wc.clim.time.gg
  
  ggsave(neus.clim.time.gg, filename=here("results","sst_climatologies_time_neus.png"), height=5, width=4, dpi=160)
  ggsave(wc.clim.time.gg, filename=here("results","sst_climatologies_time_wc.png"), height=5, width=4, dpi=160)
  
  # order months by mean temperature for plotting colors
  neus.month.labels <- neus_oisst_df_clim %>%
    group_by(month) %>%
    summarise(overall_clim = mean(sst_month_clim)) %>%
    arrange(-overall_clim) %>%
    mutate(month = as.factor(as.character(month)))
  
  wc.month.labels <- wc_oisst_df_clim %>%
    group_by(month) %>%
    summarise(overall_clim = mean(sst_month_clim)) %>%
    arrange(-overall_clim) %>%
    mutate(month = as.factor(as.character(month)))
  
  # hacky way to expand the RColorBrewer red-to-blue palette to 12 distinct colors 
  month.palette <- colorRampPalette(colors = brewer.pal(10, "RdYlBu"))(12)
  
  # monthly anomalies in both datasets 
  neus.anom.time.gg <- neus_hadisst_df_clim %>% 
    select(-sst_month_clim) %>%
    bind_rows(neus_oisst_df_clim %>% select(-sst_month_clim)) %>% # get anomalies only from full time-series of both datasets 
    group_by(date_join, month, dataset) %>%
    summarise(reg_month_anom = mean(sst_month_anom)) %>% # get mean across all cells in the region for that date & that dataset
    ungroup() %>%
    mutate(month = factor(as.character(month), levels=neus.month.labels$month)) %>% # reorder months by mean climatology for visual effect in plot
    ggplot(aes(x=date_join, y=reg_month_anom, color=month, group=month)) +
    labs(x="Time", y="SST Anomaly (°C)", title="Northeast", color="Month") +
    geom_line(size=0.8) +
    geom_hline(yintercept=0, color="black", linetype="dashed", size=0.8) +
    facet_wrap(~dataset) +
    theme_bw() + 
    scale_color_manual(values=month.palette) + 
    scale_x_date(breaks = lubridate::ymd(seq(1968, 2018, 10), truncated = 2L), date_labels = "%Y") +
    scale_y_continuous(breaks=seq(-3, 3, 1)) +
    theme( 
      legend.position=c(0.15, 0.88),
      legend.direction = "horizontal",
      legend.background = element_rect("transparent"))
  neus.anom.time.gg
  ggsave(neus.anom.time.gg, filename=here("results","sst_anomalies_time_neus.png"), height=5, width=10, dpi=160)
  
  wc.anom.time.gg <- wc_hadisst_df_clim %>% 
    select(-sst_month_clim) %>%
    bind_rows(wc_oisst_df_clim %>% select(-sst_month_clim)) %>%  
    group_by(date_join, month, dataset) %>%
    summarise(reg_month_anom = mean(sst_month_anom)) %>%  
    ungroup() %>%
    mutate(month = factor(as.character(month), levels=neus.month.labels$month)) %>% 
    ggplot(aes(x=date_join, y=reg_month_anom, color=month, group=month)) +
    labs(x="Time", y="SST Anomaly (°C)", title="West Coast", color="Month") +
    geom_line(size=0.8) +
    geom_hline(yintercept=0, color="black", linetype="dashed", size=0.8) +
    facet_wrap(~dataset) +
    theme_bw() + 
    scale_color_manual(values=month.palette) + 
    scale_x_date(breaks = lubridate::ymd(seq(1978, 2018, 10), truncated = 2L), date_labels = "%Y") +
    scale_y_continuous(breaks=seq(-3, 3, 1)) +
    theme( 
      #  legend.position=c(0.15, 0.88),
      #   legend.direction = "horizontal",
      #  legend.background = element_rect("transparent")
      
      legend.position="none") # get rid of legend which plots over a peak--will be next to NEUS plot anyway
  wc.anom.time.gg
  ggsave(wc.anom.time.gg, filename=here("results","sst_anomalies_time_wc.png"), height=5, width=10, dpi=160)
  
  # density plots of all cell-specific anomalies in each dataset
  neus.anom.pdf.gg <- neus_hadisst_df_clim %>% 
    select(-sst_month_clim) %>%
    bind_rows(neus_oisst_df_clim %>% select(-sst_month_clim)) %>% # get anomalies only from full time-series of both datasets 
    ggplot() +
    geom_density(aes(x=sst_month_anom, group=dataset, fill=dataset), color="black", alpha=0.4) +
    scale_fill_manual(values=c("#0D592C","mediumblue")) + 
    labs(x="SST Anomaly (°C)", y="Density", title="Northeast", fill="Dataset") +
    theme_bw() + 
    coord_cartesian(xlim=c(-5, 5), ylim=c(0, 0.5)) +
    theme( 
      legend.position=c(0.2, 0.88),
      legend.direction = "vertical",
      legend.background = element_rect("transparent"))
  neus.anom.pdf.gg
  ggsave(neus.anom.pdf.gg, filename=here("results","sst_anomalies_density_neus.png"), height=5, width=5, dpi=160)
  
  wc.anom.pdf.gg <- wc_hadisst_df_clim %>% 
    select(-sst_month_clim) %>%
    bind_rows(wc_oisst_df_clim %>% select(-sst_month_clim)) %>% # get anomalies only from full time-series of both datasets 
    ggplot() +
    geom_density(aes(x=sst_month_anom, group=dataset, fill=dataset), color="black", alpha=0.4) +
    scale_fill_manual(values=c("#0D592C","mediumblue")) + 
    labs(x="SST Anomaly (°C)", y="Density", title="West Coast", fill="Dataset") +
    theme_bw() + 
    coord_cartesian(xlim=c(-5, 5), ylim=c(0, 0.5)) +
    theme( 
      legend.position=c(0.2, 0.88),
      legend.direction = "vertical",
      legend.background = element_rect("transparent"))
  wc.anom.pdf.gg
  ggsave(wc.anom.pdf.gg, filename=here("results","sst_anomalies_density_wc.png"), height=5, width=5, dpi=160)
  
  # make maps of anomalies to demonstrate they are randomly distributed in space 
  
  map_dates <- as_date(c("1985-10-16", "1995-07-16", "2005-04-16","2015-01-16"))
  
  neus.hadisst.anom.map <- neus_hadisst_df_clim %>% 
    filter(date_join %in% map_dates) %>%
    ggplot() + 
    geom_sf(data=usoutline, color="#999999") +
    geom_point(aes(x=x, y=y, fill=sst_month_anom, color=sst_month_anom)) +
    scale_x_continuous(limits=neus_lonrange) +
    scale_y_continuous(limits=neus_latrange) +
    labs(x=element_blank(), y=element_blank(), title="Northeast HadISST", fill="Anomaly (°C)", color="Anomaly (°C)") +
    theme_bw() +
    ggplot2::scale_color_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    ggplot2::scale_fill_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    theme(legend.position = "right",
          axis.text.x = element_text(angle=45)) +
    facet_wrap(~date_join, ncol=4) +
    NULL
  
  neus.oisst.anom.map <- neus_oisst_df_clim %>% 
    filter(date_join %in% map_dates) %>%
    ggplot() + 
    geom_sf(data=usoutline, color="#999999") +
    geom_point(aes(x=x, y=y, fill=sst_month_anom, color=sst_month_anom)) +
    scale_x_continuous(limits=neus_lonrange) +
    scale_y_continuous(limits=neus_latrange) +
    labs(x=element_blank(), y=element_blank(), title="Northeast OISST", fill="Anomaly (°C)", color="Anomaly (°C)") +
    theme_bw() +
    ggplot2::scale_color_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    ggplot2::scale_fill_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    theme(legend.position = "right",
          axis.text.x = element_text(angle=45)) +
    facet_wrap(~date_join, ncol=4) +
    NULL
  
  wc.hadisst.anom.map <- wc_hadisst_df_clim %>% 
    filter(date_join %in% map_dates) %>%
    ggplot() + 
    geom_sf(data=usoutline, color="#999999") +
    geom_point(aes(x=x, y=y, fill=sst_month_anom, color=sst_month_anom)) +
    scale_x_continuous(limits=wc_lonrange) +
    scale_y_continuous(limits=wc_latrange) +
    labs(x=element_blank(), y=element_blank(), title="West Coast HadISST", fill="Anomaly (°C)", color="Anomaly (°C)") +
    theme_bw() +
    ggplot2::scale_color_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    ggplot2::scale_fill_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle=45)) +
    facet_wrap(~date_join, ncol=4) +
    NULL
  
  
  wc.oisst.anom.map <- wc_oisst_df_clim %>% 
    filter(date_join %in% map_dates) %>%
    ggplot() + 
    geom_sf(data=usoutline, color="#999999") +
    geom_point(aes(x=x, y=y, fill=sst_month_anom, color=sst_month_anom)) +
    scale_x_continuous(limits=wc_lonrange) +
    scale_y_continuous(limits=wc_latrange) +
    labs(x=element_blank(), y=element_blank(), title="West Coast OISST", fill="Anomaly (°C)", color="Anomaly (°C)") +
    theme_bw() +
    ggplot2::scale_color_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    ggplot2::scale_fill_gradient2(low="#042C5B", mid="white", high="#5B0404") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle=45)) +
    facet_wrap(~date_join, ncol=4) +
    NULL
  
  ggsave(neus.hadisst.anom.map, dpi=160, filename=here("results","sst_anomalies_map_neus_hadisst.png"), width=6.5, height=2.5, scale=1.8)
  ggsave(neus.oisst.anom.map, dpi=160, filename=here("results","sst_anomalies_map_neus_oisst.png"), width=6.5, height=2.5, scale=1.8)
  ggsave(wc.hadisst.anom.map, dpi=160, filename=here("results","sst_anomalies_map_wc_hadisst.png"), width=6.5, height=3.5, scale=1.5)
  ggsave(wc.oisst.anom.map, dpi=160, filename=here("results","sst_anomalies_map_wc_oisst.png"), width=6.5, height=3.5, scale=1.5)
}
