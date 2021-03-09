# this script uses maps of the US to generate a smoothed coastline for NEUS, WC, and a Middle Domain axis for EBS 
# in each case, a smoother function was applied incrementally until unwanted coastal features (bays, estuaries, etc.) vanished from the smoothed coastline

# https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html

library(here)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(raster)
library(smoothr)
library(ggplot2)
library(gridExtra)

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

wc.ymin <- 34
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


# EBS

# set line endpoints 
ebs.lon <- c(-176.5,-161)
ebs.lat <- c(62,56)

# draw line
ebs.line <- data.frame(ebs.lon, ebs.lat) %>%
  rename('x'=ebs.lon,'y'=ebs.lat) %>%
  st_as_sf(coords=c('x','y')) %>%
  summarise() %>%
  st_cast("LINESTRING") %>% 
  smoothr::densify(n=99) %>% # choose number of line segments: here it's 99 for 100 points
  st_cast("MULTIPOINT") 

ebs.points <- data.frame(st_coordinates(ebs.line)) %>%
  rename("x"=X, "y"=Y)

ebs.dists <- pointDistance(ebs.points[-nrow(ebs.points), c("x", "y")], ebs.points[-1, c("x", "y")], lonlat=TRUE)

ebs.axisdistdat <- data.frame(ebs.points[, c('x','y')], seglength=c(0, ebs.dists))
ebs.axisdistdat$lengthfromhere <- rev(cumsum(rev(ebs.axisdistdat[,"seglength"])))

saveRDS(ebs.axisdistdat, here("processed-data","ebs_axisdistdat.rds"))

# make plots of each region's smoothing for Supplement

usoutline <- rnaturalearth::ne_states("united states of america", returnclass = "sf") %>% 
  st_sf()

# just need one plot for WC
wc.smoothmap <- ggplot() +
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=wcmap %>% 
            smoothr::smooth(method="ksmooth", smoothness=1), color="darkblue", lwd=2) + 
  scale_x_continuous(limits=c(wc.xmin, wc.xmax)) +
  scale_y_continuous(limits=c(wc.ymin, wc.ymax)) +
  theme_bw() +
  labs(title="Smoothness=1") +
  theme(axis.text.x=element_text(angle=45))
wc.smoothmap

# iterate over NEUS smoothness values that we used 
neus_smoothplot <- function(smoothval){
  out <- ggplot() +
    geom_sf(data=usoutline, color="#999999") +
    geom_sf(data=neusmap %>% 
              smoothr::smooth(method="ksmooth", smoothness=smoothval), color="darkblue", lwd=1.5) + 
    scale_x_continuous(limits=c(neus.xmin, neus.xmax)) +
    scale_y_continuous(limits=c(neus.ymin, neus.ymax)) +
    theme_bw() +
    labs(title=paste0("Smoothness=", smoothval)) +
    theme(axis.text.x=element_text(angle=45)) 
}

neus.smoothmap.list <- lapply(seq(1, 8, 1), neus_smoothplot)

neus.smoothmap <- do.call("grid.arrange", c(neus.smoothmap.list, ncol=3))

ggsave(neus.smoothmap, dpi=160, height=7, width=7, file=here("results","smoothing_map_neus.png"))

ggsave(wc.smoothmap, dpi=160, height=7, width=3, filename=here("results","smoothing_map_wc.png"))
#rm(list=ls())

# make an EBS map too

ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

ebs.linemap <- ggplot() + 
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=ebs.line %>% st_set_crs(st_crs(usoutline))) +
  scale_x_continuous(limits=ebs_lonrange) +
  scale_y_continuous(limits=ebs_latrange) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45)) 

ggsave(ebs.linemap, dpi=160, height=7, width=3, filename=here("results","axis_map_ebs.png"))
