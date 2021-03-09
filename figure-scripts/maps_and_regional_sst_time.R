# make plot of the three regions and temperature trends 
library(tidyverse)
library(here)
library(sf)
library(purrr)
library(broom)
library(knitr)
here <- here::here

# load all temperature datasets 
neus_sst <- read_rds(here("processed-data","neus_sst_corrected.rds"))
wc_sst <- read_rds(here("processed-data","wc_sst_corrected.rds"))
ebs_sst <- read_rds(here("processed-data","ebs_sst_corrected.rds"))

# note that temperature is calculated by: 1. getting monthly means for each cell 2. getting the means and extremes for a year for each cell (hottest, coldest, and average months) 3. calculating the average hottest, coldest, and average months across all cells in the region every year. this is designed to capture extreme temperature variation without being too sensitive to daily or cell-specific variation. in other words, it is measuring if hot months got hotter, cold months got less cold, or average months are creeping up in SST *on average across the region*. 

neus_sst_summary <- neus_sst %>%
  group_by(x, y, year) %>%
  mutate(cell.min = min(sst),
         cell.mean = mean(sst),
         cell.max = max(sst)) %>%
  ungroup() %>%
  select(x, y, year, cell.min, cell.mean, cell.max) %>%
  distinct()%>%
  group_by(year) %>%
  mutate(region_mean_min_sst = mean(cell.min),
         region_mean_mean_sst = mean(cell.mean),
         region_mean_max_sst = mean(cell.max)) %>%
  ungroup() %>%
  select(year, region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst) %>%
  distinct() %>%
  mutate(region="neus")

wc_sst_summary <- wc_sst %>%
  group_by(x, y, year) %>%
  mutate(cell.min = min(sst),
         cell.mean = mean(sst),
         cell.max = max(sst)) %>%
  ungroup() %>%
  select(x, y, year, cell.min, cell.mean, cell.max) %>%
  distinct()%>%
  group_by(year) %>%
  mutate(region_mean_min_sst = mean(cell.min),
         region_mean_mean_sst = mean(cell.mean),
         region_mean_max_sst = mean(cell.max)) %>%
  ungroup() %>%
  select(year, region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst) %>%
  distinct() %>%
  mutate(region="wc")

ebs_sst_summary <- ebs_sst %>%
  group_by(x, y, year) %>%
  mutate(cell.min = min(sst),
         cell.mean = mean(sst),
         cell.max = max(sst)) %>%
  ungroup() %>%
  select(x, y, year, cell.min, cell.mean, cell.max) %>%
  distinct()%>%
  group_by(year) %>%
  mutate(region_mean_min_sst = mean(cell.min),
         region_mean_mean_sst = mean(cell.mean),
         region_mean_max_sst = mean(cell.max)) %>%
  ungroup() %>%
  select(year, region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst) %>%
  distinct() %>%
  mutate(region="ebs")

neusgg <- neus_sst_summary %>%
  pivot_longer(cols=c(region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst), names_to = "sstvar", values_to = "sstvalue") %>%
  ggplot(aes(x=year, y=sstvalue, group=sstvar, color=sstvar)) +
  geom_line(size=1.1) +
  scale_color_manual(values=c("#DF2301","#FF8E02","#3A4ED0")) +
  scale_x_continuous(limits = c(1967, 2018), breaks=seq(1968, 2017, 4), expand = c(0, 0)) +
  labs(x="Year", y="SST (°C)") +
  theme_bw() +
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.position="none",
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title=element_text(family="sans",size=12,color="black"),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank()
        ) + 
  NULL

wcgg <- wc_sst_summary %>%
  pivot_longer(cols=c(region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst), names_to = "sstvar", values_to = "sstvalue") %>%
  ggplot(aes(x=year, y=sstvalue, group=sstvar, color=sstvar)) +
  geom_line(size=1.1) +
  scale_color_manual(values=c("#DF2301","#FF8E02","#3A4ED0")) +
  scale_x_continuous(limits = c(1976, 2018), breaks=seq(1976, 2017, 4), expand = c(0, 0)) +
  labs(x="Year", y="SST (°C)") +
  theme_bw() +
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.position="none",
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title=element_text(family="sans",size=12,color="black"),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank()
  ) + 
  NULL

ebsgg <- ebs_sst_summary %>%
  pivot_longer(cols=c(region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst), names_to = "sstvar", values_to = "sstvalue") %>%
  ggplot(aes(x=year, y=sstvalue, group=sstvar, color=sstvar)) +
  geom_line(size=1.1) +
  scale_color_manual(values=c("#DF2301","#FF8E02","#3A4ED0")) +
  scale_x_continuous(limits = c(1988, 2018), breaks=seq(1988, 2017, 4), expand = c(0, 0)) +
  labs(x="Year", y="SST (°C)") +
  theme_bw() +
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.position="none",
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title=element_text(family="sans",size=12,color="black"),
        strip.text.x = element_blank(),
        panel.grid.minor = element_blank()
    ) + 
  NULL

ggsave(neusgg, filename=here("results","neus_sst_inset_plot.png"), height=25, width=31, units="mm", scale=1.7, dpi=600)
ggsave(wcgg, filename=here("results","wc_sst_inset_plot.png"), height=25, width=31, scale=1.7,units="mm",dpi=600)
ggsave(ebsgg, filename=here("results","ebs_sst_inset_plot.png"), height=25, width=31, scale=1.7,units="mm", dpi=600)

# get rid of whitespace
plot_crop(here("results","neus_sst_inset_plot.png"))
plot_crop(here("results","wc_sst_inset_plot.png"))
plot_crop(here("results","ebs_sst_inset_plot.png"))

# make maps of each region using the bathymetric masks generated in prep_sst.R 

neus_bathy <- st_read(here("processed-data","neus_bathy_mask.shp"))
wc_bathy <- st_read(here("processed-data","wc_bathy_mask.shp"))
ebs_bathy <- st_read(here("processed-data","ebs_bathy_mask.shp"))
neus_coast <- st_read(here("processed-data","neus_smoothed_coastline.shp"))
wc_coast <- st_read(here("processed-data","wc_smoothed_coastline.shp"))

usoutline <- rnaturalearth::ne_states("united states of america", returnclass = "sf") %>% 
  st_sf()

# set line endpoints (see get_axes_of_measurement.R)
ebs.lon <- c(-176.5,-161)
ebs.lat <- c(62,56)

# draw line
ebs.line <- data.frame(ebs.lon, ebs.lat) %>%
  rename('x'=ebs.lon,'y'=ebs.lat) %>%
  st_as_sf(coords=c('x','y')) %>%
  summarise() %>%
  st_cast("LINESTRING") %>%
  st_set_crs(st_crs(ebs_bathy))

# set bounding boxes
neus_latrange <- c(34, 46)
neus_lonrange <- c(-78, -66) 
wc_latrange <- c(31, 50)
wc_lonrange <- c(-126, -116)
ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

# make objects to plot as origin points 
neus_coastdistdat <- readRDS(here("processed-data","neus_coastdistdat.rds"))
wc_coastdistdat <- readRDS(here("processed-data","wc_coastdistdat.rds"))

neus.origin <- neus_coastdistdat %>% 
  filter(lengthfromhere == min(lengthfromhere))%>%
  st_as_sf(coords=c('x','y')) %>%
  st_cast("POINT") %>%
  st_set_crs(st_crs(neus_bathy))

wc.origin <- wc_coastdistdat %>% 
  filter(lengthfromhere == min(lengthfromhere))%>%
  st_as_sf(coords=c('x','y')) %>%
  st_cast("POINT") %>%
  st_set_crs(st_crs(wc_bathy))

ebs.origin <- data.frame(ebs.lon[2], ebs.lat[2]) %>%
  rename('x'=ebs.lon.2.,'y'=ebs.lat.2.) %>%
  st_as_sf(coords=c('x','y')) %>%
  summarise() %>%
  st_cast("POINT") %>%
  st_set_crs(st_crs(ebs_bathy))

# make waypoints
neus.coastdistrefs <- neus_coastdistdat %>% 
  filter(seglength > 0,
         seglength < 100000) %>% # get rid of weird points
  mutate(coastdist = lengthfromhere/1000,
         coastdistround = round(coastdist, digits = -2), # round to nearest hundred 
         coastdistdiff = abs(coastdist-coastdistround)) %>% 
  group_by(coastdistround) %>% 
  filter(coastdistdiff == min(coastdistdiff), # get points closest to rounded waypoint value
         coastdistround <= 1400) %>% 
  ungroup() %>% 
  dplyr::select(x, y, coastdistround)

wc.coastdistrefs <- wc_coastdistdat %>% 
  mutate(coastdist = lengthfromhere/1000,
         coastdistround = round(coastdist, digits = -2), # round to nearest hundred 
         coastdistdiff = abs(coastdist-coastdistround)) %>% 
  group_by(coastdistround) %>% 
  filter(coastdistdiff == min(coastdistdiff), 
         coastdistround > 0) %>% # don't need to label the origin
  ungroup() %>% 
  dplyr::select(x, y, coastdistround) %>% 
  filter(coastdistround %in% seq(200, 2000, 200)) # pare down to every 200 for plotting 

ebs.linedistrefs <- readRDS(here("processed-data","ebs_axisdistdat.rds"))%>% 
  mutate(linedist = lengthfromhere/1000,
         linedistround = round(linedist, digits = -2), # round to nearest hundred 
         linedistdiff = abs(linedist-linedistround)) %>% 
  group_by(linedistround) %>% 
  filter(linedistdiff == min(linedistdiff), 
         linedistround > 0) %>% # don't need to label the origin
  ungroup() %>% 
  dplyr::select(x, y, linedistround)

# make maps 
neusmap <- ggplot() + 
  geom_sf(data=neus_bathy, color="skyblue3") +
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=neus_coast, color="black", linetype="dashed", lwd=1.5) +
  geom_sf(data=neus.origin, color="black", fill="transparent",shape=4, size=4, stroke=4)+
  geom_point(data=neus.coastdistrefs, aes(x=x, y=y), color="white") +
  geom_text(data=neus.coastdistrefs, aes(x=x, y=y, label=coastdistround),hjust=0, nudge_x = 0.5, fontface="bold", size=3) +
  scale_x_continuous(limits = c(neus_lonrange[1], neus_lonrange[2]+1), expand = c(0, 0)) + # add some extra space for labels
  scale_y_continuous(limits = neus_latrange, expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none",
        text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_blank()) +
  NULL

wcmap <- ggplot() + 
  geom_sf(data=wc_bathy, color="skyblue3") +
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=wc_coast, color="black", linetype="dashed", lwd=1.5) +
  geom_sf(data=wc.origin, color="black", fill="transparent",shape=4, size=4, stroke=4)+
  geom_point(data=wc.coastdistrefs, aes(x=x, y=y), color="white") +
  geom_text(data=wc.coastdistrefs, aes(x=x, y=y, label=coastdistround),hjust=0, nudge_x = 0.5, vjust=-0.5, fontface="bold", size=3) + 
  scale_x_continuous(limits = c(wc_lonrange[1], wc_lonrange[2]+1), expand = c(0, 0)) + # add extra space for labels
  scale_y_continuous(limits = wc_latrange, expand = c(0, 0), breaks = seq(32, 50, 2)) +
  theme_bw() +
  theme(legend.position = "none",
        text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_blank()) +
  NULL

ebsmap <- ggplot() + 
  geom_sf(data=ebs_bathy, color="skyblue3") +
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=ebs.line, color="black", linetype="dashed", lwd=1.5) +
  geom_sf(data=ebs.origin, color="black", fill="transparent",shape=4, size=4, stroke=4)+
  geom_point(data=ebs.linedistrefs, aes(x=x, y=y), color="white") +
  geom_text(data=ebs.linedistrefs, aes(x=x, y=y, label=linedistround), nudge_x = -1.85, fontface="bold", size=3) + 
  scale_x_continuous(limits = ebs_lonrange, expand = c(0, 0), breaks=seq(-180, -154, 2)) +
  scale_y_continuous(limits = ebs_latrange, expand = c(0, 0), breaks=seq(54, 66, 2)) +
  theme_bw() +
  theme(legend.position = "none",
        text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.title=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)) +
  NULL

ggsave(neusmap, filename=here("results","region_map_neus.png"), width=35, units="mm", dpi=600, scale=1.6)
ggsave(wcmap, filename=here("results","region_map_wc.png"), width=30, units="mm", dpi=600, scale=1.6)
ggsave(ebsmap, filename=here("results","region_map_ebs.png"), dpi=600, width=40, units="mm", scale=1.6)

# get rid of unwanted whitespace
plot_crop(here("results","region_map_ebs.png"))
plot_crop(here("results","region_map_wc.png"))
plot_crop(here("results","region_map_neus.png"))
