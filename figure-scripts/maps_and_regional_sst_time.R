# make plot of the three regions and temperature trends 
library(tidyverse)
library(here)
library(sf)
library(purrr)
library(broom)
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

wcgg <- wc_sst_summary %>%
  pivot_longer(cols=c(region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst), names_to = "sstvar", values_to = "sstvalue") %>%
  ggplot(aes(x=year, y=sstvalue, group=sstvar, color=sstvar)) +
  geom_line(size=1.1) +
  scale_color_manual(values=c("#DF2301","#FF8E02","#3A4ED0")) +
  scale_x_continuous(limits = c(1976, 2018), breaks=seq(1976, 2017, 4), expand = c(0, 0)) +
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

ebsgg <- ebs_sst_summary %>%
  pivot_longer(cols=c(region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst), names_to = "sstvar", values_to = "sstvalue") %>%
  ggplot(aes(x=year, y=sstvalue, group=sstvar, color=sstvar)) +
  geom_line(size=1.1) +
  scale_color_manual(values=c("#DF2301","#FF8E02","#3A4ED0")) +
  scale_x_continuous(limits = c(1988, 2018), breaks=seq(1988, 2017, 4), expand = c(0, 0)) +
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

ggsave(neusgg, filename=here("results","neus_sst_inset_plot.png"), height=2, width=2, dpi=160)
ggsave(wcgg, filename=here("results","wc_sst_inset_plot.png"), height=2, width=2,dpi=160)
ggsave(ebsgg, filename=here("results","ebs_sst_inset_plot.png"), height=2,width=2, dpi=160)

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

neusmap <- ggplot() + 
  geom_sf(data=neus_bathy, color="skyblue3") +
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=neus_coast, color="black", linetype="dashed", lwd=1.5) +
  scale_x_continuous(limits = neus_lonrange, expand = c(0, 0)) +
  scale_y_continuous(limits = neus_latrange, expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none",
        text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.title=element_blank(),
        plot.margin=margin(t = 5, r = 0, b = 15, l = 5, unit = "pt")) +
  NULL
neusmap

wcmap <- ggplot() + 
  geom_sf(data=wc_bathy, color="skyblue3") +
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=wc_coast, color="black", linetype="dashed", lwd=1.5) +
  scale_x_continuous(limits = wc_lonrange, expand = c(0, 0)) +
  scale_y_continuous(limits = wc_latrange, expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none",
        text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_blank(),
        plot.margin=margin(t = 5, r = 0, b = 15, l = 5, unit = "pt")) +
  NULL

ebsmap <- ggplot() + 
  geom_sf(data=ebs_bathy, color="skyblue3") +
  geom_sf(data=usoutline, color="#999999") +
  geom_sf(data=ebs.line, color="black", linetype="dashed", lwd=1.5) +
#  geom_segment(aes(x=-158, xend=-170, y=54, yend=66),color="black", linetype="dashed", lwd=1.5) + 
  scale_x_continuous(limits = ebs_lonrange, expand = c(0, 0), breaks=seq(-180, -154, 2)) +
  scale_y_continuous(limits = ebs_latrange, expand = c(0, 0), breaks=seq(54, 66, 2)) +
  theme_bw() +
  theme(legend.position = "none",
        text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.title=element_blank(),
        plot.margin=margin(t = 5, r = 0, b = 15, l = 5, unit = "pt"),
        axis.text.x = element_text(angle=45, hjust=1)) +
  NULL
ebsmap

ggsave(neusmap, filename=here("results","region_map_neus.png"), height=4, dpi=160)
ggsave(wcmap, filename=here("results","region_map_wc.png"), height=4, dpi=160)
ggsave(ebsmap, filename=here("results","region_map_ebs.png"), height=4, width=8, dpi=160)
