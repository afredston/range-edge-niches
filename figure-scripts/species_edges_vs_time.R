library(tidyverse)
library(here)
library(stringr)

edge.spp.dat <- readRDS(here("processed-data","all_edge_spp_df.rds"))%>%
  ungroup() %>% # undo rowwise nature
  mutate(axis = as.character(axis)) # convert from factor

#####
# make example plots for methods schematic 
#####
ex.spp.ebs <- "limanda proboscidea" # longhead dab
ex.spp.neus <- "homarus americanus" # lobster
ex.spp.wc <- "sebastes pinniger" # canary rockfish

# make time-series plots for example figure:
ex.gg.ebs <- edge.spp.dat %>%
  filter(region=="ebs",
         axis=="line_km",
         species==ex.spp.ebs) %>% 
  ggplot(aes(x=year, y=Estimate)) +
  geom_point(color="grey", fill="grey") +
  geom_line(color="grey") +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error), color="grey") +
  theme(legend.position = "none") +
  labs(x="Year", y="Axis Distance (km)") +
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  NULL

ex.gg.neus <- edge.spp.dat %>%
  filter(region=="neus",
         axis=="coast_km",
         species==ex.spp.neus) %>% 
  ggplot(aes(x=year, y=Estimate)) +
  geom_point(color="grey", fill="grey") +
  geom_line(color="grey") +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error), color="grey") +
  theme(legend.position = "none") +
  labs(x="Year", y="Coastal Distance (km)") +
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  NULL

ex.gg.wc <- edge.spp.dat %>%
  filter(region=="wc",
         axis=="coast_km",
         species==ex.spp.wc) %>% 
  ggplot(aes(x=year, y=Estimate)) +
  geom_point(color="grey", fill="grey") +
  geom_line(color="grey") +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error), color="grey") +
  theme(legend.position = "none") +
  labs(x="Year", y="Coastal Distance (km)") +
  scale_x_continuous(limits=c(1976, 2018), breaks=seq(1976, 2018, 5))+
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  NULL

ggsave(ex.gg.ebs, dpi=600, width=1.5, height=1.5, filename=here("results",paste0("example_edge_",ex.spp.ebs,".png")), scale=1.5)
ggsave(ex.gg.neus, dpi=600, width=1.5, height=1.5, filename=here("results",paste0("example_edge_",ex.spp.neus,".png")), scale=1.5)
ggsave(ex.gg.wc, dpi=600, width=1.5, height=1.5, filename=here("results",paste0("example_edge_",ex.spp.wc,".png")), scale=1.5)


#####
# plot cold edges over time, both 0.01 and 0.05 quantiles (for supplement)
#####

# need to pull in the original dfs because edge.spp.dat only has 0.99 and 0.01 quantiles 
# some code to harmonize name format because this is the raw VAST output

# get edge positions with comparison quantile 
neus05 <- readRDS(here("processed-data","neus_relative_SE_vast_edge_df.rds")) %>% 
  filter(quantile=="quantile_0.05",axis=="coast_km") %>%  
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) %>% 
  filter(species %in% edge.spp.dat[edge.spp.dat$region=="neus" & edge.spp.dat$quantile=="quantile_0.01",]$species)  

# plot the edge positions we used against the other metric for warm edges 
neus.eq.gg <- edge.spp.dat %>% 
  filter(region=="neus", quantile=="quantile_0.01",axis=="coast_km") %>% 
  full_join(neus05) %>% 
  mutate(quantile = recode(quantile, "quantile_0.01"="Quantile=0.01",
                           "quantile_0.05"="Quantile=0.05"),
         species = str_to_sentence(species)) %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  scale_color_manual(values=c("black","grey")) +
  facet_wrap(~species, ncol=4)+
  labs(x="Year", y="Coastal Distance (km)",title="Northeast Equatorward Edges") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.title = element_blank())+
  NULL
neus.eq.gg

neus95 <- readRDS(here("processed-data","neus_relative_SE_vast_edge_df.rds")) %>% 
  filter(quantile=="quantile_0.95",axis=="coast_km") %>%  
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) %>% 
  filter(species %in% edge.spp.dat[edge.spp.dat$region=="neus" & edge.spp.dat$quantile=="quantile_0.99",]$species)  

# same for cold edges
neus.pol.gg <- edge.spp.dat %>% 
  filter(region=="neus", quantile=="quantile_0.99",axis=="coast_km") %>% 
  full_join(neus95) %>% 
  mutate(quantile = recode(quantile, "quantile_0.99"="Quantile=0.99",
                           "quantile_0.95"="Quantile=0.95"),
         species = str_to_sentence(species)) %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  scale_color_manual(values=c("grey","black")) +
  facet_wrap(~species, ncol=4)+
  labs(x="Year", y="Coastal Distance (km)", title="Northeast Poleward Edges") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.title = element_blank())+
  NULL
neus.pol.gg

wc05 <- readRDS(here("processed-data","wc_relative_SE_vast_edge_df.rds")) %>% 
  filter(quantile=="quantile_0.05",axis=="coast_km") %>%  
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) %>% 
  filter(species %in% edge.spp.dat[edge.spp.dat$region=="wc" & edge.spp.dat$quantile=="quantile_0.01",]$species)  

wc.eq.gg <- edge.spp.dat %>% 
  filter(region=="wc", quantile=="quantile_0.01",axis=="coast_km") %>% 
  full_join(wc05) %>% 
  mutate(quantile = recode(quantile, "quantile_0.01"="Quantile=0.01",
                           "quantile_0.05"="Quantile=0.05"),
         species = str_to_sentence(species)) %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  scale_color_manual(values=c("black","grey")) +
  facet_wrap(~species, ncol=3)+
  labs(x="Year", y="Coastal Distance (km)", title="West Coast Equatorward Edges") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.title = element_blank())+
  NULL
wc.eq.gg

wc95 <- readRDS(here("processed-data","wc_relative_SE_vast_edge_df.rds")) %>% 
  filter(quantile=="quantile_0.95",axis=="coast_km") %>%  
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) %>% 
  filter(species %in% edge.spp.dat[edge.spp.dat$region=="wc" & edge.spp.dat$quantile=="quantile_0.99",]$species)  

wc.pol.gg <- edge.spp.dat %>% 
  filter(region=="wc", quantile=="quantile_0.99",axis=="coast_km") %>% 
  full_join(wc95) %>% 
  mutate(quantile = recode(quantile, "quantile_0.99"="Quantile=0.99",
                           "quantile_0.95"="Quantile=0.95"),
         species = str_to_sentence(species)) %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  scale_color_manual(values=c("grey","black")) +
  facet_wrap(~species, ncol=3)+
  labs(x="Year", y="Coastal Distance (km)", title="West Coast Poleward Edges") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.title = element_blank())+
  NULL
wc.pol.gg

ebs05 <- readRDS(here("processed-data","ebs_relative_SE_vast_edge_df.rds")) %>% 
  filter(quantile=="quantile_0.05",axis=="line_km") %>%  
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) %>% 
  filter(species %in% edge.spp.dat[edge.spp.dat$region=="ebs" & edge.spp.dat$quantile=="quantile_0.01",]$species)  

ebs.eq.gg <- edge.spp.dat %>% 
  filter(region=="ebs", quantile=="quantile_0.01",axis=="line_km") %>% 
  full_join(ebs05) %>% 
  mutate(quantile = recode(quantile, "quantile_0.01"="Quantile=0.01",
                           "quantile_0.05"="Quantile=0.05"),
         species = str_to_sentence(species)) %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  scale_color_manual(values=c("black","grey")) +
  facet_wrap(~species, ncol=4)+
  labs(x="Year", y="Middle Domain Axis (km)", title="Eastern Bering Sea Equatorward Edges") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.title = element_blank())+
  NULL
ebs.eq.gg

ebs95 <- readRDS(here("processed-data","ebs_relative_SE_vast_edge_df.rds")) %>% 
  filter(quantile=="quantile_0.95",axis=="line_km") %>%  
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) %>% 
  filter(species %in% edge.spp.dat[edge.spp.dat$region=="ebs" & edge.spp.dat$quantile=="quantile_0.99",]$species)  

ebs.pol.gg <- edge.spp.dat %>% 
  filter(region=="ebs", quantile=="quantile_0.99",axis=="line_km") %>% 
  full_join(ebs95) %>% 
  mutate(quantile = recode(quantile, "quantile_0.99"="Quantile=0.99",
                           "quantile_0.95"="Quantile=0.95"),
         species = str_to_sentence(species)) %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  scale_color_manual(values=c("grey","black")) +
  facet_wrap(~species, ncol=4)+
  labs(x="Year", y="Middle Domain Axis (km)", title="Eastern Bering Sea Poleward Edges") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom", legend.title = element_blank())+
  NULL
ebs.pol.gg

ggsave(neus.eq.gg, filename=here("results","VAST_equatorward_edges_neus.png"),height=10.5,width=8.25,dpi=160, scale=1.1)
ggsave(neus.pol.gg, filename=here("results","VAST_poleward_edges_neus.png"),height=10.5,width=8.25,dpi=160, scale=1.1)
ggsave(wc.eq.gg, filename=here("results","VAST_equatorward_edges_wc.png"),height=10.5,width=8.25,dpi=160, scale=1.1)
ggsave(wc.pol.gg, filename=here("results","VAST_poleward_edges_wc.png"),height=10.5,width=8.25,dpi=160, scale=1.1)
ggsave(ebs.eq.gg, filename=here("results","VAST_equatorward_edges_ebs.png"),height=10.5,width=8.25,dpi=160, scale=1.1)
ggsave(ebs.pol.gg, filename=here("results","VAST_poleward_edges_ebs.png"),height=10.5,width=8.25,dpi=160, scale=1.1)
