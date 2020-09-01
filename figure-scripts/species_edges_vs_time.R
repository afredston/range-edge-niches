library(tidyverse)
library(here)
library(stringr)

edge.spp.dat <- readRDS(here("processed-data","all_edge_spp_df.rds"))%>%
  ungroup() %>% # undo rowwise nature
  mutate(axis = as.character(axis)) # convert from factor

# make example plots for methods schematic 
ex.spp1 <- "chionoecetes opilio" # snow crab 
ex.spp2 <- "centropristis striata" # black seabass
# ex.spp2 <- "gadus morhua" # atlantic cod
ex.spp3 <- "sebastes pinniger" # canary rockfish

# make time-series plots for example figure:
ex1gg <- edge.spp.dat %>%
  filter(region=="ebs",
         axis=="line_km",
         species==ex.spp1) %>% 
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

ex2gg <- edge.spp.dat %>%
  filter(region=="neus",
         axis=="coast_km",
         species==ex.spp2) %>% 
  ggplot(aes(x=year, y=Estimate)) +
  geom_point(color="grey", fill="grey") +
  geom_line(color="grey") +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error), color="grey") +
  theme(legend.position = "none") +
  labs(x="Year", y="Coastal Distance (km)") +
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 4)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  NULL

ex3gg <- edge.spp.dat %>%
  filter(region=="wc",
         axis=="coast_km",
         species==ex.spp3) %>% 
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

ggsave(ex1gg, dpi=600, width=4, height=4, filename=here("results",paste0("example_edge_",ex.spp1,".png")), scale=1.2)
ggsave(ex2gg, dpi=600, width=4, height=4, filename=here("results",paste0("example_edge_",ex.spp2,".png")))
ggsave(ex3gg, dpi=600, width=1.5, height=1.5, filename=here("results",paste0("example_edge_",ex.spp3,".png")), scale=1.5)

neus.gg <- edge.spp.dat %>%
  filter(region=="neus",
         axis=="coast_km") %>%
  mutate(quantile = recode(quantile, "quantile_0.99"="Poleward Edge",
                           "quantile_0.01"="Equatorward Edge"),
         species = str_to_sentence(species)) %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  scale_color_manual(values=c("grey","black")) +
  facet_wrap(~species, ncol=5)+
  theme(legend.position = "none") +
  labs(x="Year", y="Coastal Distance (km)") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")+
  NULL
neus.gg
ggsave(neus.gg, filename=here("results","VAST_edges_neus.png"),height=10.5,width=8.25,dpi=160, scale=1.1)

ebs.gg <- edge.spp.dat %>%
  filter(region=="ebs",
         axis=="line_km") %>%
  mutate(quantile = recode(quantile, "quantile_0.99"="Poleward Edge",
                           "quantile_0.01"="Equatorward Edge"),
         species = str_to_sentence(species)) %>% 
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() + 
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) + 
  scale_color_manual(values=c("grey","black")) +
  facet_wrap(~species, ncol=6)+
  theme(legend.position = "none") +
  labs(x="Year", y="Middle Domain Axis (km)") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")+
  NULL
ggsave(ebs.gg, filename=here("results","VAST_edges_ebs.png"),height=10.5,width=8.25,dpi=160, scale=1.1)

wc.gg <- edge.spp.dat %>%
  filter(region=="wc",
         axis=="coast_km")  %>%
  mutate(quantile = recode(quantile, "quantile_0.99"="Poleward Edge",
                           "quantile_0.01"="Equatorward Edge"),
         species = str_to_sentence(species)) %>% 
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() + 
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) + 
  scale_color_manual(values=c("grey","black")) +
  facet_wrap(~species, ncol=5)+
  theme(legend.position = "none") +
  labs(x="Year", y="Coastal Distance (km)") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")+
  NULL
ggsave(wc.gg, filename=here("results","VAST_edges_wc.png"),height=10.5,width=8.25,dpi=160, scale=1.1)
rm(list=ls())
