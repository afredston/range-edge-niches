library(tidyverse)
library(here)

edge.spp.dat <- readRDS(here("processed-data","all_edge_spp_df.rds"))

neus.gg <- edge.spp.dat %>%
  filter(region=="neus",
         axis=="coast_km") %>%
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) +
  facet_wrap(~species, ncol=5)+
  theme(legend.position = "none") +
  labs(x="Year", y="Coastal Distance (km)") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")+
  NULL
ggsave(neus.gg, filename=here("results","VAST_edges_neus_no_SE_filter.png"),height=12,width=7,dpi=160, scale=1.1)

ebs.gg <- edge.spp.dat %>%
  filter(region=="ebs",
         axis=="NW_km") %>% 
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() + 
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) + 
  facet_wrap(~species, ncol=5)+
  theme(legend.position = "none") +
  labs(x="Year", y="Northwesting (km)") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")+
  NULL
ggsave(ebs.gg, filename=here("results","VAST_edges_ebs_no_SE_filter.png"),height=12,width=7,dpi=160, scale=1.1)

wc.gg <- edge.spp.dat %>%
  filter(region=="wc",
         axis=="coast_km") %>% 
  ggplot(aes(x=year, y=Estimate, group=quantile, color=quantile)) +
  geom_point() +
  geom_line() + 
  geom_errorbar(aes(ymin=Estimate-Std.Error, ymax=Estimate+Std.Error, group=quantile, color=quantile)) + 
  facet_wrap(~species, ncol=5)+
  theme(legend.position = "none") +
  labs(x="Year", y="Coastal Distance (km)") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1), legend.position="bottom")+
  NULL
ggsave(wc.gg, filename=here("results","VAST_edges_wc_no_SE_filter.png"),height=12,width=7,dpi=160, scale=1.1)
rm(list=ls())
