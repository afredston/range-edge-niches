library(tidyverse)
library(here)

spp.bayes.niche.lm.stats <- read_csv( here("results","species_bayes_niche_lm_summary.csv"))

dat.predict.niche <- read_csv(here("processed-data","species_thermal_niche_v_time.csv"))

# make example plots for methods schematic 
ex.spp.ebs <- "paralithodes camtschaticus" # red king crab
ex.spp.neus <- "gadus morhua" # atlantic cod
ex.spp.wc <- "sebastes pinniger" # canary rockfish

ex.niche.gg.ebs <- dat.predict.niche %>%
  filter(region=="ebs",quantile=="quantile_0.99", species==ex.spp.ebs) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="SST at Edge (°C)", color=NULL) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))+
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+ 
  NULL

ex.niche.gg.neus <- dat.predict.niche %>%
  filter(region=="neus",quantile=="quantile_0.01", species==ex.spp.neus) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="SST at Edge (°C)", color=NULL) +
  theme(legend.position="none",  
        axis.text.x=element_text(angle=45, hjust=1))+
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+
  NULL

ex.niche.gg.wc <- dat.predict.niche %>%
  filter(region=="wc",quantile=="quantile_0.01", species==ex.spp.wc) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="SST at Edge (°C)", color=NULL) +
  theme(legend.position="none", 
        axis.text.x=element_text(angle=45, hjust=1))+
  scale_x_continuous(limits=c(1976, 2018), breaks=seq(1976, 2018, 5))+
  NULL

ggsave(ex.niche.gg.ebs, dpi=600, width=1.5, height=1.5, filename=here("results",paste0("example_niche_",ex.spp.ebs,".png")), scale=1.5)
ggsave(ex.niche.gg.neus, dpi=600, width=1.5, height=1.5, filename=here("results",paste0("example_niche_",ex.spp.neus,".png")), scale=1.5)
ggsave(ex.niche.gg.wc, dpi=600, width=1.5, height=1.5, filename=here("results",paste0("example_niche_",ex.spp.wc,".png")), scale=1.5)


# make regional plots for supplement 
neus.cold.niche.time.gg <- dat.predict.niche %>%
  filter(region=="neus", quantile=="quantile_0.99") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", title="Northeast Poleward Edges", color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+ 
  facet_wrap(~species, ncol=5) +
  NULL
neus.cold.niche.time.gg

neus.warm.niche.time.gg <- dat.predict.niche %>%
  filter(region=="neus", quantile=="quantile_0.01") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)",title="Northeast Equatorward Edges",  color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+ 
  facet_wrap(~species, ncol=5) +
  NULL
neus.warm.niche.time.gg

wc.cold.niche.time.gg <- dat.predict.niche %>%
  filter(region=="wc", quantile=="quantile_0.99") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", title="West Coast Poleward Edges", color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1976, 2018), breaks=seq(1976, 2018, 4))+ 
  facet_wrap(~species, ncol=3) +
  NULL
wc.cold.niche.time.gg

wc.warm.niche.time.gg <- dat.predict.niche %>%
  filter(region=="wc", quantile=="quantile_0.01") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)",title="West Coast Equatorward Edges",  color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1976, 2018), breaks=seq(1976, 2018, 4))+ 
  facet_wrap(~species, ncol=5) +
  NULL
wc.warm.niche.time.gg

ebs.cold.niche.time.gg <- dat.predict.niche %>%
  filter(region=="ebs", quantile=="quantile_0.99") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", title="Eastern Bering Sea Poleward Edges", color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+ 
  facet_wrap(~species, ncol=5) +
  NULL
ebs.cold.niche.time.gg

ebs.warm.niche.time.gg <- dat.predict.niche %>%
  filter(region=="ebs", quantile=="quantile_0.01") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot(aes(x=year_match, y=sst, ymin=sst-sstSE, ymax=sst+sstSE, color=predicted.var)) +
  geom_point() +
  geom_line() +
  geom_errorbar() +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)",title="Eastern Bering Sea Equatorward Edges",  color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+ 
  
  facet_wrap(~species, ncol=5) +
  NULL
ebs.warm.niche.time.gg


ggsave(neus.cold.niche.time.gg, filename=here("results","neus_cold_edge_niche_shifts_v_time.png"), height=10.5,width=8.25, dpi=160, scale=1.2)
ggsave(neus.warm.niche.time.gg, filename=here("results","neus_warm_edge_niche_shifts_v_time.png"), height=10.5,width=8.25,dpi=160, scale=1.2)

ggsave(wc.cold.niche.time.gg, filename=here("results","wc_cold_edge_niche_shifts_v_time.png"), height=10.5,width=8.25,dpi=160, scale=1.2)
ggsave(wc.warm.niche.time.gg, filename=here("results","wc_warm_edge_niche_shifts_v_time.png"), height=10.5,width=8.25, dpi=160, scale=1.2)

ggsave(ebs.cold.niche.time.gg, filename=here("results","ebs_cold_edge_niche_shifts_v_time.png"), height=10.5,width=8.25,dpi=160, scale=1.2)
ggsave(ebs.warm.niche.time.gg, filename=here("results","ebs_warm_edge_niche_shifts_v_time.png"), height=10.5,width=8.25, dpi=160, scale=1.2)
