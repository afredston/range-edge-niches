library(tidyverse)
library(here)

spp.bayes.niche.lm.stats <- read_csv( here("results","species_bayes_niche_lm_summary.csv"))

dat.predict.niche <- read_csv(here("processed-data","species_thermal_niche_v_time.csv"))

neus.cold.niche.time.gg <- dat.predict.niche %>%
  filter(region=="neus", quantile=="quantile_0.99") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", title="Northeast Cold Edges", color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+ 
  facet_wrap(~species, ncol=5) +
  NULL
neus.cold.niche.time.gg

neus.warm.niche.time.gg <- dat.predict.niche %>%
  filter(region=="neus", quantile=="quantile_0.01") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)",title="Northeast Warm Edges",  color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+ 
  facet_wrap(~species, ncol=5) +
  NULL
neus.warm.niche.time.gg

wc.cold.niche.time.gg <- dat.predict.niche %>%
  filter(region=="wc", quantile=="quantile_0.99") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", title="West Coast Cold Edges", color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1976, 2018), breaks=seq(1976, 2018, 4))+ 
  facet_wrap(~species, ncol=5) +
  NULL
wc.cold.niche.time.gg

wc.warm.niche.time.gg <- dat.predict.niche %>%
  filter(region=="wc", quantile=="quantile_0.01") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)",title="West Coast Warm Edges",  color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1976, 2018), breaks=seq(1976, 2018, 4))+ 
  facet_wrap(~species, ncol=5) +
  NULL
wc.warm.niche.time.gg

ebs.cold.niche.time.gg <- dat.predict.niche %>%
  filter(region=="ebs", quantile=="quantile_0.99") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", title="Eastern Bering Sea Cold Edges", color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+ 
  facet_wrap(~species, ncol=5) +
  NULL
ebs.cold.niche.time.gg

ebs.warm.niche.time.gg <- dat.predict.niche %>%
  filter(region=="ebs", quantile=="quantile_0.01") %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)",title="Eastern Bering Sea Warm Edges",  color=NULL) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90))+
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+ 
  
  facet_wrap(~species, ncol=5) +
  NULL
ebs.warm.niche.time.gg

