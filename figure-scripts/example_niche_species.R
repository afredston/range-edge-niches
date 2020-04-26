library(here)
library(tidyverse)

# make example plots for methods schematic 
ex.spp1 <- "gadus macrocephalus" # good tracker, cold edge, EBS
ex.spp2 <- "sebastes pinniger" # non tracker, warm edge, WC
ex.spp3 <- "paralichthys oblongus" # lagged tracker, cold edge, NEUS

ex.spp.bayes.gg1 <- spp.bayes.niche.filter %>%
  filter(species==ex.spp1) %>%
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>%
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
  theme(legend.position="bottom") +
  NULL
ex.spp.bayes.gg1

ex.spp.bayes.gg2 <- spp.bayes.niche.filter %>%
  filter(species==ex.spp2) %>%
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>%
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
  theme(legend.position="bottom") +
  NULL
ex.spp.bayes.gg2

ex.spp.bayes.gg3 <- spp.bayes.niche.filter %>%
  filter(species==ex.spp3) %>%
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>%
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
  theme(legend.position="bottom") +
  NULL
ex.spp.bayes.gg3

ex.spp.time.gg1 <- dat.predict.niche %>%
  filter(species==ex.spp1) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
  theme(legend.position="bottom")+
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+
  NULL
ex.spp.time.gg1

ex.spp.time.gg2 <- dat.predict.niche %>%
  filter(species==ex.spp2) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
  theme(legend.position="bottom")+
  scale_x_continuous(limits=c(1982, 2018), breaks=seq(1982, 2018, 4))+
  NULL
ex.spp.time.gg2

ex.spp.time.gg3 <- dat.predict.niche %>%
  filter(species==ex.spp3) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
  theme(legend.position="bottom")+
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+
  NULL
ex.spp.time.gg3

ggsave(ex.spp.bayes.gg1, dpi=160, width=4, height=4, filename=here("results","example_1_posterior.png"))
ggsave(ex.spp.bayes.gg2, dpi=160, width=4, height=4, filename=here("results","example_2_posterior.png"))
ggsave(ex.spp.bayes.gg3, dpi=160, width=4, height=4, filename=here("results","example_3_posterior.png"))
ggsave(ex.spp.time.gg1, dpi=160, width=4, height=4, filename=here("results","example_1_niche.png"))
ggsave(ex.spp.time.gg2, dpi=160, width=4, height=4, filename=here("results","example_2_niche.png"))
ggsave(ex.spp.time.gg3, dpi=160, width=4, height=4, filename=here("results","example_3_niche.png"))