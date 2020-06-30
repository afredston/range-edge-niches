library(here)
library(tidyverse)

spp.bayes.niche.groups <- read_csv(here("results","species_by_thermal_niche_group.csv"))
spp.bayes.lm.df.summary <- read_csv(here("results","species_edge_shifts_vs_time.csv"))

edge.shift.df <- spp.bayes.niche.groups %>%
  left_join(spp.bayes.lm.df.summary) %>%
  arrange(desc(species)) %>%
  mutate(quantile = recode(quantile, "quantile_0.01"="Warm Edge","quantile_0.99"="Cold Edge"),
         niche.group = recode(niche.group, "good_tracker" = "TNH", "partial_tracker" = "PTH", "non_tracker" = "TIH"),
    species = str_to_sentence(species),
         species = factor(species, unique(species)),
         quantile = as.factor(quantile),
    niche.group = factor(niche.group, levels=c("TNH","PTH","TIH"))) 

neus.edge.gg <- edge.shift.df %>%
  filter(region=="neus") %>%
  ggplot(aes(y=species, x=median, xmin=lower, xmax=upper, color=quantile, fill=quantile, shape=niche.group)) +
  geom_errorbarh() +
  geom_point() +
  scale_color_manual(values=c(`Warm Edge`='#C7361D', `Cold Edge`='#3A4ED0')) +
  geom_vline(xintercept=0, color="black") +
  labs(y=NULL, x="Edge Shift (km/yr)", title="Northeast", fill="", color="", shape="") +
  scale_x_continuous(breaks=seq(-60, 60, 10), limits=c(-68, 68)) + 
  theme_bw() +
  theme(legend.position = c(0.2, 0.9),
  #      legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="sans",size=10,color="black"),
        legend.text = element_text(size=10),
        axis.title=element_text(family="sans",size=12,color="black"),
        axis.text=element_text(family="sans",size=10,color="black"),
        plot.margin=unit(c(5.5, 15, 5.5, 5.5), "points")) +
  NULL
neus.edge.gg

wc.edge.gg <- edge.shift.df %>%
  filter(region=="wc") %>%
  ggplot(aes(y=species, x=median, xmin=lower, xmax=upper, color=quantile, fill=quantile, shape=niche.group)) +
  geom_errorbarh() +
  geom_point() +
  scale_color_manual(values=c(`Warm Edge`='#C7361D', `Cold Edge`='#3A4ED0')) +
  geom_vline(xintercept=0, color="black") +
  labs(y=NULL, x="Edge Shift (km/yr)", title="West Coast", fill="", color="", shape="") +
  scale_x_continuous(breaks=seq(-60, 60, 10), limits=c(-68, 68)) + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="sans",size=10,color="black"),
        legend.text = element_text(size=10),
        axis.title=element_text(family="sans",size=12,color="black"),
        axis.text=element_text(family="sans",size=10,color="black"),
        plot.margin=unit(c(5.5, 15, 5.5, 5.5), "points")) +
  NULL
wc.edge.gg

ebs.edge.gg <- edge.shift.df %>%
  filter(region=="ebs") %>%
  ggplot(aes(y=species, x=median, xmin=lower, xmax=upper, color=quantile, fill=quantile, shape=niche.group)) +
  geom_errorbarh() +
  geom_point() +
  scale_color_manual(values=c(`Warm Edge`='#C7361D', `Cold Edge`='#3A4ED0')) +
  geom_vline(xintercept=0, color="black") +
  labs(y=NULL, x="Edge Shift (km/yr)", title="Eastern Bering Sea", fill="", color="", shape="") +
  scale_x_continuous(breaks=seq(-60, 60, 10), limits=c(-75, 75)) + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="sans",size=10,color="black"),
        legend.text = element_text(size=10),
        axis.title=element_text(family="sans",size=12,color="black"),
        axis.text=element_text(family="sans",size=10,color="black"),
        plot.margin=unit(c(5.5, 15, 5.5, 5.5), "points")) +
  NULL
ebs.edge.gg

ggsave(neus.edge.gg, filename=here("results","neus_spp_edge_shifts_forest_plot.png"), height=10, width=5, dpi=160, scale=1.2)
ggsave(wc.edge.gg, filename=here("results","wc_spp_edge_shifts_forest_plot.png"), height=10, width=5, dpi=160, scale=1.2)
ggsave(ebs.edge.gg, filename=here("results","ebs_spp_edge_shifts_forest_plot.png"), height=10, width=5, dpi=160, scale=1.2)
