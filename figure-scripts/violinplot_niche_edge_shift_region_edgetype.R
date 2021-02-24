library(tidyverse)
library(here)
library(stringr)

spp.bayes.edge.lm.df.summary <- read_csv(here("results","species_edge_shifts_vs_time.csv")) # edge shifts over time 
spp.bayes.niche.results <- read_csv(here("results","edge_thermal_extreme_tracked_summary.csv")) # thermal extreme tracked by each edge (warm, cold, both, or neither)

# note that the final version of this figure had annotated outliers. these were added manually in Inkscape by sorting the spp.bayes.edge.lm.df.summary dataframe by the median column to identify the most extreme range shifts.

gg.neus.violin <- spp.bayes.edge.lm.df.summary %>%
  left_join(spp.bayes.niche.results, by=c("species","region","quantile")) %>%
  filter(region=="neus") %>%
  mutate(quantile = recode(quantile, "quantile_0.99"="Poleward Edge",
                           "quantile_0.01"="Equatorward Edge"),
         quantile=as.factor(quantile),
         varTracked = recode(varTracked, 
                             "both" = "Both",
                             "predict.sstmin" = "Winter",
                             "predict.sstmax"="Summer",
                             "none" = "Neither"),
         varTracked = factor(varTracked, levels=c("Both","Winter","Summer","Neither"))) %>% 
  ggplot(aes(median, factor(varTracked), fill=varTracked)) + 
  geom_violin() +
  geom_vline(aes(xintercept=0),linetype="dashed",color="black")+
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(values=c("#9900cc","#3A4ED0","#DF2301","grey"), drop=FALSE) +
  coord_cartesian(xlim=c(-25, 25)) +
  scale_x_continuous(labels=seq(-25, 25, 5), breaks=seq(-25, 25, 5)) +  theme_bw() +
  theme_bw() +
  facet_wrap(~quantile, ncol=1) + 
  labs(x=NULL, y="Temperature Extreme Tracked", title="Northeast") +
  theme(legend.position = "none") +
  # theme(legend.position = "none", axis.text.x = element_blank()) +
  NULL
gg.neus.violin
ggsave(gg.neus.violin, filename=here("results","neus_edge_vs_niche_shifts.png"), width=110, units="mm", height=60, dpi=600, scale=1.5)

gg.wc.violin <- spp.bayes.edge.lm.df.summary %>%
  left_join(spp.bayes.niche.results, by=c("species","region","quantile")) %>%
  filter(region=="wc") %>%
  mutate(quantile = recode(quantile, "quantile_0.99"="Poleward Edge",
                           "quantile_0.01"="Equatorward Edge"),
         quantile=as.factor(quantile),
         varTracked = recode(varTracked, 
                             "both" = "Both",
                             "predict.sstmin" = "Winter",
                             "predict.sstmax"="Summer",
                             "none" = "Neither"),
         varTracked = factor(varTracked, levels=c("Both","Winter","Summer","Neither"))) %>% 
  ggplot(aes(median, factor(varTracked), fill=varTracked)) + 
  geom_violin() +
  geom_vline(aes(xintercept=0),linetype="dashed",color="black")+
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(values=c("#9900cc","#3A4ED0","#DF2301","grey"), drop=FALSE) +
  coord_cartesian(xlim=c(-25, 25)) +
  scale_x_continuous(labels=seq(-25, 25, 5), breaks=seq(-25, 25, 5)) +  theme_bw() +
  facet_wrap(~quantile, ncol=1) + 
  #  labs(x=NULL, y="Temperature Extreme Tracked", title="West Coast") +
  labs(x=NULL, y=NULL, title="West Coast") +
  theme(legend.position = "none") +
  # theme(legend.position = "none", axis.text.x = element_blank()) +
  NULL
gg.wc.violin
ggsave(gg.wc.violin, filename=here("results","wc_edge_vs_niche_shifts.png"), width=100, units="mm", height=35, dpi=600, scale=1.5)

gg.ebs.violin <- spp.bayes.edge.lm.df.summary %>%
  left_join(spp.bayes.niche.results, by=c("species","region","quantile")) %>%
  filter(region=="ebs") %>%
  mutate(quantile = recode(quantile, "quantile_0.99"="Poleward Edge",
                           "quantile_0.01"="Equatorward Edge"),
         quantile=as.factor(quantile),
         varTracked = recode(varTracked, 
                             "both" = "Both",
                             "predict.sstmin" = "Winter",
                             "predict.sstmax"="Summer",
                             "none" = "Neither"),
         varTracked = factor(varTracked, levels=c("Both","Winter","Summer","Neither"))) %>% 
  ggplot(aes(median, factor(varTracked), fill=varTracked)) + 
  geom_violin() +
  geom_vline(aes(xintercept=0),linetype="dashed",color="black")+
  geom_jitter(height = 0, width = 0.1) +
  scale_fill_manual(values=c("#9900cc","#3A4ED0","#DF2301","grey"), drop=FALSE) +
  coord_cartesian(xlim=c(-25, 25)) +
  scale_x_continuous(labels=seq(-25, 25, 5), breaks=seq(-25, 25, 5)) +
  theme_bw() +
  facet_wrap(~quantile, ncol=1) + 
  labs(x="Edge Shift (km/year)", y=NULL, title="Eastern Bering Sea") +
  theme(legend.position = "none") +
  NULL
gg.ebs.violin
ggsave(gg.ebs.violin, filename=here("results","ebs_edge_vs_niche_shifts.png"), width=100, units="mm", height=45, dpi=600, scale=1.5)
