library(here)
library(tidyverse)


spp.bayes.edge.lm.df.summary <- read_csv(here("results","species_edge_shifts_vs_time.csv")) # edge shifts over time

spp.edge.niche.summary <- read_csv(here("results","species_bayes_niche_lm_summary.csv")) %>% # niche shifts over time
  rename('mean.niche'=mean,
         'median.niche' = median,
         'lower.niche' = lower,
         'upper.niche' = upper) %>% 
  left_join(spp.bayes.edge.lm.df.summary, by=c('species','region','quantile')) %>% 
  mutate(quantile = recode(quantile, "quantile_0.01"="Equatorward Edge","quantile_0.99"="Poleward Edge"),
         predicted.var = recode(predicted.var, "predict.sstmax" = "Summer Temperature", "predict.sstmin" = "Winter Temperature"),
         region = recode(region, "ebs" = "Eastern Bering Sea","wc" = 'West Coast','neus' = 'Northeast')
         )

gg.edge.niche.scatter <- spp.edge.niche.summary %>% 
  filter(!species=='merluccius albidus') %>% # eliminating one outlier that shifted >20 km/yr into summer temperatures that were 0.2 C/yr hotter!
  ggplot(aes(x=mean, y=mean.niche, xmin=lower, xmax=upper, ymin=lower.niche, ymax=upper.niche, color=quantile, fill=quantile, shape=region)) +
  scale_color_manual(values=c('#00441b','#238b45')) +
  scale_fill_manual(values=c('#00441b','#238b45')) +
  geom_errorbar() +
  geom_errorbarh() + 
  geom_hline(yintercept=0, color="black", linetype="dashed") +
  geom_vline(xintercept=0, color="black", linetype="dashed") +
  geom_point(size=2.5, alpha=0.5) +
  scale_x_continuous(limits=c(-33, 33), breaks=seq(-30, 30, 10)) +
  scale_y_continuous(limits = c(-.12, .12), breaks=seq(-.10, .10, .05)) + 
  facet_wrap(~predicted.var) +
  theme_bw() + 
  labs(x='Edge Shift (km/yr)', y = 'Niche Shift (°C/yr)') +
  theme(
    legend.position = "bottom",
    legend.spacing.x = unit(1.0, 'mm'),
    legend.title = element_blank(),
    strip.background =element_rect(fill="white")
  ) +
  NULL

ggsave(gg.edge.niche.scatter, filename=here("results","edge_vs_niche_shift_scatterplot.png"), width=110, units="mm", height=60, dpi=600, scale=1.5)
