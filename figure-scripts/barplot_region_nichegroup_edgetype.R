library(here)
library(tidyverse)

spp.bayes.niche.groups <- read_csv(here("results","species_by_thermal_niche_group.csv"))
spp.bayes.niche.lm.stats <- read_csv( here("results","species_bayes_niche_lm_summary.csv"))

# make barplot summarizing number of edges by: region, edge type, and niche group
gg.niche.barplot <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(quantile=factor(quantile, levels=c('quantile_0.01','quantile_0.99')),
         quantile=recode(quantile,
                         quantile_0.01="Warm Edge",
                         quantile_0.99="Cold Edge"),
         region=factor(region, levels=c('neus','wc','ebs')),
         region=recode(region,
                       ebs="Eastern Bering Sea",
                       neus="Northeast",
                       wc="West Coast"),
         niche.group=factor(niche.group, levels=c('good_tracker','partial_tracker','non_tracker')),
         niche.group=recode(niche.group,
                            good_tracker="TNH",
                            partial_tracker="PTH",
                            non_tracker="TIH")) %>%
  group_by(quantile, region, niche.group) %>%
  mutate(count = length(species)) %>%
  select(count, region, niche.group, quantile) %>%
  ungroup() %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(y=forcats::fct_rev(niche.group), x=count, fill=quantile), stat="identity", color="black") +
  scale_fill_manual(values=c("gray50","black")) +
  theme_bw() +
  labs(x=NULL, y=NULL)+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  facet_wrap(~region, ncol=1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x="Number of Range Edges") +
  NULL
gg.niche.barplot
ggsave(gg.niche.barplot, filename=here("results","barplot_all_groups.png"),dpi=160, width=6, height=4)
