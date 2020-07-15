library(tidyverse)
library(here)

# DON'T FORGET TO RE-RUN VALIDATE_EDGE_SPP AFTER RE-RUNNING VAST!
dat.models <- readRDS(here("processed-data","all_edge_spp_df.rds")) %>%
  filter(axis %in% c('coast_km','NW_km')) %>%
  ungroup()

# how many species?

dat.summary <- dat.models %>% 
  group_by(species, quantile, region) %>%
  summarise() %>% 
  ungroup() %>%
  mutate(species=str_to_sentence(species),
         region=recode(region,
                       ebs="Eastern Bering Sea",
                       neus="Northeast",
                       wc="West Coast"),
         quantile=recode(quantile,
                         quantile_0.01="Warm Limit",
                         quantile_0.99="Cold Limit"))
write_csv(dat.summary, here("results","edge_species.csv"))

dat.models.groups <- dat.models %>%
  select(species, edgetype, region, taxongroup) %>%
  distinct() %>%
  mutate(species = as.factor(species), 
         edgetype=as.factor(edgetype),
         region=as.factor(region),
         taxongroup=as.factor(taxongroup))

# read in results from Bayesian models 

spp.bayes.edge.lm.df.summary <- read_csv(here("results","species_edge_shifts_vs_time.csv")) # edge shifts over time
dat.predict.niche <- read_csv(here("processed-data","species_thermal_niche_v_time.csv")) # niche data
spp.bayes.niche.lm.stats <- read_csv(here("results","species_bayes_niche_lm_summary.csv")) # niche shifts over time 


#######################
### classify results by niche hypothesis 
#######################

# identify which hypothesis the posterior beta is consistent with, and write out 
spp.bayes.niche.groups <- spp.bayes.niche.lm.stats %>%
  left_join(dat.models.groups %>% select(-edgetype), by=c("region","species")) %>% # add fish/invert column 
  rowwise() %>%
  mutate(crosses0 = ifelse(lower<0 & upper>0, TRUE, FALSE)) %>%
  group_by(region, quantile, species) %>% 
  mutate(niche.group = ifelse(min(abs(mean)) < 0.01, "good_tracker", # if ONE temp extreme has zero change -> good tracker 
                              ifelse(TRUE %in% crosses0, "partial_tracker", "non_tracker")) # if ONE temp extreme crosses zero -> lagged tracker 
  ) %>%
  select(region, quantile, species, niche.group) %>% 
  distinct()
write_csv(spp.bayes.niche.groups, here("results","species_by_thermal_niche_group.csv"))

# calculate stats
spp.bayes.niche.groups %>% left_join(dat.models.groups) %>% group_by(niche.group) %>% summarise(n=n()) # break down by niche grouping
spp.bayes.niche.groups %>% left_join(dat.models.groups) %>% group_by(niche.group, taxongroup) %>% summarise(n=n()) # break down by fish/invert and niche grouping
spp.bayes.niche.groups %>% left_join(dat.models.groups) %>% group_by(region, taxongroup) %>% summarise(n=n()) # by region and fish/invert
spp.bayes.niche.groups %>% group_by(niche.group, region) %>% summarise(n=n()) # by region and niche grouping
spp.bayes.niche.groups %>% group_by(niche.group, quantile) %>% summarise(n=n()) # by edge type and niche grouping
spp.bayes.niche.groups %>% group_by(region, quantile) %>% summarise(n=n()) # by region and edge type
spp.bayes.niche.groups %>% group_by(region, quantile, niche.group) %>% summarise(n=n()) # by region, niche grouping, and edge type


# of the temperature-independent species, which had a significant edge shift over time? 
spp.tih <- spp.bayes.niche.groups %>%
  filter(niche.group=="non_tracker") %>%
  left_join(spp.bayes.edge.lm.df.summary) %>%
  rowwise() %>%
  mutate(crosses0 = ifelse(lower<0 & upper>0, TRUE, FALSE)) %>%
  mutate(signif.shift = ifelse(TRUE %in% crosses0, "N", "Y"))

spp.tnh <- spp.bayes.niche.groups %>%
  filter(niche.group=="good_tracker") %>%
  left_join(spp.bayes.edge.lm.df.summary) %>%
  rowwise() %>%
  mutate(crosses0 = ifelse(lower<0 & upper>0, TRUE, FALSE)) %>%
  mutate(signif.shift = ifelse(TRUE %in% crosses0, "N", "Y"),
         sign.shift = ifelse(mean>0, "positive", "negative")) %>%
  group_by(signif.shift, sign.shift) %>% 
  summarise(n=n())
