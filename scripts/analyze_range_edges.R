# this script calculates assorted statistics reported in the manuscript, and saves some results summaries for reference 
library(tidyverse)
library(here)

################
### load data
################

# DON'T FORGET TO RE-RUN VALIDATE_EDGE_SPP AFTER RE-RUNNING VAST!
dat.models <- readRDS(here("processed-data","all_edge_spp_df.rds")) %>%
  ungroup() %>% # undo rowwise nature
  mutate(axis = as.character(axis)) %>% # convert from factor
  filter(axis %in% c('coast_km','line_km')) 

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
                         quantile_0.01="Warm Edge",
                         quantile_0.99="Cold Edge"))
# write out all the edge species used in the analysis, with species names, regions, and edge type 
write_csv(dat.summary, here("results","edge_species.csv"))

# read in results from Bayesian models 
spp.bayes.edge.lm.df.summary <- read_csv(here("results","species_edge_shifts_vs_time.csv")) # edge shifts over time
dat.predict.niche <- read_csv(here("processed-data","species_thermal_niche_v_time.csv")) # niche data
spp.bayes.niche.lm.stats <- read_csv(here("results","species_bayes_niche_lm_summary.csv")) # niche shifts over time 

################
### summarize niche shifts 
################ 

# did edge thermal niche shift over time? (using 90% Bayesian credible intervals)
spp.bayes.niche.results <- spp.bayes.niche.lm.stats %>%
  rowwise() %>% 
  mutate(crosses0 = ifelse(lower<0 & upper>0, TRUE, FALSE)) %>% # for every row, identify credible intervals that include zero
  group_by(species, region, quantile) %>%
  mutate(varTracked = ifelse(!TRUE %in% crosses0, "none", 
                             ifelse(!FALSE %in% crosses0, "both",
                             "one"))) %>% # classify range edges by whether they tracked both temperature extremes, neither, or just one--if just one, still need to paste through which one it is 
  rowwise() %>%
  mutate(varTracked = ifelse(varTracked=="one" & crosses0==TRUE, paste0(predicted.var),
                             ifelse(varTracked=="one" & crosses0==FALSE, NA, varTracked))) %>% # edit column to copy over which "predicted.var" has a credible interval that crossed 0, meaning the range edge remained at the same values of that temperature metric over time 
  group_by(species, region, quantile) %>%
  arrange(varTracked) %>%
  fill(varTracked) %>% # fill in NAs so each row shows the tracked temperature 
  ungroup() %>%
  group_by(species, region, quantile, varTracked) %>% 
  summarise()
# save list of range edges with which thermal extreme they tracked (cold, warm, both, or neither)
write_csv(spp.bayes.niche.results, here("results","edge_thermal_extreme_tracked_summary.csv"))
