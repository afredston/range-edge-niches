# this script collates all statistics reported in-text in the manuscript. it is ordered to approximately parallel the order in which these statistics are described. 

library(tidyverse)
library(here)
library(purrr)
library(broom)

#################
# axis length 
#################
neus.coastdistdat <- readRDS(here("processed-data","neus_coastdistdat.rds"))
wc.coastdistdat <- readRDS(here("processed-data","wc_coastdistdat.rds"))
ebs.axisdistdat <- readRDS(here("processed-data","ebs_axisdistdat.rds"))

neus.len <- (max(neus.coastdistdat$lengthfromhere)-min(neus.coastdistdat$lengthfromhere))/1000 # subtract largest distance from smallest distance (need the latter because origin is not at exactly zero); change from m to km
wc.len <- (max(wc.coastdistdat$lengthfromhere)-min(wc.coastdistdat$lengthfromhere))/1000 
ebs.len <- (max(ebs.axisdistdat$lengthfromhere)-min(ebs.axisdistdat$lengthfromhere))/1000

#################
# range edge filters 
#################
# all species with models fitted by VAST: 

dat.models <- readRDS(here("processed-data","all_edge_spp_df.rds")) # just using this name to be consistent with the rest of the scripts 

#################
# summary of study species
#################

# how many species converged in VAST? 
vast.spp <- read_csv(here("processed-data","spp_taxonomy.csv"))

# breakdown by region
vast.spp %>% group_by(region) %>% summarise(n=n())

# how many found in >1 region? 
vast.spp %>% group_by(query) %>% summarise(n=n()) %>% filter(n>1) 

# how many unique species?
length(unique(vast.spp$query))

# how many actually had range edges?
dat.summary <- readRDS(here("processed-data","all_edge_spp_df.rds")) %>%
  ungroup() %>% # undo rowwise nature
  mutate(axis = as.character(axis)) %>% # convert from factor
  filter(axis %in% c('coast_km','line_km'))  %>%
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
                         quantile_0.99="Cold Edge")) %>%
  left_join(vast.spp %>% select(Species, Class) %>% distinct(), by=c("species"="Species"))

# break down by region (and compare to same test above to see how many spp were eliminated for not having edges)
dat.summary %>% group_by(region) %>% summarise(n=n())

# write out for appendix
write_csv(dat.summary, here("results","species_summary_for_supplement.csv"))

#################
# change in temperature over time
#################

# load all temperature datasets 
neus_sst <- read_rds(here("processed-data","neus_sst_corrected.rds"))
wc_sst <- read_rds(here("processed-data","wc_sst_corrected.rds"))
ebs_sst <- read_rds(here("processed-data","ebs_sst_corrected.rds"))

# note that temperature is calculated by: 1. getting monthly means for each cell 2. getting the means and extremes for a year for each cell (hottest, coldest, and average months) 3. calculating the average hottest, coldest, and average months across all cells in the region every year. this is designed to capture extreme temperature variation without being too sensitive to daily or cell-specific variation. in other words, it is measuring if hot months got hotter, cold months got less cold, or average months are creeping up in SST *on average across the region*. 

neus_sst_summary <- neus_sst %>%
  drop_na() %>% # tiny number of cells have NAs for some reason, just drop them here 
  group_by(x, y, year) %>%
  mutate(cell.min = min(sst),
         cell.mean = mean(sst),
         cell.max = max(sst)) %>%
  ungroup() %>%
  select(x, y, year, cell.min, cell.mean, cell.max) %>%
  distinct()%>%
  group_by(year) %>%
  mutate(region_mean_min_sst = mean(cell.min),
         region_mean_mean_sst = mean(cell.mean),
         region_mean_max_sst = mean(cell.max)) %>%
  ungroup() %>%
  select(year, region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst) %>%
  distinct() %>%
  mutate(region="neus")

wc_sst_summary <- wc_sst %>%
  drop_na() %>%
  group_by(x, y, year) %>%
  mutate(cell.min = min(sst),
         cell.mean = mean(sst),
         cell.max = max(sst)) %>%
  ungroup() %>%
  select(x, y, year, cell.min, cell.mean, cell.max) %>%
  distinct()%>%
  group_by(year) %>%
  mutate(region_mean_min_sst = mean(cell.min),
         region_mean_mean_sst = mean(cell.mean),
         region_mean_max_sst = mean(cell.max)) %>%
  ungroup() %>%
  select(year, region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst) %>%
  distinct() %>%
  mutate(region="wc")

ebs_sst_summary <- ebs_sst %>%
  drop_na() %>%
  group_by(x, y, year) %>%
  summarise(cell.min = min(sst),
            cell.mean = mean(sst),
            cell.max = max(sst)) %>%
  group_by(year) %>%
  mutate(region_mean_min_sst = mean(cell.min),
         region_mean_mean_sst = mean(cell.mean),
         region_mean_max_sst = mean(cell.max)) %>%
  ungroup() %>%
  select(year, region_mean_min_sst, region_mean_mean_sst, region_mean_max_sst) %>%
  distinct() %>%
  mutate(region="ebs")

# estimate trends in each temperature value with simple linear regressions

lm_min_sst <- bind_rows(neus_sst_summary, wc_sst_summary, ebs_sst_summary) %>% 
  group_by(region) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~lm(region_mean_min_sst ~ year, data=.x)),
    tidymodel = purrr::map(model, tidy)
  ) %>%
  unnest(tidymodel) %>%
  select(-data, -model) 
lm_min_sst 

lm_mean_sst <- bind_rows(neus_sst_summary, wc_sst_summary, ebs_sst_summary) %>% 
  group_by(region) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~lm(region_mean_mean_sst ~ year, data=.x)),
    tidymodel = purrr::map(model, tidy)
  ) %>%
  unnest(tidymodel) %>%
  select(-data, -model) 
lm_mean_sst  

lm_max_sst <- bind_rows(neus_sst_summary, wc_sst_summary, ebs_sst_summary) %>% 
  group_by(region) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~lm(region_mean_max_sst ~ year, data=.x)),
    tidymodel = purrr::map(model, tidy)
  ) %>%
  unnest(tidymodel) %>%
  select(-data, -model) 
lm_max_sst 

#################
# range shift stats
#################

spp.bayes.lm.df.summary <- read_csv(here("results","species_edge_shifts_vs_time.csv")) 

# most extreme shifts 
spp.bayes.lm.df.summary[spp.bayes.lm.df.summary$mean==max(spp.bayes.lm.df.summary$mean),]
spp.bayes.lm.df.summary[spp.bayes.lm.df.summary$mean==min(spp.bayes.lm.df.summary$mean),]

#################
# edge thermal niche stats
#################

spp.bayes.niche.results <- read_csv(here("results","edge_thermal_extreme_tracked_summary.csv"))

# add in taxonomy

dat.models <- readRDS(here("processed-data","all_edge_spp_df.rds")) %>%
  ungroup() %>% # undo rowwise nature
  mutate(axis = as.character(axis)) %>% # convert from factor
  filter(axis %in% c('coast_km','line_km')) %>%
  select(species, region, Class, taxongroup) %>%
  distinct()

spp.bayes.niche.groups <- spp.bayes.niche.results %>% 
  left_join(dat.models)

# calculate stats
spp.bayes.niche.groups %>% group_by(varTracked) %>% summarise(n=n()) # break down by niche grouping
spp.bayes.niche.groups %>% group_by(varTracked, region) %>% summarise(n=n())
