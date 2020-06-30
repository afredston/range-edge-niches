library(tidyverse)
library(here)
library(purrr)
library(broom)

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
lm_min_sst # no significant increases anywhere 

lm_mean_sst <- bind_rows(neus_sst_summary, wc_sst_summary, ebs_sst_summary) %>% 
  group_by(region) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~lm(region_mean_mean_sst ~ year, data=.x)),
    tidymodel = purrr::map(model, tidy)
  ) %>%
  unnest(tidymodel) %>%
  select(-data, -model) 
lm_mean_sst # significant increase in EBS and (especially) NEUS

lm_max_sst <- bind_rows(neus_sst_summary, wc_sst_summary, ebs_sst_summary) %>% 
  group_by(region) %>%
  nest() %>%
  mutate(
    model = purrr::map(data, ~lm(region_mean_max_sst ~ year, data=.x)),
    tidymodel = purrr::map(model, tidy)
  ) %>%
  unnest(tidymodel) %>%
  select(-data, -model) 
lm_max_sst # significant increase in NEUS, EBS
