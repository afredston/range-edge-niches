library(tidyverse)
library(here)
dat.predict.lag <- read_csv(here("processed-data","species_thermal_niche_v_time_with_lags.csv")) 

lagdf <- dat.predict.lag %>% 
  group_by(species, region, quantile, predicted.var) %>% 
  filter(year_match > min(year_match)) %>% 
  arrange(year_match) %>% 
  mutate(sst_lag = lag(sst), # get SST last year at last year's edge 
         sst_diff = sst_at_last_years_edge - sst_lag, # calculate change in SST at last year's edge: temperature there this year - temperature there last year 
         edge_position_diff = edge_position - edge_position_last_year # calculate shift 
  ) %>% 
  ungroup()

# plot of all species
gglag <- lagdf %>% 
  mutate(
    region=recode(region,
                  ebs="Eastern Bering Sea",
                  neus="Northeast",
                  wc="West Coast"),
    predicted.var=recode(predicted.var,
                         predict.sstmax="Summer SST",
                         predict.sstmin="Winter SST")) %>% 
  ggplot(aes(x=sst_diff, y=edge_position_diff)) + 
  geom_point() + 
  theme_bw() +
  facet_grid(predicted.var~region) +
  labs(x="SST Change at Last Year's Edge (°C)", y="One-Year Change in Edge Position (km)") +
  NULL

ggsave(gglag, width=6, height=3.5, filename=here("results","year_over_year_edge_vs_sst.png"), dpi=160, scale=1.3)

# plot of some example species

exspp <- c('dipturus laevis','stenotomus chrysops','sebastes pinniger','sebastes semicinctus','chionoecetes opilio','limanda aspera')

ggexspp <- lagdf %>% 
  filter(species %in% exspp) %>% 
  mutate(
    region=recode(region,
                  ebs="Eastern Bering Sea",
                  neus="Northeast",
                  wc="West Coast"),
    predicted.var=recode(predicted.var,
                         predict.sstmax="Summer SST",
                         predict.sstmin="Winter SST"),
    species = str_to_sentence(species)) %>% 
  ggplot(aes(x=sst_diff, y=edge_position_diff)) + 
  geom_point() + 
  theme_bw() +
  facet_wrap(~species) +
  labs(x="SST Change at Last Year's Edge (°C)", y="One-Year Change in Edge Position (km)") +
  NULL

ggsave(ggexspp, width=6, height=3.5, filename=here("results","year_over_year_edge_vs_sst_examples.png"), dpi=160, scale=1.3)
