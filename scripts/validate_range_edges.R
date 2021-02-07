# IDENTIFY WHICH SPECIES HAVE RANGE EDGES IN EACH REGION, AND ADD TRAITS
# currently just using a distance buffer that is X% of the total measurement axis length (edge must fall in middle 80%)
library(tidyverse)
library(magrittr)
library(here)

buffer_amount <- 0.10 # what fraction of the total axis length is "too close" to the region edge on either side? if 0.05, range edges need to be in the middle 90% of the region 
SEoptions <- c("relative","absolute")
SEtype <- SEoptions[1]
quantiles_to_use <- c('quantile_0.01','quantile_0.99')

##################
# import and harmonize species dfs
##################

# get taxonomy to bind to edge df at the end
spp.taxonomy <- read_csv(here("processed-data","spp_taxonomy.csv")) %>%
  select(query, Class) %>%
  distinct()

# reformat VAST output 
ebs.vast <- if(SEtype=="relative"){readRDS(here("processed-data","ebs_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","ebs_absolute_SE_vast_edge_df.rds"))} 

ebs.df <- ebs.vast %>%
  filter(quantile %in% quantiles_to_use) %>% # get rid of centroid, other quantiles
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>% # harmonize name format 
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`)  

neus.vast <- if(SEtype=="relative"){readRDS(here("processed-data","neus_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","neus_absolute_SE_vast_edge_df.rds"))}

neus.df <- neus.vast %>%
  filter(quantile %in% quantiles_to_use) %>%
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) 
                                                          
wc.vast <- if(SEtype=="relative"){readRDS(here("processed-data","wc_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","wc_absolute_SE_vast_edge_df.rds"))} 

wc.df <- wc.vast %>%
  filter(quantile %in% quantiles_to_use) %>%
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) 

# get bounds of each region
ebs.maxnw <- max(ebs.df[ebs.df$axis=="line_km",]$Estimate)
ebs.minnw <- min(ebs.df[ebs.df$axis=="line_km",]$Estimate)

ebs.axislen = ebs.maxnw-ebs.minnw
ebs.buffer = ebs.axislen*buffer_amount

neus.maxUp <- max(neus.df[neus.df$axis=="coast_km",]$Estimate)
neus.minUp <- min(neus.df[neus.df$axis=="coast_km",]$Estimate)

neus.axislen = neus.maxUp-neus.minUp
neus.buffer = neus.axislen*buffer_amount

wc.maxUp <- max(wc.df[wc.df$axis=="coast_km",]$Estimate)
wc.minUp <- min(wc.df[wc.df$axis=="coast_km",]$Estimate)

wc.axislen = wc.maxUp-wc.minUp
wc.buffer = wc.axislen*buffer_amount

# identify range edge species along relevant axis 
ebs.nw.edges <- ebs.df %>%
  filter(axis=="line_km",
         quantile %in% quantiles_to_use) %>%
  group_by(species, quantile) %>%
  mutate(meanedge = mean(Estimate)) %>%
  ungroup() %>%
  select(meanedge, species, quantile) %>%
  distinct() %>%
  rowwise() %>%
  mutate(edgecheck = between(meanedge, ebs.minnw+ebs.buffer, ebs.maxnw-ebs.buffer)) %>%
  filter(edgecheck==TRUE)

neus.coast.edges <- neus.df %>%
  filter(axis=="coast_km",
         quantile %in% quantiles_to_use) %>%
  group_by(species, quantile) %>%
  mutate(meanedge = mean(Estimate)) %>%
  ungroup() %>%
  select(meanedge, species, quantile) %>%
  distinct() %>%
  rowwise() %>%
  mutate(edgecheck = between(meanedge, neus.minUp+neus.buffer, neus.maxUp-neus.buffer)) %>%
  filter(edgecheck==TRUE)

wc.coast.edges <- wc.df %>%
  filter(axis=="coast_km",
         quantile %in% quantiles_to_use) %>%
  group_by(species, quantile) %>%
  mutate(meanedge = mean(Estimate)) %>%
  ungroup() %>%
  select(meanedge, species, quantile) %>%
  distinct() %>%
  rowwise() %>%
  mutate(edgecheck = between(meanedge, wc.minUp+wc.buffer, wc.maxUp-wc.buffer)) %>%
  filter(edgecheck==TRUE)

# which species have which edges?
wc.coast.pol.spp <- wc.coast.edges %>%
  filter(quantile==quantiles_to_use[2]) %>%
  pull(unique(species)) 

wc.coast.eq.spp <- wc.coast.edges %>%
  filter(quantile==quantiles_to_use[1]) %>%
  pull(unique(species))

neus.coast.pol.spp <- neus.coast.edges %>%
  filter(quantile==quantiles_to_use[2]) %>%
  pull(unique(species))

neus.coast.eq.spp <- neus.coast.edges %>%
  filter(quantile==quantiles_to_use[1]) %>%
  pull(unique(species))

ebs.nw.pol.spp <- ebs.nw.edges %>%
  filter(quantile==quantiles_to_use[2]) %>%
  pull(unique(species))

ebs.nw.eq.spp <- ebs.nw.edges %>%
  filter(quantile==quantiles_to_use[1]) %>%
  pull(unique(species))

# label edges by edge type, and keep only estimates of the relevant quantile (e.g., throw out 0.05 quantile for cold edges)
ebs.prep <- ebs.df %>% 
  mutate(region="ebs") %>% # add columns for binding later 
  rowwise() %>%
  mutate(edgetype=ifelse(species %in% ebs.nw.pol.spp & species %in% ebs.nw.eq.spp, "both", ifelse(species %in% ebs.nw.pol.spp, "coldedge", ifelse(species%in% ebs.nw.eq.spp, "warmedge", "neither")))) %>%
  filter(!edgetype=="neither") %>% # keep only rows with estimates from the range edge, whichever it is
  mutate(ref = ifelse(edgetype=="coldedge" & quantile==quantiles_to_use[2], "keep", ifelse(edgetype=="warmedge" & quantile==quantiles_to_use[1],"keep",ifelse(edgetype=="both","keep","drop")))) %>% # match the edge type with the correct quantile--we did calculate 0.01 and 0.99 quantiles for each range, but only want to keep 0.01 for warm edges and 0.99 for cold edges 
  filter(ref=="keep") %>%
  select(-ref)

neus.prep <- neus.df %>% 
  mutate(region="neus") %>%
  rowwise() %>%
  mutate(edgetype=ifelse(species %in% neus.coast.pol.spp & species %in% neus.coast.eq.spp, "both", ifelse(species %in% neus.coast.pol.spp, "coldedge", ifelse(species%in% neus.coast.eq.spp, "warmedge", "neither")))) %>%
  filter(!edgetype=="neither") %>% 
  mutate(ref = ifelse(edgetype=="coldedge" & quantile==quantiles_to_use[2], "keep", ifelse(edgetype=="warmedge" & quantile==quantiles_to_use[1],"keep",ifelse(edgetype=="both","keep","drop")))) %>%
  filter(ref=="keep") %>%
  select(-ref)

wc.prep <- wc.df %>% 
  mutate(region="wc")%>%
  rowwise() %>%
  mutate(edgetype=ifelse(species %in% wc.coast.pol.spp & species %in% wc.coast.eq.spp, "both", ifelse(species %in% wc.coast.pol.spp, "coldedge", ifelse(species%in% wc.coast.eq.spp, "warmedge", "neither")))) %>%
  filter(!edgetype=="neither") %>% 
  mutate(ref = ifelse(edgetype=="coldedge" & quantile==quantiles_to_use[2], "keep", ifelse(edgetype=="warmedge" & quantile==quantiles_to_use[1],"keep",ifelse(edgetype=="both","keep","drop")))) %>%
  filter(ref=="keep") %>%
  select(-ref)

dat <- rbind(wc.prep, neus.prep, ebs.prep) %>%
  left_join(spp.taxonomy, by=c("species"="query")) %>%
  mutate(taxongroup = ifelse(Class %in% c("Actinopterygii","Elasmobranchii"),"fish","invertebrates")) 

# pull in fish traits from Beukhof et al. 2019
# data: https://doi.org/10.1594/PANGAEA.900866
# I downloaded the Excel file and saved the first sheet as a .csv for import into R

traits_raw <- read_csv(here("raw-data","beukhof_2019_traits.csv")) %>% 
  mutate(taxon = str_to_lower(taxon)) %>% # for matching with our dataset 
  select(taxon, FAO, habitat, feeding.mode, tl, offspring.size, spawning.type, fecundity, age.maturity, length.infinity, age.max) %>% # keep only traits of interest 
  distinct()

# join by species name and FAO region 
traits_dat <- dat %>% 
  select(species, region, edgetype) %>% 
  distinct() %>% 
  mutate(FAO = ifelse(region == "neus", 21, # some species from region 31 could also end up in NEUS, but I checked the Beukhof data and all of those species are also in their list from 21
                      ifelse(region == "wc", 77, 67))) %>% 
  inner_join(traits_raw, by=c("FAO","species"="taxon")) %>% 
  group_by(species, FAO, region, edgetype, habitat, feeding.mode, spawning.type) %>% 
  summarise(across(tl:age.max, mean))# you should check at this stage that there are no longer multiple rows for species*FAO combos

dat <- dat %>% 
  left_join(traits_dat, by=c('species','edgetype','region'))


# try Wald test (Jim's idea)
# need to expand all combinations of years 
# dat_complete <- dat %>%
#   select(-taxongroup, -edgetype, -Class) %>%
#   group_by(quantile, species, region, axis) %>%
#   mutate(year_compare = year,
#          Estimate_compare = Estimate,
#          Std.Error_compare = Std.Error
#   ) %>%
#   complete(year, year_compare) %>%
#   ungroup() %>% # this produces rows with the full factorial combinations of all years, but doesn't fill in the estimate and SE columns 
# group_by(quantile, species, region, axis, year)  %>%
#   fill(c(Estimate, Std.Error), .direction="downup") %>% # the direction command tells it to look anywhere in the grouped df for the value--could be above or below the row in question
#   ungroup() %>%
#   group_by(quantile, species, region, axis, year_compare)  %>%
#   fill(c(Estimate_compare, Std.Error_compare), .direction="downup") %>%  
#   ungroup() %>%
#   filter(!year==year_compare) %>% # get rid of rows where we're comparing the same years 
#   rowwise() %>%
#   mutate(mu.diff = Estimate - Estimate_compare, 
#          se.diff = sqrt(Std.Error^2+Std.Error_compare^2),
#          Zscore = mu.diff / se.diff,
#          p.val = 2*(1-pnorm(abs(Zscore)))
#          ) %>% # use Wald test to get a p-value for each set of pairs 
#   group_by(quantile, species, region, axis) %>%
#   mutate(min.pval = min(p.val),
#     wald.signif = ifelse(min(p.val)<=0.1, "Y","N")) %>%
#   ungroup() %>%
#   select(quantile, species, region, axis, wald.signif, min.pval) %>%
#   distinct()
# 
# dat %<>% left_join(dat_complete)

saveRDS(dat, here("processed-data","all_edge_spp_df.rds"))
rm(list=ls())
