# IDENTIFY WHICH SPECIES HAVE RANGE EDGES IN EACH REGION
# currently just using a distance buffer 
library(tidyverse)
library(magrittr)
library(here)

# choose buffer size and SE type from VAST output
km.buffer <- 150
# se.buffer <- 200

SEoptions <- c("relative","absolute")
SEtype <- SEoptions[1]

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
  filter(!quantile=="quantile_0.5") %>% # get rid of centroid 
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>% # harmonize name format 
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`)  

neus.vast <- if(SEtype=="relative"){readRDS(here("processed-data","neus_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","neus_absolute_SE_vast_edge_df.rds"))}

neus.df <- neus.vast %>%
  filter(!quantile=="quantile_0.5") %>%
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) 
                                                          
wc.vast <- if(SEtype=="relative"){readRDS(here("processed-data","wc_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","wc_absolute_SE_vast_edge_df.rds"))} 

wc.df <- wc.vast %>%
  filter(!quantile=="quantile_0.5") %>%
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) 

# get bounds of each region
ebs.maxnw <- max(ebs.df[ebs.df$axis=="NW_km",]$Estimate)
ebs.minnw <- min(ebs.df[ebs.df$axis=="NW_km",]$Estimate)

neus.maxUp <- max(neus.df[neus.df$axis=="coast_km",]$Estimate)
neus.minUp <- min(neus.df[neus.df$axis=="coast_km",]$Estimate)

wc.maxUp <- max(wc.df[wc.df$axis=="coast_km",]$Estimate)
wc.minUp <- min(wc.df[wc.df$axis=="coast_km",]$Estimate)

# identify range edge species along relevant axis 
ebs.nw.edges <- ebs.df %>%
  filter(axis=="NW_km",
         quantile%in%c("quantile_0.95","quantile_0.05")) %>%
  group_by(species, quantile) %>%
  mutate(meanedge = mean(Estimate)) %>%
  ungroup() %>%
  select(meanedge, species, quantile) %>%
  distinct() %>%
  rowwise() %>%
  mutate(edgecheck = between(meanedge, ebs.minnw+km.buffer, ebs.maxnw-km.buffer)) %>%
  filter(edgecheck==TRUE)

neus.coast.edges <- neus.df %>%
  filter(axis=="coast_km",
         quantile%in%c("quantile_0.95","quantile_0.05")) %>%
  group_by(species, quantile) %>%
  mutate(meanedge = mean(Estimate)) %>%
  ungroup() %>%
  select(meanedge, species, quantile) %>%
  distinct() %>%
  rowwise() %>%
  mutate(edgecheck = between(meanedge, neus.minUp+km.buffer, neus.maxUp-km.buffer)) %>%
  filter(edgecheck==TRUE)

wc.coast.edges <- wc.df %>%
  filter(axis=="coast_km",
         quantile%in%c("quantile_0.95","quantile_0.05")) %>%
  group_by(species, quantile) %>%
  mutate(meanedge = mean(Estimate)) %>%
  ungroup() %>%
  select(meanedge, species, quantile) %>%
  distinct() %>%
  rowwise() %>%
  mutate(edgecheck = between(meanedge, wc.minUp+km.buffer, wc.maxUp-km.buffer)) %>%
  filter(edgecheck==TRUE)

# which species have which edges?
wc.coast.pol.spp <- wc.coast.edges %>%
  filter(quantile=="quantile_0.95") %>%
  pull(unique(species)) 

wc.coast.eq.spp <- wc.coast.edges %>%
  filter(quantile=="quantile_0.05") %>%
  pull(unique(species))

neus.coast.pol.spp <- neus.coast.edges %>%
  filter(quantile=="quantile_0.95") %>%
  pull(unique(species))

neus.coast.eq.spp <- neus.coast.edges %>%
  filter(quantile=="quantile_0.05") %>%
  pull(unique(species))

ebs.nw.pol.spp <- ebs.nw.edges %>%
  filter(quantile=="quantile_0.95") %>%
  pull(unique(species))

ebs.nw.eq.spp <- ebs.nw.edges %>%
  filter(quantile=="quantile_0.05") %>%
  pull(unique(species))

# label edges by edge type, and keep only estimates of the relevant quantile (e.g., throw out 0.05 quantile for cold edges)
ebs.prep <- ebs.df %>% 
  mutate(region="ebs") %>% # add columns for binding later 
  filter(!quantile=="quantile_0.5") %>%
  rowwise() %>%
  mutate(edgetype=ifelse(species %in% ebs.nw.pol.spp & species %in% ebs.nw.eq.spp, "both", ifelse(species %in% ebs.nw.pol.spp, "coldedge", ifelse(species%in% ebs.nw.eq.spp, "warmedge", "neither")))) %>%
  filter(!edgetype=="neither") %>% # keep only rows with estimates from the range edge, whichever it is
  mutate(ref = ifelse(edgetype=="coldedge" & quantile=="quantile_0.95", "keep", ifelse(edgetype=="warmedge" & quantile=="quantile_0.05","keep",ifelse(edgetype=="both","keep","drop")))) %>%
  filter(ref=="keep") %>%
  select(-ref)

neus.prep <- neus.df %>% 
  mutate(region="neus") %>%
  filter(!quantile=="quantile_0.5") %>%
  rowwise() %>%
  mutate(edgetype=ifelse(species %in% neus.coast.pol.spp & species %in% neus.coast.eq.spp, "both", ifelse(species %in% neus.coast.pol.spp, "coldedge", ifelse(species%in% neus.coast.eq.spp, "warmedge", "neither")))) %>%
  filter(!edgetype=="neither") %>% 
  mutate(ref = ifelse(edgetype=="coldedge" & quantile=="quantile_0.95", "keep", ifelse(edgetype=="warmedge" & quantile=="quantile_0.05","keep",ifelse(edgetype=="both","keep","drop")))) %>% 
  filter(ref=="keep") %>%
  select(-ref)

wc.prep <- wc.df %>% 
  mutate(region="wc")%>%
  filter(!quantile=="quantile_0.5") %>%
  rowwise() %>%
  mutate(edgetype=ifelse(species %in% wc.coast.pol.spp & species %in% wc.coast.eq.spp, "both", ifelse(species %in% wc.coast.pol.spp, "coldedge", ifelse(species%in% wc.coast.eq.spp, "warmedge", "neither")))) %>%
  filter(!edgetype=="neither") %>% 
  mutate(ref = ifelse(edgetype=="coldedge" & quantile=="quantile_0.95", "keep", ifelse(edgetype=="warmedge" & quantile=="quantile_0.05","keep",ifelse(edgetype=="both","keep","drop")))) %>%
  filter(ref=="keep") %>%
  select(-ref)

dat <- rbind(wc.prep, neus.prep, ebs.prep) %>%
  left_join(spp.taxonomy, by=c("species"="query")) %>%
  mutate(taxongroup = ifelse(Class %in% c("Actinopterygii","Elasmobranchii"),"fish","invertebrates"))

# try Wald test (Jim's idea)
# need to expand all combinations of years 
dat_complete <- dat %>%
  select(-taxongroup, -edgetype, -Class) %>%
  group_by(quantile, species, region, axis) %>%
  mutate(year_compare = year,
         Estimate_compare = Estimate,
         Std.Error_compare = Std.Error
  ) %>%
  complete(year, year_compare) %>%
  ungroup() %>% # this produces rows with the full factorial combinations of all years, but doesn't fill in the estimate and SE columns 
group_by(quantile, species, region, axis, year)  %>%
  fill(c(Estimate, Std.Error), .direction="downup") %>% # the direction command tells it to look anywhere in the grouped df for the value--could be above or below the row in question
  ungroup() %>%
  group_by(quantile, species, region, axis, year_compare)  %>%
  fill(c(Estimate_compare, Std.Error_compare), .direction="downup") %>%  
  ungroup() %>%
  filter(!year==year_compare) %>% # get rid of rows where we're comparing the same years 
  rowwise() %>%
  mutate(mu.diff = Estimate - Estimate_compare, 
         se.diff = sqrt(Std.Error^2+Std.Error_compare^2),
         Zscore = mu.diff / se.diff,
         p.val = 2*(1-pnorm(abs(Zscore)))
         ) %>% # use Wald test to get a p-value for each set of pairs 
  group_by(quantile, species, region, axis) %>%
  mutate(min.pval = min(p.val),
    wald.signif = ifelse(min(p.val)<=0.1, "Y","N")) %>%
  ungroup() %>%
  select(quantile, species, region, axis, wald.signif, min.pval) %>%
  distinct()

dat %<>% left_join(dat_complete)

saveRDS(dat, here("processed-data","all_edge_spp_df.rds"))
rm(list=ls())
