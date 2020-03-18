# IDENTIFY WHICH SPECIES HAVE RANGE EDGES IN EACH REGION
# currently just using a distance buffer 
library(tidyverse)
library(here)

# choose buffer size and SE type from VAST output
km.buffer <- 150

SEoptions <- c("relative","absolute")
SEtype <- SEoptions[1]

##################
# import and harmonize species dfs
##################

# spp.classes <- read_csv(here("processed-data","spp_taxonomy.csv")) %>%
#   select(query, Class) %>%
#   distinct()

# reformat VAST output 
ebs.df <- if(SEtype=="relative"){readRDS(here("processed-data","ebs_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","ebs_absolute_SE_vast_edge_df.rds"))} %>%
  filter(!quantile=="quantile_0.5") %>% # get rid of centroid 
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>% # harmonize name format 
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`)

neus.df <- if(SEtype=="relative"){readRDS(here("processed-data","neus_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","neus_absolute_SE_vast_edge_df.rds"))} %>%
  filter(!quantile=="quantile_0.5") %>%
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`)

testplot1 <- neus.df %>% 
  filter(species==tmp$species[1]) %>%
  mutate(axis=as.character(axis)) %>%
  filter(axis=="coast_km") %>% 
  ggplot(aes(x=year, y=Estimate))+ 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(x=year, ymin=Estimate-Std.Error, ymax=Estimate+Std.Error))
                                                          
wc.df <- readRDS(here("processed-data","wc_vast_edge_df.rds")) %>%
  filter(!quantile=="quantile_0.5") %>%
  mutate(species = gsub("_", " ", species),
         species=tolower(species)) %>%
  pivot_wider(names_from=quantity, values_from=value ) %>%
  rename("Std.Error"=`Std. Error`) 

# identify range edge species
ebs.maxnw <- max(ebs.df[ebs.df$axis=="NW_km",]$Estimate)
ebs.minnw <- min(ebs.df[ebs.df$axis=="NW_km",]$Estimate)

neus.maxUp <- max(neus.df[neus.df$axis=="coast_km",]$Estimate)
neus.minUp <- min(neus.df[neus.df$axis=="coast_km",]$Estimate)

wc.maxUp <- max(wc.df[wc.df$axis=="coast_km",]$Estimate)
wc.minUp <- min(wc.df[wc.df$axis=="coast_km",]$Estimate)

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
  filter(!edgetype=="neither") %>% # keep only rows with estimates from the range edge, whichever it is
  mutate(ref = ifelse(edgetype=="coldedge" & quantile=="quantile_0.95", "keep", ifelse(edgetype=="warmedge" & quantile=="quantile_0.05","keep",ifelse(edgetype=="both","keep","drop")))) %>% # keep only quantiles matching the relevant edges 
  filter(ref=="keep") %>%
  select(-ref)

wc.prep <- wc.df %>% 
  mutate(region="wc")%>%
  filter(!quantile=="quantile_0.5") %>%
  rowwise() %>%
  mutate(edgetype=ifelse(species %in% wc.coast.pol.spp & species %in% wc.coast.eq.spp, "both", ifelse(species %in% wc.coast.pol.spp, "coldedge", ifelse(species%in% wc.coast.eq.spp, "warmedge", "neither")))) %>%
  filter(!edgetype=="neither") %>% # keep only rows with estimates from the range edge, whichever it is
  mutate(ref = ifelse(edgetype=="coldedge" & quantile=="quantile_0.95", "keep", ifelse(edgetype=="warmedge" & quantile=="quantile_0.05","keep",ifelse(edgetype=="both","keep","drop")))) %>%
  filter(ref=="keep") %>%
  select(-ref)

dat <- rbind(wc.prep, neus.prep, ebs.prep) %>%
  left_join(spp.taxonomy, by=c("species"="query")) %>%
  mutate(taxongroup = ifelse(class %in% c("Teleostei","Chondrichthyes"),"fish","invertebrates"))

saveRDS(dat, here("processed-data","all_edge_spp_df.rds"))
rm(list=ls())
