# get higher taxonomy of species found in trawl surveys
# because of how taxize works with WORMS, you can't just hit run on this script--it needs manual input when there are multiple taxonomic records
# I just select the one labeled "accepted" when prompted

## load packages and data 
library(tidyverse)
library(here)
library(taxize)

SEoptions <- c("relative","absolute")
SEtype <- SEoptions[1]

ebs.vast <- if(SEtype=="relative"){readRDS(here("processed-data","ebs_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","ebs_absolute_SE_vast_edge_df.rds"))} 

neus.vast <- if(SEtype=="relative"){readRDS(here("processed-data","neus_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","neus_absolute_SE_vast_edge_df.rds"))}

wc.vast <- if(SEtype=="relative"){readRDS(here("processed-data","wc_relative_SE_vast_edge_df.rds"))
}else{readRDS(here("processed-data","wc_absolute_SE_vast_edge_df.rds"))} 

neus.spplist <- c(unique(neus.vast$species)) %>% tolower() %>% str_replace_all("_", " ")
wc.spplist <- c(unique(wc.vast$species)) %>% tolower() %>% str_replace_all("_", " ")
ebs.spplist <- c(unique(ebs.vast$species)) %>% tolower() %>% str_replace_all("_", " ")

# get higher taxonomy with taxize

neus.tax <- rbind(classification(neus.spplist, db="worms")) %>%
  filter(rank %in% c('Phylum','Class','Order','Family','Genus','Species')) %>%
  dplyr::select(-id) %>%
  group_by(query) %>%
  spread(key=rank, value=name) %>%
  ungroup()%>%
  mutate(region="Northeast")

wc.tax <- rbind(classification(wc.spplist, db="worms")) %>%
  filter(rank %in% c('Phylum','Class','Order','Family','Genus','Species')) %>%
  dplyr::select(-id) %>%
  group_by(query) %>%
  spread(key=rank, value=name) %>%
  ungroup()%>%
  mutate(region="West Coast")

ebs.tax <- rbind(classification(ebs.spplist, db="worms")) %>%
  filter(rank %in% c('Phylum','Class','Order','Family','Genus','Species')) %>%
  dplyr::select(-id) %>%
   group_by(query) %>%
  spread(key=rank, value=name) %>%
  ungroup() %>%
  mutate(region="Eastern Bering Sea")

spp.taxonomy <- rbind(neus.tax, ebs.tax, wc.tax)
write_csv(spp.taxonomy, here("processed-data","spp_taxonomy.csv"))
rm(list=ls())
