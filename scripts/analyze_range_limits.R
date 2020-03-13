#######################
### load packages 
#######################

library(data.table)
library(tidyverse)
library(here)
library(purrr)
library(broom)
library(lmerTest)
library(mgcv)
library(rstanarm)
library(tidybayes)
#library(pwr)

#######################
### load data 
#######################

# NEED TO RE-RUN VALIDATE_EDGE_SPP AFTER RE-RUNNING VAST!
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
                         quantile_0.05="Warm Limit",
                         quantile_0.95="Cold Limit"))
write_csv(dat.summary, here("results","edge_species.csv"))

dat.models.groups <- dat.models %>%
  select(species, edgetype, region, taxongroup) %>%
  distinct() %>%
  mutate(species = as.factor(species), 
         edgetype=as.factor(edgetype),
         region=as.factor(region),
         taxongroup=as.factor(taxongroup))

# prep data for models

ebs.oisst.summary <- read_rds(here("processed-data","ebs_oisst_rotated.rds")) %>% 
  group_by(year_match, month, lat, lon) %>%
  mutate(sst.month.mean = mean(sst)) %>% # calculate monthly means
  ungroup() %>%
  group_by(year_match, lat, lon) %>%
  mutate(cell.annual.sst = mean(sst.month.mean)) %>% # calculate annual mean of monthly SSTs for each grid cell
  ungroup() %>%
  group_by(year_match) %>%
  mutate(mean.annual.sst = mean(cell.annual.sst)) %>% # calculate annual mean of all grid cells
  select(year_match, mean.annual.sst) %>%
  distinct() %>%
  mutate(region="ebs")

wc.oisst.summary <- read_rds(here("processed-data","wc.oisst.coastdist.rds")) %>% 
  rename(lon="x",lat="y") %>%
  group_by(year_match, month, lat, lon) %>%
  mutate(sst.month.mean = mean(sst)) %>% # calculate monthly means
  ungroup() %>%
  group_by(year_match, lat, lon) %>%
  mutate(cell.annual.sst = mean(sst.month.mean)) %>% # calculate annual mean of monthly SSTs for each grid cell
  ungroup() %>%
  group_by(year_match) %>%
  mutate(mean.annual.sst = mean(cell.annual.sst)) %>% # calculate annual mean of all grid cells
  select(year_match, mean.annual.sst) %>%
  distinct() %>%
  mutate(region="wc")

neus.hadisst.summary <- read_rds(here("processed-data","neus.hadisst.coastdist.rds")) %>% 
  rename(lon="x",lat="y") %>%
  group_by(year_match, lat, lon) %>%
  mutate(cell.annual.sst = mean(sst)) %>% # calculate annual mean of monthly SSTs for each grid cell
  ungroup() %>%
  group_by(year_match) %>%
  mutate(mean.annual.sst = mean(cell.annual.sst)) %>% # calculate annual mean of all grid cells
  select(year_match, mean.annual.sst) %>%
  distinct() %>%
  mutate(region="neus")

dat.sst <- rbind(neus.hadisst.summary, wc.oisst.summary, ebs.oisst.summary)

dat.edges.sst <- dat.models %>%
  left_join(dat.sst, by=c("region","year"="year_match"))

#######################
### single-species Bayesian models
#######################
# spp.lms <- dat.edges.sst %>%
#   filter(axis %in% c('coast_km','NW_km')) %>%
#   group_by(species, region, quantile, edgetype) %>%
#   nest() %>%
#   mutate(
#     model = purrr::map(data, ~lm(Estimate ~ mean.annual.sst, data = .x)), 
#     tidymodel = purrr::map(model, tidy)
#   ) %>% 
#   unnest(tidymodel) %>% 
#   select(-data, -model) %>%
#   filter(!term=="(Intercept)")

# species shifts vs T 
spp.bayes.sst.lm.df <- NULL
for(i in unique(dat.edges.sst$species)) {
  dfprep1 <- dat.edges.sst[dat.edges.sst$species==i,] # subdivide by species 
  for(j in unique(dfprep1$region)) {
    dfprep2 <- dat.edges.sst[dat.edges.sst$species==i & dat.edges.sst$region==j,] # subdivide by region, for the few spp found in multiple regions 
    for(k in unique(dfprep2$quantile)) {
      df <- dat.edges.sst[dat.edges.sst$species==i & dat.edges.sst$region==j & dat.edges.sst$quantile==k,] # subdivide by quantile, for the few spp with both range edges
      spp.bayes.lm <- try(stan_glm(Estimate ~ mean.annual.sst, 
                                   data=df, 
                                   family=gaussian(), 
                                   iter = 12000,
                                   warmup = 2000,
                                   adapt_delta = 0.95,
                                   chains = 4,
                                   cores = 1,
                                   prior = normal()))
      if(!class(spp.bayes.lm)[1] == "try-error") { # adding try() here because some edges are so invariant that the model fails 
        spp.bayes.lm.tidy <- tidy_draws(spp.bayes.lm) %>%
          mutate(species = paste0(i),
                 region = paste0(j),
                 quantile = paste0(k))
        spp.bayes.sst.lm.df <- rbind(spp.bayes.sst.lm.df, spp.bayes.lm.tidy)
      }
    }
  }
}

bayes.lm.sst.edgetype.gg <- spp.bayes.sst.lm.df  %>%
  group_by(.draw, quantile, region) %>%
  mutate(mean.param = mean(mean.annual.sst) ) %>%
  ungroup() %>% 
  select(.draw, mean.param, quantile, region) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=quantile), color="black", alpha=0.5) +
  scale_fill_brewer(type="seq", palette="YlGnBu") +
  labs(x="Coefficient of Edge vs Temperature (km/C)", y="Density", fill="Range Limit Type") +
  theme(legend.position=c(0.2, 0.8)) +
  facet_wrap(~region) +
  NULL
bayes.lm.sst.edgetype.gg
ggsave(bayes.lm.sst.edgetype.gg, filename=here("results","edge_coefficients_by_reg.png"), dpi=160, width=8, height=5)

bayes.lm.sst.gg <- spp.bayes.sst.lm.df  %>%
  group_by(.draw, quantile) %>%
  mutate(mean.param = mean(mean.annual.sst) ) %>%
  ungroup() %>% 
  select(.draw, mean.param, quantile) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=quantile), color="black", alpha=0.5) +
  scale_fill_brewer(type="seq", palette="YlGnBu") +
  labs(x="Coefficient of Edge vs Temperature (km/C)", y="Density", fill="Range Limit Type") +
  theme(legend.position=c(0.2, 0.8)) +
  NULL
bayes.lm.sst.gg
ggsave(bayes.lm.sst.gg, filename=here("results","edge_coefficients.png"), dpi=160, width=5, height=5)

# species shifts vs time

spp.bayes.lm.df <- NULL
for(i in unique(dat.models$species)) {
  dfprep1 <- dat.models[dat.models$species==i,] # subdivide by species 
  for(j in unique(dfprep1$region)) {
    dfprep2 <- dat.models[dat.models$species==i & dat.models$region==j,] # subdivide by region, for the few spp found in multiple regions 
    for(k in unique(dfprep2$quantile)) {
      df <- dat.models[dat.models$species==i & dat.models$region==j & dat.models$quantile==k,] # subdivide by quantile, for the few spp with both range edges
      spp.bayes.lm <- try(stan_glm(Estimate ~ year, 
                                   data=df, 
                                   family=gaussian(), 
                                   iter = 12000,
                                   warmup = 2000,
                                   adapt_delta = 0.95,
                                   chains = 4,
                                   cores = 1,
                                   prior = normal()))
      if(!class(spp.bayes.lm)[1] == "try-error") { # adding try() here because some edges are so invariant that the model fails 
        spp.bayes.lm.tidy <- tidy_draws(spp.bayes.lm) %>%
          mutate(species = paste0(i),
                 region = paste0(j),
                 quantile = paste0(k))
        spp.bayes.lm.df <- rbind(spp.bayes.lm.df, spp.bayes.lm.tidy)
      }
    }
  }
}
bayes.lm.time.gg <- spp.bayes.lm.df  %>%
  group_by(.draw, region) %>%
  mutate(mean.param = mean(year) ) %>%
  ungroup() %>% 
  select(.draw, mean.param, region) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=region), color="black", alpha=0.5) +
  scale_fill_brewer(type="seq", palette="YlGnBu", labels=c("Eastern Bering Sea","Northeast","West Coast")) +
  labs(x="Coefficient of Edge vs Time (km/yr)", y="Density", fill="Region") +
  theme(legend.position="bottom") +
  NULL
bayes.lm.time.gg
ggsave(bayes.lm.time.gg, width=3, height=4, dpi=160, filename=here("results","edge_coefficients_time.png"), scale=1.6)

bayes.lm.time.edgetype.gg <- spp.bayes.lm.df  %>%
  group_by(.draw, region, quantile) %>%
  mutate(mean.param = mean(year) ) %>%
  ungroup() %>% 
  select(.draw, mean.param, region, quantile) %>%
  distinct() %>%
  mutate(quantile=recode(quantile,
                         quantile_0.05="Warm Limit",
                         quantile_0.95="Cold Limit")) %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=region), color="black", alpha=0.5) +
  scale_fill_brewer(type="seq", palette="YlGnBu", labels=c("Eastern Bering Sea","Northeast","West Coast")) +
  labs(x="Coefficient of Edge vs Time (km/yr)", y="Density", fill="Region") +
  theme(legend.position="bottom") +
  facet_wrap(~quantile) +
  NULL
bayes.lm.time.edgetype.gg
ggsave(bayes.lm.time.edgetype.gg, width=6, height=4, dpi=160, filename=here("results","edge_coefficients_time_edgetype.png"), scale=1.6)

bayes.lm.time.stats <- spp.bayes.lm.df %>%
  group_by(.draw, region) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(region) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))

bayes.lm.time.edgetype.stats <- spp.bayes.lm.df %>%
  group_by(.draw, region, quantile) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(region, quantile) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))

#######################
### predict temperatures at edges 
#######################

# split by region because temp datasets are different

# set up temp data
ebs.oisst.rotated <- readRDS(here("processed-data","ebs_oisst_rotated.rds"))
wc.oisst.coastdist <- readRDS(here("processed-data","wc.oisst.coastdist.rds"))
neus.hadisst.coastdist <- readRDS(here("processed-data","neus.hadisst.coastdist.rds"))

# note that daily datasets are converted into monthly means for comparability among regions 
ebs.oisst.NWprep <- ebs.oisst.rotated %>% 
  group_by(N_km, E_km, month) %>%
  mutate(cell.month.mean=mean(sst)) %>%
  ungroup() %>%
  group_by(NW_km, year_match) %>%
  mutate(sstmean = mean(cell.month.mean),
         sstmax = max(cell.month.mean),
         sstmin = min(cell.month.mean)) %>%
  ungroup() %>%
  dplyr::select(year_match, NW_km, sstmean, sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

wc.oisst.prepgam <- wc.oisst.coastdist %>% 
  group_by(x, y, month) %>%
  mutate(cell.month.mean = mean(sst)) %>%
  group_by(coast_km, year_match) %>%
  mutate(sstmean = mean(cell.month.mean),
         sstmax = max(cell.month.mean),
         sstmin = min(cell.month.mean)) %>%
  ungroup() %>%
  dplyr::select(year_match, coast_km, sstmean,sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

neus.hadisst.prepgam <- neus.hadisst.coastdist %>% 
  group_by(coast_km, year_match) %>%
  mutate(sstmean = mean(sst),
         sstmax = max(sst),
         sstmin = min(sst)) %>%
  ungroup() %>%
  dplyr::select(year_match, coast_km, sstmean,sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

# set up GAMs
ebs.oisst.nw.temp.gam.mean <- gam(sstmean ~ year_match + s(NW_km, by=year_match), data=ebs.oisst.NWprep)
ebs.oisst.nw.temp.gam.99 <- gam(sstmax ~ year_match + s(NW_km, by=year_match), data=ebs.oisst.NWprep)
ebs.oisst.nw.temp.gam.01 <- gam(sstmin ~ year_match + s(NW_km, by=year_match), data=ebs.oisst.NWprep)

wc.oisst.temp.gam.mean <- gam(sstmean ~ year_match + s(coast_km, by=year_match), data=wc.oisst.prepgam)
wc.oisst.temp.gam.99 <- gam(sstmax ~ year_match + s(coast_km, by=year_match), data=wc.oisst.prepgam)
wc.oisst.temp.gam.01 <- gam(sstmin ~ year_match + s(coast_km, by=year_match), data=wc.oisst.prepgam)

neus.hadisst.temp.gam.mean <- gam(sstmean ~ year_match + s(coast_km, by=year_match), data=neus.hadisst.prepgam)
neus.hadisst.temp.gam.99 <- gam(sstmax ~ year_match + s(coast_km, by=year_match), data=neus.hadisst.prepgam)
neus.hadisst.temp.gam.01 <- gam(sstmin ~ year_match + s(coast_km, by=year_match), data=neus.hadisst.prepgam)

# predict temp from edge position--prep datasets
ebs.minyear <- min(ebs.oisst.rotated$year_match)
wc.minyear <- min(wc.oisst.coastdist$year_match)
neus.minyear <- min(neus.hadisst.coastdist$year_match)

ebs.pred <- dat.models %>%
  filter(region=="ebs",
         year >= ebs.minyear) %>%
  dplyr::select( -axis) %>%
  rename(NW_km=Estimate,
         year_match=year) 

wc.pred <- dat.models %>%
  filter(region=="wc",
         year >= wc.minyear) %>%
  dplyr::select( -axis) %>%
  rename(coast_km=Estimate,
         year_match=year) 

neus.pred <- dat.models %>%
  filter(region=="neus",
         year >= neus.minyear) %>%
  dplyr::select( -axis) %>%
  rename(coast_km=Estimate,
         year_match=year) 

# add columns with predicted temperature at edge every year 
ebs.pred$predict.sstmean <- predict.gam(ebs.oisst.nw.temp.gam.mean, ebs.pred)
ebs.pred$predict.sstmean.se <- predict.gam(ebs.oisst.nw.temp.gam.mean, ebs.pred, se.fit=TRUE)$se.fit
ebs.pred$predict.sstmax <- predict.gam(ebs.oisst.nw.temp.gam.99, ebs.pred)
ebs.pred$predict.sstmax.se <- predict.gam(ebs.oisst.nw.temp.gam.99, ebs.pred,se.fit=TRUE)$se.fit
ebs.pred$predict.sstmin <- predict.gam(ebs.oisst.nw.temp.gam.01, ebs.pred)
ebs.pred$predict.sstmin.se <- predict.gam(ebs.oisst.nw.temp.gam.01, ebs.pred,se.fit=TRUE)$se.fit

wc.pred$predict.sstmean <- predict.gam(wc.oisst.temp.gam.mean, wc.pred)
wc.pred$predict.sstmean.se <- predict.gam(wc.oisst.temp.gam.mean, wc.pred, se.fit=TRUE)$se.fit
wc.pred$predict.sstmax <- predict.gam(wc.oisst.temp.gam.99, wc.pred)
wc.pred$predict.sstmax.se <- predict.gam(wc.oisst.temp.gam.99, wc.pred,se.fit=TRUE)$se.fit
wc.pred$predict.sstmin <- predict.gam(wc.oisst.temp.gam.01, wc.pred)
wc.pred$predict.sstmin.se <- predict.gam(wc.oisst.temp.gam.01, wc.pred, se.fit=TRUE)$se.fit

neus.pred$predict.sstmean <- predict.gam(neus.hadisst.temp.gam.mean, neus.pred)
neus.pred$predict.sstmean.se <- predict.gam(neus.hadisst.temp.gam.mean, neus.pred, se.fit=TRUE)$se.fit
neus.pred$predict.sstmax <- predict.gam(neus.hadisst.temp.gam.99, neus.pred)
neus.pred$predict.sstmax.se <- predict.gam(neus.hadisst.temp.gam.99, neus.pred, se.fit=TRUE)$se.fit
neus.pred$predict.sstmin <- predict.gam(neus.hadisst.temp.gam.01, neus.pred)
neus.pred$predict.sstmin.se <- predict.gam(neus.hadisst.temp.gam.01, neus.pred, se.fit=TRUE)$se.fit

neus.pred <- rename(neus.pred, edge_position=coast_km)
neus.pred$axis <- "coast_km"

wc.pred <- rename(wc.pred, edge_position=coast_km)
wc.pred$axis <- "coast_km"

ebs.pred <- rename(ebs.pred, edge_position=NW_km)
ebs.pred$axis <- "NW_km"

dat.predict1 <- rbind(neus.pred, wc.pred, ebs.pred)%>%
  select(-predict.sstmean.se, -predict.sstmax.se, -predict.sstmin.se) %>%
  pivot_longer(cols=c(predict.sstmean, predict.sstmax, predict.sstmin), names_to="predicted.var",values_to="sst") 


dat.predict <- rbind(neus.pred, wc.pred, ebs.pred)%>%
  select(-predict.sstmean, -predict.sstmax, -predict.sstmin) %>%
  pivot_longer(cols=c(predict.sstmean.se, predict.sstmax.se, predict.sstmin.se), names_to="predicted.var",values_to="sstSE") %>%
  mutate(predicted.var=str_replace(predicted.var, ".se","")) %>%
  inner_join(dat.predict1)

# plot thermal niches over time
neus.thermal.niche.gg <- dat.predict %>%
  filter(region=="neus") %>%
  ggplot() +
  geom_line(aes(x=year_match, y=sst, group=predicted.var, color=predicted.var)) +
  geom_point(aes(x=year_match, y=sst, group=predicted.var, color=predicted.var)) +
  geom_errorbar(aes(x=year_match, ymax=sst+sstSE, ymin=sst-sstSE, group=predicted.var, color=predicted.var))+
  facet_wrap(quantile ~ species) +
  theme_bw() +
  theme(legend.position = "bottom")
neus.thermal.niche.gg
ggsave(neus.thermal.niche.gg, filename=here("results","neus_edge_niches_time.png"),dpi=160, width=12,height=11,scale=1.3)

wc.thermal.niche.gg <- dat.predict %>%
  filter(region=="wc")%>%
  ggplot(aes(x=year_match, y=sst, group=predicted.var, color=predicted.var)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(x=year_match, ymax=sst+sstSE, ymin=sst-sstSE, group=predicted.var, color=predicted.var))+
  facet_wrap(quantile ~ species)+
  theme_bw() +
  theme(legend.position = "bottom")
wc.thermal.niche.gg
ggsave(wc.thermal.niche.gg, filename=here("results","wc_edge_niches_time.png"),dpi=160, width=12,height=11,scale=1.3)

ebs.thermal.niche.gg <- dat.predict %>%
  filter(region=="ebs")%>%
  ggplot(aes(x=year_match, y=sst, group=predicted.var, color=predicted.var)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(x=year_match, ymax=sst+sstSE, ymin=sst-sstSE, group=predicted.var, color=predicted.var))+
  facet_wrap(quantile ~ species)+
  theme_bw() +
  theme(legend.position = "bottom")
ebs.thermal.niche.gg
ggsave(ebs.thermal.niche.gg, filename=here("results","ebs_edge_niches_time.png"),dpi=160, width=12,height=11,scale=1.3)

# Bayesian test for edge thermal niche change over time 

dat.predict.niche <- dat.predict %>%
  filter(!predicted.var=="predict.sstmean")

spp.bayes.niche.lm.df <- NULL
for(i in unique(dat.predict.niche$species)) {
  dfprep1 <- dat.predict.niche[dat.predict.niche$species==i,] # subdivide by species 
  for(j in unique(dfprep1$region)) {
    dfprep2 <- dat.predict.niche[dat.predict.niche$species==i & dat.predict.niche$region==j,] # subdivide by region, for the few spp found in multiple regions 
    for(k in unique(dfprep2$quantile)) {
      dfprep3 <- dat.predict.niche[dat.predict.niche$species==i & dat.predict.niche$region==j & dat.predict.niche$quantile==k,] # subdivide by quantile, for the few spp with both range edges
      for(l in unique(dfprep3$predicted.var)){
        df <- dat.predict.niche[dat.predict.niche$species==i & dat.predict.niche$region==j & dat.predict.niche$quantile==k & dat.predict.niche$predicted.var==l,] # split by the two temperature extremes that we want to analyze separately
        spp.bayes.lm <- try(stan_glm(sst ~ year_match, 
                                     data=df, 
                                     family=gaussian(), 
                                     iter = 40000,
                                     warmup = 10000,
                                     adapt_delta = 0.99,
                                     chains = 4,
                                     cores = 1,
                                     prior = normal(),
                                     control = list(max_treedepth = 20)))
        if(!class(spp.bayes.lm)[1] == "try-error") { # adding try() here because some edges are so invariant that the model fails 
          spp.bayes.lm.tidy <- tidy_draws(spp.bayes.lm) %>%
            mutate(species = paste0(i),
                   region = paste0(j),
                   quantile = paste0(k),
                   predicted.var = paste0(l),
                   intercept.rhat = summary(spp.bayes.lm)[,"Rhat"][1],
                   year_match.rhat = summary(spp.bayes.lm)[,"Rhat"][2],
                   sigma.rhat = summary(spp.bayes.lm)[,"Rhat"][3])
          spp.bayes.niche.lm.df <- rbind(spp.bayes.niche.lm.df, spp.bayes.lm.tidy)
        }
      }
    }
  }
}
quantile(spp.bayes.niche.lm.df$intercept.rhat)
quantile(spp.bayes.niche.lm.df$year_match.rhat) # this one is not causing estimation problems 
quantile(spp.bayes.niche.lm.df$sigma.rhat) # this is the problem child 

spp.bayes.niche.filter <- spp.bayes.niche.lm.df %>% 
  group_by(region, species, quantile) %>%
  mutate(max.rhat = max(intercept.rhat, year_match.rhat, sigma.rhat)) %>%
  filter(max.rhat <= 1.1) # get rid of spp*region*edge combos where one of the SST extreme models didn't converge

setdiff(spp.bayes.niche.lm.df %>% select(region, species, quantile) %>% distinct(), spp.bayes.niche.filter %>% select(region, species, quantile) %>% distinct()) # just one row removed: neptunea lyrata, ebs, cold edge

spp.bayes.niche.lm.stats <- spp.bayes.niche.filter %>% 
  group_by(.draw, species, region, quantile, predicted.var) %>%
  summarise(beta.mean = mean(year_match)) %>%
  group_by(species, region, quantile, predicted.var) %>%
  summarise(mean=mean(beta.mean),
            median=median(beta.mean),
            lower=quantile(beta.mean, 0.05),
            upper=quantile(beta.mean, 0.95)) %>%
  select(species, region, quantile, predicted.var, mean, median, lower, upper) %>%
  distinct()

quantile(spp.bayes.niche.lm.stats$mean)
spp.bayes.niche.lm.stats %>% 
  ggplot() +
  geom_histogram(aes(x=mean))

spp.bayes.niche.groups <- spp.bayes.niche.lm.stats %>%
  left_join(dat.models.groups %>% select(-edgetype), by=c("region","species")) %>%
  rowwise() %>%
  mutate(crosses0 = ifelse(lower<0 & upper>0, TRUE, FALSE)) %>%
  group_by(region, quantile, species) %>% 
  mutate(niche.group = ifelse(min(abs(mean)) < 0.01, "good_tracker", 
                              ifelse(TRUE %in% crosses0, "lag_tracker", "non_tracker"))) %>%
  select(region, quantile, species, niche.group) %>% 
  distinct()
write_csv(spp.bayes.niche.groups, here("results","species_by_thermal_niche_group.csv"))

# calculate stats
spp.bayes.niche.groups %>% left_join(dat.models.groups) %>% group_by(niche.group, taxongroup) %>% summarise(n=n())
spp.bayes.niche.groups %>% left_join(dat.models.groups) %>% group_by(region, taxongroup) %>% summarise(n=n())
spp.bayes.niche.groups %>% group_by(niche.group, region) %>% summarise(n=n())
spp.bayes.niche.groups %>% group_by(niche.group, quantile) %>% summarise(n=n())
spp.bayes.niche.groups %>% group_by(region, quantile) %>% summarise(n=n())

dat.thermal.gg <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(region=factor(region, levels=c('neus','wc','ebs')),
         region=recode(region,
                       ebs="Eastern Bering Sea",
                       neus="Northeast",
                       wc="West Coast"),
         niche.group=factor(niche.group, levels=c('good_tracker','lag_tracker','non_tracker')),
         niche.group=recode(niche.group,
                            good_tracker="TNH",
                            lag_tracker="TLH",
                            non_tracker="TIH")) %>%
  ggplot(aes(x=mean_predict.sstmin, y=mean_predict.sstmax, group=niche.group)) +
  geom_point(aes(shape=niche.group, color=niche.group, fill=niche.group)) +
  scale_color_manual(values = c("TIH"="#fb8072","TNH"="deepskyblue4","TLH"="mediumpurple3")) +
  geom_errorbar(aes(x=mean_predict.sstmin, ymin=lower_predict.sstmax, ymax=upper_predict.sstmax, color=niche.group)) +
  geom_errorbarh(aes(y=mean_predict.sstmax, xmin=lower_predict.sstmin, xmax=upper_predict.sstmin, color=niche.group)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() +
  labs(x="Change in Cold Extreme at Edge (°C/year)", y="Change in Warm Extreme at Edge (°C/year)") +
  facet_wrap(~region) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  NULL
dat.thermal.gg
ggsave(dat.thermal.gg, filename=here("results","thermal_niche_time_coefficients.png"),dpi=160, width=10, height=4)

gg.niche.regions <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(region=factor(region, levels=c('neus','wc','ebs')),
         region=recode(region,
                       ebs="Eastern Bering Sea",
                       neus="Northeast",
                       wc="West Coast"),
         niche.group=factor(niche.group, levels=c('good_tracker','lag_tracker','non_tracker')),
         niche.group=recode(niche.group,
                            good_tracker="TNH",
                            lag_tracker="TLH",
                            non_tracker="TIH")) %>%
  ggplot(aes(region)) +
  geom_bar(aes(fill=niche.group), color="black") +
  scale_fill_manual(values=c("white","darkgrey","black")) +
  theme_bw() +
  labs(x=NULL, y=NULL)+
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  NULL
gg.niche.regions

gg.niche.quantile <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(quantile=factor(quantile, levels=c('quantile_0.05','quantile_0.95')),
         quantile=recode(quantile,
                         quantile_0.05="Warm Edge",
                         quantile_0.95="Cold Edge"),
         niche.group=factor(niche.group, levels=c('good_tracker','lag_tracker','non_tracker')),
         niche.group=recode(niche.group,
                            good_tracker="TNH",
                            lag_tracker="TLH",
                            non_tracker="TIH")) %>%
  ggplot(aes(quantile)) +
  geom_bar(aes(fill=niche.group), color="black") +
  scale_fill_manual(values=c("white","darkgrey","black")) +
  theme_bw() +
  labs(x=NULL, y=NULL)+
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  NULL
gg.niche.quantile

gg.niche.taxa <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  left_join(dat.models.groups) %>% 
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(
    niche.group=factor(niche.group, levels=c('good_tracker','lag_tracker','non_tracker')),
    niche.group=recode(niche.group,
                       good_tracker="TNH",
                       lag_tracker="TLH",
                       non_tracker="TIH")) %>%
  ggplot(aes(taxongroup)) +
  geom_bar(aes(fill=niche.group), color="black") +
  scale_fill_manual(values=c("white","darkgrey","black")) +
  theme_bw() +
  labs(x=NULL, y=NULL)+
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  NULL
gg.niche.taxa
ggsave(gg.niche.regions, filename=here("results","niche_barplot_region.png"), dpi=300, width=4, height=7, scale=0.8)
ggsave(gg.niche.quantile, filename=here("results","niche_barplot_quantile.png"), dpi=300, width=3, height=7, scale=0.8)
ggsave(gg.niche.taxa, filename=here("results","niche_barplot_taxa.png"), dpi=300, width=3, height=7, scale=0.8)

# make example plots
ex.spp1 <- "gadus macrocephalus" # good tracker, cold edge, EBS
ex.spp2 <- "sebastes pinniger" # non tracker, warm edge, WC
ex.spp3 <- "paralichthys oblongus" # lagged tracker, cold edge, NEUS

ex.spp.bayes.gg1 <- spp.bayes.niche.filter %>% 
  filter(species==ex.spp1) %>% 
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>% 
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
  theme(legend.position="bottom") +
  NULL
ex.spp.bayes.gg1

ex.spp.bayes.gg2 <- spp.bayes.niche.filter %>% 
  filter(species==ex.spp2) %>% 
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>% 
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
  theme(legend.position="bottom") +
  NULL
ex.spp.bayes.gg2

ex.spp.bayes.gg3 <- spp.bayes.niche.filter %>% 
  filter(species==ex.spp3) %>% 
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>% 
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
  theme(legend.position="bottom") +
  NULL
ex.spp.bayes.gg3

ex.spp.time.gg1 <- dat.predict.niche %>% 
  filter(species==ex.spp1) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
  theme(legend.position="bottom")+
  scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+
  NULL
ex.spp.time.gg1

ex.spp.time.gg2 <- dat.predict.niche %>% 
  filter(species==ex.spp2) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
  theme(legend.position="bottom")+
  scale_x_continuous(limits=c(1982, 2018), breaks=seq(1982, 2018, 4))+
  NULL
ex.spp.time.gg2

ex.spp.time.gg3 <- dat.predict.niche %>% 
  filter(species==ex.spp3) %>%
  mutate(species = str_to_sentence(species)) %>%
  ggplot() +
  geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
  geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
  scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  theme_bw() +
  labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
  theme(legend.position="bottom")+
  scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+
  NULL
ex.spp.time.gg3

ggsave(ex.spp.bayes.gg1, dpi=160, width=4, height=4, filename=here("results","example_1_posterior.png"))
ggsave(ex.spp.bayes.gg2, dpi=160, width=4, height=4, filename=here("results","example_2_posterior.png"))
ggsave(ex.spp.bayes.gg3, dpi=160, width=4, height=4, filename=here("results","example_3_posterior.png"))
ggsave(ex.spp.time.gg1, dpi=160, width=4, height=4, filename=here("results","example_1_niche.png"))
ggsave(ex.spp.time.gg2, dpi=160, width=4, height=4, filename=here("results","example_2_niche.png"))
ggsave(ex.spp.time.gg3, dpi=160, width=4, height=4, filename=here("results","example_3_niche.png"))

############### 
# test for edge thermal niche change over time 

# dat.predict.lms <- dat.predict %>%
#   group_by(species, region, quantile, predicted.var) %>%
#   filter(!predicted.var=="predict.sstmean") %>%
#   nest() %>%
#   mutate(model=purrr::map(data, ~lm(sst~year_match, data=.x)),
#          tidymodel=purrr::map(model, tidy)) %>%
#   unnest(tidymodel) %>% 
#   ungroup() %>%
#   filter(!term=="(Intercept)") %>%
#   select(-data, -model)
#  
# clusters <- kmeans(dat.predict.lms[dat.predict.lms$p.value>0.05,]$std.error, centers=2)
# 
# zero.slopes <- dat.predict.lms %>%
#   filter(p.value > 0.05) %>%
#   add_column(clusterID=clusters$cluster) %>%
#   group_by(clusterID) %>%
#   mutate(clustermean = mean(std.error)) %>%
#   ungroup() %>%
#   mutate(group = ifelse(clusterID==2, "zeroSlope_lowSE","zeroSlope_highSE")) %>%
#   select(species, region, predicted.var, quantile, clustermean, group)
# 
# # gives 2 rows to each species--one for each thermal extreme 
# dat.thermal.niche <- dat.predict.lms %>%
#   left_join(zero.slopes, by=c("species","region","predicted.var","quantile")) %>%
#   mutate(group=replace_na(group, "nonzeroSlope"))
# 
# dat.thermal.niche.grouped <- dat.thermal.niche %>% 
#   group_by(species, quantile, region) %>%
#   mutate(spp.biogeo.class = ifelse("zeroSlope_lowSE" %in% group, "Fast Temperature Tracker", ifelse("zeroSlope_highSE" %in% group, "Lagged Temperature Tracker","Non Temperature Tracker"))) %>% # classify by range edge type, not by temperature extreme
#   ungroup() %>% 
#   select(-term, -clustermean, -statistic, -group) %>%
#   pivot_wider(values_from=c(estimate, std.error, p.value), names_from = predicted.var)
# 
# thermal.niche.summary <- dat.thermal.niche.grouped %>% 
#   group_by(spp.biogeo.class) %>%
#   summarise(n=n())
# 
# dat.thermal.gg <- dat.thermal.niche.grouped %>%
#   mutate(region=factor(region, levels=c('neus','wc','ebs'),),
#          region=recode(region,
#                        ebs="Eastern Bering Sea",
#                        neus="Northeast",
#                        wc="West Coast")) %>%
#   ggplot(aes(x=estimate_predict.sstmin, y=estimate_predict.sstmax, group=spp.biogeo.class)) +
#   geom_point(aes(shape=spp.biogeo.class, color=spp.biogeo.class, fill=spp.biogeo.class)) +
#   scale_color_manual(values = c("Non Temperature Tracker"="#fb8072","Fast Temperature Tracker"="deepskyblue4","Lagged Temperature Tracker"="mediumpurple3")) +
#   geom_errorbar(aes(x=estimate_predict.sstmin, ymin=estimate_predict.sstmax-std.error_predict.sstmax, ymax=estimate_predict.sstmax+std.error_predict.sstmax, color=spp.biogeo.class)) +
#   geom_errorbarh(aes(y=estimate_predict.sstmax, xmin=estimate_predict.sstmin-std.error_predict.sstmin, xmax=estimate_predict.sstmin+std.error_predict.sstmin, color=spp.biogeo.class)) +
#   geom_vline(xintercept=0, linetype="dashed") +
#   geom_hline(yintercept=0, linetype="dashed") +
#   theme_bw() +
#   labs(x="Change in Cold Extreme at Edge (°C/year)", y="Change in Warm Extreme at Edge (°C/year)") +
# facet_wrap(~region) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank()) +
#   NULL
# dat.thermal.gg
# ggsave(dat.thermal.gg, filename=here("results","thermal_niche_time_coefficients.png"),dpi=160, width=10, height=4)
# 
# #######################
# ### are range limits also niche limits? 
# #######################
# # many species did have a significant change in max SST over time, but with a very small effect size, so I think it's fine to use time-series means. 
# wc.niche <- readRDS(here("processed-data","wc_niche.rds")) %>%
#   pivot_longer(cols=c('spp.sst.mean','spp.sst.max','spp.sst.min'), names_to="sstvar", values_to="sst") %>%
#   group_by(species, sstvar) %>% 
#   mutate(niche.sst.mean = mean(sst),
#          niche.sst.sd = sd(sst)) %>%
#   ungroup() %>%
#   select(species, niche.sst.mean, niche.sst.sd, sstvar) %>%
#   distinct() %>%
#   mutate(region="wc",
#          sstvar=recode(sstvar, 
#          spp.sst.mean="sst.mean",
#          spp.sst.max="sst.max",
#          spp.sst.min="sst.min"),
#          species = tolower(gsub("_"," ",species)))
# 
# neus.niche <- readRDS(here("processed-data","neus_niche.rds")) %>%
#   pivot_longer(cols=c('spp.sst.mean','spp.sst.max','spp.sst.min'), names_to="sstvar", values_to="sst") %>%
#   group_by(species, sstvar) %>% 
#   mutate(niche.sst.mean = mean(sst, na.rm=TRUE),
#          niche.sst.sd = sd(sst, na.rm=TRUE)) %>%
#   ungroup() %>%
#   select(species, niche.sst.mean, niche.sst.sd, sstvar) %>%
#   distinct() %>%
#   mutate(region="neus",
#          sstvar=recode(sstvar, 
#                        spp.sst.mean="sst.mean",
#                        spp.sst.max="sst.max",
#                        spp.sst.min="sst.min"),
#          species = tolower(gsub("_"," ",species)))
# 
# ebs.niche <- readRDS(here("processed-data","ebs_niche.rds")) %>%
#   pivot_longer(cols=c('spp.sst.mean','spp.sst.max','spp.sst.min'), names_to="sstvar", values_to="sst") %>%
#   group_by(species, sstvar) %>% 
#   mutate(niche.sst.mean = mean(sst, na.rm=TRUE),
#          niche.sst.sd = sd(sst, na.rm=TRUE)) %>%
#   ungroup() %>%
#   select(species, niche.sst.mean, niche.sst.sd, sstvar) %>%
#   distinct() %>%
#   mutate(region="ebs",
#          sstvar=recode(sstvar, 
#                        spp.sst.mean="sst.mean",
#                        spp.sst.max="sst.max",
#                        spp.sst.min="sst.min"),
#          species = tolower(gsub("_"," ",species)))
# 
# dat.niche <- rbind(wc.niche, neus.niche, ebs.niche)
# 
# # only use species with a range edge in the region based on dat.predict 
# dat.predict.means <- dat.predict %>%
#   rename("sstvar"=predicted.var) %>%
#   group_by(species, region, quantile, sstvar) %>%
#   mutate(edge.sst.mean = mean(sst),
#          edge.sst.sd = sd(sst)) %>%
#   ungroup() %>%
#   select(species, region, quantile, sstvar, edge.sst.mean, edge.sst.sd) %>%
#   distinct() %>%
#   mutate(
#          sstvar=recode(sstvar, 
#                        predict.sstmean="sst.mean",
#                        predict.sstmax="sst.max",
#                        predict.sstmin="sst.min")) %>%
#   left_join(dat.niche, by=c("region","species","sstvar"))
# 
# ebs.niche.gg <- dat.predict %>%
#   filter(region=="ebs") %>%
#   ggplot(aes(x=year_match, y=sst, group=predicted.var))+
#   geom_line(aes(color=predicted.var)) + 
# facet_wrap(~species) 
# 
# wc.niche.gg <- dat.predict %>%
#   filter(region=="wc") %>%
#   ggplot(aes(x=year_match, y=sst, group=predicted.var))+
#   geom_line(aes(color=predicted.var)) + 
#   facet_wrap(~species) 
# 
# neus.niche.gg <- dat.predict %>%
#   filter(region=="neus") %>%
#   ggplot(aes(x=year_match, y=sst, group=predicted.var))+
#   geom_line(aes(color=predicted.var)) + 
#   facet_wrap(~species) 
# 
# cold.niche.r2 <- summary(lm(edge.sst.mean~niche.sst.mean, 
#                             data=dat.predict.means %>%
#                filter(quantile=="quantile_0.95",
#                       sstvar=="sst.min")))$r.squared
# cold.niche.r2.exp <- paste("R^2 == ", round(cold.niche.r2,2))
# 
# cold.niche.gg <- dat.predict.means %>%
#   filter(quantile=="quantile_0.95",
#          sstvar=="sst.min") %>%
#   ggplot(aes(x=niche.sst.mean, y=edge.sst.mean)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin=niche.sst.mean-niche.sst.sd, xmax=niche.sst.mean+niche.sst.sd)) +
#   geom_errorbar(aes(ymin=edge.sst.mean-edge.sst.sd, ymax=edge.sst.mean+edge.sst.sd)) +
#   annotate(geom="text",label=cold.niche.r2.exp, x=11, y=13.5, parse=TRUE, size=5) +
#   geom_abline(slope=1, intercept=0,linetype="dashed") +
#   scale_x_continuous(limits=c(-2, 15)) + 
#   scale_y_continuous(limits=c(-2, 15)) + 
#   labs(x="Range Cold Extreme (°C)", y="Edge Cold Extreme (°C)", title="Cold Range Limits") +
#   theme_bw() +
#   theme(text=element_text(family="sans",size=12,color="black"),
#         legend.position="none",
#         axis.text=element_text(family="sans",size=8,color="black"), 
#         axis.title=element_text(family="sans",size=12,color="black"),
#         panel.grid.minor = element_blank()
#   ) + 
#   NULL
# cold.niche.gg
# ggsave(cold.niche.gg, filename=here("results","edge_vs_range_niche_cold.png"), dpi=160, width=5, height=5)
# 
# warm.niche.r2 <- summary(lm(edge.sst.mean~niche.sst.mean, 
#                             data=dat.predict.means %>%
#                               filter(quantile=="quantile_0.05",
#                                      sstvar=="sst.max")))$r.squared
# warm.niche.r2.exp <- paste("R^2 == ", round(warm.niche.r2,2))
# 
# warm.niche.gg <- dat.predict.means %>%
#   filter(quantile=="quantile_0.05",
#          sstvar=="sst.max") %>%
#   ggplot(aes(x=niche.sst.mean, y=edge.sst.mean)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin=niche.sst.mean-niche.sst.sd, xmax=niche.sst.mean+niche.sst.sd)) +
#   geom_errorbar(aes(ymin=edge.sst.mean-edge.sst.sd, ymax=edge.sst.mean+edge.sst.sd)) +
#   annotate(geom="text",label=warm.niche.r2.exp, x=24, y=28, parse=TRUE, size=5) +
#   geom_abline(slope=1, intercept=0,linetype="dashed") +
#   scale_x_continuous(limits=c(8, 30)) + 
#   scale_y_continuous(limits=c(8, 30)) + 
#   labs(x="Range Warm Extreme (°C)", y="Edge Warm Extreme (°C)", title="Warm Range Limits") +
#   theme_bw() +
#   theme(text=element_text(family="sans",size=12,color="black"),
#         legend.position="none",
#         axis.text=element_text(family="sans",size=8,color="black"), 
#         axis.title=element_text(family="sans",size=12,color="black"),
#         panel.grid.minor = element_blank()
#   ) + 
#   NULL
# warm.niche.gg
# ggsave(warm.niche.gg, filename=here("results","edge_vs_range_niche_warm.png"), dpi=160, width=5, height=5)
# 
# alpha=0.05
# power.lvl=0.95
# nsamples=dim(dat.summary)[1]
# powerTOSTone(alpha=alpha, statistical_power = power.lvl, N=nsamples)# here using N as num. samples, and then for each model, N will become num. years....
# testlmdf <- dat.predict.lms[214,]
# TOSTone(m=testlmdf$estimate, mu=0, sd=testlmdf$std.error, n=34, low_eqbound_d = -0.35, high_eqbound_d = 0.35)

# corresponds to Cohen's d of +/- 0.35 

# edge.niche.gg <- dat.predict.means %>%
#   mutate(sstvar=recode(sstvar, 
#                        sst.max="Maximum Monthly SST",
#                        sst.min="Minimum Monthly SST", 
#                        sst.mean="Mean Monthly SST"
#          ),
#          region=recode(region,
#                        ebs="Eastern Bering Sea",
#                        neus="Northeast",
#                        wc="West Coast")) %>%
#   ggplot(aes(x=niche.sst.mean, y=edge.sst.mean)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin=niche.sst.mean-niche.sst.sd, xmax=niche.sst.mean+niche.sst.sd)) +
#   geom_errorbar(aes(ymin=edge.sst.mean-edge.sst.sd, ymax=edge.sst.mean+edge.sst.sd)) +
#   facet_grid(region~sstvar)+
#   geom_abline(slope=1, intercept=0,linetype="dashed") +
#   scale_x_continuous(limits=c(-2, 30)) + 
#   scale_y_continuous(limits=c(-2, 30)) + 
#   labs(x="Range Thermal Niche (°C)", y="Edge Thermal Niche (°C)") +
#   theme_bw() +
#   theme(text=element_text(family="sans",size=12,color="black"),
#         legend.position="none",
#         axis.text=element_text(family="sans",size=8,color="black"), 
#        # axis.text.x = element_text(angle = 90, hjust = 1), 
#         axis.title=element_text(family="sans",size=12,color="black"),
#  #       strip.text.x = element_blank(),
#    #     panel.grid.major = element_blank() 
#         panel.grid.minor = element_blank()
#  ) + 
#   NULL
# edge.niche.gg
# ggsave(edge.niche.gg, filename=here("results","edge_vs_range_niche.png"), dpi=160, width=6, height=6)
