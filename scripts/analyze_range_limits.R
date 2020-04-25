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

# prep explanatory variables for models using region-wide change in T 

ebs.sst.summary <- read_rds(here("processed-data","ebs_sst_rotated.rds")) %>% 
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

# ebs.sst.scalemean = mean(ebs.sst.summary$mean.annual.sst)
# ebs.sst.scalesd = sd(ebs.sst.summary$mean.annual.sst)
# ebs.sst.summary$mean.annual.sst.scale <- (ebs.sst.summary$mean.annual.sst-ebs.sst.scalemean)/ebs.sst.scalesd

wc.sst.summary <- read_rds(here("processed-data","wc_sst_coastdist.rds")) %>% 
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

# wc.sst.scalemean = mean(wc.sst.summary$mean.annual.sst)
# wc.sst.scalesd = sd(wc.sst.summary$mean.annual.sst)
# wc.sst.summary$mean.annual.sst.scale <- (wc.sst.summary$mean.annual.sst-wc.sst.scalemean)/wc.sst.scalesd

neus.sst.summary <- read_rds(here("processed-data","neus_sst_coastdist.rds")) %>% 
  rename(lon="x",lat="y") %>%
  group_by(year_match, lat, lon) %>%
  mutate(cell.annual.sst = mean(sst)) %>% # calculate annual mean of monthly SSTs for each grid cell
  ungroup() %>%
  group_by(year_match) %>%
  mutate(mean.annual.sst = mean(cell.annual.sst)) %>% # calculate annual mean of all grid cells
  select(year_match, mean.annual.sst) %>%
  distinct() %>%
  mutate(region="neus")

dat.sst.summary <- rbind(neus.sst.summary, wc.sst.summary, ebs.sst.summary)

dat.edges.sst.summary <- dat.models %>%
  left_join(dat.sst.summary, by=c("region","year"="year_match"))

#######################
### single-species Bayesian models of edge position ~ regional T or time 
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
for(i in unique(dat.edges.sst.summary$species)) {
  dfprep1 <- dat.edges.sst.summary[dat.edges.sst.summary$species==i,] # subdivide by species 
  for(j in unique(dfprep1$region)) {
    dfprep2 <- dat.edges.sst.summary[dat.edges.sst.summary$species==i & dat.edges.sst.summary$region==j,] # subdivide by region, for the few spp found in multiple regions 
    for(k in unique(dfprep2$quantile)) {
      df <- dat.edges.sst.summary[dat.edges.sst.summary$species==i & dat.edges.sst.summary$region==j & dat.edges.sst.summary$quantile==k,] # subdivide by quantile, for the few spp with both range edges
      spp.bayes.lm <- try(stan_glm(Estimate ~ mean.annual.sst, 
                                   data=df, 
                                   family=gaussian(), 
                                   iter = 12000,
                                   warmup = 2000,
                                   adapt_delta = 0.95,
                                   chains = 4,
                                   cores = 1,
                                   prior = normal(0, 1000),# Burrows paper reports temperature spatial gradients up to 0.02 degrees C / km. taking the inverse, this could range from 50km shift for degree C up to 1000 km or higher. centering on 0 because relationship is not necessarily positive for all species. 
                                   weights = 1/(Std.Error)^2)) 
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
                                   prior = normal(0, 100),
                                   weights = 1/(Std.Error)^2
                                   )) # centered on 0 to allow negative coefficients for some species. SD of 100 makes it pretty flat but still bounds to a reasonable order of magnitude in km/yr
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
                         quantile_0.01="Warm Limit",
                         quantile_0.99="Cold Limit")) %>%
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
ebs.sst <- readRDS(here("processed-data","ebs_sst_rotated.rds"))
wc.sst <- readRDS(here("processed-data","wc_sst_coastdist.rds"))
neus.sst <- readRDS(here("processed-data","neus_sst_coastdist.rds"))

# note that datasets have both monthly readings (pre-1982, from hadisst) and daily (post-1982, from oisst) that are first converted into monthly means for comparability among regions/years
ebs.sst.prepgam <- ebs.sst %>% 
  group_by(N_km, E_km, year_match, month) %>%
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

wc.sst.prepgam <- wc.sst %>% 
  group_by(x, y, year_match, month) %>%
  mutate(cell.month.mean = mean(sst)) %>%
  group_by(coast_km, year_match) %>%
  mutate(sstmean = mean(cell.month.mean),
         sstmax = max(cell.month.mean),
         sstmin = min(cell.month.mean)) %>%
  ungroup() %>%
  dplyr::select(year_match, coast_km, sstmean,sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

neus.sst.prepgam <- neus.sst %>% 
  group_by(x, y, month, year_match) %>%
  mutate(cell.month.mean = mean(sst)) %>%
  group_by(coast_km, year_match) %>%
  mutate(sstmean = mean(cell.month.mean),
         sstmax = max(cell.month.mean),
         sstmin = min(cell.month.mean)) %>%
  ungroup() %>%
  dplyr::select(year_match, coast_km, sstmean,sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

# set up GAMs for coastal and NW axes 
ebs.sst.temp.gam.mean <- gam(sstmean ~ year_match + s(NW_km, by=year_match), data=ebs.sst.prepgam)
ebs.sst.temp.gam.99 <- gam(sstmax ~ year_match + s(NW_km, by=year_match), data=ebs.sst.prepgam)
ebs.sst.temp.gam.01 <- gam(sstmin ~ year_match + s(NW_km, by=year_match), data=ebs.sst.prepgam)

wc.sst.temp.gam.mean <- gam(sstmean ~ year_match + s(coast_km, by=year_match), data=wc.sst.prepgam)
wc.sst.temp.gam.99 <- gam(sstmax ~ year_match + s(coast_km, by=year_match), data=wc.sst.prepgam)
wc.sst.temp.gam.01 <- gam(sstmin ~ year_match + s(coast_km, by=year_match), data=wc.sst.prepgam)

neus.sst.temp.gam.mean <- gam(sstmean ~ year_match + s(coast_km, by=year_match), data=neus.sst.prepgam)
neus.sst.temp.gam.99 <- gam(sstmax ~ year_match + s(coast_km, by=year_match), data=neus.sst.prepgam)
neus.sst.temp.gam.01 <- gam(sstmin ~ year_match + s(coast_km, by=year_match), data=neus.sst.prepgam)

# predict temp from edge position--prep datasets

ebs.pred <- dat.models %>%
  filter(region=="ebs") %>%
  dplyr::select( -axis) %>%
  rename(NW_km=Estimate,
         year_match=year) 

wc.pred <- dat.models %>%
  filter(region=="wc") %>%
  dplyr::select( -axis) %>%
  rename(coast_km=Estimate,
         year_match=year) 

neus.pred <- dat.models %>%
  filter(region=="neus") %>%
  dplyr::select( -axis) %>%
  rename(coast_km=Estimate,
         year_match=year) 

# add columns with predicted temperature at edge every year 
ebs.pred$predict.sstmean <- predict.gam(ebs.sst.temp.gam.mean, ebs.pred)
ebs.pred$predict.sstmean.se <- predict.gam(ebs.sst.temp.gam.mean, ebs.pred, se.fit=TRUE)$se.fit
ebs.pred$predict.sstmax <- predict.gam(ebs.sst.temp.gam.99, ebs.pred)
ebs.pred$predict.sstmax.se <- predict.gam(ebs.sst.temp.gam.99, ebs.pred,se.fit=TRUE)$se.fit
ebs.pred$predict.sstmin <- predict.gam(ebs.sst.temp.gam.01, ebs.pred)
ebs.pred$predict.sstmin.se <- predict.gam(ebs.sst.temp.gam.01, ebs.pred,se.fit=TRUE)$se.fit

wc.pred$predict.sstmean <- predict.gam(wc.sst.temp.gam.mean, wc.pred)
wc.pred$predict.sstmean.se <- predict.gam(wc.sst.temp.gam.mean, wc.pred, se.fit=TRUE)$se.fit
wc.pred$predict.sstmax <- predict.gam(wc.sst.temp.gam.99, wc.pred)
wc.pred$predict.sstmax.se <- predict.gam(wc.sst.temp.gam.99, wc.pred,se.fit=TRUE)$se.fit
wc.pred$predict.sstmin <- predict.gam(wc.sst.temp.gam.01, wc.pred)
wc.pred$predict.sstmin.se <- predict.gam(wc.sst.temp.gam.01, wc.pred, se.fit=TRUE)$se.fit

neus.pred$predict.sstmean <- predict.gam(neus.sst.temp.gam.mean, neus.pred)
neus.pred$predict.sstmean.se <- predict.gam(neus.sst.temp.gam.mean, neus.pred, se.fit=TRUE)$se.fit
neus.pred$predict.sstmax <- predict.gam(neus.sst.temp.gam.99, neus.pred)
neus.pred$predict.sstmax.se <- predict.gam(neus.sst.temp.gam.99, neus.pred, se.fit=TRUE)$se.fit
neus.pred$predict.sstmin <- predict.gam(neus.sst.temp.gam.01, neus.pred)
neus.pred$predict.sstmin.se <- predict.gam(neus.sst.temp.gam.01, neus.pred, se.fit=TRUE)$se.fit

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
                                     prior = normal(0, 2), # relatively flat but constrained, SST vs year is often order 0.1
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
quantile(spp.bayes.niche.lm.df$year_match.rhat)  
quantile(spp.bayes.niche.lm.df$sigma.rhat) # this previously caused estimation problems, doesn't appear to be doing that anymore  

spp.bayes.niche.filter <- spp.bayes.niche.lm.df %>% 
  group_by(region, species, quantile) %>%
  mutate(max.rhat = max(intercept.rhat, year_match.rhat, sigma.rhat)) %>%
  filter(max.rhat <= 1.1) # get rid of spp*region*edge combos where one of the SST extreme models didn't converge

setdiff(spp.bayes.niche.lm.df %>% select(region, species, quantile) %>% distinct(), spp.bayes.niche.filter %>% select(region, species, quantile) %>% distinct()) # 0 at present 

# SLOW
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

quantile(spp.bayes.niche.lm.stats$mean) # should be distributed in the neighborhood of zero
spp.bayes.niche.lm.stats %>% 
  ggplot() +
  geom_histogram(aes(x=mean))

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
spp.bayes.niche.groups %>% left_join(dat.models.groups) %>% group_by(niche.group, taxongroup) %>% summarise(n=n())
spp.bayes.niche.groups %>% left_join(dat.models.groups) %>% group_by(region, taxongroup) %>% summarise(n=n())
spp.bayes.niche.groups %>% group_by(niche.group, region) %>% summarise(n=n())
spp.bayes.niche.groups %>% group_by(niche.group, quantile) %>% summarise(n=n())
spp.bayes.niche.groups %>% group_by(region, quantile) %>% summarise(n=n())

# mega barplot 
gg.niche.barplot <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(quantile=factor(quantile, levels=c('quantile_0.01','quantile_0.99')),
         quantile=recode(quantile,
                         quantile_0.01="Warm Edge",
                         quantile_0.99="Cold Edge"),
         region=factor(region, levels=c('neus','wc','ebs')),
         region=recode(region,
                       ebs="Eastern Bering Sea",
                       neus="Northeast",
                       wc="West Coast"),
         niche.group=factor(niche.group, levels=c('good_tracker','partial_tracker','non_tracker')),
         niche.group=recode(niche.group,
                            good_tracker="TNH",
                            partial_tracker="PTH",
                            non_tracker="TIH")) %>%
  group_by(quantile, region, niche.group) %>%
  mutate(count = length(species)) %>%
  select(count, region, niche.group, quantile) %>%
  ungroup() %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(y=forcats::fct_rev(niche.group), x=count, fill=quantile), stat="identity", color="black") +
  scale_fill_manual(values=c("white","black")) +
  theme_bw() +
  labs(x=NULL, y=NULL)+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())+
  facet_wrap(~region, ncol=1) +
   guides(fill = guide_legend(reverse = TRUE)) +
  labs(x="Number of Range Edges") +
NULL
gg.niche.barplot
ggsave(gg.niche.barplot, filename=here("results","barplot_all_groups.png"),dpi=160, width=6, height=4)

# scatterplot 
# dat.thermal.gg <- spp.bayes.niche.lm.stats %>%
#   left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
#   ungroup() %>%
#   pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
#   mutate(region=factor(region, levels=c('neus','wc','ebs')),
#          region=recode(region,
#                        ebs="Eastern Bering Sea",
#                        neus="Northeast",
#                        wc="West Coast"),
#          niche.group=factor(niche.group, levels=c('good_tracker','partial_tracker','non_tracker')),
#          niche.group=recode(niche.group,
#                             good_tracker="TNH",
#                             partial_tracker="PTH",
#                             non_tracker="TIH")) %>%
#   ggplot(aes(x=mean_predict.sstmin, y=mean_predict.sstmax, group=niche.group)) +
#   geom_point(aes(shape=niche.group, color=niche.group, fill=niche.group)) +
#   scale_color_manual(values = c("TIH"="#fb8072","TNH"="deepskyblue4","PTH"="mediumpurple3")) +
#   geom_errorbar(aes(x=mean_predict.sstmin, ymin=lower_predict.sstmax, ymax=upper_predict.sstmax, color=niche.group)) +
#   geom_errorbarh(aes(y=mean_predict.sstmax, xmin=lower_predict.sstmin, xmax=upper_predict.sstmin, color=niche.group)) +
#   geom_vline(xintercept=0, linetype="dashed") +
#   geom_hline(yintercept=0, linetype="dashed") +
#   theme_bw() +
#   labs(x="Change in Cold Extreme at Edge (°C/year)", y="Change in Warm Extreme at Edge (°C/year)") +
#   facet_wrap(~region) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank()) +
#   NULL
# dat.thermal.gg
# ggsave(dat.thermal.gg, filename=here("results","thermal_niche_time_coefficients.png"),dpi=160, width=10, height=4)

# barplot by region 
gg.niche.regions <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(region=factor(region, levels=c('neus','wc','ebs')),
         region=recode(region,
                       ebs="Eastern Bering Sea",
                       neus="Northeast",
                       wc="West Coast"),
         niche.group=factor(niche.group, levels=c('good_tracker','partial_tracker','non_tracker')),
         niche.group=recode(niche.group,
                            good_tracker="TNH",
                            partial_tracker="PTH",
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

# barplot by edge type 
gg.niche.quantile <- spp.bayes.niche.lm.stats %>%
  left_join(spp.bayes.niche.groups, by=c("region","quantile","species")) %>%
  ungroup() %>%
  pivot_wider(names_from=predicted.var, values_from=c(mean,median,lower,upper)) %>%
  mutate(quantile=factor(quantile, levels=c('quantile_0.01','quantile_0.99')),
         quantile=recode(quantile,
                         quantile_0.01="Warm Edge",
                         quantile_0.99="Cold Edge"),
         niche.group=factor(niche.group, levels=c('good_tracker','partial_tracker','non_tracker')),
         niche.group=recode(niche.group,
                            good_tracker="TNH",
                            partial_tracker="PTH",
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
    niche.group=factor(niche.group, levels=c('good_tracker','partial_tracker','non_tracker')),
    niche.group=recode(niche.group,
                       good_tracker="TNH",
                       partial_tracker="PTH",
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

# # make example plots for talks--species choices might be a little out of date 
# ex.spp1 <- "gadus macrocephalus" # good tracker, cold edge, EBS
# ex.spp2 <- "sebastes pinniger" # non tracker, warm edge, WC
# ex.spp3 <- "paralichthys oblongus" # lagged tracker, cold edge, NEUS
# 
# ex.spp.bayes.gg1 <- spp.bayes.niche.filter %>% 
#   filter(species==ex.spp1) %>% 
#   group_by(.draw, predicted.var) %>%
#   mutate(mean.param = mean(year_match) ) %>%
#   ungroup() %>% 
#   select(.draw, mean.param, predicted.var) %>%
#   distinct() %>%
#   ggplot() +
#   theme_bw() +
#   geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
#   scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
#   labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
#   theme(legend.position="bottom") +
#   NULL
# ex.spp.bayes.gg1
# 
# ex.spp.bayes.gg2 <- spp.bayes.niche.filter %>% 
#   filter(species==ex.spp2) %>% 
#   group_by(.draw, predicted.var) %>%
#   mutate(mean.param = mean(year_match) ) %>%
#   ungroup() %>% 
#   select(.draw, mean.param, predicted.var) %>%
#   distinct() %>%
#   ggplot() +
#   theme_bw() +
#   geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
#   scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
#   labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
#   theme(legend.position="bottom") +
#   NULL
# ex.spp.bayes.gg2
# 
# ex.spp.bayes.gg3 <- spp.bayes.niche.filter %>% 
#   filter(species==ex.spp3) %>% 
#   group_by(.draw, predicted.var) %>%
#   mutate(mean.param = mean(year_match) ) %>%
#   ungroup() %>% 
#   select(.draw, mean.param, predicted.var) %>%
#   distinct() %>%
#   ggplot() +
#   theme_bw() +
#   geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
#   scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
#   labs(x="Posterior Distribution of Coefficient (°C/year)",y="Density", fill=NULL) +
#   theme(legend.position="bottom") +
#   NULL
# ex.spp.bayes.gg3
# 
# ex.spp.time.gg1 <- dat.predict.niche %>% 
#   filter(species==ex.spp1) %>%
#   mutate(species = str_to_sentence(species)) %>%
#   ggplot() +
#   geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
#   geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
#   scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
#   theme_bw() +
#   labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
#   theme(legend.position="bottom")+
#   scale_x_continuous(limits=c(1988, 2018), breaks=seq(1988, 2018, 4))+
#   NULL
# ex.spp.time.gg1
# 
# ex.spp.time.gg2 <- dat.predict.niche %>% 
#   filter(species==ex.spp2) %>%
#   mutate(species = str_to_sentence(species)) %>%
#   ggplot() +
#   geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
#   geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
#   scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
#   theme_bw() +
#   labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
#   theme(legend.position="bottom")+
#   scale_x_continuous(limits=c(1982, 2018), breaks=seq(1982, 2018, 4))+
#   NULL
# ex.spp.time.gg2
# 
# ex.spp.time.gg3 <- dat.predict.niche %>% 
#   filter(species==ex.spp3) %>%
#   mutate(species = str_to_sentence(species)) %>%
#   ggplot() +
#   geom_point(aes(x=year_match, y=sst, color=predicted.var)) +
#   geom_line(aes(x=year_match, y=sst, color=predicted.var)) +
#   scale_color_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
#   theme_bw() +
#   labs(x="Year",y="Sea Surface Temperature at Edge (°C)", color=NULL) +
#   theme(legend.position="bottom")+
#   scale_x_continuous(limits=c(1968, 2018), breaks=seq(1968, 2018, 5))+
#   NULL
# ex.spp.time.gg3
# 
# ggsave(ex.spp.bayes.gg1, dpi=160, width=4, height=4, filename=here("results","example_1_posterior.png"))
# ggsave(ex.spp.bayes.gg2, dpi=160, width=4, height=4, filename=here("results","example_2_posterior.png"))
# ggsave(ex.spp.bayes.gg3, dpi=160, width=4, height=4, filename=here("results","example_3_posterior.png"))
# ggsave(ex.spp.time.gg1, dpi=160, width=4, height=4, filename=here("results","example_1_niche.png"))
# ggsave(ex.spp.time.gg2, dpi=160, width=4, height=4, filename=here("results","example_2_niche.png"))
# ggsave(ex.spp.time.gg3, dpi=160, width=4, height=4, filename=here("results","example_3_niche.png"))
