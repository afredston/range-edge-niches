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
  ungroup() %>% # undo rowwise nature
  mutate(axis = as.character(axis)) %>% # convert from factor
  filter(axis %in% c('coast_km','line_km')) 

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
                         quantile_0.01="Equatorward Edge",
                         quantile_0.99="Poleward Edge"))


dat.summary.traits <- dat.models %>% 
  select(-Estimate, -Std.Error, -axis, -year) %>% # remove all cols that vary for a given edge
  distinct()

dat.models.groups <- dat.models %>%
  select(species, edgetype, region, taxongroup) %>%
  distinct() %>%
  mutate(species = as.factor(species), 
         edgetype=as.factor(edgetype),
         region=as.factor(region),
         taxongroup=as.factor(taxongroup))

#######################
### species shifts vs time
#######################

spp.bayes.edge.lm.df <- NULL
for(i in unique(dat.models$species)) {
  dfprep1 <- dat.models[dat.models$species==i,] # subdivide by species 
  for(j in unique(dfprep1$region)) {
    dfprep2 <- dat.models[dat.models$species==i & dat.models$region==j,] # subdivide by region, for the few spp found in multiple regions 
    for(k in unique(dfprep2$quantile)) {
      df <- dat.models[dat.models$species==i & dat.models$region==j & dat.models$quantile==k,] # subdivide by quantile, for the few spp with both range edges
      spp.bayes.lm <- try(stan_glm(Estimate ~ year, 
                                   data=df, 
                                   family=gaussian(), 
                                   iter = 40000,
                                   warmup = 10000,
                                   adapt_delta = 0.99,
                                   chains = 4,
                                   cores = 1,
                                   prior = normal(0, 50),
                                   weights = 1/(Std.Error^2)
      )) # centered on 0 to allow negative coefficients for some species. SD of 50 intended to exceed upper range of marine climate velocities (around 200 km/dec) in Burrows et al 2011
      if(!class(spp.bayes.lm)[1] == "try-error") { # adding try() here because some edges are so invariant that the model fails 
        spp.bayes.lm.tidy <- tidy_draws(spp.bayes.lm) %>%
          mutate(species = paste0(i),
                 region = paste0(j),
                 quantile = paste0(k),
                 intercept.rhat = summary(spp.bayes.lm)[,"Rhat"][1],
                 year_match.rhat = summary(spp.bayes.lm)[,"Rhat"][2],
                 sigma.rhat = summary(spp.bayes.lm)[,"Rhat"][3])
        spp.bayes.edge.lm.df <- rbind(spp.bayes.edge.lm.df, spp.bayes.lm.tidy)
      }
    }
  }
}

# check for convergence 
quantile(spp.bayes.edge.lm.df$intercept.rhat)
quantile(spp.bayes.edge.lm.df$year_match.rhat)  
quantile(spp.bayes.edge.lm.df$sigma.rhat) # this previously caused estimation problems, doesn't appear to be doing that anymore  

spp.bayes.edge.filter <- spp.bayes.edge.lm.df %>% 
  group_by(region, species, quantile) %>%
  mutate(max.rhat = max(intercept.rhat, year_match.rhat, sigma.rhat)) %>%
  filter(max.rhat <= 1.1) # get rid of spp*region*edge combos where one of the edge ~ time models didn't converge--just a check--may not get rid of any 

setdiff(spp.bayes.edge.lm.df %>% select(region, species, quantile) %>% distinct(), spp.bayes.edge.filter %>% select(region, species, quantile) %>% distinct()) # 0 at present, all models converged

rm(spp.bayes.edge.filter) # if they all converged this is the same as spp.bayes.edge.lm.df

# plot posteriors 
# bayes.lm.time.gg <- spp.bayes.edge.lm.df  %>%
#   group_by(.draw, region) %>%
#   mutate(mean.param = mean(year) ) %>%
#   ungroup() %>%
#   select(.draw, mean.param, region) %>%
#   distinct() %>%
#   ggplot() +
#   theme_bw() +
#   geom_density(aes(x=mean.param, fill=region), color="black", alpha=0.5) +
#   scale_fill_brewer(type="seq", palette="YlGnBu", labels=c("Eastern Bering Sea","Northeast","West Coast")) +
#   labs(x="Coefficient of Edge Position on Time (km/yr)", y="Density", fill="Region") +
#   theme(legend.position="bottom") +
#   NULL
# bayes.lm.time.gg
# ggsave(bayes.lm.time.gg, width=3, height=4, dpi=160, filename=here("results","edge_coefficients_time.png"), scale=1.6)

bayes.lm.time.edgetype.gg <- spp.bayes.edge.lm.df  %>%
  group_by(.draw, region, quantile) %>%
  mutate(mean.param = mean(year) ) %>%
  ungroup() %>%
  select(.draw, mean.param, region, quantile) %>%
  distinct() %>%
  mutate(quantile=recode(quantile,
                         quantile_0.01="Equatorward Edge",
                         quantile_0.99="Poleward Edge"),
         region=factor(region, levels=c('neus','wc','ebs'))) %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=region), color="black", alpha=0.5) +
  scale_fill_brewer(type="seq", palette="YlGnBu", labels=c("Northeast","West Coast","Eastern Bering Sea")) +
  labs(x="Coefficient of Edge Position on Time (km/yr)", y="Density", fill="Region") +
  theme(legend.position="bottom") +
  facet_wrap(~quantile) +
  NULL
bayes.lm.time.edgetype.gg
ggsave(bayes.lm.time.edgetype.gg, width=168, height=100,
       units="mm", dpi=600, filename=here("results","edge_coefficients_time_edgetype.png"))

# generate posteriors grouped different ways, some reported in-text--need full model output for this 

# edge shift stats pooled by region
spp.bayes.edge.lm.df %>%
  group_by(.draw, region) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(region) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))


# edge shift stats pooled by edgetype
spp.bayes.edge.lm.df %>%
  group_by(.draw, quantile) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(quantile) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))

# edge shift stats pooled by region and edgetype
spp.bayes.edge.lm.df %>%
  group_by(.draw, region, quantile) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(region, quantile) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))

# edge shift stats pooled by taxon group 
spp.bayes.edge.lm.df %>%
  left_join(dat.models.groups) %>%
  group_by(.draw, taxongroup) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(taxongroup) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))

# pooled by taxon group and edge type
spp.bayes.edge.lm.df %>%
  left_join(dat.models.groups) %>%
  group_by(.draw, taxongroup, quantile) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(taxongroup, quantile) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))

# species-specific results
spp.bayes.edge.lm.df.summary <- spp.bayes.edge.lm.df %>%
  group_by(.draw, species, region, quantile) %>%
  summarise(mean.year = mean(year)) %>%
  group_by(species, region, quantile) %>%
  summarise(mean=mean(mean.year),
            median=median(mean.year),
            lower=quantile(mean.year, 0.05),
            upper=quantile(mean.year, 0.95))

write_csv(spp.bayes.edge.lm.df.summary, here("results","species_edge_shifts_vs_time.csv"))

#######################
### predict temperatures at edges 
#######################

# split by region because temp datasets are different

# set up temp data
ebs.sst <- readRDS(here("processed-data","ebs_sst_linedist.rds"))
wc.sst <- readRDS(here("processed-data","wc_sst_coastdist.rds"))
neus.sst <- readRDS(here("processed-data","neus_sst_coastdist.rds"))

# prepare df of all SST values at all edge positions 
ebs.sst.prepgam <- ebs.sst %>% 
  group_by(line_km, year_match) %>%
  mutate(sstmean = mean(sst),
         sstmax = max(sst),
         sstmin = min(sst)) %>%
  ungroup() %>%
  dplyr::select(year_match, line_km, sstmean, sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

wc.sst.prepgam <- wc.sst %>% 
  group_by(coast_km, year_match) %>%
  mutate(sstmean = mean(sst),
         sstmax = max(sst),
         sstmin = min(sst)) %>%
  ungroup() %>%
  dplyr::select(year_match, coast_km, sstmean,sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

neus.sst.prepgam <- neus.sst %>% 
  group_by(coast_km, year_match) %>%
  mutate(sstmean = mean(sst),
         sstmax = max(sst),
         sstmin = min(sst)) %>%
  ungroup() %>%
  dplyr::select(year_match, coast_km, sstmean,sstmax, sstmin) %>%
  distinct() %>%
  mutate(year_match = as.factor(year_match)) 

# set up GAMs to predict temperature at position of edge this year 
ebs.sst.temp.gam.max <- gam(sstmax ~ year_match + s(line_km, by=year_match), data=ebs.sst.prepgam)
ebs.sst.temp.gam.min <- gam(sstmin ~ year_match + s(line_km, by=year_match), data=ebs.sst.prepgam)

wc.sst.temp.gam.max <- gam(sstmax ~ year_match + s(coast_km, by=year_match), data=wc.sst.prepgam)
wc.sst.temp.gam.min <- gam(sstmin ~ year_match + s(coast_km, by=year_match), data=wc.sst.prepgam)

neus.sst.temp.gam.max <- gam(sstmax ~ year_match + s(coast_km, by=year_match), data=neus.sst.prepgam)
neus.sst.temp.gam.min <- gam(sstmin ~ year_match + s(coast_km, by=year_match), data=neus.sst.prepgam)

# predict temp from edge position--prep datasets

ebs.pred <- dat.models %>%
  filter(region=="ebs") %>%
  dplyr::select( -axis) %>%
  rename(line_km=Estimate,
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

# set up dfs with same column names but lagged edge positions to predict temperature this year at the position where a species' edge was last year
ebs.pred.lag <- ebs.pred %>% 
  arrange(year_match) %>% 
  group_by(species, quantile) %>% 
  mutate(line_km_lag = lag(line_km)) %>% 
  ungroup() %>% 
  select(-line_km) %>% 
  rename(line_km = line_km_lag) %>% 
  filter(year_match > min(year_match))

wc.pred.lag <- wc.pred %>% 
  arrange(year_match) %>% 
  group_by(species, quantile) %>% 
  mutate(coast_km_lag = lag(coast_km)) %>% 
  ungroup() %>% 
  select(-coast_km) %>% 
  rename(coast_km = coast_km_lag) %>% 
  filter(year_match > min(year_match))

neus.pred.lag <- neus.pred %>% 
  arrange(year_match) %>% 
  group_by(species, quantile) %>% 
  mutate(coast_km_lag = lag(coast_km)) %>% 
  ungroup() %>% 
  select(-coast_km) %>% 
  rename(coast_km = coast_km_lag) %>% 
  filter(year_match > min(year_match))

# add columns with predicted temperature at edge every year 
ebs.pred$predict.sstmax <- predict.gam(ebs.sst.temp.gam.max, ebs.pred)
ebs.pred$predict.sstmax.se <- predict.gam(ebs.sst.temp.gam.max, ebs.pred,se.fit=TRUE)$se.fit
ebs.pred$predict.sstmin <- predict.gam(ebs.sst.temp.gam.min, ebs.pred)
ebs.pred$predict.sstmin.se <- predict.gam(ebs.sst.temp.gam.min, ebs.pred,se.fit=TRUE)$se.fit

wc.pred$predict.sstmax <- predict.gam(wc.sst.temp.gam.max, wc.pred)
wc.pred$predict.sstmax.se <- predict.gam(wc.sst.temp.gam.max, wc.pred,se.fit=TRUE)$se.fit
wc.pred$predict.sstmin <- predict.gam(wc.sst.temp.gam.min, wc.pred)
wc.pred$predict.sstmin.se <- predict.gam(wc.sst.temp.gam.min, wc.pred, se.fit=TRUE)$se.fit

neus.pred$predict.sstmax <- predict.gam(neus.sst.temp.gam.max, neus.pred)
neus.pred$predict.sstmax.se <- predict.gam(neus.sst.temp.gam.max, neus.pred, se.fit=TRUE)$se.fit
neus.pred$predict.sstmin <- predict.gam(neus.sst.temp.gam.min, neus.pred)
neus.pred$predict.sstmin.se <- predict.gam(neus.sst.temp.gam.min, neus.pred, se.fit=TRUE)$se.fit

neus.pred <- rename(neus.pred, edge_position=coast_km)
neus.pred$axis <- "coast_km"

wc.pred <- rename(wc.pred, edge_position=coast_km)
wc.pred$axis <- "coast_km"

ebs.pred <- rename(ebs.pred, edge_position=line_km)
ebs.pred$axis <- "line_km"

# tidy columns and combine
dat.predict1 <- rbind(neus.pred, wc.pred, ebs.pred)%>%
  select(-predict.sstmax.se, -predict.sstmin.se) %>%
  pivot_longer(cols=c(predict.sstmax, predict.sstmin), names_to="predicted.var",values_to="sst") 

dat.predict <- rbind(neus.pred, wc.pred, ebs.pred)%>%
  select(-predict.sstmax, -predict.sstmin) %>%
  pivot_longer(cols=c(predict.sstmax.se, predict.sstmin.se), names_to="predicted.var",values_to="sstSE") %>%
  mutate(predicted.var=str_replace(predicted.var, ".se","")) %>%
  inner_join(dat.predict1)

write_csv(dat.predict, here("processed-data","species_thermal_niche_v_time.csv"))

# exploring isotherm tracking: what is the predicted temperature this year at last year's edge position?
ebs.pred.lag$predict.sstmax <- predict.gam(ebs.sst.temp.gam.max, ebs.pred.lag)
ebs.pred.lag$predict.sstmax.se <- predict.gam(ebs.sst.temp.gam.max, ebs.pred.lag,se.fit=TRUE)$se.fit
ebs.pred.lag$predict.sstmin <- predict.gam(ebs.sst.temp.gam.min, ebs.pred.lag)
ebs.pred.lag$predict.sstmin.se <- predict.gam(ebs.sst.temp.gam.min, ebs.pred.lag,se.fit=TRUE)$se.fit

wc.pred.lag$predict.sstmax <- predict.gam(wc.sst.temp.gam.max, wc.pred.lag)
wc.pred.lag$predict.sstmax.se <- predict.gam(wc.sst.temp.gam.max, wc.pred.lag,se.fit=TRUE)$se.fit
wc.pred.lag$predict.sstmin <- predict.gam(wc.sst.temp.gam.min, wc.pred.lag)
wc.pred.lag$predict.sstmin.se <- predict.gam(wc.sst.temp.gam.min, wc.pred.lag, se.fit=TRUE)$se.fit

neus.pred.lag$predict.sstmax <- predict.gam(neus.sst.temp.gam.max, neus.pred.lag)
neus.pred.lag$predict.sstmax.se <- predict.gam(neus.sst.temp.gam.max, neus.pred.lag, se.fit=TRUE)$se.fit
neus.pred.lag$predict.sstmin <- predict.gam(neus.sst.temp.gam.min, neus.pred.lag)
neus.pred.lag$predict.sstmin.se <- predict.gam(neus.sst.temp.gam.min, neus.pred.lag, se.fit=TRUE)$se.fit

# tidy results and label lag columns appropriately
neus.pred.lag <- rename(neus.pred.lag, edge_position_last_year=coast_km)
neus.pred.lag$axis <- "coast_km"

wc.pred.lag <- rename(wc.pred.lag, edge_position_last_year=coast_km)
wc.pred.lag$axis <- "coast_km"

ebs.pred.lag <- rename(ebs.pred.lag, edge_position_last_year=line_km)
ebs.pred.lag$axis <- "line_km"

# first tidy means
dat.predict1.lag <- rbind(neus.pred.lag, wc.pred.lag, ebs.pred.lag)%>%
  select(-predict.sstmax.se, -predict.sstmin.se) %>%
  pivot_longer(cols=c(predict.sstmax, predict.sstmin), names_to="predicted.var",values_to="sst_at_last_years_edge")

# then tidy SEs and join 
dat.predict.lag.prep <- rbind(neus.pred.lag, wc.pred.lag, ebs.pred.lag)%>%
  select(-predict.sstmax, -predict.sstmin) %>%
  pivot_longer(cols=c(predict.sstmax.se, predict.sstmin.se), names_to="predicted.var",values_to="sst_SE_at_last_years_edge") %>%
  mutate(predicted.var=str_replace(predicted.var, ".se","")) %>%
  inner_join(dat.predict1.lag) %>% 
  select(region, species, quantile, year_match, predicted.var, axis, sst_at_last_years_edge, sst_SE_at_last_years_edge, edge_position_last_year) %>% 
  distinct()

dat.predict.lag <- dat.predict %>% 
  left_join(dat.predict.lag.prep)

write_csv(dat.predict.lag, here("processed-data","species_thermal_niche_v_time_with_lags.csv"))

#######################
### estimate change in edge thermal niche over time
#######################

# Bayesian test for edge thermal niche change over time 

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
                                     prior = normal(0, 0.1), # exceeds highest rates of warming we found in the paper which were around 0.04 C/yr
                                     control = list(max_treedepth = 20),
                                     weights = 1/(sstSE^2)
        ))
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
  filter(max.rhat <= 1.1) # get rid of spp*region*edge combos where one of the SST extreme models didn't converge--just a check--may not get rid of any 

setdiff(spp.bayes.niche.lm.df %>% select(region, species, quantile) %>% distinct(), spp.bayes.niche.filter %>% select(region, species, quantile) %>% distinct()) # 0 at present 

# SLOW

# summarize posterior distributions and write out 
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
write_csv(spp.bayes.niche.lm.stats, here("results","species_bayes_niche_lm_summary.csv"))

quantile(spp.bayes.niche.lm.stats$mean) # should be distributed in the neighborhood of zero
spp.bayes.niche.lm.stats %>% 
  ggplot() +
  geom_histogram(aes(x=mean))

##########################
# Figure 1 example plots
##########################

# while most plots are generated in figure-scripts, these require the full STAN output to generate posteriors 
# the other small plots in the methods figure are generated in the respective figure-scripts files (e.g., time-series of range edges are generated in figure-scripts/species_edges_vs_time.R)

# non tracker - lobster - neus 

# make example plots for methods schematic 
ex.spp.ebs <- "paralithodes camtschaticus" # red king crab
ex.spp.neus <- "gadus morhua" # atlantic cod
ex.spp.wc <- "sebastes pinniger" # canary rockfish

summary.spp.ebs <- spp.bayes.niche.lm.stats %>% 
  filter(species==ex.spp.ebs)

summary.spp.neus <- spp.bayes.niche.lm.stats %>% 
  filter(species==ex.spp.neus)

summary.spp.wc <- spp.bayes.niche.lm.stats %>% 
  filter(species==ex.spp.wc, region=="wc")

# if species are updated, be sure to change year limits in the time-series figures below 

ex.spp.bayes.gg.ebs <- spp.bayes.niche.filter %>%
  filter(species==ex.spp.ebs) %>%
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>%
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  geom_vline(aes(xintercept=0), color="black", linetype="dotted") + 
  geom_segment(aes(x=summary.spp.ebs[summary.spp.ebs$predicted.var=="predict.sstmax",]$lower, xend=summary.spp.ebs[summary.spp.ebs$predicted.var=="predict.sstmax",]$upper, y=-1, yend=-1), color="#DF2301", lwd=2) + # add 90% credible interval 
  geom_segment(aes(x=summary.spp.ebs[summary.spp.ebs$predicted.var=="predict.sstmin",]$lower, xend=summary.spp.ebs[summary.spp.ebs$predicted.var=="predict.sstmin",]$upper, y=0, yend=0), color="#3A4ED0", lwd=2) +
  scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Coefficient (°C/year)",y="Density", fill=NULL) +
  scale_x_continuous(breaks=c(-0.1, -0.05, 0, 0.05)) +
  theme(legend.position="none") +
  NULL
ex.spp.bayes.gg.ebs

ex.spp.bayes.gg.neus <- spp.bayes.niche.filter %>%
  filter(species==ex.spp.neus) %>%
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>%
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  geom_vline(aes(xintercept=0), color="black", linetype="dotted") + 
  geom_segment(aes(x=summary.spp.neus[summary.spp.neus$predicted.var=="predict.sstmax",]$lower, xend=summary.spp.neus[summary.spp.neus$predicted.var=="predict.sstmax",]$upper, y=-1, yend=-1), color="#DF2301", lwd=2) + # add 90% credible interval 
  geom_segment(aes(x=summary.spp.neus[summary.spp.neus$predicted.var=="predict.sstmin",]$lower, xend=summary.spp.neus[summary.spp.neus$predicted.var=="predict.sstmin",]$upper, y=0, yend=0), color="#3A4ED0", lwd=2) +scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Coefficient (°C/year)",y="Density", fill=NULL) +
  scale_x_continuous(breaks=c(-0.1, -0.05, 0, 0.05)) +
  theme(legend.position="none") +
  NULL
ex.spp.bayes.gg.neus

ex.spp.bayes.gg.wc <- spp.bayes.niche.filter %>%
  filter(species==ex.spp.wc) %>%
  group_by(.draw, predicted.var) %>%
  mutate(mean.param = mean(year_match) ) %>%
  ungroup() %>%
  select(.draw, mean.param, predicted.var) %>%
  distinct() %>%
  ggplot() +
  theme_bw() +
  geom_density(aes(x=mean.param, fill=predicted.var), color="black", alpha=0.5) +
  geom_vline(aes(xintercept=0), color="black", linetype="dotted") + 
  geom_segment(aes(x=summary.spp.wc[summary.spp.wc$predicted.var=="predict.sstmax",]$lower, xend=summary.spp.wc[summary.spp.wc$predicted.var=="predict.sstmax",]$upper, y=-1, yend=-1), color="#DF2301", lwd=2) + # add 90% credible interval 
  geom_segment(aes(x=summary.spp.wc[summary.spp.wc$predicted.var=="predict.sstmin",]$lower, xend=summary.spp.wc[summary.spp.wc$predicted.var=="predict.sstmin",]$upper, y=0, yend=0), color="#3A4ED0", lwd=2) +scale_fill_manual(values=c("#DF2301","#3A4ED0"), labels=c("Warm Extreme","Cold Extreme")) +
  labs(x="Coefficient (°C/year)",y="Density", fill=NULL) +
  scale_x_continuous(breaks=c(-0.1, -0.05, 0, 0.05)) +
  theme(legend.position="none") +
  NULL
ex.spp.bayes.gg.wc

ggsave(ex.spp.bayes.gg.ebs, dpi=600, width=1.5, height=1.4, filename=here("results",paste0("example_posterior_",ex.spp.ebs, ".png")),scale=1.5)
ggsave(ex.spp.bayes.gg.neus, dpi=600, width=1.5, height=1.4, filename=here("results",paste0("example_posterior_",ex.spp.neus,".png")),scale=1.5)
ggsave(ex.spp.bayes.gg.wc, dpi=600, width=1.5, height=1.4, filename=here("results",paste0("example_posterior_",ex.spp.wc,".png")),scale=1.5)


##########################
# Traits analysis for appendix
##########################

## EDGE SHIFTS ~ TRAITS
# spp.bayes.edge.traits <- spp.bayes.edge.lm.df %>% 
#   left_join(dat.summary.traits)
# 
# gg.edge.habitat <- spp.bayes.edge.traits %>% 
#   filter(!is.na(habitat)) %>%  # filter rows with NAs -- leaves only fishes
#   group_by(.draw, region, habitat) %>% 
#   mutate(mean.param = mean(year)) %>% 
#   ungroup() %>% 
#   select(.draw, region, habitat, mean.param) %>% 
#   distinct() %>%
#   mutate(
#     region=recode(region,
#                   ebs="Eastern Bering Sea",
#                   neus="Northeast",
#                   wc="West Coast")) %>% 
#   ggplot() +
#   theme_bw() +
#   facet_wrap(~region) +
#   geom_density(aes(x=mean.param, fill=habitat), color="black", alpha=0.5) +
#   labs(x="Edge Shift (km/yr)", y="Density") +
#   theme(legend.position="bottom") +
#   NULL
# ggsave(gg.edge.habitat, width=6, height=2, filename=here("results","edge_shifts_v_habitat.png"), dpi=160, scale=1.2)
# 
# gg.edge.spawning.type <- spp.bayes.edge.traits %>% 
#   filter(!is.na(spawning.type)) %>%  # filter rows with NAs -- leaves only fishes
#   group_by(.draw, region, spawning.type) %>% 
#   mutate(mean.param = mean(year)) %>% 
#   ungroup() %>% 
#   select(.draw, region, spawning.type, mean.param) %>% 
#   distinct() %>%
#   mutate(
#     region=recode(region,
#                   ebs="Eastern Bering Sea",
#                   neus="Northeast",
#                   wc="West Coast")) %>% 
#   ggplot() +
#   theme_bw() +
#   facet_grid(~ region) +
#   geom_density(aes(x=mean.param, fill=spawning.type), color="black", alpha=0.5) +
#   labs(x="Edge Shift (km/yr)", y="Density") +
#   theme(legend.position="bottom") +
#   NULL
# ggsave(gg.edge.spawning.type, width=6, height=2, filename=here("results","edge_shifts_v_spawn_type.png"), dpi=160, scale=1.2)
# 
# gg.edge.fecundity <-  spp.bayes.edge.traits %>% 
#   filter(!is.na(fecundity)) %>%
#   mutate(fecundity.q = ifelse(fecundity < quantile(dat.summary.traits$fecundity, 0.25, na.rm=TRUE), "0-0.25",ifelse(fecundity < quantile(dat.summary.traits$fecundity, 0.5, na.rm=TRUE), "0.25-0.5", ifelse(fecundity < quantile(dat.summary.traits$fecundity, 0.75, na.rm=TRUE),"0.5-0.75", "0.75-1")))) %>%  
#   group_by(.draw, region, fecundity.q) %>% 
#   mutate(mean.param = mean(year)) %>% 
#   ungroup() %>% 
#   select(.draw, region, fecundity.q, mean.param) %>% 
#   distinct() %>%
#   mutate(
#     region=recode(region,
#                   ebs="Eastern Bering Sea",
#                   neus="Northeast",
#                   wc="West Coast")) %>% 
#   ggplot() +
#   theme_bw() +
#  facet_grid(~ region) +
#   geom_density(aes(x=mean.param, fill=fecundity.q), color="black", alpha=0.5) +
#   labs(x="Edge Shift (km/yr)", y="Density") +
#   theme(legend.position="bottom") +
#   NULL
# ggsave(gg.edge.fecundity, width=6, height=2, filename=here("results","edge_shifts_v_fecundity.png"), dpi=160, scale=1.2)
# 
# gg.edge.tl <- spp.bayes.edge.lm.df.summary %>% 
#   left_join(dat.summary.traits)  %>%
#   mutate(
#     region=recode(region,
#                   ebs="Eastern Bering Sea",
#                   neus="Northeast",
#                   wc="West Coast")) %>%  
#   ggplot() +
#   theme_bw() +
#   geom_point(aes(x=tl, y=mean)) +
#   geom_errorbar(aes(x=tl, y=mean, ymin=lower, ymax=upper)) +
#   facet_wrap(~region) +
#   labs(x="Trophic Level", y="Edge Shift (km/yr)") +
#   NULL
# ggsave(gg.edge.tl, width=6, height=2, filename=here("results","edge_shifts_v_trophic_level.png"), dpi=160, scale=1.2)

# NICHE SHIFTS ~ TRAITS

spp.bayes.niche.traits <- spp.bayes.niche.filter %>% 
  left_join(dat.summary.traits)

gg.niche.habitat <- spp.bayes.niche.traits %>% 
  filter(!is.na(habitat)) %>%  
  group_by(.draw, region, predicted.var, habitat) %>% 
  mutate(mean.param = mean(year_match)) %>% 
  ungroup() %>% 
  select(.draw, region, predicted.var, habitat, mean.param) %>% 
  distinct() %>%
  mutate(
    region=recode(region,
                  ebs="Eastern Bering Sea",
                  neus="Northeast",
                  wc="West Coast"),
    predicted.var=recode(predicted.var,
                         predict.sstmax="Summer SST",
                         predict.sstmin="Winter SST")) %>% 
  ggplot() +
  theme_bw() +
  facet_grid(predicted.var~region) +
  geom_density(aes(x=mean.param, fill=habitat), color="black", alpha=0.5) +
  labs(x="Niche Shift (°C/yr)", y="Density") +
  theme(legend.position="bottom") +
  NULL
ggsave(gg.niche.habitat, width=6, height=3.5, filename=here("results","niche_shifts_v_habitat.png"), dpi=160, scale=1.3)


gg.niche.spawning.type <- spp.bayes.niche.traits %>% 
  filter(!is.na(spawning.type)) %>%  
  group_by(.draw, region, predicted.var, spawning.type) %>% 
  mutate(mean.param = mean(year_match)) %>% 
  ungroup() %>% 
  select(.draw, region, predicted.var, spawning.type, mean.param) %>% 
  distinct() %>%
  mutate(
    region=recode(region,
                  ebs="Eastern Bering Sea",
                  neus="Northeast",
                  wc="West Coast"),
    predicted.var=recode(predicted.var,
                         predict.sstmax="Summer SST",
                         predict.sstmin="Winter SST")) %>% 
  ggplot() +
  theme_bw() +
  facet_grid(predicted.var~region) +
  geom_density(aes(x=mean.param, fill=spawning.type), color="black", alpha=0.5) +
  labs(x="Niche Shift (°C/yr)", y="Density") +
  theme(legend.position="bottom") +
  NULL
ggsave(gg.niche.spawning.type, width=6, height=3.5, filename=here("results","niche_shifts_v_spawn_type.png"), dpi=160, scale=1.3)

gg.niche.fecundity <-  spp.bayes.niche.traits %>% 
  filter(!is.na(fecundity)) %>%
  mutate(fecundity.q = ifelse(fecundity < quantile(dat.summary.traits$fecundity, 0.25, na.rm=TRUE), "0-0.25",ifelse(fecundity < quantile(dat.summary.traits$fecundity, 0.5, na.rm=TRUE), "0.25-0.5", ifelse(fecundity < quantile(dat.summary.traits$fecundity, 0.75, na.rm=TRUE),"0.5-0.75", "0.75-1")))) %>%  
  group_by(.draw, region, predicted.var, fecundity.q) %>% 
  mutate(mean.param = mean(year_match)) %>% 
  ungroup() %>% 
  select(.draw, region, predicted.var, fecundity.q, mean.param) %>% 
  distinct() %>%
  mutate(
    region=recode(region,
                  ebs="Eastern Bering Sea",
                  neus="Northeast",
                  wc="West Coast"),
    predicted.var=recode(predicted.var,
                         predict.sstmax="Summer SST",
                         predict.sstmin="Winter SST")) %>% 
  rename("Fecundity Quantile"=fecundity.q) %>% 
  ggplot() +
  theme_bw() +
  facet_grid(predicted.var~ region) +
  geom_density(aes(x=mean.param, fill=`Fecundity Quantile`), color="black", alpha=0.5) +
  labs(x="Niche Shift (°C/yr)", y="Density") +
  theme(legend.position="bottom") +
  NULL
ggsave(gg.niche.fecundity, width=6, height=3.5, filename=here("results","niche_shifts_v_fecundity.png"), dpi=160, scale=1.3)

gg.niche.tl <- spp.bayes.niche.lm.stats %>% 
  left_join(dat.summary.traits)  %>%
  mutate(
    region=recode(region,
                  ebs="Eastern Bering Sea",
                  neus="Northeast",
                  wc="West Coast"),
    predicted.var=recode(predicted.var,
                         predict.sstmax="Summer SST",
                         predict.sstmin="Winter SST")) %>%  
  ggplot() +
  theme_bw() +
  geom_point(aes(x=tl, y=mean)) +
  geom_errorbar(aes(x=tl, y=mean, ymin=lower, ymax=upper)) +
  facet_grid(predicted.var~region) +
  labs(x="Trophic Level", y="Niche Shift (°C/yr)") +
  NULL
ggsave(gg.niche.tl, width=6, height=3.5, filename=here("results","niche_shifts_v_trophic_level.png"), dpi=160, scale=1.3)


gg.niche.len <- spp.bayes.niche.lm.stats %>% 
  left_join(dat.summary.traits)  %>%
  mutate(
    region=recode(region,
                  ebs="Eastern Bering Sea",
                  neus="Northeast",
                  wc="West Coast"),
    predicted.var=recode(predicted.var,
                         predict.sstmax="Summer SST",
                         predict.sstmin="Winter SST")) %>%  
  ggplot() +
  theme_bw() +
  geom_point(aes(x=length.infinity, y=mean)) +
  geom_errorbar(aes(x=length.infinity, y=mean, ymin=lower, ymax=upper)) +
  facet_grid(predicted.var~region) +
  labs(x="Infinity Length (cm)", y="Niche Shift (°C/yr)") +
  NULL
ggsave(gg.niche.len, width=6, height=3.5, filename=here("results","niche_shifts_v_length.png"), dpi=160, scale=1.3)


# temporarily saving the model outputs too, for further traits analysis 
#saveRDS(spp.bayes.edge.lm.df, paste0(getwd(), "/processed-data/full_output_edges_v_time.rds"))
#saveRDS(spp.bayes.niche.filter, paste0(getwd(), "/processed-data/full_output_niches_v_time.rds"))
