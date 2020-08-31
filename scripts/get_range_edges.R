# fit single-species models and get range limits in all study regions using VAST 

########################
### load packages, functions 
########################

library(VAST)
library(FishStatsUtils)
library(TMB)
library(FishData)
Version = get_latest_version( package="VAST" )
source(paste0(getwd(),'/functions/get_range_edge.R'))
source(paste0(getwd(),'/functions/get_density.R'))
library(doParallel)
library(foreach)

########################
### model settings
########################

# see Methods and VAST documentation for more background
quantiles_of_interest = c(0.01, 0.5, 0.99)

n_x = 100 
fine_scale = TRUE

FieldConfig = c("Omega1"=0, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=1)

RhoConfig = c("Beta1"=4, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(2,1)

Options =  c("Calculate_Range"=1)

strata.limits <- data.frame('STRATA'="All_areas")

Species_set = 100

Regions = c("Northwest_Atlantic","California_current","Eastern_Bering_Sea")
regs = c("neus","wc","ebs")

########################
### start for loop for each region
########################

for(i in 1:length(Regions)){
  Region = Regions[i]
  reg = regs[i]
  
  # delete spatial file; this needs to be re-generated for each new region
  unlink(paste0(getwd(),"/Kmeans-100.Rdata")) 
  
  # create file to store output 
  RegionFile = paste0(getwd(),'/VAST_',reg,'/')
  if (!dir.exists(RegionFile)){
    dir.create(RegionFile)
  } else {
    print(paste0(RegionFile, " already exists"))
  }
  
  
  ########################
  ### get catch data, region-specific
  ########################
  
  # gets rid of taxa that are not IDed to species, and harmonizes some taxonomic names 
  if(reg=="neus"){
    
    Species_omit <- c("Crustacea shrimp","Myctophidae spp.")
    
    Catch_rates <- as.data.frame(readRDS(paste0(getwd(),"/processed-data/neus_catch_rates.rds")))
    
    # fix incorrect taxa 
    Catch_rates$Sci <- gsub("Peprilus tracanthus", "Peprilus triacanthus", Catch_rates$Sci)
    Catch_rates$Sci <- gsub("Geryon quinquedens", "Chaceon quinquedens", Catch_rates$Sci)
    Catch_rates$Year <- as.numeric(Catch_rates$Year)
    
    Catch_rates <- Catch_rates[!Catch_rates$Sci %in% Species_omit,]
    
    Species_list <- unique(Catch_rates$Sci) 
    
  }
  
  if(reg=="wc"){
    Species_omit <- c("Actiniaria","Scyphozoa","Thaliacea","Myctophidae")
    
    Catch_rates_wcg <- download_catch_rates(survey='West_coast_groundfish_bottom_trawl_survey', species_set=Species_set, measurement_type="biomass", add_zeros = TRUE, error_tol=0.1)   
    Catch_rates_wcg$AreaSwept_ha <- NULL
    Catch_rates_wcg$Survey <- "annual"
    Catch_rates_wct <- download_catch_rates(survey='West_coast_triennial', species_set=Species_set, measurement_type="biomass", add_zeros = TRUE) 
    Catch_rates_wct$AreaSwept_ha <- NULL
    Catch_rates_wct$Survey <- "triennial"
    
    Species_list <- intersect(unique(Catch_rates_wcg$Sci), unique(Catch_rates_wct$Sci)) # keep only species found in both surveys
    
    Species_list <- setdiff(Species_list, Species_omit)
    
    Catch_rates_wcg <- Catch_rates_wcg[Catch_rates_wcg$Sci %in% Species_list,]
    Catch_rates_wct <- Catch_rates_wct[Catch_rates_wct$Sci %in% Species_list,]
    
    Catch_rates <- rbind(Catch_rates_wct, Catch_rates_wcg)
    
    # fix incorrect taxa
    Catch_rates$Sci <- gsub("Doryteuthis_opalescens", "Loligo_opalescens", Catch_rates$Sci)
    Catch_rates$Sci <- gsub("Bathyraja_kincaidii", "Bathyraja_interrupta", Catch_rates$Sci)
    Catch_rates$Sci <- gsub("Cancer_magister", "Metacarcinus_magister", Catch_rates$Sci)
    
    Species_list <- unique(Catch_rates$Sci)
  }
  
  if(reg=="ebs"){
    Species_omit <- c("Nudibranchia","Gastropoda","Bryozoa","Rajidae","Porifera","Actiniaria","Scyphozoa","Paguridae","gastropod_egg","Henricia_sp.","Metridium_sp.","Volutopsius_sp.","Aplidium_sp._A_(Clark_2006)","Argis_sp.","Buccinum_sp.","Crangon_sp.","Gersemia_sp.","Lepidopsetta_sp.","Chionoecetes_hybrid") 
    
    Catch_rates <- download_catch_rates(survey="Eastern_Bering_Sea", species_set=Species_set) 
    Catch_rates$Sci <- as.character(Catch_rates$Sci)
    
    # filter Catch_rates 
    Catch_rates <- Catch_rates[!Catch_rates$Sci %in% Species_omit,] 
    Catch_rates <- Catch_rates[Catch_rates$Year>=1989,] # eliminating early years for now because of sampling design change
    Species_list <- unique(Catch_rates$Sci) 
  }
  
  ########################
  ### other model set-up
  ########################
  
  settings = make_settings( "n_x"=n_x, "Region"=Region, purpose="index", bias.correct=FALSE, use_anisotropy=FALSE,
                            "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "OverdispersionConfig"=OverdispersionConfig,
                            "Options"=Options,"ObsModel"=ObsModel  )
  
  # derived objects that are constant among all models in a region
  Year_Set <- sort(unique(Catch_rates$Year))
  
  # check if axis conversion dataframe exists, set up code to generate if not (has to be derived from the VAST model further down)
  # this dataframe will provide conversions between VAST defaults (easting/northing), lat/lon, and the regional coastal distance / NW axis
  Z_gmFile <- paste0(RegionFile,'Z_gm.rds')
  if (!file.exists(Z_gmFile) & reg %in% c('neus','wc')){
    # load in coastal distance data
    coastdistdat <- readRDS(paste0(getwd(),'/processed-data/',reg,'_coastdistdat.rds'))
    
    # function to get coastal length of lat/lon coords
    get_length <- function(lon, lat, distdf) {
      tmp <- distdf %>%
        mutate(abs.diff.x2 = abs(x-lon)^2,
               abs.diff.y2 = abs(y-lat)^2,
               abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
               )) %>%
        filter(abs.diff.xy == min(abs.diff.xy)) %>%
        dplyr::select(lengthfromhere) %>%
        pull()
      return(tmp)
    }
    # separate for EBS which doesn't use coastal distance
    
  } else if (!file.exists(Z_gmFile) & reg=="ebs"){
    # load in axis distance data
    axisdistdat <- readRDS(paste0(getwd(),'/processed-data/',reg,'_axisdistdat.rds'))
    
    # function to get coastal length of lat/lon coords
    get_length <- function(lon, lat, distdf) {
      tmp <- distdf %>%
        mutate(abs.diff.x2 = abs(x-lon)^2,
               abs.diff.y2 = abs(y-lat)^2,
               abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
               )) %>%
        filter(abs.diff.xy == min(abs.diff.xy)) %>%
        dplyr::select(lengthfromhere) %>%
        pull()
      return(tmp)
    }
  } 
  
  ########################
  ###  start single-species models in parallel
  ########################
  
  detectCores()
  registerDoParallel(18) # UPDATE FOR YOUR OWN MACHINE
  
  foreach(j = Species_list) %dopar%{
    library(VAST)
    library(TMB)
    
    ########################
    ###  prep catch data 
    ########################
    
    DF <- Catch_rates[Catch_rates$Sci==j,]
    
    if(reg %in% c('neus','ebs')){
      Data_Geostat = data.frame( "spp"=DF[,"Sci"], "Year"=DF[,"Year"], "Catch_KG"=DF[,"Wt"], "AreaSwept_km2"=0.01, "Vessel"=0, "Lat"=DF[,"Lat"], "Lon"=DF[,"Long"] )}
    
    # west coast needs Survey column for Q_ik catchability covariate 
    if(reg=="wc"){
      Data_Geostat = data.frame( "spp"=DF[,"Sci"], "Year"=DF[,"Year"], "Catch_KG"=DF[,"Wt"], "AreaSwept_km2"=0.01, "Vessel"=0, "Lat"=DF[,"Lat"], "Lon"=DF[,"Long"], "Survey"=DF[,"Survey"] )}
    
    # below I'm getting rid of years where the species was observed 100% of the time 
    Data_Geostat$nonzero = ifelse(Data_Geostat$Catch_KG > 0, "Y","N")
    
    Data_Geostat$prop_nonzero = ave(Data_Geostat$nonzero, Data_Geostat$Year, FUN=function(x) mean(x=='Y'))
    
    Data_Geostat <- Data_Geostat[Data_Geostat$prop_nonzero<1,]
    Data_Geostat <- Data_Geostat[Data_Geostat$prop_nonzero>0,] # also drop years where species was never observed
    
    Years2Include <- sort(unique(Data_Geostat$Year))
    
    # Pad first years
    Data_Placeholder = Data_Geostat[1:5,c('Year','Catch_KG','AreaSwept_km2','Lat','Lon')]
    Data_Placeholder[,'Catch_KG'] = NA
    Data_Placeholder[,'Year'] = min(Data_Placeholder[,'Year']) - 5:1
    
    # separate WC again to preserve Survey column for Q_ik
    if(reg %in% c('neus','ebs')){
      Data_Geostat <- rbind( Data_Placeholder, Data_Geostat[c('Year','Catch_KG','AreaSwept_km2','Lat','Lon')] )}
    if(reg=='wc'){
      Data_Placeholder[,'Survey'] = "triennial"
      Data_Geostat <- rbind( Data_Placeholder, Data_Geostat[c('Year','Catch_KG','AreaSwept_km2','Lat','Lon','Survey')] )}
    
    
    ########################
    ###  Build model without running, and extract objects
    ########################
    
    fit = 
      fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
                "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                "a_i"=Data_Geostat[,'AreaSwept_km2'],
                fine_scale=fine_scale,
                anisotropy=FALSE,
                run_model=FALSE
      )
    
    # WEST COAST ONLY--set up catchability covariate 
    if(reg=="wc"){
      Q_ik <- matrix( ifelse(Data_Geostat$Survey=='annual', 1, 0), ncol=1 ) 
    }
    
    ########################
    ###  write out Z_gm if necessary 
    ########################
    
    # if Z_gm matrix does not exist yet, generate it and save it 
    if (!file.exists(Z_gmFile) & reg %in% c('neus','wc')){
      library(tidyverse)
      # edit Z_gm to add new axes:
      Z_gm = fit$spatial_list$loc_g
      # convert northings/eastings to lat/lon
      tmpUTM = cbind('PID'=1,'POS'=1:nrow(Z_gm),'X'=Z_gm[,'E_km'],'Y'=Z_gm[,'N_km'])
      attr(tmpUTM,"projection") = "UTM"
      attr(tmpUTM,"zone") = fit$extrapolation_list$zone - ifelse( fit$extrapolation_list$flip_around_dateline==TRUE, 30, 0 )
      latlon_g = PBSmapping::convUL(tmpUTM)                                                         #$
      latlon_g = cbind( 'Lat'=latlon_g[,"Y"], 'Lon'=latlon_g[,"X"])
      latlon_g <- as.data.frame(latlon_g)
      # use lat/lon to get coastal distance
      coast_km = NULL
      for(k in 1:nrow(latlon_g)){
        out = get_length(lon=latlon_g$Lon[k], lat=latlon_g$Lat[k], distdf = coastdistdat)/1000 # works better in a for loop than in apply where it's hard to get the rowwise nature to work--revisit sometime?
        coast_km <- c(coast_km, out)
      }
      Z_gm = cbind( Z_gm, "coast_km"=coast_km )
      Z_gm_axes = colnames(Z_gm)
      saveRDS(Z_gm, Z_gmFile)
    } else if (reg %in% c('neus','wc')) {
      Z_gm = readRDS(Z_gmFile)
      Z_gm_axes <- colnames(Z_gm)
      print(paste0(Z_gmFile, " already exists"))
    }
    
    # different approach to coordinate conversion for EBS, which uses a rotated NW axis not coastal distance
    if (!file.exists(Z_gmFile) & reg=="ebs"){
      library(tidyverse)
      Z_gm = fit$spatial_list$loc_g
      tmpUTM = cbind('PID'=1,'POS'=1:nrow(Z_gm),'X'=Z_gm[,'E_km'],'Y'=Z_gm[,'N_km'])
      attr(tmpUTM,"projection") = "UTM"
      attr(tmpUTM,"zone") = 2  # DANGER!! THIS SHOULD BE fit$extrapolation_list$zone - ifelse( fit$extrapolation_list$flip_around_dateline==TRUE, 30, 0 ) BUT IT RESULTS IN A NEGATIVE NUMBER SO I'M OVERWRITING IT HERE
      latlon_g = PBSmapping::convUL(tmpUTM)                                                         #$
      latlon_g = cbind( 'Lat'=latlon_g[,"Y"], 'Lon'=latlon_g[,"X"])
      latlon_g <- as.data.frame(latlon_g)
      
      line_km=NULL
      for(k in 1:nrow(latlon_g)){
        out = get_length(lon=latlon_g$Lon[k], lat=latlon_g$Lat[k], distdf = axisdistdat)/1000 
        line_km <- c(line_km, out)
      }
      Z_gm = cbind( Z_gm, "line_km"=line_km )
      Z_gm_axes = colnames(Z_gm)
      saveRDS(Z_gm, Z_gmFile)
    } else if (reg=="ebs") {
      Z_gm = readRDS(Z_gmFile)
      Z_gm_axes <- colnames(Z_gm)
      print(paste0(Z_gmFile, " already exists"))
    }
    
    # write out coordinate conversion df if it doesn't exist already 
    
    refdfFile <- paste0(getwd(),'/processed-data/',reg,'_coords_conversion.rds')
    if(!file.exists(refdfFile)) {
      Z_gm1 = Z_gm
      # revert to lat/lon using code from FishStatsUtils
      tmpUTM = cbind('PID'=1,'POS'=1:nrow(Z_gm1),'X'=Z_gm1[,'E_km'],'Y'=Z_gm1[,'N_km'])
      attr(tmpUTM,"projection") = "UTM"
      attr(tmpUTM,"zone") = fit$extrapolation_list$zone - ifelse( fit$extrapolation_list$flip_around_dateline==TRUE, 30, 0 )
      latlon_g = PBSmapping::convUL(tmpUTM)                                                         
      latlon_g = cbind( 'Lat'=latlon_g[,"Y"], 'Lon'=latlon_g[,"X"])
      latlon_g <- as.data.frame(latlon_g)
      refdf <- cbind(Z_gm1, latlon_g) 
      saveRDS(refdf,refdfFile)
    }
    
    ########################
    ###  re-build and run model 
    ########################
    
    # this script includes progressive checks for model fitting issues, and re-runs fit_model with different specifications as necessary.  
    # WC is separated because of the extra catchability coefficient for their survey change 
    
    # attempt 0: estimate parameters without estimating standard errors, just so we have them if needed later. this should run fine for all species. 
    if(reg %in% c('ebs','neus')){
      fit0 =try(
        fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
                  "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                  "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                  lower=-Inf, upper=Inf,
                  test_fit = FALSE,
                  fine_scale=fine_scale,
                  anisotropy=FALSE,
                  Use_REML=TRUE,
                  Z_gm = Z_gm,
                  getsd=FALSE,
                  newtonsteps=0
        ))
      
      # attempt 1: re-run model with SE estimation. if this model has no convergence issues as tested below, it will be the last model we run for this species. 
      fit = try(
        fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
                  "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                  "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                  lower=-Inf, upper=Inf,
                  test_fit = FALSE,
                  fine_scale=fine_scale,
                  anisotropy=FALSE,
                  Use_REML=TRUE,
                  Z_gm = Z_gm,
                  getsd=TRUE,
                  newtonsteps=1
        ))
      
      if(!class(fit)=="try-error"){fit$parameter_estimates$adjustments_for_convergence <- "none"}
      
      
      # check that all maximum gradients have an absolute value below 0.01 and the model didn't throw an error 
      if(class(fit)=="try-error" || (!class(fit)=="try-error" &  max(abs(fit$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
        
        # attempt 2: if fit didn't work above, try re-running it with the parameter estimates from fit0
        fit = try(
          fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
                    "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                    "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                    "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                    "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                    lower=-Inf, upper=Inf,
                    test_fit = FALSE,
                    fine_scale=fine_scale,
                    anisotropy=FALSE,
                    Use_REML=TRUE,
                    Z_gm = Z_gm,
                    getsd=TRUE,
                    newtonsteps=1,
                    parameters=fit0$ParHat
          ))
        
        if(!class(fit)=="try-error"){fit$parameter_estimates$adjustments_for_convergence <- "used fit0 parameters"}
        
      }
      
      
      # check that all maximum gradients have an absolute value below 0.01 and the model didn't throw an error 
      
      if(class(fit)=="try-error" || (!class(fit)=="try-error" &  max(abs(fit$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
        
        # attempt 3: if fit still failed / didn't converge, it could be because L_beta1_cf or L_beta2_ct is approaching zero; use a different RhoConfig 
        
        # if it's a model, but didn't converge... 
        if(!class(fit)=="try-error") {
          if("L_beta1_cf" %in% fit$parameter_estimates$diagnostics$Param & "L_beta2_ct" %in% fit$parameter_estimates$diagnostics$Param){ # AND if these parameters exist in the model... 
            if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001) {
              newRhoConfig = c("Beta1"=3, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
            }
            if(abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
              newRhoConfig = c("Beta1"=4, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
            }
            if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001 & 
               abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
              newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
            }}
        }
        
        # if it's an error...
        if(class(fit)=="try-error"){
          newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
        }
        
        fit = try(
          fit_model("settings"=settings, 
                    "RhoConfig"= newRhoConfig,
                    "Lat_i"=Data_Geostat[,'Lat'],
                    "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                    "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                    "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                    "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                    lower=-Inf, upper=Inf,
                    test_fit = FALSE,
                    fine_scale=fine_scale,
                    anisotropy=FALSE,
                    Use_REML=TRUE,
                    Z_gm = Z_gm,
                    getsd=TRUE,
                    newtonsteps=1
                    
          ))
        if(!class(fit)=="try-error"){fit$parameter_estimates$adjustments_for_convergence <- "changed RhoConfig"}
        
      }
    }
    
    if(reg=="wc"){
      
      fit0 =try(
        fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
                  "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                  "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                  lower=-Inf, upper=Inf,
                  test_fit = FALSE,
                  fine_scale=fine_scale,
                  anisotropy=FALSE,
                  Use_REML=TRUE,
                  Z_gm = Z_gm,
                  getsd=FALSE,
                  newtonsteps=0,
                  Q_ik=Q_ik # ONLY FOR WEST COAST, adding relative catchability factor the survey the data is from
        ))
      
      fit = try(
        fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
                  "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                  "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                  "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                  "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                  lower=-Inf, upper=Inf,
                  test_fit = FALSE,
                  fine_scale=fine_scale,
                  anisotropy=FALSE,
                  Use_REML=TRUE,
                  Z_gm = Z_gm,
                  getsd=TRUE,
                  newtonsteps=1,
                  Q_ik=Q_ik
        ))
      if(!class(fit)=="try-error"){fit$parameter_estimates$adjustments_for_convergence <- "none"}
      
      if(class(fit)=="try-error" || (!class(fit)=="try-error" &  max(abs(fit$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
        
        fit = try(
          fit_model("settings"=settings, "Lat_i"=Data_Geostat[,'Lat'],
                    "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                    "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                    "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                    "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                    lower=-Inf, upper=Inf,
                    test_fit = FALSE,
                    fine_scale=fine_scale,
                    anisotropy=FALSE,
                    Use_REML=TRUE,
                    Z_gm = Z_gm,
                    getsd=TRUE,
                    newtonsteps=1,
                    parameters=fit0$ParHat,
                    Q_ik=Q_ik
          ))
        if(!class(fit)=="try-error"){fit$parameter_estimates$adjustments_for_convergence <- "used fit0 parameters"}
        
      }
      
      if(class(fit)=="try-error" || (!class(fit)=="try-error" &  max(abs(fit$parameter_estimates$diagnostics$final_gradient)) > 0.01)){
        
        if(!class(fit)=="try-error"){
          if("L_beta1_cf" %in% fit$parameter_estimates$diagnostics$Param & "L_beta2_ct" %in% fit$parameter_estimates$diagnostics$Param){
            if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001) {
              newRhoConfig = c("Beta1"=3, "Beta2"=4, "Epsilon1"=4, "Epsilon2"=4)
            }
            if(abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
              newRhoConfig = c("Beta1"=4, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
            }
            if(abs(fit$parameter_estimates$par["L_beta1_cf"]) < 0.001 & 
               abs(fit$parameter_estimates$par["L_beta2_ct"]) < 0.001) {
              newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
            }}
        } 
        
        if(class(fit)=="try-error"){
          newRhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=4)
        }
        
        fit = try(
          fit_model("settings"=settings, 
                    "RhoConfig"= newRhoConfig,
                    "Lat_i"=Data_Geostat[,'Lat'],
                    "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'],
                    "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'],
                    "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                    "getReportCovariance"=FALSE, "getJointPrecision"=TRUE,
                    lower=-Inf, upper=Inf,
                    test_fit = FALSE,
                    fine_scale=fine_scale,
                    anisotropy=FALSE,
                    Use_REML=TRUE,
                    Z_gm = Z_gm,
                    getsd=TRUE,
                    newtonsteps=1,
                    Q_ik=Q_ik
                    
          ))
        if(!class(fit)=="try-error"){fit$parameter_estimates$adjustments_for_convergence <- "changed RhoConfig"}
        
      }
      
    }
    
    
    # save species-specific parameter_estimates 
    if(!class(fit)=="try-error"){
      capture.output( fit$parameter_estimates, file=file.path(paste0(RegionFile,j,"_parameter_estimates",".txt") ))
    }else{
      capture.output( fit0$parameter_estimates, file=file.path(paste0(RegionFile,j,"_parameter_estimates_run0",".txt") ))
      
    }
    
    
    # calculate range edges with two methods of estimating SE (only one used in the paper)
    if(!class(fit)=="try-error"){
      # out_absolute <- try(
      #   get_range_edge( fit.model=fit,
      #                   working_dir=paste0(getwd(),"/"), 
      #                   Year_Set=Year_Set, 
      #                   Years2Include=Years2Include, # drops years with no data (or 100% data)
      #                   n_samples=100, 
      #                   quantiles=quantiles_of_interest,
      #                   "Z_gm_axes"=Z_gm_axes,
      #                   "calculate_relative_to_average" = FALSE)
      # )
      
      source(paste0(getwd(),'/functions/get_range_edge.R'))
      source(paste0(getwd(),'/functions/get_density.R'))    
      
      out_relative <- try(
        get_range_edge( fit.model=fit,
                        working_dir=paste0(getwd(),"/"), 
                        Year_Set=Year_Set, 
                        Years2Include=Years2Include, # drops years with no data (or 100% data)
                        n_samples=100, 
                        quantiles=quantiles_of_interest,
                        "Z_gm_axes"=Z_gm_axes,
                        "calculate_relative_to_average" = TRUE)
      )
      
      # these files are huge, which is why the code to save them is commented out below. density isn't used in the paper, just useful to see and I wanted to keep the method for doing this visible 
      out_density <- try(
        get_density(fit.model=fit,
                    Years2Include=Years2Include)
      )
    }
    
    # if(!class(out_absolute)=='try-error'){
    #   out_absolute$species <- paste0(j) 
    #   write.csv(out_absolute, file.path(paste0(RegionFile,j,"_absolute_SE_edges.csv")))}
    
    if(!class(out_relative)=='try-error'){
      out_relative$species <- paste0(j) 
      write.csv(out_relative, file.path(paste0(RegionFile,j,"_relative_SE_edges.csv")))}
    
    # if(!class(out_density)=='try-error'){
    #   out_density$species <- paste0(j) 
    #   write.csv(out_density, file.path(paste0(RegionFile,j,"_density.csv")))}        
    # not saving this right now because each df for each species is ~1M rows, but keeping the code so anyone running line-by-line can easily access the density df
    
    
  }
  
  ########################
  ###  save outputs
  ########################
  
  # not making a collated df for density because they're individually huge 
  
  rel_edge_files <- list.files(path=RegionFile, pattern="relative_SE_edges.csv", full.names=TRUE)
  rel_edge_df <- dplyr::bind_rows(lapply(rel_edge_files, read.csv))
  rel_edge_df$X <- NULL
  saveRDS(rel_edge_df, paste0(getwd(),"/processed-data/",reg,'_relative_SE_vast_edge_df.rds'))
  
  # abs_edge_files <- list.files(path=RegionFile, pattern="absolute_SE_edges.csv", full.names=TRUE)
  # abs_edge_df <- dplyr::bind_rows(lapply(abs_edge_files, read.csv))
  # abs_edge_df$X <- NULL
  # saveRDS(abs_edge_df, paste0(getwd(),"/processed-data/",reg,'_absolute_SE_vast_edge_df.rds'))
  
  capture.output( settings, file=file.path(RegionFile,'settings.txt'))
  
  spp_not_converged <- setdiff(Species_list, unique(rel_edge_df$species)) 
  saveRDS(spp_not_converged, paste0(RegionFile, "spp_not_converged.rds"))
  
  
}# end of parallel


