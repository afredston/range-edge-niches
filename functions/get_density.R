get_density <- function(fit.model, Years2Include){
  
  # Get object
  obj<- fit.model$tmb_list[["Obj"]]
  
  # Extract
  report<- obj$report()
  
  # Slots worth describing -- density at each knot? Can project this after to get east/northing, or could pull east/northing out of the fit.model$extrapolation_list$Data_Extrap object...
  dens.df<- data.frame("lon" = fit.model$spatial_list$latlon_g[,"Lon"],
                       "lat" = fit.model$spatial_list$latlon_g[,"Lat"],
                       "density" = report$D_gcy)
  
  dens.df <- reshape2::melt(dens.df, id.vars=c('lon','lat'))
  colnames(dens.df) <- c('lon','lat','year',
                         'density')
  
  dens.df$year <- as.numeric(gsub("density.","", dens.df$year))
  dens.df$year <- fit$year_labels[dens.df$year] # convert to real years 
  
  dens.df <- dens.df[dens.df$year %in% Years2Include,] # keep only years with real data 
  
  return(dens.df)
  }