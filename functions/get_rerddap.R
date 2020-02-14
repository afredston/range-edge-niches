# this function uses the NOAA ERDDAP server to access climate data
# while it will fetch the data, users still need to select: 
# bounding box: latrange=c(minlat, maxlat), lonrange=c(minlon, maxlon)
# the dataset from https://coastwatch.pfeg.noaa.gov/erddap/index.html
# from the dataset, users should specify the datasetID; the fields requested (in addition to the "dimension variables", i.e., x, y, and z if applicable, which come by default), based on the "grid variables" names (can put multiple e.g. c('salt','temp'), ); and the years, not to exceed the date range in ERDDAP 

library(rerddap) 
library(progress)

get_rerddap <- function(datasetID, latrange, lonrange, startyear, endyear, fields, verbose) {
  rerddap.df <- NULL
  
  # set up progress bar 
  if(verbose==TRUE) {
  pb <- progress_bar$new(
    format = "  downloading [:bar] :percent in :elapsed",
    total = length(startyear:endyear), clear = FALSE, width= 60)
  
  for(i in startyear:endyear) {
    pb$tick()
    daterange = c(paste0(i,'-01-01'),paste0(i,'-12-31'))
    out = griddap(datasetID, latitude = latrange, longitude = lonrange, time = daterange, fields = fields)$data
    rerddap.df <- rbind(rerddap.df, out)
  }
  return(rerddap.df)
  }
  else {
    for(i in startyear:endyear) {
      daterange = c(paste0(i,'-01-01'),paste0(i,'-12-31'))
      out = griddap(datasetID, latitude = latrange, longitude = lonrange, time = daterange, fields = fields)$data
      rerddap.df <- rbind(rerddap.df, out)
    }
    return(rerddap.df)
  }
}
