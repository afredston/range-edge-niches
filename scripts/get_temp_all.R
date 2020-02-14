# THIS SCRIPT TAKES HOURS TO RUN, so don't re-run the whole thing if you don't have to! the OISST datasets in particular are slow.

# get all temperature datasets from ERDDAP for all regions:
# https://coastwatch.pfeg.noaa.gov/erddap/index.html
library(here)
source(here("functions","get_rerddap.R"))

# select bounding boxes for all regions 
neus_latrange <- c(35, 45)
neus_lonrange <- c(-78, -66) 
wc_latrange <- c(30, 50)
wc_lonrange <- c(-126, -117)
ebs_latrange <- c(54, 66)
ebs_lonrange <- c(-179.5, -154)

# program call for each dataset; get from the ERDDAP pages and update end years if more data is added
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdHadISST.html
# https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst2Agg_LonPM180.html

hadisst <- "erdHadISST"
hadisst_years <- c(1967, 2018)
hadisst_fields <- "sst"
oisst <- "ncdcOisst2Agg_LonPM180"
oisst_years <- c(1982, 2018)
oisst_fields <- "sst"

neus_hadisst <- get_rerddap(datasetID = hadisst, latrange = neus_latrange, lonrange = neus_lonrange, startyear=hadisst_years[1], endyear=hadisst_years[2], fields=hadisst_fields, verbose=TRUE)
neus_oisst <- get_rerddap(datasetID = oisst, latrange = neus_latrange, lonrange = neus_lonrange, startyear=oisst_years[1], endyear=oisst_years[2], fields=oisst_fields, verbose=TRUE)

wc_hadisst <- get_rerddap(datasetID = hadisst, latrange = wc_latrange, lonrange = wc_lonrange, startyear=hadisst_years[1], endyear=hadisst_years[2], fields=hadisst_fields, verbose=TRUE)
wc_oisst <- get_rerddap(datasetID = oisst, latrange = wc_latrange, lonrange = wc_lonrange, startyear=oisst_years[1], endyear=oisst_years[2], fields=oisst_fields, verbose=TRUE)

ebs_hadisst <- get_rerddap(datasetID = hadisst, latrange = ebs_latrange, lonrange = ebs_lonrange, startyear=hadisst_years[1], endyear=hadisst_years[2], fields=hadisst_fields, verbose=TRUE)
ebs_oisst <- get_rerddap(datasetID = oisst, latrange = ebs_latrange, lonrange = ebs_lonrange, startyear=oisst_years[1], endyear=oisst_years[2], fields=oisst_fields, verbose=TRUE)

# if there is an error that says "Error in R_nc4_open: NetCDF: Unknown file format"
# find the directory C:\Users\alexafh\AppData\Local\Cache/R/rerddap/ and delete its contents 

dfs <- list(neus_oisst, neus_hadisst, wc_hadisst, wc_oisst,ebs_hadisst, ebs_oisst)
names <- c("neus_oisst", "neus_hadisst", "wc_hadisst", "wc_oisst","ebs_hadisst", "ebs_oisst")
names(dfs) <- names

for(i in 1:length(dfs)){
  saveRDS(dfs[i], here("processed-data",paste0(names[i],"_raw.rds")))
}
rm(list=ls())
