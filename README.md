# Realized thermal niche tracking at range limits of North American marine species

### A. Fredston-Hermann, M. Pinsky, B. Selden, C. Szuwalski, J. T. Thorson, S. D. Gaines, B. S. Halpern 

To reference the data or methods here, please cite the manuscript. Contact A. Fredston-Hermann with questions at fredston@rutgers.edu. 

## Overview 

This repository contains code to:

* Estimate the annual position of range limits from NOAA trawl survey data using [VAST](https://github.com/James-Thorson-NOAA/VAST) 
* Use hindcast SST datasets to reconstruct the extreme temperatures found at species' range limits each year 
* Conduct statistical analyses of range limit shifts and edge thermal niche shifts over time 

To fully reproduce the analysis, users will need to install VAST. 

The repository is organized as follows:

* `raw-data` contains raw datasets obtained from other sources that have not been changed in any way. Where possible, data are downloaded within scripts, so future users do not need to have data files. Where that was not possible, scripts using the data files contain instructions to download them. 
* `scripts` contains all code to analyze or transform data. 
* `processed-data` contains outputs of scripts that filter, clean, summarize, or analyze data, e.g., a dataframe with the higher taxonomy of species used in the final analysis. 
* `functions` contains homemade functions called in `scripts`.
* `figures` contains **code** to produce figures in the manuscript (the outputs are in `results`).
* `results` contains figures, tables, and other outputs that are used in the manuscript. 

There are some additional directories for model outputs that are not version controlled. Large files (e.g., raw temperature datasets), images, and PDFs are not version controlled. 

## Instructions

Scripts should be run in the following order:

1. `get_neus_data.R` imports Northeast data (the only dataset requiring manual download) and reformats it to match those downloaded with [FishData](https://github.com/James-Thorson/FishData)
1. `get_neus_wc_coastlines.R` creates a coastal distance axis for use with VAST for those two regions
1. `get_range_limits.R` uses VAST to calculate range edges for all three regions. It runs in parallel (don't forget to update the number of cores based on your machine!) but even so may take several days for all species and regions. Note that the output data frames are also in the repository
1. `get_taxonomy.R` fetches higher taxonomy of study species from [WORMS](http://marinespecies.org/aphia.php?p=search) using [taxize](https://github.com/ropensci/taxize/)
1. `validate_range_limits.R` filters the VAST output for only range limits that actually fall in the study region, based on passing certain filters. This needs to be re-run every time VAST is re-run in `get_range_limits.R` 
1. `prep_sst.R` fetches historical SST data from the [NOAA ERDDAP server](https://coastwatch.pfeg.noaa.gov/erddap/index.html) for each region, crops rasters to the extent of a bathymetric mask for each region that is also created in this script, performs a mean bias correction to combine different SST datasets, and writes the SST data out as dataframes
1. `match_sst_to_axis.R` matches SST values to points along the axis of range limit measurement for each study region, and combines SST datasets where necessary. This doesn't need to be re-run if new edges are generated, because it takes as input the VAST coordinates, not the range edge positions.  
1. `analyze_range_limits.R` executes all models to analyze range edges. 
