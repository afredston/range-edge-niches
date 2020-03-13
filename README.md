# Realized thermal niche tracking at range limits of North American marine species

### A. Fredston-Hermann, M. Pinsky, B. Selden, C. Szuwalski, J. T. Thorson, S. D. Gaines, B. S. Halpern 

To reference the data or methods here, please cite the manuscript. Contact A. Fredston-Hermann with questions. 

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


Large files (e.g., raw temperature datasets), images, and PDFs are not version controlled. 

## Instructions

Scripts should be run in the following order:

1. `get_neus_data.R` imports Northeast data (the only dataset requiring manual download) and reformats it to match those downloaded with [FishData](https://github.com/James-Thorson/FishData)
1. `get_neus_wc_coastlines.R` creates a coastal distance axis for use with VAST for those two regions 
1. `get_range_limits.R` uses VAST to calculate range edges for all three regions. It runs in parallel but even so may take several days for all species and regions. Note that the output data frames are also in the repository
1. `get_taxonomy.R` fetches higher taxonomy of study species from [WORMS](http://marinespecies.org/aphia.php?p=search) using [taxize](https://github.com/ropensci/taxize/)
1. `validate_range_limits.R` filters the VAST output for only range limits that actually fall in the study region, based on passing certain filters. This needs to be re-run every time VAST is re-run in `get_range_limits.R` 
1. `get_temp_all.R` fetches historical SST data from the [NOAA ERDDAP server](https://coastwatch.pfeg.noaa.gov/erddap/index.html) for each region. `get_temp_all.R` and `crop_temp_all.R` only need to be run once
1. `crop_temp_all.R` takes the SST datasets and crops them to within the US EEZ and a depth cutoff
1. `match_sst_to_axis.R` matches SST values to points along the axis of range limit measurement for each study region, and combines SST datasets where necessary 
1. `analyze_range_limits.R` 
