# Realized thermal niche tracking at range limits of North American marine species

### A. Fredston, M. Pinsky, B. Selden, C. Szuwalski, J. T. Thorson, S. D. Gaines, B. S. Halpern 

To reference the data or methods here, please cite the manuscript. Contact A. Fredston with questions at fredston@rutgers.edu. 

## Overview 

This repository contains code to:

* Estimate the annual position of range edges from NOAA trawl survey data using [VAST](https://github.com/James-Thorson-NOAA/VAST) 
* Use hindcast SST datasets to reconstruct the extreme temperatures found at species' range edges each year 
* Conduct statistical analyses of range edge shifts and edge thermal niche shifts over time 

To fully reproduce the analysis, users will need to install VAST. However, this repository contains outputs from VAST and subsequent statistical analyses for users who just want to examine the results or reproduce figures. 

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
1. `get_axes_of_measurement.R` creates a coastal distance axis for use with VAST for those two regions
1. `get_range_edges.R` uses VAST to calculate range edges for all three regions. It runs in parallel (don't forget to update the number of cores based on your machine!) but even so may take several days for all species and regions. Note that the output data frames are also in the repository
1. `validate_range_edges.R` filters the VAST output for only range edges that actually fall in the study region, based on passing certain filters. This needs to be re-run every time VAST is re-run in `get_range_edges.R` 
1. `get_taxonomy.R` fetches higher taxonomy of study species from [WORMS](http://marinespecies.org/aphia.php?p=search) using [taxize](https://github.com/ropensci/taxize/)
1. `prep_sst.R` fetches historical SST data from the [NOAA ERDDAP server](https://coastwatch.pfeg.noaa.gov/erddap/index.html) for each region, crops rasters to the extent of a bathymetric mask for each region that is also created in this script, performs a mean bias correction to combine different SST datasets, and writes the SST data out as dataframes. This does not need to be repeated unless the source data is updated
1. `match_sst_to_axis.R` matches SST values to points along the axis of range limit measurement for each study region, and combines SST datasets where necessary. This doesn't need to be re-run if new edges are generated, because it takes as input the VAST coordinates, not the range edge positions.
1. `calculate_edge_thermal_niches.R` fits Bayesian models to range edges and edge thermal niches. This script was developed on a remote server and may overwhelm personal computers' processing capacity.  
1. `analyze_range_edges.R` conducts the main analyses in the paper. 
1. `paper_stats.R` calculates miscellaneous statistical results reported in the manuscript. 

Code to generate figures in the manuscript and supplementary materials can typically be found in the `figure-scripts` folder. A few exceptions: 

* Figure 1 requires the full output from the Bayesian models, which is quite large; those figures are generated at the end of `calculate_edge_thermal_niches.R`. 
* The figures in Appendix 2 showing how the axes of measurement were developed are generated in `get_axes_of_measurement.R`

## Computational requirements 

Analyses were conducted in R 3.6.0 on a machine with the following specifications: 
Windows Server 2019
2-Intel 6154 Xeon    3.0Ghz 18cores(36 threads)each
512 GB memory
6TB- SSD storage
18TB- HDD storage
NVIDIA P4000 8Gb
10Gb ethernet

On this machine, the VAST models (`get_range_edges.R`) for all regions were run in parallel on 18 cores overnight. Bayesian models in `calculate_edge_thermal_niches.R` are not currently scripted to run in parallel; on this machine, each set of models takes several hours to run. The Bayesian models, especially the second set fitting edge thermal niche change over time, are highly memory-intensive and may crash systems with low memory capacity. 

## Use, problems, and feedback

If you encounter any errors or problems, please create an Issue here. Likewise, please consider starting an Issue for any questions so that others can view conversations about the analysis and code. Again, don't hesitate to contact the lead author, A. Fredston, at fredston@rutgers.edu. 
