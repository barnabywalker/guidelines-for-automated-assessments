# Evidence-based guidelines for developing automated conservation assessment methods

This is the R code used to perform the analysis described in the paper *Evidence-based guidelines for developing automated conservation assessment methods*. This includes functions and scripts used to gather data and perform the analysis, summarise the results, plot figures for the paper, and write the paper. Some manual reformatting was done on the paper after compiling it, so running the manuscript files will not reproduce the paper exactly. 

## Running the scripts

You can either re-run the whole analysis from end to end, or download the archived data outputs and just run the bit that you are interested in.

In both cases, you'll need to first download or clone this repository. After that, you'll need to download either the original data sources or the archived data outputs.

### Downloading the data sources

To re-run the whole analysis, the original datasources need to be downloaded and saved into a folder called `data`.

The data sources are:
* [A human footprint index raster from SEDAC CIESIN](https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint/data-download) - save in a folder named `data/rasters`.
* [A human population density raster from SEDAC CIESIN](https://sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse) - save in a folder named `data/rasters`.
* [Elevation SRTM raster tiles](https://drive.google.com/drive/folders/0B_J08t5spvd8RWRmYmtFa2puZEE) - save and extract the file `SRTMv4.1.zip` in a folder named `data/rasters/elevation`.
* [Myrtaceae (family containing *Myrcia*) occurrences from GBIF](https://doi.org/10.15468/dl.fyf5g2) - download and extract in the `data` folder and edit `myrcia_gbif_path` in `analysis/03_process_occurrences.R`.
* [Fabaceae (legumes) occurrences from GBIF](https://doi.org/10.15468/dl.nm4p3y) - download and extract in the `data` folder and edit `legume_gbif_path` in `analysis/03_process_occurrences.R`.
* [Orchidaceae (orchids) occurrences from GBIF](https://doi.org/10.15468/dl.wsvw3m) - download and extract in the `data` folder and edit `orchid_gbif_path` in `analysis/03_process_occurrences.R`.
* [All IUCN Red List assessments for vascular plants](https://www.iucnredlist.org/search?dl=true&permalink=bec1e3e1-6aea-4f4e-9fbd-2a34c6d0270f) - download and extract in the `data` folder and edit the `rl_assessments` and `rl_taxonomy` variables in `analysis/02_collate_species.R`. This is a saved search accessed at IUCN Red List version 2021-1 and will be different once the Red List is updated, but the assessments can be filtered based on the version they were added under to recreate the list used here.

Some data sources have not been published yet, and so have not been included in the data archive. These are:
* Occurrence records from the database prepared for a monographic treatment of *Myrcia*, which has not been published yet.
* Unpublished assessments for species in the genus *Myrcia*, which may change category during review.
* The list of legume species assessed as part of the Sampled Red List Index.

While the scripts processing these data sources will not run, outputs from these scripts are preserved in the data archive.

### Downloading the archived data

To run scripts that use the archived output data, you will need to download the outputs from the archive [here](INSERT ARCHIVE LINK) and save them into a folder called `output`.

Files containing occurrence records from the monographic database of *Myrcia* species have been removed, but predictors calculated from these files are preserved. As such, steps 7 to 11 of the analysis should re-run fine from the archived data.

## Software

All analysis was carried out in R version 4.0.5.

The packages and versions used were:
- [here](https://here.r-lib.org/) - 0.1
- [tidyverse](https://www.tidyverse.org/) - 1.3.0
- [tidymodels](https://www.tidymodels.org/) - 0.1.1
- [vroom](https://vroom.r-lib.org/) - 1.4.0
- [sf](https://r-spatial.github.io/sf/) - 0.9.6
- [raster](https://cran.r-project.org/web/packages/raster/raster.pdf) - 3.4.5
- [exactextractr](https://isciences.gitlab.io/exactextractr/) - 0.5.1
- [CoordinateCleaner](https://docs.ropensci.org/CoordinateCleaner/) - 2.0.18
- [rCAT](https://cran.r-project.org/web/packages/rCAT/rCAT.pdf) - 0.1.6
- [shapper](https://modeloriented.github.io/shapper/) - 0.1.3
- [kewr](https://barnabywalker.github.io/kewr/) - 0.4.0
- [patchwork](https://patchwork.data-imaginist.com/) - 1.1.0
- [ggridges](https://wilkelab.org/ggridges/) - 0.5.2
- [ggforce](https://ggforce.data-imaginist.com/) - 0.3.2
- [furrr](https://furrr.futureverse.org/) - 0.2.0
- [randomForest](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf) - 4.6.14
- [progress](https://github.com/r-lib/progress) - 1.2.2