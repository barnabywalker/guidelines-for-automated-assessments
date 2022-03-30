[![DOI](https://zenodo.org/badge/347164961.svg)](https://zenodo.org/badge/latestdoi/347164961)

# Evidence-based guidelines for developing automated conservation assessment methods

This is the R code used to perform the analysis described in the paper *Evidence-based guidelines for developing automated conservation assessment methods*. This includes functions and scripts used to gather data and perform the analysis, summarise the results, plot figures for the paper, and write the paper. Some manual reformatting was done on the paper after compiling it, so running the manuscript files will not reproduce the paper exactly. 

The outputs of these scripts are archived on Zenodo, at https://doi.org/10.5281/zenodo.4899925.

These scripts use the [`tidymodels`](https://www.tidymodels.org/) framework for building and evaluating the models. We wrote a package to more easily incorporate some of the automated assessment methods into this framework, which can be found at [`tidyassessments`](https://www.github.com/barnabywalker/tidyassessments).

## Method specification

All automated assessment methods are specified in their own script in the `methods` folder. This includes the data preprocessing steps, model setup, the cross-validation or bootstrapping scheme, any hyperparameter values to tune over, and some variables to tell the evaluation script what quantities to calculate.

You can add more methods by adapting these scripts to the model or method you want to use, as long as it can be used with tidymodels.

## Running the scripts

All scripts have been written so that you can `source` them or run them via the command line using `RScript analysis/SCRIPT_NAME.R`. Check each script for the arguments you need to feed into them.

You can either re-run the whole analysis from end to end (using `run_analysis.R`), or download the archived data outputs and just run the bit that you are interested in.

In both cases, you'll need to first download or clone this repository. After that, you'll need to download either the original data sources or the archived data outputs.

### Downloading the data sources

To re-run the whole analysis, the original datasources need to be downloaded and saved into a folder called `data`.

The data sources are:
* [Myrtaceae (family containing *Myrcia*) occurrences from GBIF](https://doi.org/10.15468/dl.ehgw5u) - download and extract in the `data` folder and use as the `occurrence_file` input to `analysis/03_process_occurrences.R`.
* [Fabaceae (legumes) occurrences from GBIF](https://doi.org/10.15468/dl.ryr37h) - download and extract in the `data` folder and use as the `occurrence_file` input to `analysis/03_process_occurrences.R`.
* [Orchidaceae (orchids) occurrences from GBIF](https://doi.org/10.15468/dl.u5mkrs) - download and extract in the `data` folder and use as the `occurrence_file` input to `analysis/03_process_occurrences.R`.
* [All IUCN Red List assessments for vascular plants](https://www.iucnredlist.org/search?dl=true&permalink=bec1e3e1-6aea-4f4e-9fbd-2a34c6d0270f) - download and extract in the `data` folder and using the resulting folder as the `redlist_dir` input to `analysis/02_collate_species.R`. This is a saved search accessed at IUCN Red List version 2021-3 and will be different once the Red List is updated, but the assessments can be filtered based on the version they were added under to recreate the list used here.
* [Shape file of WGSRPD botanical countries]() - download and extract into the `data` folder and use the resulting level 3 shape file as the `shape_file` input to `analysis/05_clean_occurrences.R`.

The `analysis/01_compile_rasters.R` script can be used to download all raster datasources.

Some data sources have not been published yet, and so have not been included in the data archive. These are:
* Occurrence records from the database prepared for a monographic treatment of *Myrcia*, which has not been published yet.
* Unpublished assessments for species in the genus *Myrcia*, which may change category during review, but will eventually be published on the IUCN Red List.
* The list of legume species assessed as part of the Sampled Red List Index. The species included in the Sampled Red List Index is updated frequently and the latest version can be requested from the IUCN.

While the scripts processing these data sources will not run, outputs from these scripts are archived on Zenodo, at https://doi.org/10.5281/zenodo.4899925.

### Downloading the archived data

To run scripts that use the archived output data, you will need to download the outputs from the archive [here](INSERT ARCHIVE LINK) and save them into a folder called `output`.

Files containing occurrence records from the monographic database of *Myrcia* species have been removed, but predictors calculated from these files are preserved. As such, steps 7 to 10 of the analysis should re-run fine from the archived data.

## Software

All analysis was carried out in R version 4.1.2.

The packages and versions used were:
- [here](https://here.r-lib.org/) - 1.0.1
- [cli](https://cli.r-lib.org/) - 3.1.1
- [tidyverse](https://www.tidyverse.org/) - 1.3.1
- [tidymodels](https://www.tidymodels.org/) - 0.1.4
- [vroom](https://vroom.r-lib.org/) - 1.5.7
- [sf](https://r-spatial.github.io/sf/) - 1.0.6
- [raster](https://cran.r-project.org/web/packages/raster/raster.pdf) - 3.5.15
- [exactextractr](https://isciences.gitlab.io/exactextractr/) - 0.7.2
- [CoordinateCleaner](https://docs.ropensci.org/CoordinateCleaner/) - 2.0.20
- [rCAT](https://cran.r-project.org/web/packages/rCAT/rCAT.pdf) - 0.1.6
- [ConR]() - 1.3.0
- [fastshap](https://bgreenwell.github.io/fastshap/) - 0.0.7
- [vip](https://koalaverse.github.io/vip/) - 0.3.2
- [kewr](https://barnabywalker.github.io/kewr/) - 0.6.0
- [patchwork](https://patchwork.data-imaginist.com/) - 1.1.1
- [ggridges](https://wilkelab.org/ggridges/) - 0.5.3
- [ggforce](https://ggforce.data-imaginist.com/) - 0.3.3
- [ggdist](https://mjskay.github.io/ggdist/) - 3.1.1
- [furrr](https://furrr.futureverse.org/) - 0.2.3
- [randomForest](https://cran.r-project.org/web/packages/randomForest/randomForest.pdf) - 4.6.14
- [keras](https://keras.rstudio.com/) - 2.8.0
- [reticulate](https://rstudio.github.io/reticulate/) - 1.24

The `analysis/01_compile_rasters.R` script also uses the package [velox](https://github.com/hunzikp/velox) (version 0.2.1), which has been removed from CRAN because it is not actively maintained. We use it for fast aggregation of raster tiles, specifically to make the forest loss raster. You can still install the package from the GitHub page using:

```
remotes::install_github("hunzikp/velox")
```
