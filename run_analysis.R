#' Run analysis for the paper "Evidence-based guidelines for developing automated
#'  conservation assessment methods".
#' 
#' This will:
#' 1. download rasters (if needed) and process them to the desired resolution
#' 2. collate lists of accepted species and their IUCN Red list assessments for
#'    the 3 study groups (Myrica, legumes, orchids)
#' 3. process occurrence data from GBIF and reconcile names to WCVP taxonomy
#' 4. extract environmental and threat-related values at each occurrence record
#' 5. generate sets of occurrences at different levels of cleaning
#' 6. summarise per-species predictor values for each of these sets
#' 7. run automated assessment methods on these predictor sets and evaluate their performance
#' 8. summarise the results for easier comparison and analysis
#' 9. make the figures for the paper
#' 10. make figures and tables for the supplementary materials
#' 
#' All scripts should load data from the `data` folder and save outputs to the `output` folder.
#' 

# libraries ----
library(here)    # handles file paths
library(glue)    # string interpolation
library(stringr) # manipulate strings
library(vroom)   # fast text file reading/writing

# 1. download and process rasters ----
overwrite <- TRUE
download_dir <- "data/raster-downloads"
output_dir <- "data/rasters"
source(here("analysis/01_compile_rasters.R"))

# 2. make species lists and match to assessments ----
redlist_dir <- "data/redlist_plants_2021_3"
distributions <- TRUE
output_dir <- "output/species-lists"
redlist_dir <- "data/redlist_plants_2021_3"
wcvp_dir <- "data/wcvp"

name <- "orchids"
family <- "Orchidaceae"
source(here("analysis/02_collate_species.R"))

name <- "legumes"
family <- "Fabaceae"
srli_file <- "data/srli_assessed_legumes.xlsx"
source(here("analysis/02_collate_species.R"))

srli_file <- NULL
family <- NULL
name <- "myrcia"
genus <- "Myrcia,Gomidesia,Marlierea,Mitranthes,Calyptranthes"
section_file <- "data/myrcia_sections.csv"
working_set <- "data/myrcias_from_working_set.csv"
source(here("analysis/02_collate_species.R"))

# 3. match gbif occurrences to accepted names ----
output_dir <- "output/occurrence-records"
wcvp_dir <- "data/wcvp"

file_list <- list(
  list(
    occurrences="data/myrcia-db_raw-occurrences_2022-02-09.csv", 
    species="output/species-lists/myrcia_species-list.csv"
  ),
  list(
    occurrences="data/myrcia-db_raw-occurrences_2022-02-09.csv", 
    species="output/species-lists/myrcia_species-list.csv"
  ),
  list(
    occurrences="data/myrcia-db_raw-occurrences_2022-02-09.csv", 
    species="output/species-lists/myrcia_species-list.csv"
  ),
  list(
    occurrences="data/myrcia-db_raw-occurrences_2022-02-09.csv", 
    species="output/species-lists/myrcia_species-list.csv"
  )
)

for (item in file_list) {
  occurrence_file <- item[1]
  species_list <- item[2]
  source(here("analysis/03_process_occurrences.R")) 
}

# 4. annotate occurrences with raster values ----
raster_dir <- "data/rasters"
output_dir <- "output/labelled-occurrences"

files <- list.files("output/occurrence-records/", pattern="_occurrences.csv",
                    full.names=TRUE)
for (occurrence_file in files) {
  source(here("analysis/04_annotate_points.R"))
}

# 5. generate cleaned occurrences sets ----
output_dir <- "output/cleaned_occurrences"
distribution_dir <- "output/species-lists"
shape_file <- "data/wgsrpd/level3/level3.shp"

files <- list.files("output/labelled-occurrences", pattern="_labelled-occurrences.csv",
                    full.names=TRUE)
for (occurrence_file in files) {
  # don't want to run automated cleaning on the monographic db
  already_cleaned <- str_detect(occurrence_file, "myrcia-db")
  source(here("analysis/05_clean_occurrences.R"))
}

# 6. prepare species-level predictors ----
output_dir <- "output/predictors"
occurrence_dir <- "output/labelled-occurrences"
assessment_dir <- "output/species-lists"

files <- list.files("output/cleaned-occurrences")
for (records_file in files) {
  source(here("analysis/06_prepare_predictors.R"))
}

# create predictor sets for all three groups together
for (cmbn in step_cmbns) {
  in_pattern <- glue("[a-z]+-gbif_filter-{cmbn$filter}_clean-{cmbn$clean}.csv")
  out_name <- glue("all_filter-{cmbn$filter}_clean-{cmbn$clean}.csv")
  predictor_files <- list.files(here("output/predictors/"), pattern=in_pattern,
                                 full.names=TRUE)
  predictors <- vroom(predictor_files)
  vroom_write(predictors, here("output/predictors/", out_name))
}

# 7. Evaluate AA methods ----
output_dir <- "output/method-results"
method_dir <- "methods"
model_dir <- "output/trained-models"
downsample <- FALSE
random_seed <- 1989
srli <- FALSE

files <- list.files("output/predictors", pattern=".csv")
methods <- c("eoo-threshold", "conr", "decision-stump", "decision-tree",
             "random-forest", "iucnn")

for (method in methods) {
  for (predictor_file in files) {
    source(here("analysis/07_evaluate_methods.R")) 
  }
}

# run each method on the SRLI legume dataset
legume_files <- list.files("output/predictors", pattern="legumes.*.csv")
srli <- TRUE
for (method in methods) {
  for (predictor_file in legume_files) {
    source(here("analysis/07_evaluate_methods.R")) 
  }
}

# downsample just the ML methods
methods <- c("decision-stump", "decision-tree", "random-foret", "iucnn")
srli <- FALSE
downsample <- TRUE
for (method in methods) {
  for (predictor_file in files) {
    source(here("analysis/07_evaluate_methods.R")) 
  }
}

srli <- TRUE
for (method in methods) {
  for (predictor_file in legume_files) {
    source(here("analysis/07_evaluate_methods.R")) 
  }
}

# use grouped cv on the ML methods
srli <- FALSE
downsample <- FALSE
for (method in methods) {
  for (predictor_file in files) {
    if (str_detect(predictor_file, "myrcia")) {
      group <- "section"
    } else if (str_detect(predictor_file, "all_")) {
      group <- "group"
    } else {
      group <- "genus"
    }
    source(here("analysis/07_evaluate_methods.R")) 
  }
}

downsample <- TRUE
for (method in methods) {
  for (predictor_file in files) {
    if (str_detect(predictor_file, "myrcia")) {
      group <- "section"
    } else if (str_detect(predictor_file, "all_")) {
      group <- "group"
    } else {
      group <- "genus"
    }
    source(here("analysis/07_evaluate_methods.R")) 
  }
}

# 8. summarise results for the manuscript ----
species_dir <- "output/species-lists"
occurrence_dir <- "output/occurrence-records"
predictor_dir <- "output/predictors"
results_dir <- "output/method-results"
output_dir <- "output/summarised-results"
model_dir <- "output/trained-models"

source(here("analysis/08_summarise_results.R"))

# 9. make figures for the manuscript ----
results_dir <- "output/summarised-results"
figure_dir <- "figures"

source(here("analysis/09_make_figures.R"))

# 10. make figures and tables for supplementary ----
results_dir <- "output/summarised-results"
figure_dir <- "figures"

source(here("analysis/10_supplementary_materials.R"))
