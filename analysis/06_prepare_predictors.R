#' Prepare species-level predictors from labelled occurrence records.
#' 
#' This script takes a list of occurrence records to include (the output of script `05_clean_occurrences.R`) 
#' as input and summarises the predictor values in the corresponding labelled points file 
#' (the output of script `04_labels_points.R`) for each species.
#' 
#' The values that are calculated are:
#'  - Extent of Occurrence (EOO) (km^2)
#'  - latitude of the centroid of the occurrence records (decimal degrees)
#'  - mean human footprint index
#'  - minimum human population density (persons/km^2)
#'  - mean forest loss (proportion of pixels showing forest loss)
#'  - mean annual average temperature (degrees C)
#'  - mean precipitation in the driest quarter (mm)
#'  - maximum elevation (m)
#'  - number of occurrence records
#'  - EOO, as calculated by ConR
#'  - Area of occupancy (km^2), as calculated by ConR
#'  - Number of locations, as calculated by ConR
#' 
#' The script then joins these predictors to the prediction target - IUCN RL assessment status.
#' Species without an assessment will be retained for prediction and the category will be NA.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save predictors to
#'  - `occurrence_dir`: path to a directory containing files of labelled occurrences,
#'        possibly the outputs from `analysis/04_label_points.R`
#'  - `assessment_dir`: path to a directory containing files with assessments for accepted species,
#'        possibly the outputs from `analysis/02_collate_species.R`
#'  - `records_file`: path to file containing list of records to include,
#'        possibly an output from `analysis/05_clean_coordinates.R`
#' 
#' EXAMPLE CLI:
#'  Rscript analysis/06_prepare_predictors.R --records_file=output/cleaned-occurrences/myrcia-gbif_filter-1_clean-A.csv
#' 
#' EXAMPLE SOURCE:
#'  output_dir="output/predictors"
#'  occurrence_dir <- "output/labelled-occurrences"
#'  assessment_dir <- "output/species-lists"
#'  records_file <- "output/cleaned-occurrences/myrcia-gbif_filter-1_clean-A.csv"
#'  source("analysis/06_prepare_predictors.R")
#'  

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(here)         # handle file paths
shhlibrary(cli)          # nice command line interface
shhlibrary(multidplyr)   # parallelise data manipulation
shhlibrary(dplyr)        # manipulate data
shhlibrary(sf)           # handle spatial data
shhlibrary(vroom)        # fast reading/writing for text files
shhlibrary(stringr)      # handle strings
shhlibrary(glue)         # format strings

source(here("R/predictor_functions.R"))

# CLI ----
cli_h1("Preparing species-level predictors")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output/predictors",
    occurrence_dir="output/labelled-occurrences",
    assessment_dir="output/species-lists"
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  output_dir <- args$output_dir
  occurrence_dir <- args$occurrence_dir
  records_file <- args$records_file
  assessment_dir <- args$assessment_dir
}

if (! exists("records_file", mode="character")) {
  cli_abort(c(
    "no path to records provided",
    "x"="You must provide a path to a file of occurrences to include as {.var records_file}."
  ))
}

if (! exists("occurrence_dir", mode="character")) {
  cli_abort(c(
    "no path to occurrences provided",
    "x"="You must provide a path to a directory containing labelled occurrences as {.var occurrence_dir}."
  ))
}

if (! exists("assessment_dir", mode="character")) {
  cli_abort(c(
    "no path to assessments provided",
    "x"="You must provide a path to a directory containing assessments as {.var assessment_dir}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save labelled records provided",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

dir.create(output_dir, showWarnings=FALSE)
name <- str_extract(basename(records_file), "^[\\w-]+?(?=_)")

occurrence_file <- file.path(occurrence_dir, glue("{name}_labelled-occurrences.csv"))

taxa_name <- str_extract(name, "^\\w+")
assessment_file <- file.path(assessment_dir, glue("{taxa_name}_species-list.csv"))

cli_alert_info("Preparing predictors from from {.file {occurrence_file}}")
cli_alert_info("Using occurrences records specified in {.file {records_file}}")
cli_alert_info("Using assessments from {.file {assessment_file}}")
cli_alert_info("Saving predictors to {.file {output_dir}}")

# load occurrences ----
included_records <- vroom(records_file,
                          col_types=cols(wcvp_id=col_character(),
                                         specimen_id=col_character()),
                          show_col_types=FALSE, progress=FALSE)

labelled_points <- vroom(occurrence_file,
                         col_types=cols(specimen_id=col_character()),
                         show_col_types=FALSE, progress=FALSE)

# only use occurrences kept by cleaning process in calculations
occurrences <- tidylog::inner_join(included_records, 
                          labelled_points, 
                          by=c("specimen_id", "wcvp_id", "wcvp_name"))

# set up parallel computing ----
ncores <- parallelly::availableCores()
cluster <- new_cluster(ncores)
cli_alert_info("Using {.strong {ncores}} cores")

cluster_assign(
  cluster, 
  calculate_eoo=calculate_eoo,
  calculate_conr_aoo=calculate_conr_aoo,
  calculate_conr_eoo=calculate_conr_eoo,
  calculate_conr_locations=calculate_conr_locations,
  calculate_conr_unique=calculate_conr_unique,
  project_conr_points=project_conr_points
)

cluster_library(cluster, packages=c("dplyr", "rCAT", "sf"))

# summarise predictors ----

cli_alert_info("Calculating geographic and raster predictors.")
predictors <-
  occurrences %>%
  filter(! is.na(decimalLatitude), ! is.na(decimalLongitude)) %>%
  st_as_sf(coords=c("decimalLongitude", "decimalLatitude"),
           crs=4326,
           remove=FALSE) %>%
  st_transform("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs") %>%
  group_by(wcvp_id) %>%
  partition(cluster=cluster) %>%
  summarise(eoo=calculate_eoo(decimalLongitude, decimalLatitude),
            human_footprint=mean(human_footprint, na.rm=TRUE),
            human_population=min(human_population, na.rm=TRUE),
            forest_loss=mean(forest_loss, na.rm=TRUE),
            annual_temperature=mean(annual_temperature, na.rm=TRUE),
            precipitation_driest=mean(precipitation_driest, na.rm=TRUE),
            elevation=max(elevation, na.rm=TRUE),
            conr_eoo=calculate_conr_eoo(decimalLongitude, decimalLatitude),
            conr_aoo=calculate_conr_aoo(decimalLongitude, decimalLatitude),
            conr_locs=calculate_conr_locations(decimalLongitude, decimalLatitude),
            conr_obs=calculate_conr_unique(decimalLongitude, decimalLatitude),
            n_specimens=n(),
            geometry=st_union(geometry)
  ) %>%
  collect() %>%
  mutate(conr_eoo=pmax(conr_eoo, conr_aoo)) %>%
  # infinite values can break things and are harder to replace than NAs
  mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), NA_real_, .x)))

# calculate latitude of range centroid
predictors <- 
  predictors %>%
  mutate(centroid=st_centroid(geometry)) %>%
  mutate(centroid=st_transform(centroid, 4326)) %>%
  mutate(centroid_latitude=st_coordinates(centroid)[,2]) %>% 
  select(-centroid, -geometry)

cli_alert_success("Calculated predictors for {.strong {nrow(predictors)}} species")

# join to target ----
rl_assessments <- 
  assessment_file %>%
  vroom(show_col_types=FALSE, progress=FALSE) %>%
  select(one_of("id", "name", "category", "genus", "section", "srli"))

predictors <- left_join(predictors, rl_assessments, by=c("wcvp_id"="id"))
cli_alert_success("Added assessments for {.strong {sum(!is.na(predictors$category))}} species")

# save predictors to file ----
output_file <- file.path(output_dir, basename(records_file))
vroom_write(predictors, output_file)

