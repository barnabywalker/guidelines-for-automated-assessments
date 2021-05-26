#' Prepare species-level predictors from labelled occurrence records.
#' 
#' This script takes a list of occurrence records labelled with predictor values (the output
#' of script `04_annotate_points.R`) and a list of occurrence records to include (the output
#' of script `05_clean_occurrences.R`) as input and summarises the predictor values for each species.
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
#' 
#' The script then joins these predictors to the prediction target - IUCN RL assessment status.
#' Species without an assessment will be retained for prediction and the category will be NA.
#' 
#' EXPECTED INPUTS:
#'  - `records_file`: path to file containing list of records to include (output of cleaning script)
#'  - `occurrence_file`: path to file containing labelled occurrence records
#'  - `output_dir`: path to a directory to save predictors (must exist)

# libraries ----
library(here)         # handle file paths
library(multidplyr)   # parallelise data manipulation
library(dplyr)        # manipulate data
library(sf)           # handle spatial data
library(vroom)        # fast reading/writing for text files

source(here("R/helper_functions.R"))

# load occurrences ----
included_records <- vroom(records_file,
                          col_types=cols(wcvp_id=col_character(),
                                         specimen_id=col_character()))

labelled_points <- vroom(occurrence_file,
                         col_types=cols(specimen_id=col_character()))

# only use occurrences kept by cleaning process in calculations
occurrences <- inner_join(included_records, 
                          labelled_points, 
                          by=c("specimen_id"))

# set up parallel computing ----
cluster <- new_cluster(parallel::detectCores())
cluster_assign(cluster, calculate_eoo=calculate_eoo)
cluster_library(cluster, packages=c("dplyr", "rCAT", "sf"))

# summarise predictors ----
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
            hfi=mean(hfi, na.rm=TRUE),
            hpd=min(hpd, na.rm=TRUE),
            forest_loss=mean(forest_loss, na.rm=TRUE),
            temperature_annual=mean(temperature_annual, na.rm=TRUE),
            precipitation_driest=mean(precipitation_driest, na.rm=TRUE),
            elevation=max(elevation, na.rm=TRUE),
            n_specimens=n(),
            geometry=st_union(geometry)
  ) %>%
  collect()

# calculate latitude of range centroid
predictors <- 
  predictors %>%
  mutate(centroid=st_centroid(geometry)) %>%
  mutate(centroid=st_transform(centroid, 4326)) %>%
  mutate(centroid_latitude=st_coordinates(centroid)[,2]) %>% 
  select(-centroid, -geometry)

# join to target ----
rl_assessments <- vroom(assessment_file,
                        col_select=c(id, name, category))

predictors <- left_join(predictors, rl_assessments, by=c("wcvp_id"="id"))

# save predictors to file ----
vroom_write(predictors, output_file)