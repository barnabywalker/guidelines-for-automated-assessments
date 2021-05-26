#' Apply cleaning and filtering steps to occurrence files.
#' 
#' Part of this study is to evaluate the effects of different
#' automated cleaning steps on the performance of automated
#' assessment methods across different species groups. The cleaning
#' is split into two groups: filtering and coordinate cleaning.
#' 
#' The filtering steps are:
#' 1. No filtering
#' 2. Keeping only occurrences based on preserved specimens
#' 3. Removing occurrences at exactly the same coordinates for each species
#' 4. Both 2. and 3. applied together
#' 
#' The coordinate cleaning steps are:
#' A. No coordinate cleaning.
#' B. Removing occurrences at (0,0) coordinates
#' C. Removing occurrences in the sea, at institution locations,
#'    and at country centroids
#' D. Removing occurrences that are not in the species native range,
#'    as recorded in POWO.
#'    
#' Each coordinate cleaning step is applied on top of the previous ones.
#' 
#' This script saves a new file for each step.
#' 
#' EXPECTED INPUTS:
#'  - `occurrence_file`: path to a file with occurrences to be cleaned
#'  - `distribution_file`: path to file containing distributions for species in the occurrence file
#'  - `output_dir`: path to a directory to save in (which must exist)

# libraries ----
library(here)               # handle file paths
library(dplyr)              # manipulate data
library(vroom)              # fast reading/writing for text files
library(glue)               # string interpolation
library(CoordinateCleaner)  # utilities for cleaning occurrence coordinates
library(sf)                 # handle spatial data

# load data ----
output_cols <- c("specimen_id", "wcvp_name", "wcvp_id")
occurrences <- vroom(occurrence_file)

distributions <- vroom(distribution_file)
wgsrpd3 <- st_read(here("data/wgsrpd/level3/level3.shp"))

# loop through filter steps and clean ----
for (step in 1:4) {
  filtered <- as_tibble(occurrences)
  
  # applies to steps 2 and 4 (only keep specimens)
  if (step %% 2 == 0) {
    filtered <- 
      filtered %>%
      filter(basisOfRecord %in% "PRESERVED_SPECIMEN")
  }
  
  # applies to steps 3 and 4 (remove duplicates)
  if (step > 2) {
    filtered <- 
      filtered %>%
      distinct(wcvp_id, decimalLongitude, decimalLatitude, .keep_all=TRUE)
  }
  
  ## cleaning step A ----
  cleaned <- 
    filtered %>%
    filter(hasCoordinate, ! hasGeospatialIssue)
  
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-A.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
  
  ## cleaning step B ----
  cleaned <-
    cleaned %>%
    cc_zero(
      lon="decimalLongitude", 
      lat="decimalLatitude",
      verbose=FALSE,
      value="clean"
    )
  
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-B.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
  
  ## cleaning step C ----
  cleaned <- 
    cleaned %>%
    cc_equ(
      lat="decimalLatitude", 
      lon="decimalLongitude",
      verbose=FALSE,
      value="clean"
    ) %>%
    cc_sea(
      lat="decimalLatitude", 
      lon="decimalLongitude",
      verbose=FALSE,
      value="clean"
    ) %>%
    clean_coordinates(
      lat="decimalLatitude", 
      lon="decimalLongitude",
      species="wcvp_id",
      tests=c("capitals", "institutions", "gbif", "centroids"),
      value="clean")
  
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-C.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
  
  ## cleaning step D ----
  cleaned <-
    cleaned %>%
    st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), 
             crs=st_crs("EPSG:4326")) %>%
    st_join(
      wgsrpd3 %>% select(distribution=LEVEL3_COD)
    ) %>%
    tidylog::inner_join(
      distributions,
      by=c("wcvp_id"="id", "distribution")
    ) %>%
    st_drop_geometry() %>%
    distinct(specimen_id, .keep_all=TRUE)
  
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-D.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
}