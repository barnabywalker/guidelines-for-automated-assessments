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

# libraries ----
library(here)
library(dplyr)
library(vroom)
library(glue)
library(CoordinateCleaner)
library(sf)

# load occurrence file ----
output_cols <- c("specimen_id", "wcvp_name", "wcvp_id")
occurrences <- vroom(occurrence_file)

distributions <- vroom(distribution_file)
wgsrpd3 <- st_read(here("data/wgsrpd/level3/level3.shp"))

# loop through filter steps and clean ----
for (step in 1:4) {
  filtered <- as_tibble(occurrences)
  
  if (step %% 2 == 0) {
    filtered <- 
      filtered %>%
      filter(basisOfRecord %in% "PRESERVED_SPECIMEN")
  }
  
  if (step > 2) {
    filtered <- 
      filtered %>%
      distinct(wcvp_id, decimalLongitude, decimalLatitude, .keep_all=TRUE)
  }
  
  cleaned <- 
    filtered %>%
    filter(hasCoordinate, ! hasGeospatialIssue)
  
  # save cleaning step A
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-A.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
  
  cleaned <-
    cleaned %>%
    cc_zero(
      lon="decimalLongitude", 
      lat="decimalLatitude",
      verbose=FALSE,
      value="clean"
    )
  
  # save cleaning step B
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-B.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
  
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
  
  # save cleaning step C
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-C.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
  
  # TODO: CHECK THIS!!!
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
  
  # save cleaning step C
  output_path <- glue(output_dir, "/{group}_filter-{step}_clean-D.csv")
  cleaned %>%
    select(all_of(output_cols)) %>%
    vroom_write(output_path)
}