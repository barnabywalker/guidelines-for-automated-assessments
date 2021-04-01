#' Annotate occurrence records that have coordinates with values
#' that will be used to calculated predictors for species.
#' 
#' Most of the predictors used, apart from EOO and latitude of range centroid,
#' use values extracted from rasters and aggregated across occurrence records
#' for a species. This script extracts these values for every occurrence record
#' that it can. Before extracting raster values, each occurrence record is
#' buffered by 5 km.
#' 
#' 

# libraries ----
library(here)
library(vroom)
library(dplyr)
library(sf)
library(rgee)

source(here("R/helper_functions.R"))

# set up earth engine ----
ee_Initialize(drive=TRUE)

# set buffer size ----
buffer_size <- 5000

# load occurrence data ----

if (! exists("occurrence_file")) {
  stop("You must provide the path to a file of occurrence records")
}

keep_cols <- c("specimen_id", 
               "decimalLatitude", "decimalLongitude", 
               "hasCoordinate", "hasGeospatialIssue")

occurrences <- vroom(occurrence_file,
                     col_select=all_of(keep_cols),
                     col_types=cols(specimen_id=col_character(),
                                    wcvp_id=col_character()))

# remove missing occurrences ----
missing_points <- filter(occurrences, ! hasCoordinate | hasGeospatialIssue)
occurrences <- filter(occurrences, hasCoordinate, ! hasGeospatialIssue)

# make occurrences spatial object ----
occurrences_sf <-
  occurrences %>%
  st_as_sf(coords=c("decimalLongitude", "decimalLatitude"),
           crs="EPSG:4326",
           remove=FALSE)
  
# label with raster values ----

# TODO: Make this work with EE rasters for elevation and forest loss

# buffering takes a lot of memory, so we buffer and extract in chunks
# make a list of all the rasters we want to use
raster_files <- list(
  # elevation
  elevation=here("output/rasters/elevation_merged.tif"),
  # hansen forest loss
  forest_loss=here("output/rasters/forest_loss.tif"),
  # human footprint index
  hfi=here("output/rasters/human_footprint.tif"),
  # human population density
  hpd=here("output/rasters/human_population.tif"),
  # environmental rasters from worldclim
  temperature_annual=here("output/rasters/annual_temperature.tif"),
  precipitation_driest=here("output/rasters/precipitation_driest.tif")
)

# define the aggregation functions to use when extracting from the rasters
fcns <- list(
  # elevation
  elevation="mean",
  # hansen forest loss
  forest_loss=function(values, coverage_fraction) mean(values > 0, na.rm=TRUE),
  # human footprint index
  hfi="mean",
  # human population density
  hpd="mean",
  # environmental rasters from worldclim
  temperature_annual="mean",
  precipitation_driest="mean"
)

extracted_values <- extract_buffered_chunks(occurrences_sf, raster_files, fcns,
                                            chunk_size=100000)

# join up all the chunks into a single data frame
labelled_occurrences <-
  occurrences_sf %>%
  st_drop_geometry() %>%
  bind_cols(extracted_values)

# save annotated points ----
labelled_occurrences <- bind_rows(labelled_occurrences, missing_points)

vroom_write(labelled_occurrences, output_occurrence_file)
