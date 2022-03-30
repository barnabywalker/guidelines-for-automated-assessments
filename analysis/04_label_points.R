#' Annotate occurrence records that have coordinates with values
#' that will be used to calculated predictors for species.
#' 
#' Most of the predictors used (apart from EOO, latitude of range centroid, and the ConR parameters)
#' use values extracted from rasters and aggregated across occurrence records
#' for a species. This script extracts these values for every occurrence record
#' that it can. Before extracting raster values, each occurrence record is
#' buffered by a 5 km radius.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save the labelled occurrences to
#'  - `raster_dir`: path to a directory with rasters to extract values from
#'  - `occurrence_file`: path to file containing occurrences to be labelled
#'  - `buffer`: the radius to buffer points by, in km (defaults to 5 km)
#' 
#' EXAMPLE CLI:
#'  Rscript analysis/04_label_points.R --occurrence_file=myrcia-gbif_occurrences.csv
#' 
#' EXAMPLE SOURCE:
#'  raster_dir <- "data/rasters",
#'  output_dir <- "output/labelled-occurrences"
#'  buffer <- 5
#'  occurrence_file <- myrcia-gbif_occurrences.csv
#'  source("analysis/04_label_points.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(here)   # handle file paths
shhlibrary(cli)    # nice command line interfaces
shhlibrary(vroom)  # fast reading/writing for text files
shhlibrary(dplyr)  # manipulating data
shhlibrary(sf)     # hand spatial data
shhlibrary(glue)   # string interpolation
shhlibrary(stringr) # handle strings

source(here("R/predictor_functions.R"))

# CLI ----
cli_h1("Annotate occurrence records with raster values")

if (sys.nframe() == 0L) {
  default_args <- list(
    raster_dir="data/rasters",
    output_dir="output/labelled-occurrences",
    buffer=5
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  raster_dir <- args$raster_dir
  output_dir <- args$output_dir
  buffer <- args$buffer
  occurrence_file <- args$occurrence_file
}

if (! exists("raster_dir")) {
  cli_abort(c(
    "no path to rasters provided",
    "x"="You must provide a path to a folder of rasters as {.var raster_dir}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save labelled records provided",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

dir.create(output_dir, showWarnings=FALSE)
name <- str_extract(basename(occurrence_file), "^[\\w-]+(?=_)")
savefile <- file.path(output_dir, glue("{name}_labelled-occurrences.csv"))

cli_alert_info("Labelling occurrence records from {.file {occurrence_file}}")
cli_alert_info("Loading rasters from {.file {raster_dir}}")
cli_alert_info("Saving labelled occurrence records to {.file {savefile}}")
cli_alert_info("Occurrence points will be buffered by {.strong {buffer} km}")

buffer <- buffer * 1000

# load occurrence data ----
keep_cols <- c("specimen_id", "wcvp_id", "wcvp_name",
               "decimalLatitude", "decimalLongitude", "basisOfRecord", 
               "hasCoordinate", "hasGeospatialIssue")

occurrences <- vroom(occurrence_file,
                     col_select=all_of(keep_cols),
                     col_types=cols(specimen_id=col_character(),
                                    wcvp_id=col_character()),
                     progress=FALSE)

cli_alert_success("{.strong {nrow(occurrences)}} records loaded from {.file {occurrence_file}}")

# remove missing occurrences ----
missing_points <- filter(occurrences, ! hasCoordinate | hasGeospatialIssue)
occurrences <- filter(occurrences, hasCoordinate, ! hasGeospatialIssue)

cli_alert_info("{.strong {nrow(occurrences)}} records have coordinates")

# make occurrences spatial object ----
occurrences_sf <-
  occurrences %>%
  st_as_sf(coords=c("decimalLongitude", "decimalLatitude"),
           crs="EPSG:4326",
           remove=FALSE)
  
# label with raster values ----
raster_files <- list.files(raster_dir, full.names=TRUE, pattern="\\w.tiff?$")

raster_names <- str_extract(basename(raster_files), "[a-z]+(-[a-z]+)?")
cli_alert_info("Labelling records with values from {.val {str_replace_all(raster_names, '-', ' ')}}")
raster_names <- str_replace_all(raster_names, "-", "_")
names(raster_files) <- raster_names

## extract values and join to occurrences ----
extracted_values <- extract_buffered_chunks(occurrences_sf, raster_files, chunk_size=100000)

cli_alert_success("Labelled occurrences that had coordinates")

labelled_occurrences <-
  occurrences_sf %>%
  st_drop_geometry() %>%
  bind_cols(extracted_values)

# some points get double labelled
labelled_occurrences <-
  labelled_occurrences %>%
  distinct(specimen_id, wcvp_name, wcvp_id, .keep_all=TRUE)

# save annotated points ----
labelled_occurrences <- bind_rows(labelled_occurrences, missing_points)

vroom_write(labelled_occurrences, savefile)
cli_alert_success("Labelled occurrences save to {.file {savefile}}")
