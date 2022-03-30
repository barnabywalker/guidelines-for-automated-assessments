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
#'  - `output_dir`: path to a directory to cleaned occurrence files to
#'  - `distribution_dir`: path to a directory containing species distribution files,
#'        possibly outputs from `analysis/02_collate_species.R`
#'  - `shape_file`: path to a shape file of the WGSRPD botanical countries (level 3)
#'  - `occurrence_file`: path to a file with occurrences to be cleaned, 
#'        possibly an output from `analysis/04_label_points.R`
#'  - `already_cleaned`: boolean, indicates if you just want to resave an occurrence file, 
#'        e.g. one from a monographic database
#' 
#' EXAMPLE CLI:
#'  Rscript analysis/05_clean_occurrences.R --occurrence_file=myrcia-gbif_labelled-occurrences.csv
#' 
#' EXAMPLE SOURCE:
#'  output_dir <- "output/cleaned-occurrences"
#'  distribution_dir <-  "output/species-lists"
#'  shape_file <- "data/wgsrpd/level3/level3.shp"
#'  already_cleaned <- FALSE
#'  occurrence_file <- "myrcia-gbif_labelled-occurrences.csv"
#'  source("analysis/05_clean_occurrences.R")
#'  

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(here)               # handle file paths
shhlibrary(cli)                # nice command line interface
shhlibrary(dplyr)              # manipulate data
shhlibrary(vroom)              # fast reading/writing for text files
shhlibrary(glue)               # string interpolation
shhlibrary(CoordinateCleaner)  # utilities for cleaning occurrence coordinates
shhlibrary(sf)                 # handle spatial data
shhlibrary(stringr)             # handle strings

sf::sf_use_s2(FALSE)

# CLI -----
cli_h1("Cleaning occurrence records")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output/cleaned-occurrences",
    distribution_dir="output/species-lists",
    shape_file="data/wgsrpd/level3/level3.shp",
    already_cleaned=FALSE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  output_dir <- args$output_dir
  occurrence_file <- args$occurrence_file
  distribution_dir <- args$distribution_dir
  shape_file <- args$shape_file
  already_cleaned <- args$already_cleaned
}

if (! exists("occurrence_file")) {
  cli_abort(c(
    "no path to occurrences provided",
    "x"="You must provide a path to file of occurrence records as {.var occurrence_file}."
  ))
}

if (! exists("shape_file")) {
  cli_abort(c(
    "no path to WGSRPD shape file provided",
    "x"="You must provide a path to a shape file for WGSRPD level 3 as {.var shape_file}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save labelled records provided",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

if (! exists("distribution_dir")) {
  cli_abort(c(
    "no path to directory with distribution files",
    "x"="You must provide a path to a directory with distributions as {.var distribution_dir}"
  ))
}

dir.create(output_dir, showWarnings=FALSE)
name <- str_extract(basename(occurrence_file), "^[\\w-]+(?=_)")

distribution_file <- file.path(distribution_dir, glue("{str_extract(name, '^\\\\w+(?=-)')}_distributions.csv"))

cli_alert_info("Cleaning occurrence records from {.file {occurrence_file}}")
cli_alert_info("Loading distributions from {.file {distribution_file}}")
cli_alert_info("Saving cleaned occurrence records to {.file {output_dir}}")

# load data ----
output_cols <- c("specimen_id", "wcvp_name", "wcvp_id")
occurrences <- vroom(occurrence_file, show_col_types=FALSE, progress=FALSE)

distributions <- vroom(distribution_file, show_col_types=FALSE, progress=FALSE)
wgsrpd3 <- st_read(shape_file)

# loop through filter steps and clean ----
if (already_cleaned) {
  output_path <- glue(output_dir, "/{name}_filter-{1}_clean-db.csv")
  cli_alert_info("File already cleaned, saving as is to {.file {output_path}}")
  occurrences %>%
      select(all_of(output_cols)) %>%
      vroom_write(output_path)
} else {
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
    
    output_path <- glue(output_dir, "/{name}_filter-{step}_clean-A.csv")
    cleaned %>%
      select(all_of(output_cols)) %>%
      vroom_write(output_path)
    
    cli_alert_success("Saved {.strong {nrow(cleaned)}} occurrences for cleaning step A, filter step {step}.")
    
    ## cleaning step B ----
    cleaned <-
      cleaned %>%
      cc_zero(
        lon="decimalLongitude", 
        lat="decimalLatitude",
        verbose=FALSE,
        value="clean"
      )
    
    output_path <- glue(output_dir, "/{name}_filter-{step}_clean-B.csv")
    cleaned %>%
      select(all_of(output_cols)) %>%
      vroom_write(output_path)
    cli_alert_success("Saved {.strong {nrow(cleaned)}} occurrences for cleaning step B, filter step {step}.")
    
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
        value="clean",
        verbose=FALSE)
    
    output_path <- glue(output_dir, "/{name}_filter-{step}_clean-C.csv")
    cleaned %>%
      select(all_of(output_cols)) %>%
      vroom_write(output_path)
    
    cli_alert_success("Saved {.strong {nrow(cleaned)}} occurrences for cleaning step C, filter step {step}.")
    
    ## cleaning step D ----
    cleaned <-
      cleaned %>%
      st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), 
              crs=st_crs("EPSG:4326")) %>%
      st_join(
        wgsrpd3 %>% select(distribution=LEVEL3_COD)
      ) %>%
      inner_join(
        distributions,
        by=c("wcvp_id"="id", "distribution")
      ) %>%
      st_drop_geometry() %>%
      distinct(specimen_id, .keep_all=TRUE)
    
    output_path <- glue(output_dir, "/{name}_filter-{step}_clean-D.csv")
    
    cleaned %>%
      select(all_of(output_cols)) %>%
      vroom_write(output_path)
    
    cli_alert_success("Saved {.strong {nrow(cleaned)}} occurrences for cleaning step D, filter step {step}.")
  }
}

cli_alert_success("Saved cleaned occurrences records for {.val {name}} to {.file {output_dir}}")
