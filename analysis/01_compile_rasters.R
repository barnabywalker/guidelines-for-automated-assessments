#' Compile raster data used in this study.
#' 
#' The lowest resolution rasters are 2.5 min, so all other
#' rasters will be aggregated up to this resolution. 
#' They also have different projections, but it will be easier 
#' to reproject the points before extracting values.
#' 
#' The rasters being downloaded are:
#'   - some climatic layers from Bioclim
#'   - the elevation raster used to generate the Bioclim layers
#'   - human footprint index
#'   - human population density
#'   - Hansen et al global forest loss layers that need downloading and stitching together
#' 
#' This script uses velox for fast raster aggregation, which is no longer under development.
#' If you can't install it using remotes::install_github("hunzikp/velox"), you might need to
#' swap the code out for normal raster aggregation. Velox also loads the rasters
#' into memory, so if you find you're getting memory errors, you should try and swap it out
#' for normal raster aggregation.
#' 
#' This may take a while to run, and the download of forest loss is quite big.
#' You might want to run chunks of this at a time and delete the full raster downloads when it's done.
#' 
#' EXPECTED INPUTS:
#'  - `download_dir`: path to a directory to download rasters into
#'  - `output_dir`: path to a directory to save the processed rasters to
#'  - `overwrite`: whether to overwrite the downloaded and processed rasters if they already exist
#' 
#' EXAMPLE CLI:
#'  RScript analysis/01_compile_rasters.R
#' 
#' EXAMPLE SOURCE:
#'  download_dir <- "data/raster-downloads"
#'  output_dir <- "data/rasters"
#'  overwrite <- FALSE
#'  source("analysis/01_compile_rasters.R")
#' 

# libraries ----
# you might want to remove this for debugging packages
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))

shhlibrary(here)      # handling file paths
shhlibrary(raster)    # handling rasters
shhlibrary(velox)     # fast aggregation of rasters
shhlibrary(cli)       # nice command line interface
shhlibrary(purrr)     # map functions across lists
shhlibrary(glue)      # string interpolation
shhlibrary(stringr)   # manipulating strings

source(here("R/raster_functions.R"))

# CLI ----

cli_h1("Download and process raster data")

if (sys.nframe() == 0L) {
  default_args <- list(
    download_dir="data/raster-downloads",
    output_dir="data/rasters",
    overwrite=FALSE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  download_dir <- args$download_dir
  output_dir <- args$output_dir
  overwrite <- args$overwrite
}

if (! exists("download_dir")) {
  cli_abort(c(
    "no path to download rasters provided",
    "x"="You must provide the download path as the variable {.var download_dir}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save processed rasters provided",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

if (! exists("overwrite")) {
  overwrite <- FALSE
}

dir.create(download_dir, showWarnings=FALSE)
dir.create(output_dir, showWarnings=FALSE)

cli_alert_info("Downloading rasters to {.file {download_dir}}")
cli_alert_info("Saving processed rasters to {.file {output_dir}}")
cli_alert_warning(c("Existing rasters ",
                  "{.strong will {ifelse(overwrite, '', 'not ')}be}",
                  " overwritten."))

# human footprint index ----
cli_h2("Human footprint index")
hfi_url <- "https://wcshumanfootprint.org/data/HFP2009.zip"
download_raster(hfi_url, output_dir, tempdir=download_dir, name="human footprint", overwrite=overwrite)

# human population density ----
cli_h2("Human population density")
hpd_url <- "https://sedac.ciesin.columbia.edu/downloads/data/gpw-v4/gpw-v4-population-density-rev11/gpw-v4-population-density-rev11_2020_2pt5_min_tif.zip"
download_raster(hpd_url, output_dir, tempdir=download_dir, name="human population", overwrite=overwrite)

# bioclimatic rasters ----
cli_h2("Bioclimatic variables")

bio_stack <- getData("worldclim", var="bio", res=2.5, path=download_dir)
cli_alert_success("Bioclim rasters loaded from {.file {download_dir}}.")

save_bioclim(1, output_dir, name="annual temperature", stack=bio_stack, overwrite=overwrite)
save_bioclim(4, output_dir, name="temperature seasonality", stack=bio_stack, overwrite=overwrite)
save_bioclim(12, output_dir, name="annual precipitation", stack=bio_stack, overwrite=overwrite)
save_bioclim(15, output_dir, name="precipitation seasonality", stack=bio_stack, overwrite=overwrite)
save_bioclim(17, output_dir, name="precipitation driest", stack=bio_stack, overwrite=overwrite)

# global forest loss ----
cli_h2("Forest loss")
hansen_dir <- file.path(download_dir, "hansen-forest-loss")
savefile <- file.path(output_dir, "forest-loss.tif")

if (length(list.files(hansen_dir)) < 504 & ! file.exists(savefile)) {
  download_loss_tiles(hansen_dir)  
}

if (! file.exists(savefile) | overwrite) {
  merged_raster <- merge_tiles(hansen_dir, tolerance=0.1)
  # aggregate to final desired resolution
  merged_raster <- aggregate_raster(merged_raster, desired_res=2.5/60, aggfun="mean")
  writeRaster(merged_raster, savefile, overwrite=TRUE)
  
  cli_alert_success("forest loss raster saved to {.file {savefile}}.")
} else {
  cli_alert_success("forest loss raster already exists at {.file {savefile}}.")
}

# elevation ----
cli_h2("Elevation (SRTM)")

srtm_url <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_elev.zip"
download_raster(srtm_url, output_dir, tempdir=download_dir, name="elevation", overwrite=overwrite)
