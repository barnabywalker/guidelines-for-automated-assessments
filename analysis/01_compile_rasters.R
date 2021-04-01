#' Compile raster data used in this study..
#' 
#' The lowest resolution rasters are 2.5 min, so all other
#' rasters will be aggregated up to this resolution. They also
#' have different projections, but it will be easier to reproject
#' the points before extracting values.
#' 
#' The rasters being downloaded are:
#'   - some climatic layers from Bioclim
#'   - elevation layers from CIGAR CSI that need stitching together
#'   - human footprint index
#'   - human population density
#'   - Hansen et al global forest loss layers that need downloading and stitching together
#' 
#' Before running this script you will need to download some things:
#'   - the elevation SRTM rasters using the bulk download link (https://drive.google.com/drive/folders/0B_J08t5spvd8RWRmYmtFa2puZEE) 
#'     from CIGAR (https://cgiarcsi.community/data/srtm-90m-digital-elevation-database-v4-1/).
#'   - the HFI raster from SEDAC CIESIN (https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint/data-download)
#'   - the HPD raster from SEDAC CIESIN (https://sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse)
#' 
#' This script uses velox for fast raster aggregation, which is no longer under development.
#' If you can't install it using devtools::install_github("hunzikp/velox"), you might need to
#' swap the code out for normal raster aggregation. Velox also loads the rasters
#' into memory, so if you find you're getting memory errors, you should try and swap it out
#' for normal raster aggregation.
#' 
#' This may take a while to run, and the downloads of the SRTM and forest loss are quite big.
#' You might want to run chunks of this at a time and delete the full raster downloads when it's done.
#' 

# libraries ----
library(here)      # handling file paths
library(raster)    # handling rasters
library(velox)     # fast aggregation of rasters
library(purrr)     # map functions across lists
library(glue)      # string interpolation
library(stringr)   # manipulating strings
library(progress)  # making progress bars
library(rgee)

# setup ----

message("Compiling raster data...")

data_dir <- here("data/rasters")

output_dir <- here("output/rasters")
dir.create(output_dir, showWarnings=FALSE)

# check that rasters have been downloaded ----
hfi_file <- paste(data_dir, "wildareas-v3-2009-human-footprint.tif", sep="/")
hfi_url <- "https://sedac.ciesin.columbia.edu/data/set/wildareas-v3-2009-human-footprint/data-download"
if (!file.exists(hfi_file)) {
  msg <- glue("Human footprint index raster not found.",
              " Please download from {hfi_url} to {hfi_file}.")
  stop(msg)
}

hpd_file <- paste(data_dir, "gpw_v4_population_density_rev11_2020_2pt5_min.tif", sep="/")
hpd_url <- "https://sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse"
if (!file.exists(hpd_file)) {
  msg <- glue("Human population density raster not found.",
              " Please dowpnload from {hpd_url} to {hpd_file}.")
  stop(msg)
}

# bioclimatic rasters ----
bio_stack <- getData("worldclim", var="bio", res=2.5, path=here("data/rasters"))

writeRaster(bio_stack[[1]], paste(output_dir, "annual_temperature.tif", sep="/"))
writeRaster(bio_stack[[4]], paste(output_dir, "temperature_seasonality.tif", sep="/"))
writeRaster(bio_stack[[12]], paste(output_dir, "annual_precipitation.tif", sep="/"))
writeRaster(bio_stack[[15]], paste(output_dir, "precipitation_seasonality.tif", sep="/"))
writeRaster(bio_stack[[17]], paste(output_dir, "precipitation_driest.tif", sep="/"))

message("  Downloaded bioclimatic rasters!")

# human footprint index ----
r <- raster(hfi_file)

# replace the NA placeholder with NA so it doesn't throw off aggregation
values <- values(r)
values[values == 128] <- NA
values(r) <- values

# save
writeRaster(r, paste(output_dir, "human_footprint.tif", sep="/"))
message("  Resaved human footprint raster!")

# human population density ----
r <- raster(hpd_file)
writeRaster(r, paste(output_dir, "human_population.tif", sep="/"))
message("  Resaved human population raster!")

# global forest loss ----
# using google earth engine for this!

# elevation ----
# using google earth engine for this!