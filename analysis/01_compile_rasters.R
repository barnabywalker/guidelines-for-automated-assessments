#' Compile raster data used in this study..
#' 
#' The lowest resolution rasters are 2.5 min, so all other
#' rasters will be aggregated up to this resolution (apart 
#' from the forest loss raster, which was such a small resolution 
#' aggregating to 2.5 min would take much too long). 
#' They also have different projections, but it will be easier 
#' to reproject the points before extracting values.
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
dir.create(here("data/rasters/hansen"))

# need to request the URLs to each tile to download
raster_urls <- readLines("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2018-v1.6/lossyear.txt")

pb <- progress_bar$new(total=length(raster_urls),
                       format="downloading and aggregating forest loss rasters [:bar] :current/:total (:percent)")
pb$tick(0)

# download and aggregate each tile individually
for (i in 1:length(raster_urls)) {
  pb$tick()
  
  url <- raster_urls[i]
  name <- str_extract(url, "\\d+(N|S)_\\d+(E|W)\\.tif")
  
  ras <- raster(url)
  # convert to binary loss/no loss - comes as year of loss
  values(ras) <- ifelse(values(ras) > 0, 1, 0)
  vx <- velox(ras)
  
  # aggregating with mean to get proportion of pixels in set that have loss
  vx$aggregate(factor=40, aggtype="mean")
  vx$write(path=here("data/rasters/hansen", name), overwrite=TRUE)
}

# now we load and merge the aggregated tiles
raster_files <- list.files(here("data/rasters/hansen"))
rasters <- purrr::map(raster_files, ~raster(here("data/rasters/hansen", .x)))
merged_raster <- do.call(merge, c(rasters, tolerance=0.1))

writeRaster(merged_raster, paste(output_dir, "forest_loss.tif", sep="/"),
            overwrite=TRUE)

message("  Processed and saved forest loss rasters!")

# elevation ----

# get tile filenames from download
elevation_files <- list.files(here("data/rasters/elevation/SRTMv4.1/tifs"), pattern=".+\\.tif$")


pb <- progress_bar$new(total=length(elevation_files), 
                       format="aggregating elevation rasters [:bar] :current/:total (:percent)")
pb$tick(0)
for (elevation_file in elevation_files) {
  pb$tick()
  
  tile_id <- str_extract(elevation_file, "(?<=srtm_).+(?=\\.tif)")
  out_file <- paste0("elevation_tile_", tile_id, ".tif")
  
  # aggregating to 50x resolution from 90m resolution, to match bioclim rasters
  tile <- velox(here("data/rasters/elevation/SRTMv4.1/tifs", elevation_file))
  tile$aggregate(factor=50, aggtype="mean")
  writeRaster(tile$as.RasterLayer(), paste(output_dir, "elevation", out_file, sep="/"), overwrite=TRUE)
}

# now we load in all the aggregated tiles and merge into a single raster
elevation_files <- list.files(output_dir, pattern="elevation_tile_", full.names=TRUE)
raster_list <- map(elevation_files, raster)
merged_raster <- do.call(merge, c(raster_list, tolerance=0.1))

writeRaster(merged_raster, paste(output_dir, "elevation_merged.tif", sep="/"))

message("  Processed and saved elevation rasters!")