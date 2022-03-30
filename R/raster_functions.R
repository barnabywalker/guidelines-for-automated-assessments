#' Functions for downloading and processing rasters

#' Download a raster from a URL.
#' 
#' Download a raster file that is stored at a URL in a zip file, extract the
#' raster, and save the resulting tif file to a destination.
#' 
#' @param url URL pointing to a raster file to download.
#' @param savedir directory to save the raster to.
#' @param tempdir optional directory to download raster zip file to before extracting.
#' @param name optional name to save the raster under.
#' @param overwrite whether to overwrite the file if it already exists.
#' 
download_raster <- function(url, savedir, tempdir=NULL, name=NULL, overwrite=FALSE,
                            username=NULL, password=NULL) {
  if (is.null(tempdir)) {
    tempdir <- savedir
  }
  
  if (is.null(name)) {
    name <- str_extract(url, "(?<=/).+(?=\\.zip)")
  }
  
  savename <- str_to_lower(name)
  savename <- str_replace_all(savename, " ", "-")
  
  zipfile <- file.path(tempdir, glue("{savename}-raw.zip"))
  savefile <- file.path(savedir, glue("{savename}.tiff"))
  
  if (! file.exists(savefile) | overwrite) {
    if (str_detect(url, "sedac.ciesin")) {
      res <- download_sedac(url=url, destfile=zipfile, username=username, password=password)
    } else {
      res <- download.file(url=url, destfile=zipfile, quiet=TRUE, 
                           mode="wb", cacheOK=TRUE)  
    }
    
    cli_alert_success("{name} raster downloaded from {.url {url}}.")
    
    unzipfiles <- unzip(zipfile, list=TRUE)
    unzipfile <- unzipfiles[str_detect(unzipfiles$Name, "\\.tif$"),]$Name
    unzip(zipfile, exdir=tempdir)
    removed <- file.remove(zipfile)
    
    copied <- file.copy(file.path(tempdir, unzipfile), savefile)
    
    cli_alert_success("{name} raster saved to {.file {savefile}}.")
  } else {
    cli_alert_success("{name} raster already exists at {.file {savefile}}.")
  }
  
  invisible(savefile)
}

#' Download a file from SEDAC to a destination.
#'
#' Use a username and password to log in to the NASA service
#' and then download the desired file. username and password can either be 
#' passed as arguments or set as environment variables. To sign up to the service,
#' visit https://urs.earthdata.nasa.gov/home.
#'
#' @param url URL pointing to the file from SEDAC to download.
#' @param destfile path to save the file to.
#' @param username your username to log in to the service.
#' @param password your password to log in to the service.
#'
download_sedac <- function(url, destfile, username=NULL, password=NULL) {
  if (is.null(username)) {
    username <- Sys.getenv("SEDAC_EMAIL")
  }
  
  if (is.null(password)) {
    password <- Sys.getenv("SEDAC_PASSWORD")
  }
  
  login_url <- "https://urs.earthdata.nasa.gov/home"
  sesh <- rvest::session(login_url)
  form <- rvest::html_form(sesh)[[1]]
  filled_form <- rvest::html_form_set(form, username=username, password=password)
  login_page <- rvest::session_submit(sesh, filled_form)
  download <- rvest::session_jump_to(login_page, url)
  
  writeBin(download$response$content, destfile)
  
  invisible(destfile)
}


#' Save a single bioclim raster layer to a destination.
#'
#' Save the specified raster layer from a stack of bioclim rasters to a file.
#' If no stack is provided, the bioclim rasters will be loaded from or downloaded to
#' a specified folder. The different layers in the bioclim stack are detailed here:
#' https://www.worldclim.org/data/bioclim.html.
#'
#' @param layer_idx The layer number to save from a stack of rasters.
#' @param savedir The directory to save the raster into.
#' @param tempdir An optional directory to load the bioclim rasters from or download to.
#' @param name An optional name to give the raster.
#' @param stack The stack of bioclim rasters, will be loaded/downloaded if not provided.
#' @param overwrite Whether to overwrite the save file if it already exists.
#'
save_bioclim <- function(layer_idx, savedir, tempdir=NULL, name=NULL, stack=NULL, overwrite=FALSE) {
  if (is.null(tempdir)) {
    tempdir <- savedir
  }
  
  if (is.null(stack)) {
    stack <- raster::getData("worldclim", var="bio", res=2.5, path=tempdir)
    cli::cli_alert_success("Bioclim rasters loaded from {.file {tempdir}}.")
  }
  
  if (is.null(name)) {
    name <- glue::glue("bio{layer_idx}")
  }
  
  savename <- stringr::str_to_lower(name)
  savename <- stringr::str_replace_all(savename, " ", "-")
  savefile <- file.path(savedir, glue::glue("{savename}.tiff"))
  
  if (! file.exists(savefile) | overwrite) {
    raster::writeRaster(stack[[layer_idx]], savefile, overwrite=overwrite)
    cli::cli_alert_success("Saved annual temperature raster to {.file {savefile}}")
  } else {
    cli::cli_alert_success("Annual temperature raster exists at {.file {savefile}}")
  }
  
  invisible(savefile)
}

#' Download all the tiles of the Global Forest Watch loss-year raster.
#'
#' Downloads all the tiles to the desired folder, converting the year of loss
#' to a binary loss/no loss and aggregating to the desired resolution.
#' Takes a long time.
#'
#' @param savedir the directory to save the tiles in.
#' @param resolution the desired resolution of the output rasters.
#' 
download_loss_tiles <- function(savedir, resolution=2.5) {
  dir.create(savedir, showWarnings=FALSE)
  cli_alert_info("Downloading forest lost tiles to {.file {savedir}}")
  # need to request the URLs to each tile to download
  raster_urls <- readLines("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2018-v1.6/lossyear.txt")
  
  cli_progress_bar("downloading and aggregating tiles", total=length(raster_urls))
  
  # download and aggregate each tile individually
  for (i in 1:length(raster_urls)) {
    url <- raster_urls[i]
    name <- str_extract(url, "\\d+(N|S)_\\d+(E|W)\\.tif")
    download.file(url, destfile=file.path(savedir, name), quiet=TRUE, cacheOK=TRUE)
    
    ras <- raster(file.path(savedir, name))
    # convert to binary loss/no loss - comes as year of loss
    values(ras) <- ifelse(values(ras) > 0, 1, 0)
    
    # aggregating with mean to get proportion of pixels in set that have loss
    # don't want to go too coarse otherwise the tiles won't align properly
    ras <- aggregate_raster(ras, factor=40, aggfun="mean")
    raster::writeRaster(ras, file.path(savedir, name), overwrite=TRUE)
    
    cli_progress_update()
  }
  
  invisible(savedir)
}

#' Aggregate a raster to a desired resolution.
#'
#' Use `velox` to aggregate a raster to a desired, coarser resolution.
#' Aggregation needs to be by an integer factor, so the outcome may not
#' be the exact resolution requested.
#'
#' @param ras a raster layer object for aggregation.
#' @param desired_res the desired resolution of the raster in decimal degrees.
#' @param factor an integer specifying the desired change in resolution.
#' @param aggfun a string specifying the aggregation function, e.g. 'mean' or 'max'.
#'
#' @return the aggregated raster.
#'
aggregate_raster <- function(ras, desired_res=NULL, factor=NULL, aggfun="mean") {
  if (is.null(desired_res) & is.null(factor)) {
    cli::cli_abort(c(
      "No change in resolution specified",
      x="You need to provide a value for either {.var {desired_res}} or {.var {factor}}"
    ))
  }

  if (! is.null(desired_res) & ! is.null(factor)) {
    cli::cli_alert_warning("Values provided for {.var {desired_res} and {.var {factor}, ignoring {.var {desired_res}.")
  }

  if (! is.null(factor)) {
    agg_factor <- factor
  } else {
    agg_factor <- round(desired_res / res(ras)[1]) 
  }

  vx <- velox(ras)
  vx$aggregate(factor=agg_factor, aggtype=aggfun)

  vx$as.RasterLayer()
}

#' Merge raster tiles together to make one big raster.
#'
#' @param tiledir path to a directory storing the raster tiles as individual files.
#' @param tolerance tolerance in degrees for misaligned tiles.
#'
#' @return a single raster
#'
merge_tiles <- function(tiledir, tolerance=0.1) {
  cli_alert_info("Merging tiles in {.file {tiledir}}")
  raster_files <- list.files(tiledir, full.names=TRUE)
  rasters <- purrr::map(raster_files, ~raster(.x))
  merged_raster <- do.call(merge, c(rasters, tolerance=tolerance))
}






