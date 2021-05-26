#' Functions for calculating predictor values from coordinates and rasters.
#' 

#' Extract values from a locally stored raster file.
#' 
#' Loads a raster from a file, extracts and aggregates values
#' at each point according to the specified function. This uses 
#' [exactextractr::exact_extract] to quickly extract values,
#' which expects polygons. So the points must be buffered. All 
#' commonly used aggregation functions ignore missing values in 
#' the raster. See [https://github.com/isciences/exactextractr](here) 
#' to specify your own function.
#' 
#' @param points An `sf` data frame of buffered points.
#' @param raster_file The path to a raster file.
#' @param .fun The function aggregate, either as a string or a function
#'   that accepts the arguments `value` and `coverage_fraction`.
#' 
#' @return A vector with the extracted values
#' 
extract_values <- function(points, raster_file, .fun="mean") {
  r <- raster::raster(raster_file)
  r_crs <- raster::crs(r)
  points_crs <- st_crs(points)
  
  if (r_crs != points_crs) {
    points <- st_transform(points, r_crs)
  }
  
  exactextractr::exact_extract(r, points, .fun, progress=FALSE)
}


#' Buffer points and extract values from rasters.
#' 
#' This buffers a set of points and then extracts values from 
#' a set of raster files, applying the appropriate aggregation
#' functions to the values for each point. Extracting from a set of
#' rasters hopefully saves time by only buffering the points once.
#' 
#' @param points An `sf` data frame of points.
#' @param raster_files A named list of paths to raster files.
#' @param agg_functions A named list of functions to apply to values
#'   from each raster.
#' @param buffer_size The radius of buffer to apply in meters.
#' 
#' @return A named list of vectors containing the extracted values.
#' 
extract_buffered_points_ <- function(points, raster_files, agg_functions, buffer_size=5000) {
  
  # project to a projection in meters so we can apply buffer in m
  buffered_points <- 
    points %>%
    st_transform("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs") %>%
    st_buffer(dist=buffer_size)
  
  purrr::map2(raster_files, agg_functions, ~extract_values(buffered_points, .x, .fun=.y))
}


#' Buffer points and extract values from rasters in chunks.
#' 
#' Buffering too many points can consume a lot of memory,
#' so this splits a data frame of points into chunks and does
#' it to each chunk in turn.
#' 
#' @param points An `sf` data frame of points.
#' @param raster_files A named list of raster files.
#' @param agg_functions A named list of functions to aggregate
#'   the extracted values by.
#' @param buffer_size Radius of the buffer to apply in meters.
#' @param chunk_size The size of each chunk to split the points into.
#' 
#' @return A tibble of extracted values.
#' 
extract_buffered_chunks <- function(points, raster_files, agg_functions, buffer_size=5000, chunk_size=1e6) {
  n_chunks <- ceiling(nrow(points) / chunk_size)
  
  pb <- progress::progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=n_chunks,
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(x) {
    pb$tick()
    extract_buffered_points_(x, raster_files, agg_functions, 
                             buffer_size=buffer_size)
  }
  
  points %>%
    tibble::rowid_to_column(".chunk") %>%
    mutate(.chunk=ceiling(.chunk / chunk_size)) %>%
    group_by(.chunk) %>%
    group_split() %>%
    purrr::map(f) %>%
    purrr::map_dfr(as_tibble)
}

#' Calculate the Extent of Occurrence (EOO) from points.
#' 
#' @param longitude A vector of longitudes.
#' @param latitude A vector of corresponding latitudes.
#' 
#' @return A vector of EOOs in km^2.
#' 
calculate_eoo <- function(longitude, latitude) {
  data.frame(long=longitude, lat=latitude) %>%
    filter(! is.na(long), ! is.na(lat)) %>%
    rCAT::simProjWiz(trueCOGll(.)) %>%
    rCAT::EOOarea() %>%
    "/"(1e6) %>%
    abs()
}